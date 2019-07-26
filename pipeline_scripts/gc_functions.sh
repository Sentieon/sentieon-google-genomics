#!/usr/bin/env bash

check_error()
{
    if [ "$1" = "0" ]
    then
        echo "$2 completed successfully."
        return
    fi
    echo Error peforming $2.
    exit $1
}

delta_time()
{
    start=$1
    end=$2
    dt=$((end-start))
    dh=$(echo "$dt/3600" | bc)
    dt2=$(echo "$dt-3600*$dh" | bc)
    dm=$(echo "$dt2/60" | bc)
    ds=$(echo "$dt2-60*$dm" | bc)
    echo "$dh:$dm:$ds"
}

run()
{
    cmd=$1
    what=$2
    echo "Starting to run $what: $cmd"
    start=`date +"%D %T"`
    start_s=`date +%s`
    echo "$what start time: $start"
    eval "$cmd"
    check_error $? "$what"
    end=`date +"%D %T"`
    end_s=`date +%s`
    echo "$what end time: $end"
    runtime=$(delta_time $start_s $end_s)
    echo "$what runtime: $runtime"
}

transfer()
{
    src_file=$1
    dst_file=$2
    start_s=`date +%s`
    gsutil cp "$src_file" "$dst_file"
    check_error $? "Transfer $src_file to $dst_file"
    end_s=`date +%s`
    runtime=$(delta_time $start_s $end_s)
    echo "Transfer runtime: $runtime"
}

transfer_all_sites()
{
    dst_dir=$1; shift
    src_files=("$@")
    local_sites=()
    local_str=""
    for src_file in "${src_files[@]}"; do
        # File
        local_file=$dst_dir/$(basename "$src_file")
        transfer "$src_file" "$local_file"
        local_sites+=("$local_file")
        local_str+=" -k \"$local_file\" "
        # Index
        if $(test -e "${src_file}".idx) || $(gsutil -q stat "${src_file}".idx); then
            idx="${src_file}".idx
        elif $(test -e "${src_file}".tbi) || $(gsutil -q stat "${src_file}".tbi); then
            idx="${src_file}".tbi
        else
            echo "Cannot find idx for $src_file"
            exit 1
        fi
        local_idx=$dst_dir/$(basename "$idx")
        transfer "$idx" "$local_idx"
    done
}

upload_metrics()
{
    fun_var_cmd1=$1; shift
    fun_var_cmd2=$1; shift
    fun_pid=$1; shift
    fun_metrics_files=("$@")
    eval "fun_metrics_cmd1=\$$fun_var_cmd1"
    eval "fun_metrics_cmd2=\$$fun_var_cmd2"
    if [[ -n "$fun_metrics_cmd2" && -z "$fun_metrics_cmd1" && -f "${fun_metrics_files[0]}" ]]; then
        (run "$fun_metrics_cmd2" "Plotting metrics results." &&
            gsutil cp ${fun_metrics_files[@]} "$out_metrics" &&
            rm ${fun_metrics_files[@]}) &
        eval "$fun_pid=$! "
        eval "$fun_var_cmd2=''"
    fi
}

unset_none_variables()
{
    for var in "$@"; do
        if [[ $(eval echo \$$var) == "None" ]]; then
            eval "${var}=''"
        fi
    done
}

gc_setup()
{
    ## Download the Sentieon software
    curl -L https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${SENTIEON_VERSION}.tar.gz | tar -zxf - -C /opt/sentieon
    PATH=/opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}/bin:$PATH

    ## Dirs
    work=$scratch/work
    metrics_dir=$scratch/metrics
    ref_dir=$scratch/ref
    input_dir=$scratch/inputs
    bqsr_dir=$scratch/bqsr
    realign_dir=$scratch/realign
    dbsnp_dir=$scratch/dbsnp
    license_dir=$scratch/license
    mkdir -p $work $metrics_dir $ref_dir $input_dir $bqsr_dir \
        $realign_dir $dbsnp_dir $license_dir

    out_metrics=$OUTPUT_BUCKET/metrics/
    out_variants=$OUTPUT_BUCKET/variants/
    out_bam=$OUTPUT_BUCKET/aligned_reads/

    ## Make gsutil more robust to timeouts - which may occur during streaming transfer
    if [[ ! -f ~/.boto ]]; then
        echo -e '[Boto]\nhttp_socket_timeout=300' > ~/.boto
    fi

    ## Setup license information #
    cred=$license_dir/credentials.json
    project_file=$license_dir/credentials.json.project
    python /opt/sentieon/gen_credentials.py ${EMAIL:+--email $EMAIL} $cred "$SENTIEON_KEY"
    sleep 10
    if [[ -n $SENTIEON_KEY ]]; then
        export SENTIEON_AUTH_MECH=proxy_GOOGLE
    else
        export SENTIEON_AUTH_MECH=GOOGLE
    fi
    export SENTIEON_AUTH_DATA="$cred"
    read -r SENTIEON_JOB_TAG < "$project_file"
    export SENTIEON_JOB_TAG

    ## Test DNS resolution before falling back to hard-coded IP
    if nslookup gcp.sentieon.com > /dev/null; then
        lic_srvr_addr=gcp.sentieon.com
    else
        lic_srvr_addr=35.188.47.218
    fi
    export SENTIEON_LICENSE=$lic_srvr_addr:9003
    if ! $release_dir/bin/sentieon licclnt query -s $SENTIEON_LICENSE klib; then
        echo "Error. Could not validate license."
        if ping -c 1 $lic_srvr_addr; then
            echo "License server reachable"
        else
            echo "License server unreachable"
        fi
        if $release_dir/bin/sentieon licclnt ping; then
            echo "licclnt ping success"
        else
            echo "licclnt ping failed"
        fi
        exit 1
    fi
}

download_bams()
{
    bams=$1
    dest_arr=$2
    download_input_dir=$3

    IFS=',' read -r -a bams <<< "$bams"
    tmp_bam_dest=()
    for bam in "${bams[@]}"; do
        local_bam=$download_input_dir/$(basename "$bam")
        transfer "$bam" "$local_bam"
        if $(test -e "${bam}".bai) || $(gsutil -q stat "${bam}".bai); then
            bai="${bam}".bai
        elif $(test -e "${bam%%.bam}".bai) || $(gsutil -q stat "${bam%%.bam}".bai); then
            bai="${bam%%.bam}".bai
        else
            echo "Cannot find the index file for $bam"
            exit 1
        fi
        local_bai=$download_input_dir/$(basename "$bai")
        transfer "$bai" "$local_bai"
        tmp_bam_dest+=("$local_bam")
    done

    eval "${dest_arr}=(${tmp_bam_dest})"
}

download_intervals()
{
    if [[ -n "$INTERVAL" ]]; then
        interval="$INTERVAL"
        interval_list=""
    elif [[ -n "$INTERVAL_FILE" ]]; then
        local_interval_file=$input_dir/$(basename "$INTERVAL_FILE")
        transfer "$INTERVAL_FILE" "$local_interval_file"
        interval="$local_interval_file"
        interval_list="$local_interval_file"
    else
        interval=""
        interval_list=""
    fi
}

download_reference()
{
    ref=$ref_dir/$(basename "$REF")
    transfer "$REF" "$ref"
    transfer "${REF}".fai "${ref}".fai
    if $(test -e "${REF}".dict) || $(gsutil -q stat "${REF}".dict); then
        transfer "${REF}".dict "${ref}".dict
    elif $(test -e "${REF%%.fa}".dict) || $(gsutil -q stat "${REF%%.fa}".dict); then
        transfer "${REF%%.fa}".dict "${ref%%.fa}".dict
    elif $(test -e "${REF%%.fasta}".dict) || $(gsutil -q stat "${REF%%.fasta}".dict); then
        transfer "${REF%%.fasta}".dict "${ref%%.fasta}".dict
    else
        echo "Cannot find reference dictionary"
        exit 1
    fi
    if [[ -n "$FQ1" || -n "$TUMOR_FQ1" ]]; then
        if $(test -e "${REF}".64.amb) || $(gsutil -q stat "${REF}".64.amb); then
            middle=".64"
        elif $(test -e "${REF}".amb) || $(gsutil -q stat "${REF}".amb); then
            middle=""
        else
            echo "Cannot file BWA index files"
            exit 1
        fi
        transfer "${REF}"${middle}.amb "${ref}"${middle}.amb
        transfer "${REF}"${middle}.ann "${ref}"${middle}.ann
        transfer "${REF}"${middle}.bwt "${ref}"${middle}.bwt
        transfer "${REF}"${middle}.pac "${ref}"${middle}.pac
        transfer "${REF}"${middle}.sa  "${ref}"${middle}.sa
        if $(test -e "${REF}"${middle}.alt) || $(gsutil -q stat "${REF}"${middle}.alt); then
            transfer "${REF}"${middle}.alt "${ref}"${middle}.alt
        fi
    fi
}

bwa_mem_align()
{
    fun_base=$1; shift
    fun_fq1=$1; shift
    fun_fq2=$1; shift
    fun_rgs=$1; shift
    bam_dest=$1; shift
    fun_output_ext=$1; shift
    fun_bwa_xargs=$1; shift
    fun_util_sort_xargs=$1; shift
    fun_bam_dest=()

    IFS=',' read -r -a fun_fq1 <<< "$fun_fq1"
    IFS=',' read -r -a fun_fq2 <<< "$fun_fq2"
    IFS=',' read -r -a fun_rgs <<< "$fun_rgs"

    if [[ -z "$bwt_max_mem" ]]; then
        mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
        export bwt_max_mem="$((mem_kb / 1024 / 1024 - 2))g"
    fi

    for i in $(seq 1 ${#fun_fq1[@]}); do
        i=$((i - 1))
        fq1=${fun_fq1[$i]}
        fq2=${fun_fq2[$i]}
        readgroup=${fun_rgs[$i]}
        bwa_cmd="$release_dir/bin/bwa mem ${fun_bwa_xargs} -K 10000000 -M -R \"${readgroup}\" -t $nt \"$ref\" "
        if [[ -n "$STREAM_INPUT" ]]; then
            bwa_cmd="$bwa_cmd <(gsutil cp $fq1 -) "
            if [[ -n "$fq2" ]]; then
                bwa_cmd="$bwa_cmd <(gsutil cp $fq2 -) "
            fi
        else
            local_fq1=$input_dir/$(basename "$fq1")
            transfer "$fq1" "$local_fq1"
            bwa_cmd="$bwa_cmd \"$local_fq1\""
            if [[ -n "$fq2" ]]; then
                local_fq2=$input_dir/$(basename "$fq2")
                transfer "$fq2" "$local_fq2"
                bwa_cmd="$bwa_cmd \"$local_fq2\""
            fi
        fi
        local_bam=$work/${fun_base}sorted_bam-${i}.${fun_output_ext}
        bwa_log=$work/${fun_base}_bwa.log
        bwa_cmd="$bwa_cmd 2>$bwa_log | $release_dir/bin/sentieon util sort ${fun_util_sort_xargs} --block_size 512M -o $local_bam -t $nt --sam2bam -i -"
        run "$bwa_cmd" "BWA-mem and sorting"
        gsutil cp $bwa_log "$out_bam"
        fun_bam_dest+=($local_bam)
    done
    echo "BWA ended"

    # Cleanup
    for i in $(seq 1 ${#fun_fq1[@]}); do
        i=$((i - 1))
        fq1=${fun_fq1[$i]}
        fq2=${fun_fq2[$i]}
        local_fq1=$input_dir/$(basename "$fq1")
        local_fq2=$input_dir/$(basename "$fq2")
        if [[ -f "$local_fq1" ]]; then
            rm $local_fq1 &
        fi
        if [[ -f "$local_fq2" ]]; then
            rm $local_fq2 &
        fi
    done

    eval "${bam_dest}=(${fun_bam_dest[@]})"
}

build_metrics_cmd()
{
    fun_base=$1; shift
    cmd1_dest=$1; shift
    cmd2_dest=$1; shift
    fun_metrics_files=$1; shift
    mq=$metrics_dir/${fun_base}MeanQualityByCycle_metrics.txt
    is=$metrics_dir/${fun_base}InsertSize_metrics.txt
    qd=$metrics_dir/${fun_base}QualDistribution_metrics.txt
    gc=$metrics_dir/${fun_base}GCBias_metrics.txt
    gc_summary=$metrics_dir/${fun_base}GCBias_summary.txt
    as=$metrics_dir/${fun_base}AlignmentStat_metrics.txt
    report=$metrics_dir/${fun_base}Report_metrics.pdf

    eval "${cmd1_dest}=\"--algo MeanQualityByCycle $mq --algo QualDistribution $qd --algo GCBias --summary $gc_summary $gc --algo AlignmentStat $as --algo InsertSizeMetricAlgo $is\""
    eval "${cmd2_dest}=\"$release_dir/bin/sentieon plot metrics -o $report gc=$gc qd=$qd mq=$mq isize=$is\""

    eval "$fun_metrics_files=($mq $is $qd $gc $gc_summary $as $report)"
}

run_mark_duplicates() {
    fun_base=$1; shift
    fun_dedup=$1; shift
    fun_metrics_cmd=$1; shift
    fun_bam_str=$1; shift
    fun_bam_str_dest=$1; shift
    fun_bams_dest=$1; shift
    fun_dedup_xargs=$1; shift
    fun_output_ext=$1; shift
    fun_local_bams=("$@")

    if [[ "$fun_dedup" == "nodup" || (-z "$fun_bam_str" || -z "${fun_local_bams[@]}") ]]; then
        eval "${fun_bam_str_dest}=\"$fun_bam_str\""
        eval "${fun_bams_dest}=(\"${fun_local_bams[@]}\")"
    else
        # LocusCollector
        cmd="$release_dir/bin/sentieon driver --traverse_param=200000/10000 $fun_bam_str -t $nt -r \"$ref\" --algo LocusCollector $work/${fun_base}score.txt"
        if [[ -n $(eval "echo \$$fun_metrics_cmd") ]]; then
            eval "cmd+=\" \$$fun_metrics_cmd \""
            eval "$fun_metrics_cmd=''"
        fi
        run "$cmd" "Locus collector"

        # Dedup
        dedup_bam=$work/${fun_base}dedup.${fun_output_ext}
        if [[ "$fun_dedup" != "markdup" ]]; then
            fun_dedup_xargs+=" --rmdup "
        fi
        cmd="$release_dir/bin/sentieon driver -r \"$ref\" --traverse_param=200000/10000 $fun_bam_str -t $nt --algo Dedup ${fun_dedup_xargs} --score_info $work/${fun_base}score.txt --metrics $metrics_dir/${fun_base}dedup_metrics.txt $dedup_bam"
        run "$cmd" $fun_dedup
        eval "${fun_bam_str_dest}=\" -i $dedup_bam \""
        eval "${fun_bams_dest}=(${dedup_bam})"
    fi
}

run_bqsr()
{
    fun_base=$1; shift
    fun_bam_str=$1; shift
    fun_metrics_cmd=$1; shift
    fun_bqsr2=$1; shift
    fun_bqsr3=$1; shift
    fun_bqsr4=$1; shift
    fun_table=$1; shift
    fun_plot=$1; shift

    fun_bqsr_cmd2=
    fun_bqsr_cmd3=
    fun_bqsr_cmd4=
    csv=$work/${fun_base}recal.csv
    plot=$work/${fun_base}bqsr_report.pdf
    fun_bqsr_table=
    if [[ -n "$bqsr_sites" && -n "$fun_bam_str" ]]; then
        fun_bqsr_table=$work/${fun_base}recal_data.table
        fun_bqsr_post=$work/${fun_base}recal_data.table.post
        cmd="$release_dir/bin/sentieon driver ${interval:+--interval \"$interval\"} -t $nt -r \"$ref\" $fun_bam_str --algo QualCal $bqsr_sites $fun_bqsr_table"
        if [[ -n $(eval "echo \$$fun_metrics_cmd") ]]; then
            eval "cmd+=\" \$$fun_metrics_cmd \""
            eval "$fun_metrics_cmd=''"
        fi
        run "$cmd" "BQSR"
        if [[ -z "$NO_METRICS" ]]; then
            fun_bqsr_cmd2="--algo QualCal $bqsr_sites $fun_bqsr_post"
            fun_bqsr_cmd3="$release_dir/bin/sentieon driver --algo QualCal --plot --before $fun_bqsr_table --after $fun_bqsr_post $csv"
            fun_bqsr_cmd4="$release_dir/bin/sentieon plot bqsr -o $plot $csv"
        fi
    fi

    eval "$fun_bqsr2='$fun_bqsr_cmd2'"
    eval "$fun_bqsr3='$fun_bqsr_cmd3'"
    eval "$fun_bqsr4='$fun_bqsr_cmd4'"
    eval "$fun_table='$fun_bqsr_table'"
    eval "$fun_plot='$plot'"
}

run_bqsr_post()
{
    fun_bam_str=$1; shift
    fun_bqsr2=$1; shift
    fun_bqsr3=$1; shift
    fun_bqsr4=$1; shift
    fun_table=$1; shift
    fun_plot=$1; shift
    fun_upload_pid=$1; shift

    eval "fun_bqsr_cmd2=\$$fun_bqsr2"
    eval "fun_bqsr_cmd3=\$$fun_bqsr3"
    eval "fun_bqsr_cmd4=\$$fun_bqsr4"

    if [[ -n "$fun_bqsr_cmd2" && -n "$fun_bam_str" && -f "$fun_table" ]]; then
        cmd="$release_dir/bin/sentieon driver ${interval:+--interval \"$interval\"} -t $nt -r $ref $fun_bam_str -q $fun_table $fun_bqsr_cmd2"
        run "$cmd" "BQSR post"
        run "$fun_bqsr_cmd3" "BQSR CSV"
        run "$fun_bqsr_cmd4" "BQSR plot"
        gsutil cp $fun_plot "$out_metrics" &
        eval "$fun_upload_pid=$1 "
    fi

    eval "$fun_bqsr2=\"\""
    eval "$fun_bqsr3=\"\""
    eval "$fun_bqsr4=\"\""
}

generate_nondecoy_bed()
{
    fun_reference_fai=$1; shift
    fun_output_dest=$1; shift

    grep -v "hs37d5\|chrEBV\|hs38d1\|decoy" "$fun_reference_fai" | awk 'BEGIN{OFS="\t"} {print $1,0,$2}' > "$fun_output_dest"
}

