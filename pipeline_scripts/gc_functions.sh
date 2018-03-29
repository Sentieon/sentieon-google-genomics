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
    gsutil cp $src_file $dst_file
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
        local_file=$dst_dir/$(basename $src_file)
        transfer $src_file $local_file
        local_sites+=($local_file)
        local_str+=" -k $local_file "
        # Index
        if $(test -e ${src_file}.idx) || $(gsutil -q stat ${src_file}.idx); then
            idx=${src_file}.idx
        elif $(test -e ${src_file}.tbi) || $(gsutil -q stat ${src_file}.tbi); then
            idx=${src_file}.tbi
        else
            echo "Cannot find idx for $src_file"
            exit 1
        fi
        local_idx=$dst_dir/$(basename $idx)
        transfer $idx $local_idx
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
            gsutil cp ${fun_metrics_files[@]} $out_metrics &&
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

    ## Setup license information #
    cred=$license_dir/credentials.json
    python /opt/sentieon/gen_credentials.py $cred &
    credentials_pid=$!
    sleep 10
    export SENTIEON_AUTH_MECH=GOOGLE SENTIEON_AUTH_DATA="$cred"
    export SENTIEON_LICENSE=gcp.sentieon.com:9003 # Hard-coded for now
    if ! $release_dir/libexec/licclnt query -s $SENTIEON_LICENSE klib; then
        echo "Error. Could not validate license."
        exit 1
    fi
}

download_bams()
{
    bams=$1
    dest_arr=$2
    download_input_dir=$3

    bams=($(echo $bams | tr ',' ' '))
    tmp_bam_dest=()
    for bam in ${bams[@]}; do
        local_bam=$download_input_dir/$(basename $bam)
        transfer $bam $local_bam
        if $(test -e ${bam}.bai) || $(gsutil -q stat ${bam}.bai); then
            bai=${bam}.bai
        elif $(test -e ${bam%%.bam}.bai) || $(gsutil -q stat ${bam%%.bam}.bai); then
            bai=${bam%%.bam}.bai
        else
            echo "Cannot find the index file for $bam"
            exit 1
        fi
        local_bai=$download_input_dir/$(basename $bai)
        transfer $bai $local_bai
        tmp_bam_dest+=($local_bam)
    done

    eval "${dest_arr}=(${tmp_bam_dest})"
}

download_intervals()
{
    if [[ -n "$INTERVAL" ]]; then
        interval=" --interval $INTERVAL "
        interval_list=""
    elif [[ -n "$INTERVAL_FILE" ]]; then
        local_interval_file=$input_dir/$(basename $INTERVAL_FILE)
        transfer $INTERVAL_FILE $local_interval_file
        interval=" --interval $local_interval_file "
        interval_list=" --interval_list $local_interval_file "
    else
        interval=""
        interval_list=""
    fi
}

download_reference()
{
    ref=$ref_dir/$(basename $REF)
    transfer $REF $ref
    transfer ${REF}.fai ${ref}.fai
    if $(test -e ${REF}.dict) || $(gsutil -q stat ${REF}.dict); then
        transfer ${REF}.dict ${ref}.dict
    elif $(test -e ${REF%%.fa}.dict) || $(gsutil -q stat ${REF%%.fa}.dict); then
        transfer ${REF%%.fa}.dict ${ref%%.fa}.dict
    else
        echo "Cannot find reference dictionary"
        exit 1
    fi
    if [[ -n "$FQ1" ]]; then
        transfer ${REF}.amb ${ref}.amb
        transfer ${REF}.ann ${ref}.ann
        transfer ${REF}.bwt ${ref}.bwt
        transfer ${REF}.pac ${ref}.pac
        transfer ${REF}.sa ${ref}.sa
    fi
}

bwa_mem_align()
{
    fun_base=$1; shift
    fun_fq1=$1; shift
    fun_fq2=$1; shift
    fun_rgs=$1; shift
    bam_dest=$1; shift
    fun_bam_dest=()

    fun_fq1=($(echo $fun_fq1 | tr ',' ' '))
    fun_fq2=($(echo $fun_fq2 | tr ',' ' '))
    fun_rgs=($(echo $fun_rgs | tr ',' ' '))

    for i in $(seq 1 ${#fun_fq1[@]}); do
        i=$((i - 1))
        fq1=${fun_fq1[$i]}
        fq2=${fun_fq2[$i]}
        readgroup=${fun_rgs[$i]}
        bwa_cmd="$release_dir/bin/bwa mem -M -R \"${readgroup}\" -t $nt $ref "
        if [[ -n "$STREAM_INPUT" ]]; then
            bwa_cmd="$bwa_cmd <(gsutil cp $fq1 -) "
            if [[ -n "$fq2" ]]; then
                bwa_cmd="$bwa_cmd <(gsutil cp $fq2 -) "
            fi
        else
            local_fq1=$input_dir/$(basename $fq1)
            transfer $fq1 $local_fq1
            bwa_cmd="$bwa_cmd $local_fq1"
            if [[ -n "$fq2" ]]; then
                local_fq2=$input_dir/$(basename $fq2)
                transfer $fq2 $local_fq2
                bwa_cmd="$bwa_cmd $local_fq2"
            fi
        fi
        local_bam=$work/${fun_base}sorted_bam-${i}.bam
        bwa_cmd="$bwa_cmd | $release_dir/bin/sentieon util sort -o $local_bam -t $nt --sam2bam -i -"
        run "$bwa_cmd" "BWA-mem and sorting"
        fun_bam_dest+=($local_bam)
    done
    echo "BWA ended"

    # Cleanup
    for i in $(seq 1 ${#fun_fq1[@]}); do
        i=$((i - 1))
        fq1=${fun_fq1[$i]}
        fq2=${fun_fq2[$i]}
        local_fq1=$input_dir/$(basename $fq1)
        local_fq2=$input_dir/$(basename $fq2)
        if [[ -f "$local_fq1" ]]; then
            rm $fq1 &
        fi
        if [[ -f "$local_fq2" ]]; then
            rm $fq2 &
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
    fun_local_bams=("$@")

    if [[ "$fun_dedup" == "nodup" || (-z $fun_bam_str || -z "${fun_local_bams[@]}") ]]; then
        eval "${fun_bam_str_dest}=\"$fun_bam_str\""
        eval "${fun_bams_dest}=(${fun_local_bams[@]})"
    else
        # LocusCollector
        cmd="$release_dir/bin/sentieon driver $fun_bam_str -t $nt -r $ref --algo LocusCollector $work/${fun_base}score.txt"
        if [[ -n $(eval "echo \$$fun_metrics_cmd") ]]; then
            eval "cmd+=\" \$$fun_metrics_cmd \""
            eval "$fun_metrics_cmd=''"
        fi
        run "$cmd" "Locus collector"

        # Dedup
        dedup_bam=$work/${fun_base}dedup.bam
        if [ "$fun_dedup" = "markdup" ]; then
            cmd="$release_dir/bin/sentieon driver $fun_bam_str -t $nt --algo Dedup --score_info $work/${fun_base}score.txt --metrics $metrics_dir/${fun_base}dedup_metrics.txt $dedup_bam"
        else
            cmd="$release_dir/bin/sentieon driver $fun_bam_str -t $nt --algo Dedup --score_info $work/${fun_base}score.txt --metrics $metrics_dir/${fun_base}dedup_metrics.txt --rmdup $dedup_bam"
        fi
        run "$cmd" $fun_dedup
        for bam in ${local_bams[@]}; do
            if [[ -n $bam ]]; then
                rm $bam &
            fi
            if [[ -n ${bam}.bai ]]; then
                rm ${bam}.bai &
            fi
        done
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
        cmd="$release_dir/bin/sentieon driver $interval -t $nt -r '$ref' $fun_bam_str --algo QualCal $bqsr_sites $fun_bqsr_table"
        if [[ -n $(eval "echo \$$fun_metrics_cmd") ]]; then
            eval "cmd+=\" \$$fun_metrics_cmd \""
            eval "$fun_metrics_cmd=''"
        fi
        run "$cmd" "BQSR"
        if [[ -z "$NO_METRICS" ]]; then
            bqsr_cmd2="--algo QualCal $bqsr_sites $fun_bqsr_post"
            bqsr_cmd3="$release_dir/bin/sentieon driver --algo QualCal --plot --before $fun_bqsr_table --after $fun_bqsr_post $csv"
            fun_bqsr_cmd4="$release_dir/bin/sentieon plot bqsr -o $plot $csv"
        fi
    fi

    eval "$fun_bqsr2=\"$fun_bqsr_cmd2\""
    eval "$fun_bqsr3=\"$fun_bqsr_cmd3\""
    eval "$fun_bqsr4=\"$fun_bqsr_cmd4\""
    eval "$fun_table=\"$fun_bqsr_table\""
    eval "$fun_plot=\"$plot\""
}


