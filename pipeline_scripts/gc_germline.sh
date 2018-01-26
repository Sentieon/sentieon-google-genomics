#!/usr/bin/env bash
set -xv

# TODO
#  - vqsr
#  - add DNAscope

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
    if [[ -n $metrics_cmd2 && -z $metrics_cmd1 ]]; then
        (run "$metrics_cmd2" "Plotting metrics results." &&
            gsutil cp $mq $qd $gc $gc_summary $as $is $report $out_metrics &&
            rm $mq $qd $gc $gc_summary $as $is $report) &
        upload_metrics_pid=$!
        metrics_cmd2=
    fi
}

BASEDIR=$(dirname "$0")
version="201711.01"
release_dir="/opt/sentieon/sentieon-genomics-${version}/"
platform="ILLUMINA"
scratch=/mnt/work
nt=$(nproc)

# Set "None" variables to an empty string
# collected from the YAML file
environmental_variables=(FQ1 FQ2 BAM OUTPUT_BUCKET REF SENTIEON_TOKEN \
    GOOGLE_TOKEN READGROUP DEDUP BQSR_SITES DBSNP INTERVAL \
    INTERVAL_FILE NO_METRICS NO_BAM_OUTPUT NO_HAPLOTYPER GVCF_OUTPUT \
    STREAM_INPUT)
for var in "${environmental_variables[@]}"; do
    if [[ $(eval echo \$$var) == "None" ]]; then
        eval "${var}=''"
    fi
done

readonly FQ1 FQ2 BAM OUTPUT_BUCKET REF SENTIEON_TOKEN GOOGLE_TOKEN \
    READGROUP DEDUP BQSR_SITES DBSNP INTERVAL INTERVAL_FILE NO_METRICS \
    NO_BAM_OUTPUT NO_HAPLOTYPER GVCF_OUTPUT STREAM_INPUT

# Some basic error handling #
if [[ -n "$FQ1" && -n "$BAM" ]]; then
    echo "Please supply either fastq or BAM files"
    exit 1
fi

if [[ -n "$INTERVAL" && -n "$INTERVAL_FILE" ]]; then
    echo "Please supply either an INTERVAL or INTERVAL_FILE"
    exit 1
fi

if [[ -n "$NO_BAM_OUTPUT" && -n "$NO_HAPLOTYPER" && -n "$NO_METRICS" ]]; then
    echo "Nothing to do"
    exit 1
fi

# **********************************
# 0. Setup
# **********************************
## Dirs
work=$scratch/work
metrics_dir=$scratch/metrics
ref_dir=$scratch/ref
input_dir=$scratch/inputs
bqsr_dir=$scratch/bqsr
dbsnp_dir=$scratch/dbsnp
license_dir=$scratch/license
mkdir -p $work $metrics_dir $ref_dir $input_dir $bqsr_dir \
    $dbsnp_dir $license_dir

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

## Download input files
if [[ -n "$BAM" ]]; then
    bams=($(echo $BAM | tr ',' ' '))
    if [[ -n "$STREAM_INPUT" ]]; then
        echo "Ignoring 'STREAM_INPUT' with input BAM file"
    fi
    local_bams=()
    for bam in ${bams[@]}; do
        local_bams=$input_dir/$(basename $bam)
        transfer $bam $local_bam
        if $(test -e ${bam}.bai) || $(gsutil -q stat ${bam}.bai); then
            bai=${bam}.bai
        elif $(test -e ${bam%%.bam}.bai) || $(gsutil -q $stat ${bam%%.bam}.bai); then
            bai=${bam%%.bam}.bai
        else
            echo "Cannot find the index file for $BAM"
            exit 1
        fi
        local_bai=$input_dir/$(basename $bai)
        transfer $bai $local_bai
        local_bams+=($bam)
    done
else
    local_bams=()
fi

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

## Download all reference files
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

## Handle the sites files
gs_bqsr_sites=($(echo $BQSR_SITES | tr ',' ' '))
transfer_all_sites $bqsr_dir "${gs_bqsr_sites[@]}"
local_bqsr_sites="${local_sites[@]}"
bqsr_sites="$local_str"

REALIGN_SITES=

if [[ -n "$DBSNP" ]]; then
    transfer_all_sites $dbsnp_dir $DBSNP
    dbsnp=${local_sites[0]}
fi

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
fq1_fastqs=($(echo $FQ1 | tr ',' ' '))
fq2_fastqs=($(echo $FQ2 | tr ',' ' '))
readgroups=($(echo $READGROUP | tr ',' ' '))
if [[ -n "$FQ1" ]]; then
    echo "starting bwa"
    for i in $(seq 1 ${#fq1_fastqs[@]}); do
        i=$((i - 1))
        fq1=${fq1_fastqs[$i]}
        fq2=${fq2_fastqs[$i]}
        readgroup=${readgroups[$i]}
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
        local_bam=$work/sorted_bam-${i}.bam
        bwa_cmd="$bwa_cmd | $release_dir/bin/sentieon util sort -o $local_bam -t $nt --sam2bam -i -"
        run "$bwa_cmd" "BWA-mem and sorting"
        local_bams+=($local_bam)
    done
    echo "BWA ended"
    # Cleanup
    for i in $(seq 1 ${#fq1_fastqs[@]}); do
        i=$((i - 1))
        fq1=${fq1_fastqs[$i]}
        fq2=${fq2_fastqs[$i]}
        if [[ -n "$fq1" ]]; then
            rm $fq1 &
        fi
        if [[ -n "$fq2" ]]; then
            rm $fq2 &
        fi
    done
fi

local_bams_str=""
for bam in ${local_bams[@]}; do
    local_bams_str+=" -i $bam "
done

# ******************************************
# 2. Metrics command
# ******************************************
metrics_cmd1=
metrics_cmd2=
if [[ -z "$no_metrics" ]]; then
    mq=$metrics_dir/MeanQualityByCycle_metrics.txt
    is=$metrics_dir/InsertSize_metrics.txt
    qd=$metrics_dir/QualDistribution_metrics.txt
    gc=$metrics_dir/GCBias_metrics.txt
    gc_summary=$metrics_dir/GCBias_summary.txt
    as=$metrics_dir/AlignmentStat_metrics.txt
    report=$metrics_dir/Report_metrics.pdf

    metrics_cmd1="--algo MeanQualityByCycle $mq --algo QualDistribution $qd --algo GCBias --summary $gc_summary $gc --algo AlignmentStat $as --algo InsertSizeMetricAlgo $is"
    metrics_cmd2="$release_dir/bin/sentieon plot metrics -o $report gc=$gc qd=$qd mq=$mq isize=$is"
fi

# ******************************************
# 3. Remove duplicates
# ******************************************

if [[ "$DEDUP" == "nodup" ]]; then
    dedup_bam_str=$local_bams_str
    dedup_bams=(${local_bams[@]})
else
    # LocusCollector
    cmd="$release_dir/bin/sentieon driver $local_bams_str -t $nt -r $ref --algo LocusCollector $work/score.txt"
    if [[ -n $metrics_cmd1 ]]; then
        cmd+=" $metrics_cmd1"
        metrics_cmd1=
    fi
    run "$cmd" "Locus collector"

    # Dedup
    dedup_bam=$work/dedup.bam
    if [ "$DEDUP" = "markdup" ]; then
        cmd="$release_dir/bin/sentieon driver $local_bams_str -t $nt --algo Dedup --score_info $work/score.txt --metrics $metrics_dir/dedup_metrics.txt $dedup_bam"
    else
        cmd="$release_dir/bin/sentieon driver $local_bams_str -t $nt --algo Dedup --score_info $work/score.txt --metrics $metrics_dir/dedup_metrics.txt --rmdup $dedup_bam"
    fi
    run "$cmd" $DEDUP
    if [[ -z "$no_metrics" ]]; then
        (gsutil cp $metrics_dir/dedup_metrics.txt $out_metrics &&
            rm $metrics_dir/dedup_metrics.txt) &
        upload_dedup_pid=$!
    else
        rm $metrics_dir/dedup_metrics.txt &
    fi
    for bam in ${local_bams[@]}; do
        if [[ -n $bam ]]; then
            rm $bam &
        fi
        if [[ -n ${bam}.bai ]]; then
            rm ${bam}.bai &
        fi
    done
    dedup_bam_str=" -i $dedup_bam "
    dedup_bams=($dedup_bam)
fi

if [[ -z $NO_BAM_OUTPUT ]]; then
    for bam in ${dedup_bams[@]}; do
        gsutil cp $bam ${bam}.bai $out_bam
        upload_deduped_pid=$!
    done
fi

upload_metrics

# ******************************************
# 4. Base recalibration
# ******************************************
bqsr_table=''
bqsr_cmd2=
bqsr_cmd3=
csv=$work/recal.csv
plot=$work/bqsr_report.pdf
if [[ -n "$bqsr_sites" ]]; then
    bqsr_table=$work/recal_data.table
    bqsr_post=$work/recal_data.table.post
    cmd="$release_dir/bin/sentieon driver $interval -t $nt -r '$ref' $dedup_bam_str --algo QualCal $bqsr_sites $bqsr_table"
    if [[ -n $metrics_cmd1 ]]; then
        cmd+=" $metrics_cmd1"
    fi
    run "$cmd" "BQSR"
    if [[ -z "$no_metrics" ]]; then
        bqsr_cmd2="--algo QualCal $bqsr_sites $bqsr_post"
        bqsr_cmd3="$release_dir/bin/sentieon driver --algo QualCal --plot --before $bqsr_table --after $bqsr_post $csv"
        bqsr_cmd4="$release_dir/bin/sentieon plot bqsr -o $plot $csv"
    fi
    if [[ -z $NO_BAM_OUTPUT ]]; then
        gsutil cp $bqsr_table $out_bam
        upload_bqsr_pid=$!
    fi
fi

upload_metrics

# ******************************************
# Variant Calling
# ******************************************

# ******************************************
# 5b. HC Variant caller
# ******************************************
outgvcf=$work/hc.g.vcf.gz
outvcf=$work/hc.vcf.gz

if [[ -z $NO_HAPLOTYPER ]]; then
    if [[ -n "$GVCF_OUTPUT" ]]; then
        cmd="$release_dir/bin/sentieon driver $interval -t $nt -r '$ref' $dedup_bam_str ${bqsr_table:+-q $bqsr_table} --algo Haplotyper ${dbsnp:+-d $dbsnp} --emit_mode gvcf ${outgvcf}"
        outfile=$outgvcf
    else
        cmd="$release_dir/bin/sentieon driver $interval -t $nt -r '$ref' $dedup_bam_str ${bqsr_table:+-q $bqsr_table} --algo Haplotyper ${dbsnp:+-d $dbsnp} ${outvcf}"
        outfile=$outvcf
    fi

    if [[ -n $metrics_cmd1 ]]; then
        cmd+=" $metrics_cmd1"
        metrics_cmd1=
    fi
    if [[ -n $bqsr_cmd2 ]]; then
        cmd+=" $bqsr_cmd2"
        bqsr_cmd2=
    fi

    run "$cmd" "Haplotyper variant calling"
    gsutil cp $outfile ${outfile}.tbi $out_variants &
    upload_vcf_pid=$!
fi

if [[ -n $metrics_cmd1 ]]; then
    cmd="$release_dir/bin/sentieon driver $interval -t $nt -r '$ref' $dedup_bam_str ${bqsr_table:+-q $bqsr_table}"
    cmd+=" $metrics_cmd1"
    run "$cmd" "Metrics collection"
fi

upload_metrics

if [[ -n $bqsr_cmd2 ]]; then
    cmd="$release_dir/bin/sentieon driver $interval -t $nt -r '$ref' $dedup_bam_str -q $bqsr_table $bqsr_cmd2"
    bqsr_cmd2=
    run "$cmd" "BQSR Post"
fi

if [[ -n $bqsr_cmd3 ]]; then
    run "$bqsr_cmd3" "BQSR CSV"
    run "$bqsr_cmd4" "BQSR plot"
    gsutil cp $plot $out_metrics &
    upload_bqsr_metrics_pid=$!
fi

kill $credentials_pid

# Wait for all running uploads to finish
if [[ -n $upload_metrics_pid ]]; then
    wait $upload_metrics_pid
fi
if [[ -n $upload_dedup_pid ]]; then
    wait $upload_dedup_pid
fi
if [[ -n $upload_deduped_pid ]]; then
    wait $upload_deduped_pid
fi
if [[ -n $upload_bqsr_pid ]]; then
    wait $upload_bqsr_pid
fi
if [[ -n $upload_vcf_pid ]]; then
    wait $upload_vcf_pid
fi
if [[ -n $upload_bqsr_metrics_pid ]]; then
    wait $upload_bqsr_metrics_pid
fi

exit 0
