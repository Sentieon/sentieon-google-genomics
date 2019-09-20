#!/usr/bin/env bash

set -xveo pipefail
set +H

BASEDIR=$(dirname "$0")
scratch=/mnt/work
nt=$(nproc)
source $BASEDIR/gc_functions.sh

# Set "None" variables to an empty string
environmental_variables=(FQ1 FQ2 BAM OUTPUT_BUCKET REF READGROUP BQSR_SITES \
    DBSNP INTERVAL INTERVAL_FILE NO_METRICS NO_BAM_OUTPUT NO_HAPLOTYPER \
    GVCF_OUTPUT STREAM_INPUT PIPELINE OUTPUT_CRAM_FORMAT SENTIEON_KEY EMAIL \
    SENTIEON_VERSION CALLING_ARGS DNASCOPE_MODEL CALLING_ALGO)
unset_none_variables ${environmental_variables[@]}
OUTPUT_CRAM_FORMAT="" # Not yet supported

readonly FQ1 FQ2 BAM OUTPUT_BUCKET REF READGROUP BQSR_SITES DBSNP INTERVAL \
    INTERVAL_FILE NO_METRICS NO_BAM_OUTPUT NO_HAPLOTYPER GVCF_OUTPUT \
    STREAM_INPUT PIPELINE OUTPUT_CRAM_FORMAT SENTIEON_KEY EMAIL \
    SENTIEON_VERSION CALLING_ARGS DNASCOPE_MODEL CALLING_ALGO

release_dir="/opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}/"

# Basic error handling #
if [[ -n "$FQ1" && -n "$BAM" ]]; then
    echo "Please supply either fastq or BAM files"
    exit 1
fi

if [[ -n "$INTERVAL" && -n "$INTERVAL_FILE" ]]; then
    echo "Please supply either an INTERVAL or INTERVAL_FILE, not both"
    exit 1
fi

if [[ -n "$NO_BAM_OUTPUT" && -n "$NO_HAPLOTYPER" && -n "$NO_METRICS" ]]; then
    echo "Nothing to do"
    exit 1
fi

# **********************************
# 0. Setup
# **********************************
gc_setup
export LD_PRELOAD=${release_dir}/lib/libjemalloc.so
export MALLOC_CONF=lg_dirty_mult:-1

## Download input files
if [[ -n "$BAM" ]]; then
    download_bams "$BAM" local_bams $input_dir
else
    local_bams=()
fi

download_intervals
download_reference
if [[ $CALLING_ALGO == "DNAscope" && -n "$DNASCOPE_MODEL" ]]; then
    curl -L -o ${input_dir}/dnascope.model "$DNASCOPE_MODEL"
fi

## Handle the sites files
IFS=',' read -r -a gs_bqsr_sites <<< "$BQSR_SITES"
transfer_all_sites $bqsr_dir "${gs_bqsr_sites[@]}"
local_bqsr_sites=("${local_sites[@]}")
bqsr_sites="$local_str"

if [[ -n "$DBSNP" ]]; then
    transfer_all_sites $dbsnp_dir "$DBSNP"
    dbsnp="${local_sites[0]}"
fi

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
output_ext="bam"
if [[ -n $FQ1 ]]; then
    bwa_mem_align "" "$FQ1" "$FQ2" "$READGROUP" local_bams $output_ext "-K 100000000 -Y" "$util_sort_xargs" "true"
fi

local_bams_str=""
for bam in "${local_bams[@]}"; do
    local_bams_str+=" -i \"$bam\" "
done

# ******************************************
# 2. Metrics command
# ******************************************
metrics_cmd1=
metrics_cmd2=
metrics_files=
if [[ -z "$NO_METRICS" ]]; then
    build_metrics_cmd "" metrics_cmd1 metrics_cmd2 metrics_files
fi

# ******************************************
# 3. Remove duplicates
# ******************************************
output_ext="bam"

run_mark_duplicates "" "markdup" metrics_cmd1 "$local_bams_str" dedup_bam_str dedup_bams "$dedup_xargs" $output_ext "true" "${local_bams[@]}"
if [[ -z "$NO_METRICS" ]]; then
    (gsutil cp $metrics_dir/dedup_metrics.txt "$out_metrics" &&
        rm $metrics_dir/dedup_metrics.txt) &
    upload_dedup_pid=$!
else
    rm $metrics_dir/dedup_metrics.txt &
fi
for bam in "${local_bams[@]}"; do
    if [[ -f "$bam" ]]; then
        rm "$bam" &
    fi
    if [[ -f "${bam}".bai ]]; then
        rm "${bam}".bai &
    fi
    if [[ -f "${bam}".crai ]]; then
        rm "${bam}".crai &
    fi
done

upload_metrics metrics_cmd1 metrics_cmd2 upload_metrics_pid ${metrics_files[@]}

# ******************************************
# 4. Base recalibration
# ******************************************
bqsr_intervals="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
cmd="$release_dir/bin/sentieon driver -t $nt -r \"$ref\" --interval $bqsr_intervals $dedup_bam_str --algo QualCal $bqsr_sites $work/recal_data.table"
run "$cmd" "BQSR"

cmd="$release_dir/bin/sentieon driver -t $nt -r \"$ref\" $dedup_bam_str --read_filter QualCalFilter,table=$work/recal_data.table,prior=-1.0,indel=false,levels=10/20/30,min_qual=6 --algo ReadWriter --cram_write_options version=3.0,compressor=gzip+rans $work/recalibrated.cram"
run "$cmd" "ReadWriter"

for bam in "${dedup_bams[@]}"; do
    if [[ -f "$bam" ]]; then
        rm "$bam" &
    fi
    if [[ -f "${bam}".bai ]]; then
        rm "$bam".bai &
    fi
    if [[ -f "${bam}".crai ]]; then
        rm "${bam}".crai &
    fi
done

if [[ -z "$NO_BAM_OUTPUT" ]]; then
    gsutil cp $work/recalibrated.cram $work/recalibrated.cram.crai "$out_bam" &
    upload_recal_pid=$!
fi


# ******************************************
# 5. Variant Calling
# ******************************************

## Generate a non-decoy BED file
generate_nondecoy_bed "${ref}".fai "${ref}"_nondecoy.bed
call_interval="$interval"
call_interval=${call_interval:-"${ref}"_nondecoy.bed}

outgvcf=$work/hc.g.vcf.gz
outvcf=$work/hc.vcf.gz

algo=Haplotyper
extra_vcf_args=""
if [[ $CALLING_ALGO == "DNAscope" ]]; then
    algo="DNAscope"
    if [[ -n "$DNASCOPE_MODEL" ]]; then
        extra_vcf_args="--model ${input_dir}/dnascope.model"
    else
        extra_vcf_args="--var_type snp,indel,bnd"
    fi
    outvcf=$work/dnascope.vcf.gz
    outgvcf=$work/dnascope.g.vcf.gz
    tmpvcf=$work/tmp.vcf.gz
fi

if [[ -z $NO_HAPLOTYPER ]]; then
    if [[ -n "$GVCF_OUTPUT" ]]; then
        cmd="$release_dir/bin/sentieon driver --interval \"$call_interval\" -t $nt -r \"$ref\" -i $work/recalibrated.cram --algo $algo $CALLING_ARGS ${dbsnp:+-d \"$dbsnp\"} --emit_mode gvcf ${outgvcf}"
        outfile=$outgvcf
    else
        cmd="$release_dir/bin/sentieon driver --interval \"$call_interval\" -t $nt -r \"$ref\" -i $work/recalibrated.cram --algo $algo $CALLING_ARGS $extra_vcf_args ${dbsnp:+-d \"$dbsnp\"} ${outvcf}"
        outfile=$outvcf
    fi
    run "$cmd" "Variant calling"

    # DNAscope SV calling
    if [[ $CALLING_ALGO == "DNAscope" && -z "$GVCF_OUTPUT" && -z "$DNASCOPE_MODEL" ]]; then
        mv $outvcf $tmpvcf
        mv ${outvcf}.tbi ${tmpvcf}.tbi
        cmd="$release_dir/bin/sentieon driver -t $nt -r \"$ref\" --algo SVSolver -v $tmpvcf $outvcf"
        run "$cmd" "SVSolver"
    fi

    # DNAscope model apply
    if [[ $CALLING_ALGO == "DNAscope" && -z "$GVCF_OUTPUT" && -n "$DNASCOPE_MODEL" ]]; then
        mv $outvcf $tmpvcf
        mv ${outvcf}.tbi ${tmpvcf}.tbi
        cmd="$release_dir/bin/sentieon driver -t $nt -r \"$ref\" --algo DNAModelApply --model ${input_dir}/dnascope.model -v $tmpvcf $outvcf"
        run "$cmd" "DNAscope model apply"
    fi

    gsutil cp $outfile ${outfile}.tbi "$out_variants" &
    upload_vcf_pid=$!
fi

# Wait for all running uploads to finish
set +e
if [[ -n $upload_recal_pid ]]; then
    wait $upload_recal_pid
fi
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
