#!/usr/bin/env bash

set -xveo pipefail
set +H

BASEDIR=$(dirname "$0")
scratch=/mnt/work
nt=$(nproc)
source $BASEDIR/gc_functions.sh

# Set "None" variables to an empty string
environmental_variables=(FQ1 FQ2 BAM OUTPUT_BUCKET REF READGROUP DEDUP \
    BQSR_SITES DBSNP INTERVAL INTERVAL_FILE NO_METRICS NO_BAM_OUTPUT \
    NO_HAPLOTYPER GVCF_OUTPUT STREAM_INPUT PIPELINE OUTPUT_CRAM_FORMAT \
    SENTIEON_KEY RECALIBRATED_OUTPUT EMAIL SENTIEON_VERSION)
unset_none_variables ${environmental_variables[@]}
OUTPUT_CRAM_FORMAT="" # Not yet supported

readonly FQ1 FQ2 BAM OUTPUT_BUCKET REF READGROUP DEDUP BQSR_SITES DBSNP \
    INTERVAL INTERVAL_FILE NO_METRICS NO_BAM_OUTPUT NO_HAPLOTYPER GVCF_OUTPUT \
    STREAM_INPUT PIPELINE OUTPUT_CRAM_FORMAT SENTIEON_KEY RECALIBRATED_OUTPUT \
    EMAIL SENTIEON_VERSION

release_dir="/opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}/"
export LD_PRELOAD=/opt/sentieon/sentieon-genomics-${version}/lib/libjemalloc.so
export MALLOC_CONF=lg_dirty_mult:-1

# Basic error handling #
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
gc_setup

## Download input files
if [[ -n "$BAM" ]]; then
    download_bams "$BAM" local_bams $input_dir
else
    local_bams=()
fi

download_intervals
download_reference

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
    bwa_mem_align "" "$FQ1" "$FQ2" "$READGROUP" local_bams $output_ext "$bwa_xargs" "$util_sort_xargs"
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

run_mark_duplicates "" "$DEDUP" metrics_cmd1 "$local_bams_str" dedup_bam_str dedup_bams "$dedup_xargs" $output_ext "${local_bams[@]}"
if [[ "$DEDUP" != "nodup" ]]; then
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
fi

if [[ -z "$NO_BAM_OUTPUT" && (-z "$bqsr_sites" || -z "$RECALIBRATED_OUTPUT" ) ]]; then
    upload_list=""
    for bam in "${dedup_bams[@]}"; do
        upload_list+=" \"$bam\" "
        if [[ -f "${bam}.bai" ]]; then
            upload_list+=" \"${bam}.bai\" "
        elif [[ -f "${bam}.crai" ]]; then
            upload_list+=" \"${bam}.crai\" "
        fi
    done
    eval gsutil cp $upload_list "$out_bam" &
    upload_deduped_pid=$!
fi

upload_metrics metrics_cmd1 metrics_cmd2 upload_metrics_pid ${metrics_files[@]}

# ******************************************
# 4. Base recalibration
# ******************************************
run_bqsr "" "$dedup_bam_str" metrics_cmd1 bqsr_cmd2 bqsr_cmd3 bqsr_cmd4 bqsr_table bqsr_plot

if [[ -n "$bqsr_sites" && -z "$NO_BAM_OUTPUT" && -n "$RECALIBRATED_OUTPUT" ]]; then
    outrecal=$work/recalibrated.bam
    cmd="$release_dir/bin/sentieon driver $dedup_bam_str -q $bqsr_table --algo ReadWriter $outrecal"
    (run "$cmd" "ReadWriter";
        gsutil cp $outrecal ${outrecal}.bai "$out_bam") &
    upload_recal_pid=$!
fi

if [[ -n "$bqsr_sites" && -z "$NO_BAM_OUTPUT" && -z "$RECALIBRATED_OUTPUT" ]]; then
    gsutil cp $bqsr_table "$out_bam" &
    upload_bqsr_pid=$!
fi

upload_metrics metrics_cmd1 metrics_cmd2 upload_metrics_pid ${metrics_files[@]}


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
if [[ $PIPELINE == "DNAscope" ]]; then
    algo="DNAscope"
    extra_vcf_args="--var_type snp,indel,bnd"
    outvcf=$work/dnascope.vcf.gz
    outgvcf=$work/dnascope.g.vcf.gz
    tmpvcf=$work/tmp.vcf.gz
fi

if [[ -z $NO_HAPLOTYPER ]]; then
    if [[ -n "$GVCF_OUTPUT" ]]; then
        cmd="$release_dir/bin/sentieon driver --interval \"$call_interval\" -t $nt -r \"$ref\" $dedup_bam_str ${bqsr_table:+-q $bqsr_table} --algo $algo ${dbsnp:+-d \"$dbsnp\"} --emit_mode gvcf ${outgvcf}"
        outfile=$outgvcf
    else
        cmd="$release_dir/bin/sentieon driver --interval \"$call_interval\" -t $nt -r \"$ref\" $dedup_bam_str ${bqsr_table:+-q $bqsr_table} --algo $algo $extra_vcf_args ${dbsnp:+-d \"$dbsnp\"} ${outvcf}"
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
    if [[ $PIPELINE == "DNAscope" && -z "$GVCF_OUTPUT" ]]; then
        mv $outvcf $tmpvcf
        mv ${outvcf}.tbi ${tmpvcf}.tbi
        cmd="$release_dir/bin/sentieon driver -t $nt -r \"$ref\" --algo SVSolver -v $tmpvcf $outvcf"
        run "$cmd" "SVSolver"
    fi
    gsutil cp $outfile ${outfile}.tbi "$out_variants" &
    upload_vcf_pid=$!
fi

if [[ -n $metrics_cmd1 ]]; then
    cmd="$release_dir/bin/sentieon driver ${interval:+--interval \"$interval\"} -t $nt -r \"$ref\" $dedup_bam_str ${bqsr_table:+-q $bqsr_table}"
    cmd+=" $metrics_cmd1"
    run "$cmd" "Metrics collection"
fi

upload_metrics metrics_cmd1 metrics_cmd2 upload_metrics_pid ${metrics_files[@]}


if [[ -n $bqsr_cmd2 ]]; then
    cmd="$release_dir/bin/sentieon driver ${interval:+--interval \"$interval\"} -t $nt -r \"$ref\" $dedup_bam_str -q $bqsr_table $bqsr_cmd2"
    bqsr_cmd2=
    run "$cmd" "BQSR Post"
fi

if [[ -n $bqsr_cmd3 ]]; then
    run "$bqsr_cmd3" "BQSR CSV"
    run "$bqsr_cmd4" "BQSR plot"
    gsutil cp $plot "$out_metrics" &
    upload_bqsr_metrics_pid=$!
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
