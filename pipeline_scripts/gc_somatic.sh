#!/usr/bin/env bash

set -xveo pipefail
set +H

BASEDIR=$(dirname "$0")
version="201711.05"
release_dir="/opt/sentieon/sentieon-genomics-${version}/"
scratch=/mnt/work
nt=$(nproc)
source $BASEDIR/gc_functions.sh

export LD_PRELOAD=/opt/sentieon/sentieon-genomics-${version}/lib/libjemalloc.so
export MALLOC_CONF=lg_dirty_mult:-1

# Set "None" varibles to an empty string
environmental_variables=(FQ1 FQ2 TUMOR_FQ1 TUMOR_FQ2 BAM TUMOR_BAM \
    OUTPUT_BUCKET REF READGROUP TUMOR_READGROUP DEDUP BQSR_SITES DBSNP \
    INTERVAL INTERVAL_FILE NO_METRICS NO_BAM_OUTPUT NO_VCF RUN_TNSNV \
    STREAM_INPUT PIPELINE REALIGN_SITES OUTPUT_CRAM_FORMAT SENTIEON_KEY)
unset_none_variables ${environmental_variables[@]}
OUTPUT_CRAM_FORMAT="" # Not yet supported

# Basic error handling #
if [[ -n "$TUMOR_FQ1" && -n "$TUMOR_BAM" ]]; then
    echo "Please supply either tumor fastq or tumor BAM files"
    exit 1
fi

if [[ -z "$TUMOR_FQ1" && -z "$TUMOR_BAM" ]]; then
    echo "Please supply either tumor fastq or tumor BAM files"
    exit 1
fi

if [[ -n "$FQ1" && -n "$BAM" ]]; then
    echo "Please supply either fastq or BAM files"
    exit 1
fi

if [[ -n "$RUN_TNSNV" && -z "$REALIGN_SITES" ]]; then
    echo "TNsnv requires indel realignment. Plese supply REALIGN_SITES."
    exit 1
fi

if [[ -z "$RUN_TNSNV" && -n "$REALIGN_SITES" ]]; then
    echo "WARNING: Ignoring REALIGN_SITES with TNhaplotyper or TNscope"
    REALIGN_SITES=
fi

if [[ -n "$INTERVAL" && -n "$INTERVAL_FILE" ]]; then
    echo "Please supply either an INTERVAL or INTERVAL_FILE"
    exit 1
fi

if [[ -n "$NO_BAM_OUTPUT" && -n "$NO_" && -n "$NO_VCF" ]]; then
    echo "Nothing to do"
    exit 1
fi

readonly FQ1 FQ2 TUMOR_FQ1 TUMOR_FQ2 BAM TUMOR_BAM \
    OUTPUT_BUCKET REF READGROUP TUMOR_READGROUP DEDUP BQSR_SITES DBSNP \
    INTERVAL INTERVAL_FILE NO_METRICS NO_BAM_OUTPUT NO_VCF RUN_TNSNV \
    STREAM_INPUT PIPELINE REALIGN_SITES OUTPUT_CRAM_FORMAT

# *****************************
# 0. Setup
# *****************************
gc_setup

## Download input files
if [[ -n "$BAM" ]]; then
    download_bams $BAM local_bams $input_dir
else
    local_bams=()
fi

if [[ -n "$TUMOR_BAM" ]]; then
    download_bams $TUMOR_BAM tumor_bams $input_dir
else
    tumor_bams=()
fi

download_intervals
download_reference

## Handle the sites files
gs_bqsr_sites=($(echo $BQSR_SITES | tr ',' ' '))
transfer_all_sites $bqsr_dir "${gs_bqsr_sites[@]}"
local_bqsr_sites="${local_sites[@]}"
bqsr_sites="$local_str"

gs_realign_sites=($(echo $REALIGN_SITES | tr ',' ' '))
transfer_all_sites $realign_dir "${gs_realign_sites[@]}"
local_realign_sites="${local_sites[@]}"
realign_sites="$local_str"

if [[ -n "$DBSNP" ]]; then
    transfer_all_sites $dbsnp_dir $DBSNP
    dbsnp=${local_sites[0]}
fi

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
output_ext="bam"
#if [[ -n "$NO_BAM_OUTPUT" || "$DEDUP" != "nodup" || -n "$REALIGN_SITES" ]]; then
#    util_sort_xargs="${util_sort_xargs} --bam_compression 1 "
#fi

if [[ -n "$FQ1" ]]; then
    bwa_mem_align "normal_" $FQ1 $FQ2 $READGROUP local_bams $output_ext "$bwa_xargs" "$util_sort_xargs"
fi

local_bams_str=""
for bam in ${local_bams[@]}; do
    local_bams_str+=" -i $bam "
done

if [[ -n "$TUMOR_FQ1" ]]; then
    bwa_mem_align "tumor_" $TUMOR_FQ1 $TUMOR_FQ2 $TUMOR_READGROUP tumor_bams $output_ext "$bwa_xargs" "$util_sort_xargs"
fi

tumor_bams_str=""
for bam in ${tumor_bams[@]}; do
    tumor_bams_str+=" -i $bam "
done

# Detect the tumor and normal sample names
normal_sample=""
if [[ -f ${local_bams[0]} ]]; then
    normal_sample=$(samtools view -H ${local_bams[0]} | grep "^@RG" | head -n 1 | sed 's/^.*SM:\(.*\)	.*$/\1/')
fi
tumor_sample=$(samtools view -H ${tumor_bams[0]} | grep "^@RG" | head -n 1 | sed 's/^.*SM:\(.*\)	.*$/\1/')

# ******************************************
# 2. Metrics command
# ******************************************
metrics_cmd1=
metrics_cmd2=
metrics_files=
tumor_metrics_cmd1=
tumor_metrics_cmd2=
tumor_metrics_files=
if [[ -z "$NO_METRICS" ]]; then
    if [[ -n $normal_sample ]]; then
        build_metrics_cmd "normal_" metrics_cmd1 metrics_cmd2 metrics_files
    fi
    build_metrics_cmd "tumor_" tumor_metrics_cmd1 tumor_metrics_cmd2 tumor_metrics_files
fi

# ******************************************
# 3. Remove duplicates
# ******************************************
output_ext="bam"
#if [[ -n "$NO_BAM_OUTPUT" || -n "$REALIGN_SITES" ]]; then
#    dedup_xargs=" --bam_compression 1 "
#fi

run_mark_duplicates "normal_" $DEDUP metrics_cmd1 "$local_bams_str" dedup_bam_str dedup_bams "$dedup_xargs" $output_ext ${local_bams[@]}
run_mark_duplicates "tumor_" $DEDUP tumor_metrics_cmd1 "$tumor_bams_str" tumor_dedup_bam_str tumor_dedup_bams "$dedup_xargs" $output_ext ${tumor_bams[@]}

if [[ "$DEDUP" != "nodup" ]]; then
    if [[ -z "$NO_METRICS" ]]; then
        to_upload="$metrics_dir/tumor_dedup_metrics.txt"
        if [[ -n "$dedup_bam_str" ]]; then
            to_upload+=" $metrics_dir/normal_dedup_metrics.txt"
        fi
        (gsutil cp $to_upload $out_metrics &&
            rm $metrics_dir/*_dedup_metrics.txt) &
        upload_dedup_pid=$!
    else
        rm $metrics_dir/*_dedup_metrics.txt &
    fi
    for bam in ${local_bams[@]} ${tumor_bams[@]}; do
        if [[ -f $bam ]]; then
            rm $bam &
        fi
        if [[ -f ${bam}.bai ]]; then
            rm ${bam}.bai &
        fi
        if [[ -f ${bam}.crai ]]; then
            rm ${bam}.crai &
        fi
    done
fi

if [[ -z "$NO_BAM_OUTPUT" && -z "$REALIGN_SITES" ]]; then
    upload_list=""
    for bam in ${dedup_bams[@]} ${tumor_dedup_bams[@]}; do
        upload_list+=" $bam "
        if [[ -f "${bam}.bai" ]]; then
            upload_list+=" ${bam}.bai "
        elif [[ -f "${bam}.crai" ]]; then
            upload_list+=" ${bam}.crai "
        fi
    done
    gsutil cp $upload_list $out_bam &
    upload_deduped_pid=$!
fi

upload_metrics metrics_cmd1 metrics_cmd2 upload_metrics_pid ${metrics_files[@]}
upload_metrics tumor_metrics_cmd1 tumor_metrics_cmd2 tumor_upload_metrics_pid ${tumor_metrics_files[@]}

# ******************************************
# 4. Indel Realignment
# ******************************************
output_ext="bam"
#if [[ -n "$NO_BAM_OUTPUT" || -n "$dedup_bam_str" ]]; then
#    realign_xargs=" --bam_compression 1 "
#fi

if [[ -n "$REALIGN_SITES" && -n "$RUN_TNSNV" ]]; then
    realigned_bam=$work/normal_realigned.${output_ext}
    tumor_realigned_bam=$work/tumor_realigned.${output_ext}
    if [[ -n "$dedup_bam_str" ]]; then
        cmd="$release_dir/bin/sentieon driver $dedup_bam_str -t $nt -r $ref --algo Realigner $realign_xargs $interval_list $realign_sites $realigned_bam"
        run "$cmd" "Indel Realign - Normal"
    fi
    cmd="$release_dir/bin/sentieon driver $tumor_dedup_bam_str -t $nt -r $ref --algo Realigner $realign_xargs $interval_list $realign_sites $tumor_realigned_bam"
    run "$cmd" "Indel Realign - Tumor"

    # Cleanup
    for bam in ${dedup_bams[@]} ${tumor_dedup_bams[@]}; do
        if [[ -f "$bam" ]]; then
            rm $bam &
        fi
        if [[ -f "${bam}.bai" ]]; then
            rm ${bam}.bai &
        fi
        if [[ -f "${bam}.crai" ]]; then
            rm ${bam}.crai &
        fi
    done

    if [[ -n "$dedup_bam_str" ]]; then
        realigned_bams=($realigned_bam)
        realigned_bam_str=" -i $realigned_bam "
    else
        realigned_bams=()
        realigned_bam_str=""
    fi
    tumor_realigned_bams=($tumor_realigned_bam)
    tumor_realigned_bam_str=" -i $tumor_realigned_bam "
else
    realigned_bams=${dedup_bams[@]}
    realigned_bam_str=${dedup_bam_str}
    tumor_realigned_bams=${tumor_dedup_bams[@]}
    tumor_realigned_bam_str=${tumor_dedup_bam_str}
fi


# ******************************************
# 5. Base recalibration
# ******************************************
run_bqsr "normal_" "$realigned_bam_str" metrics_cmd1 bqsr_cmd2 bqsr_cmd3 bqsr_cmd4 bqsr_table bqsr_plot
# For paired samples, BQSR-post has to be run separately
if [[ -n "$realigned_bam_str" ]]; then
    run_bqsr_post "$realigned_bam_str" bqsr_cmd2 bqsr_cmd3 bqsr_cmd4 "$bqsr_table" "$bqsr_plot" upload_normal_bqsr_plot_pid
fi
run_bqsr "tumor_" "$tumor_realigned_bam_str" tumor_metrics_cmd1 tumor_bqsr_cmd2 tumor_bqsr_cmd3 tumor_bqsr_cmd4 tumor_bqsr_table tumor_bqsr_plot
if [[ -n "$realigned_bam_str" ]]; then
    run_bqsr_post "$tumor_realigned_bam_str" tumor_bqsr_cmd2 tumor_bqsr_cmd3 tumor_bqsr_cmd4 "$tumor_bqsr_table" "$tumor_bqsr_plot" tumor_upload_normal_bqsr_plot_pid
fi

if [[ -n "$bqsr_sites" && -z "$NO_BAM_OUTPUT" ]]; then
    gsutil cp $bqsr_table $tumor_bqsr_table $out_bam &
    upload_bqsr_pid=$!
fi

upload_metrics metrics_cmd1 metrics_cmd2 upload_metrics_pid ${metrics_files[@]}
upload_metrics tumor_metrics_cmd1 tumor_metrics_cmd2 tumor_upload_metrics_pid ${tumor_metrics_files[@]}

# ******************************************
# 6. Indel corealignment
# ******************************************
output_ext="bam"
#if [[ -n "$NO_BAM_OUTPUT" ]]; then
#    corealign_xargs=" --bam_compression 1 "
#fi

if [[ -n "$REALIGN_SITES" && -n "$RUN_TNSNV" && -n "$realigned_bam_str" ]]; then
    corealigned_bam=$work/corealigned.${output_ext}
    cmd="$release_dir/bin/sentieon driver $realigned_bam_str $tumor_realigned_bam_str -t $nt -r $ref --algo Realigner $corealign_xargs $interval_list $realign_sites $corealigned_bam"
    run "$cmd" "Indel co-realignment"

    # Cleanup
    for bam in ${realigned_bams[@]} ${tumor_realigned_bams[@]}; do
        if [[ -f $bam ]]; then
            rm $bam &
        fi
        if [[ -f ${bam}.bai ]]; then
            rm ${bam}.bai &
        fi
        if [[ -f ${bam}.crai ]]; then
            rm ${bam}.crai &
        fi
    done
    for site_file in ${local_realign_sites[@]}; do
        if [[ -f $site_file ]]; then
            rm $site_file &
        fi
    done

    # Upload
    upload_list=" $corealigned_bam "
    if [[ -f "${corealigned_bam}.bai" ]]; then
        upload_list+=" ${corealigned_bam}.bai "
    elif [[ -f "${corealigned_bam}.crai" ]]; then
        upload_list+=" ${corealigned_bam}.crai "
    fi
    gsutil cp $upload_list $out_bam &
    upload_corealigned_pid=$!

    corealigned_bam_str=" -i $corealigned_bam "
elif [[ -n "$REALIGN_SITES" && -n "$RUN_TNSNV" ]]; then
    upload_list=""
    for bam in ${tumor_realigned_bams[@]}; do
        upload_list+=" $bam "
        if [[ -f "${bam}.bai" ]]; then
            upload_list+=" ${bam}.bai "
        elif [[ -f "${bam}.crai" ]]; then
            upload_list+=" ${bam}.crai "
        fi
    done
    gsutil cp $upload_list $out_bam &
    upload_corealigned_pid=$!

    corealigned_bam_str=" $tumor_realigned_bam_str "
else
    corealigned_bam_str="$realigned_bam_str $tumor_realigned_bam_str "
fi
corealigned_bqsr_str=" ${tumor_bqsr_table:+-q $tumor_bqsr_table} ${bqsr_table:+-q $bqsr_table}"

# *******************************************
# 7. Variant Calling
# *******************************************

## Generate a non-decoy BED file
generate_nondecoy_bed ${ref}.fai ${ref}_nondecoy.bed

if [[ -z "$NO_VCF" ]]; then
    extra_calling_args=""
    if [[ "$PIPELINE" == "TNscope" ]]; then
        algo="TNscope"
        vcf=$work/tnscope.vcf.gz
    elif [[ -n "$RUN_TNSNV" ]]; then
        algo="TNsnv"
        vcf=$work/tnsnv.vcf.gz
    else
        algo="TNhaplotyper"
        vcf=$work/tnhaplotyper.vcf.gz
    fi

    cmd="$release_dir/bin/sentieon driver ${interval:- --interval ${ref}_nondecoy.bed} $corealigned_bam_str $corealigned_bqsr_str -t $nt -r $ref --algo $algo ${normal_sample:+--normal_sample $normal_sample} --tumor_sample $tumor_sample ${dbsnp:+--dbsnp $dbsnp} $vcf"

    if [[ -n $tumor_metrics_cmd1 ]]; then
        cmd+=" $metrics_cmd1 "
        cmd+=" $tumor_metrics_cmd1 "
        metrics_cmd1=
        tumor_metrics_cmd1=
    fi

    if [[ -n $tumor_bqsr_cmd2 ]]; then
        # cmd+=" $bqsr_cmd2 " Only 1 output file
        cmd+=" $tumor_bqsr_cmd2 "
        bqsr_cmd2=
        tumor_bqsr_cmd2=
    fi

    run "$cmd" "Variant calling"
    gsutil cp $vcf ${vcf}.tbi $out_variants &
    upload_vcf_pid=$!
fi

if [[ -n $tumor_metrics_cmd1 ]]; then
    cmd="$release_dir/bin/sentieon driver ${interval:- --interval ${ref}_nondecoy.bed} -t $nt -r '$ref' $corealigned_bam_str $corealigned_bqsr_str"
    cmd+=" $metrics_cmd1 "
    cmd+=" $tumor_metrics_cmd1 "
    metrics_cmd1=
    tumor_metrics_cmd1=
    run "$cmd" "Metrics collection"
fi

upload_metrics metrics_cmd1 metrics_cmd2 upload_metrics_pid ${metrics_files[@]}
upload_metrics tumor_metrics_cmd1 tumor_metrics_cmd2 tumor_upload_metrics_pid ${tumor_metrics_files[@]}

if [[ -n $bqsr_cmd2 ]]; then
    cmd="$release_dir/bin/sentieon driver ${interval:- --interval ${ref}_nondecoy.bed} -t $nt -r $ref $corealigned_bam_str $corealigned_bqsr_str $bqsr_cmd2 $tumor_bqsr_cmd2"
    bqsr_cmd2=
    tumor_bqsr_cmd2=
    run "$cmd" "BQSR Post"
fi

if [[ -n $tumor_bqsr_cmd3 ]]; then
    upload_list=""
    if [[ -n $bqsr_cmd3 ]]; then
        run "$bqsr_cmd3" "BQSR CSV"
        run "$bqsr_cmd4" "BQSR plot"
        upload_list+=" $bqsr_plot "
    fi
    run "$tumor_bqsr_cmd3" "Tumor BQSR CSV"
    run "$tumor_bqsr_cmd4" "Tumor BQSR plot"
    gsutil cp $upload_list $tumor_bqsr_plot $out_metrics &
    upload_bqsr_metrics_pid=$!
fi

kill $credentials_pid

# Wait for all running uploads to finish
set +e
if [[ -n $upload_normal_bqsr_plot_pid ]]; then
    wait $upload_normal_bqsr_plot_pid
fi
if [[ -n $tumor_upload_normal_bqsr_plot_pid ]]; then
    wait $tumor_upload_normal_bqsr_plot_pid
fi
if [[ -n $upload_dedup_pid ]]; then
    wait $upload_dedup_pid
fi
if [[ -n $upload_deduped_pid ]]; then
    wait $upload_deduped_pid
fi
if [[ -n $upload_metrics_pid ]]; then
    wait $upload_metrics_pid
fi
if [[ -n $tumor_upload_metrics_pid ]]; then
    wait $tumor_upload_metrics_pid
fi
if [[ -n $upload_bqsr_pid ]]; then
    wait $upload_bqsr_pid
fi
if [[ -n $upload_corealigned_pid ]]; then
    wait $upload_corealigned_pid
fi
if [[ -n $upload_vcf_pid ]]; then
    wait $upload_vcf_pid
fi
if [[ -n $upload_bqsr_metrics_pid ]]; then
    wait $upload_bqsr_metrics_pid
fi

