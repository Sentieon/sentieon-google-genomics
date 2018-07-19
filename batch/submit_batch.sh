#!/usr/bin/env bash

set -e

n_concurrent=2
runner_script=../runner/sentieon_runner.py
polling_interval=20 # in seconds

base_json=$1; shift
batch_tsv=$1; shift

header=()
read -r header_line < "$batch_tsv"
header=($header_line) # No JSON keys should contain a space

jobs_to_run=()
b=()
job_base=$( grep -v "}" "$base_json")

while IFS='' read -r line || [[ -n "$line" ]]; do
    IFS='	' read -r -a line_arr <<< "$line"
    
    job="$job_base"
    for idx in "${!line_arr[@]}"; do
        job+=", \"${header[$idx]}\":\"${line_arr[$idx]}\""
    done
    job+=" }"
    jobs_to_run+=("$job")
done <<< "$(tail -n +2 "$batch_tsv")"

running=()
echo "${jobs_to_run[@]}"

while [[ -n "${jobs_to_run[@]}" ]]; do
    to_run="${jobs_to_run[0]}"
    unset jobs_to_run[0]
    jobs_to_run=( "${jobs_to_run[@]}" )

    # Check running jobs
    while [[ "${#running[@]}" -ge "$n_concurrent" ]]; do
        sleep "${polling_interval}"

        for idx in "${!running[@]}"; do
            if ! ps -p ${running[$idx]} > /dev/null; then
                echo "finished job with PID ${running[$idx]}"
                unset running[$idx]
                running=( "${running[@]}" )
            fi
        done
    done

    # Spawn new jobs
    echo $to_run | python2 $runner_script /dev/stdin &
    running+=( "$!" )

    # Sleep
    sleep "${polling_interval}"
done
