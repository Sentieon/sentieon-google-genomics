#!/usr/bin/env bash

while true; do
    echo -e "\n\nTesting network... - $(date)"
    getent hosts gcp.sentieon.com
    host gcp.sentieon.com 169.254.169.254
    host gcp.sentieon.com 8.8.8.8
    sleep 30
done
