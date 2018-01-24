#!/usr/bin/env python

from __future__ import print_function

r"""
Runs the Sentieon Genomics Tools workflows using Google Genomics Pipelines API
"""

import yaml
import json
import argparse
import os
import sys
import pprint
import copy
import time
import ssl

from apiclient.discovery import build
import google.auth

script_dir = os.path.dirname(os.path.realpath(__file__))
germline_yaml = script_dir + "/sentieon_germline.yaml"
tn_yaml = script_dir + "/sentieon_tn.yaml"
default_json = script_dir + "/runner_default.json"
target_url_base = "https://www.googleapis.com/compute/v1/projects/{project}/zones/{zone}/instances/{instance}"


def cloud_storage_exists(client, gs_path):
    
    try:
        bucket, blob = gs_path[5:].split('/', 1)
        bucket = client.bucket(bucket)
        blob = bucket.blob(blob)
        res = blob.exists()
    except: # Catch all exceptions
        raise ValueError("Error: Could not find {gs_path} in Google Cloud Storage".format(**locals()))
    return res


def check_inputs_exist(job_vars, credentials):
    from google.cloud import storage
    client = storage.Client(credentials=credentials)

    # The DBSNP, BQSR and Realign sites files
    sites_files = []
    sites_files += job_vars["BQSR_SITES"].split(',') if job_vars["BQSR_SITES"] else []
    sites_files += job_vars["REALIGN_SITES"].split(',') if job_vars["REALIGN_SITES"] else []
    sites_files += [job_vars["DBSNP"]] if job_vars["DBSNP"] else []
    for sites_file in sites_files:
        if not cloud_storage_exists(client, sites_file):
            sys.exit("Error: Could not find supplied file {}".format(sites_file))
        if sites_file.endswith("vcf.gz"):
            if not cloud_storage_exists(client, sites_file + ".tbi"):
                sys.exit("Error: Could not find index for file {}".format(sites_file))
        else:
            if not cloud_storage_exists(client, sites_file + ".idx"):
                sys.exit("Error: Could not find index for file {}".format(sites_file))

    # All reference files
    ref = job_vars["REF"]
    if not cloud_storage_exists(client, ref):
        sys.exit("Error: Reference file not found")
    if not cloud_storage_exists(client, ref + ".fai"):
        sys.exit("Error: Reference fai index not found")
    if (not cloud_storage_exists(client, ref + ".dict") and
        not cloud_storage_exists(client, ref[:-2] + "dict")):
        sys.exit("Error: Reference dict index not found")
    # FQ specific
    if job_vars["FQ1"]:
        for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
            if not cloud_storage_exists(client, ref + suffix):
                sys.exit("Error: Reference BWA index {} not found".format(suffix))
    # BAM specific
    else:
        if (not cloud_storage_exists(client, job_vars["BAM"] + ".bai") and
            not cloud_storage_exists(client, job_vars["BAM"][:-3] + "bai")):
            sys.exit("Error: BAM supplied but BAI not found")
    

def main(vargs=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("pipeline_config", help="The json configuration file")
    parser.add_argument("--no_check_inputs_exist", action="store_true", help="Do not check that the input files exist before running the pipeline")
    parser.add_argument("--polling_interval", type=float, default=30, help="Seconds between polling the running operation")
    args = parser.parse_args()
    polling_interval = args.polling_interval

    # Grab input arguments from the json file
    job_vars = json.load(open(default_json))
    job_vars.update(json.load(open(args.pipeline_config)))
    preemptible_tries = int(job_vars["PREEMPTIBLE_TRIES"])
    if job_vars["NONPREEMPTIBLE_TRY"]:
        non_preemptible_tries = 1
    preemptible = True if preemptible_tries > 0 else False
    credentials, project_id = google.auth.default()

    # Grab the yaml for the workflow
    pipeline_yaml = germline_yaml if job_vars["PIPELINE"] == "DNA" else tn_yaml
    pipeline_dict = yaml.load(open(pipeline_yaml))
    if job_vars["OUTPUT_BUCKET"].endswith('/'):
        job_vars["OUTPUT_BUCKET"] = job_vars["OUTPUT_BUCKET"][:-1]

    # Some basic error checking to fail early
    if not job_vars["PROJECT_ID"]:
        sys.exit("Error: Please supply a PROJECT_ID")
    if job_vars["PIPELINE"] == "DNA":
        if job_vars["FQ1"] and job_vars["BAM"]:
            sys.exit("Error: Please supply either 'FQ1' or 'BAM' (not both)")
        if not job_vars["FQ1"] and not job_vars["BAM"]:
            sys.exit("Error: Please supply either 'FQ1' or 'BAM'")
        if job_vars["INTERVAL"] and job_vars["INTERVAL_FILE"]:
            sys.exit("Error: Please supply either 'INTERVAL' or 'INTERVAL_FILE'")
        if job_vars["NO_HAPLOTYPER"] and job_vars["NO_METRICS"] and job_vars["NO_BAM_OUTPUT"]:
            sys.exit("Error: No output files requested")
        if not args.no_check_inputs_exist:
            check_inputs_exist(job_vars, credentials)

    # Construct the pipeline arguments
    args_dict = {}
    args_dict["projectId"] = job_vars["PROJECT_ID"]
    args_dict["logging"] = {"gcsPath": job_vars["OUTPUT_BUCKET"] + "/worker_logs/"}
    resources_dict = {}
    resources_dict["minimumRamGb"] = job_vars["MIN_RAM_GB"]
    resources_dict["minimumCpuCores"] = job_vars["MIN_CPU"]
    resources_dict["zones"] = job_vars["ZONES"].split(',') if job_vars["ZONES"] else []
    args_dict["resources"] = copy.copy(resources_dict)
    # Translate None back to "None"
    input_dict = {}
    for input_var in pipeline_dict["inputParameters"]:
        input_dict[input_var["name"]] = job_vars[input_var["name"]]
        if input_dict[input_var["name"]] is None:
            input_dict[input_var["name"]] = "None"
    args_dict["inputs"] = input_dict

    # Construct the pipeline object
    resources_dict = {}
    resources_dict["preemptible"] = preemptible
    # Ideally we could use mdadm to combine local SSD, but that is not possible inside the container
    disk = {}
    disk["name"] = "local-disk"
    disk["mountPoint"] = "/mnt/work"
    if int(job_vars["DISK_SIZE"]) <= 375:
        print("Disk is less than 375 GB, using a single local SSD")
        disk["type"] = "LOCAL_SSD"
    else:
        print("Disk is greater than 375 GB, using a persistant SSD")
        disk["type"] = "PERSISTENT_SSD"
        disk["sizeGb"] = int(job_vars["DISK_SIZE"])
    disks = [disk]
    resources_dict["disks"] = disks
    pipeline_dict["resources"] = resources_dict
    pipeline_dict["projectId"] = job_vars["PROJECT_ID"]
    pipeline_dict["docker"] = {
            "imageName": job_vars["DOCKER_IMAGE"],
            "cmd": "bash /opt/sentieon/" + "gc_germline.sh" if job_vars["PIPELINE"] == "DNA" else "gc_tn.sh"
    }

    # Run the pipeline #
    service = build('genomics', 'v1alpha2', credentials=credentials)
    compute_service = build("compute", "v1", credentials=credentials)
    operation = None
    counter = 0
    project = job_vars["PROJECT_ID"]
    
    while non_preemptible_tries > 0 or preemptible_tries > 0:
        if operation:
            while not operation["done"]:
                time.sleep(polling_interval)
                try:
                    operation = service.operations().get(name=operation['name']).execute()
                except ssl.SSLError:
                    print("Network error while polling running operation.")
                    sys.exit(1)
            if "error" in operation:
                pprint.pprint(operation, indent=2)
                zone = operation["metadata"]["runtimeMetadata"]["computeEngine"]["zone"]
                instance = operation["metadata"]["runtimeMetadata"]["computeEngine"]["instanceName"]
                url = target_url_base.format(**locals())
                compute_ops = compute_service.zoneOperations().list(project=project, zone=zone, filter="(targetLink eq {url}) (operationType eq compute.instances.preempted)".format(**locals())).execute()
                if "items" in compute_ops: # Failed due to preemption
                    print("Run {} failed. Retrying...".format(counter))
                else:
                    print("Run {} failed, but not due to preemption. Exit".format(counter))
                    operation = None
                    break
            else:
                pprint.pprint(operation, indent=2)
                break

        if preemptible_tries > 0:
            args_dict["resources"]["preemptible"] = True
            preemptible_tries -= 1
        else:
            args_dict["resources"]["preemptible"] = False
            non_preemptible_tries -= 1

        print("Running pipeline:")
        body = {
            "ephemeralPipeline": pipeline_dict,
            "pipelineArgs": args_dict
        }
        pprint.pprint(body, indent=2)

        operation = service.pipelines().run(body={
            "ephemeralPipeline": pipeline_dict,
            "pipelineArgs": args_dict
        }).execute()
        counter += 1

    if operation:
        while not operation["done"]:
            time.sleep(polling_interval)
            try:
                operation = service.operations().get(name=operation["name"]).execute()
            except ssl.SSLError:
                print("Network error while waiting for the final operation to finish")
                sys.exit(1)
        if "error" in operation:
            zone = operation["metadata"]["runtimeMetadata"]["computeEngine"]["zone"]
            instance = operation["metadata"]["runtimeMetadata"]["computeEngine"]["instanceName"]
            url = target_url_base.format(**locals())
            compute_ops = compute_service.zoneOperations().list(project=project, zone=zone, filter="(targetLink eq {url}) (operationType eq compute.instances.preempted)".format(**locals())).execute()
            if "items" in compute_ops:
                print("Final run failed due to preemption")
            else:
                print("Final run failed, not due to preemption")
        else:
            print("Operation succeeded")

if __name__ == "__main__":
    main()
