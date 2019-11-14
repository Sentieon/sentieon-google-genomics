#!/usr/bin/env python

from __future__ import print_function

r"""
Runs the Sentieon Genomics Tools workflows using Google Pipelines API
"""

import yaml
import json
import argparse
import os
import sys
import copy
import time
import ssl
import warnings
import logging
import random
import google.auth
import googleapiclient.errors

from apiclient.discovery import build
from pprint import pformat
from googleapiclient.errors import HttpError

script_dir = os.path.dirname(os.path.realpath(__file__))
germline_yaml = script_dir + "/germline.yaml"
somatic_yaml = script_dir + "/somatic.yaml"
ccdg_yaml = script_dir + "/ccdg.yaml"
default_json = script_dir + "/runner_default.json"
target_url_base = ("https://www.googleapis.com/compute/v1/projects/{project}/"
                   "zones/{zone}/instances/{instance}")


def cloud_storage_exists(client, gs_path):
    try:
        bucket, blob = gs_path[5:].split('/', 1)
        bucket = client.bucket(bucket)
        blob = bucket.blob(blob)
        res = blob.exists()
    except:  # Catch all exceptions
        raise ValueError("Error: Could not find {gs_path} in Google Cloud "
                         "Storage".format(**locals()))
    return res


def check_inputs_exist(job_vars, credentials):
    from google.cloud import storage
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Your application has authenticated "
                                "using end user credentials from Google Cloud "
                                "SDK")
        client = storage.Client(credentials=credentials)

    # The DBSNP, BQSR and Realign sites files
    sites_files = []
    sites_files += (job_vars["BQSR_SITES"].split(',') if
                    job_vars["BQSR_SITES"] else [])
    sites_files += (job_vars["REALIGN_SITES"].split(',') if
                    job_vars["REALIGN_SITES"] else [])
    sites_files += [job_vars["DBSNP"]] if job_vars["DBSNP"] else []
    for sites_file in sites_files:
        if not cloud_storage_exists(client, sites_file):
            logging.error("Could not find supplied file "
                          "{}".format(sites_file))
            sys.exit(-1)
        if sites_file.endswith("vcf.gz"):
            if not cloud_storage_exists(client, sites_file + ".tbi"):
                logging.error("Could not find index for file "
                              "{}".format(sites_file))
                sys.exit(-1)
        else:
            if not cloud_storage_exists(client, sites_file + ".idx"):
                logging.error("Could not find index for file "
                              "{}".format(sites_file))
                sys.exit(-1)

    # The data input files
    gs_split_files = (
            job_vars["FQ1"],
            job_vars["TUMOR_FQ1"],
            job_vars["FQ2"],
            job_vars["TUMOR_FQ2"],
            job_vars["BAM"],
            job_vars["TUMOR_BAM"])
    gs_files = ()
    for split_file in gs_split_files:
        if not split_file:
            continue
        for input_file in split_file.split(','):
            if not cloud_storage_exists(client, input_file):
                logging.error("Could not find the supplied file "
                              "{}".format(input_file))
                sys.exit(-1)
    for input_file in gs_files:
        if not cloud_storage_exists(client, input_file):
            logging.error("Could not file the supplied file "
                          "{}".format(input_file))
            sys.exit(-1)

    # All reference files
    ref = job_vars["REF"]
    ref_base = ref[:-3] if ref.endswith(".fa") else ref[:-6]
    if not cloud_storage_exists(client, ref):
        logging.error("Reference file not found")
        sys.exit(-1)
    if not cloud_storage_exists(client, ref + ".fai"):
        logging.error("Reference fai index not found")
        sys.exit(-1)
    if (not cloud_storage_exists(client, ref + ".dict") and
            not cloud_storage_exists(client, ref_base + ".dict")):
        logging.error("Reference dict index not found")
        sys.exit(-1)
    # FQ specific
    if job_vars["FQ1"] or job_vars["TUMOR_FQ1"]:
        for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
            if (not cloud_storage_exists(client, ref + suffix) and
                    not cloud_storage_exists(client, ref + ".64" + suffix)):
                logging.error("Reference BWA index {} not "
                              "found".format(suffix))
                sys.exit(-1)
    # BAM specific
    bam_vars = ("BAM", "TUMOR_BAM")
    for bam_type in bam_vars:
        if job_vars[bam_type]:
            for bam in job_vars[bam_type].split(','):
                if (not cloud_storage_exists(client, bam + ".bai") and
                        not cloud_storage_exists(client, bam + "bai")):
                    logging.error("BAM supplied but BAI not found")
                    sys.exit(-1)


def main(vargs=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("pipeline_config", help="The json configuration file")
    parser.add_argument(
            "--verbose",
            "-v",
            action="count",
            help="Increase the runner verbosity")
    parser.add_argument(
            "--no_check_inputs_exist",
            action="store_true",
            help="Do not check that the input files exist before running the "
                 "pipeline")
    parser.add_argument(
            "--polling_interval",
            type=float,
            default=30,
            help="Seconds between polling the running operation")
    args = parser.parse_args()
    polling_interval = args.polling_interval

    logging.getLogger("googleapiclient."
                      "discovery_cache").setLevel(logging.ERROR)
    log_format = "%(filename)s::%(funcName)s [%(levelname)s] %(message)s"
    if not hasattr(args, "verbose") or args.verbose is None:
        log_level = logging.WARNING
    elif args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    else:
        log_level = logging.WARNING
    logging.basicConfig(level=log_level, format=log_format)

    # Grab input arguments from the json file
    job_vars = json.load(open(default_json))
    job_vars.update(json.load(open(args.pipeline_config)))
    preemptible_tries = int(job_vars["PREEMPTIBLE_TRIES"])
    if job_vars["NONPREEMPTIBLE_TRY"]:
        non_preemptible_tries = 1
    preemptible = True if preemptible_tries > 0 else False
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Your application has authenticated "
                                "using end user credentials from Google Cloud "
                                "SDK")
        credentials, project_id = google.auth.default()

    # Warn with depreciated JSON keys
    if "MIN_RAM_GB" in job_vars or "MIN_CPU" in job_vars:
        logging.warning("'MIN_RAM_GB' and 'MIN_CPU' are now ignored. "
                        "Please use 'MACHINE_TYPE' to specify the instance "
                        "type")

    # Grab the yaml for the workflow
    pipeline = job_vars["PIPELINE"]
    if pipeline == "GERMLINE":
        pipeline_yaml = germline_yaml
    elif pipeline == "SOMATIC":
        pipeline_yaml = somatic_yaml
    elif pipeline == "CCDG":
        pipeline_yaml = ccdg_yaml
    else:
        logging.error("Pipeline '" + pipeline + "'. Valid "
                      "values are 'GERMLINE' and 'SOMATIC'")
        sys.exit(-1)
    try:
        pipeline_dict = yaml.safe_load(open(pipeline_yaml))
    except IOError:
        logging.error("No yaml \"{}\" found.".format(pipeline_yaml))
        sys.exit(-1)

    # Try not to create nearly empty directories
    while job_vars["OUTPUT_BUCKET"].endswith('/'):
        job_vars["OUTPUT_BUCKET"] = job_vars["OUTPUT_BUCKET"][:-1]

    # Some basic error checking to fail early
    if not job_vars["PROJECT_ID"]:
        logging.error("Please supply a PROJECT_ID")
        sys.exit(-1)

    # Shared errors
    if job_vars["FQ1"] and job_vars["BAM"]:
        logging.error("Please supply either 'FQ1' or 'BAM' (not both)")
        sys.exit(-1)
    if job_vars["INTERVAL"] and job_vars["INTERVAL_FILE"]:
        logging.error("Please supply either 'INTERVAL' or 'INTERVAL_FILE'")
        sys.exit(-1)
    if ((job_vars["FQ1"] and job_vars["READGROUP"]) and
            (len(job_vars["FQ1"].split(',')) !=
             len(job_vars["READGROUP"].split(',')))):
        logging.error("The number of fastq files must match the number of "
                      "supplied readgroups")
        sys.exit(-1)

    # Pipeline specific errors
    if pipeline == "GERMLINE" or pipeline == "CCDG":
        if not job_vars["FQ1"] and not job_vars["BAM"]:
            logging.error("Please supply either 'FQ1' or 'BAM'")
            sys.exit(-1)
        if (job_vars["NO_HAPLOTYPER"] and
                job_vars["NO_METRICS"] and
                job_vars["NO_BAM_OUTPUT"]):
            logging.error("No output files requested")
            sys.exit(-1)
        if job_vars["RECALIBRATED_OUTPUT"] and job_vars["BQSR_SITES"] is None:
            logging.error("Cannot output a recalibrated BAM file without "
                          "running BQSR. Please supply 'BQSR_SITES'")
            sys.exit(-1)
        valid_algos = ("Haplotyper", "DNAscope")
        if job_vars["CALLING_ALGO"] not in valid_algos:
            logging.error(job_vars["CALLING_ALGO"] + "' is not a "
                          "valid germline variant calling algo. Please set "
                          "'CALLING_ALGO' to one of " + str(valid_algos))
            sys.exit(-1)
        # Additional CCDG checks
        if pipeline == "CCDG":
            if job_vars["BQSR_SITES"] is None:
                logging.error("The CCDG pipeline requires known sites for "
                              "BQSR. Please supply 'BQSR_SITES'")
                sys.exit(-1)
    elif pipeline == "SOMATIC":
        if job_vars["TUMOR_FQ1"] and job_vars["TUMOR_BAM"]:
            logging.error("Please supply either 'TUMOR_FQ1' or 'TUMOR_BAM' "
                          "(not both)")
            sys.exit(-1)
        if (not job_vars["TUMOR_FQ1"] and
                not job_vars["TUMOR_BAM"]):
            logging.error("Please supply either 'TUMOR_FQ1' or 'TUMOR_BAM'")
            sys.exit(-1)
        if (job_vars["RUN_TNSNV"] and
                not job_vars["REALIGN_SITES"]):
            logging.error("TNsnv requires indel realignment. Please supply "
                          "'REALIGN_SITES'")
            sys.exit(-1)
        if (job_vars["NO_BAM_OUTPUT"] and
                job_vars["NO_VCF"] and job_vars["NO_METRICS"]):
            logging.error("No output files requested")
            sys.exit(-1)
        if ((job_vars["TUMOR_FQ1"] and job_vars["TUMOR_READGROUP"]) and
                (len(job_vars["TUMOR_FQ1"].split(',')) !=
                    len(job_vars["TUMOR_READGROUP"].split(',')))):
            logging.error("The number of tumor fastq files must match the "
                          "number of supplied readgroups")
            sys.exit(-1)
        valid_algos = ("TNhaplotyper", "TNhaplotyper2", "TNscope", "TNsnv")
        if job_vars["CALLING_ALGO"] not in valid_algos:
            logging.error(job_vars["CALLING_ALGO"] + "' is not a "
                          "valid somatic variant calling algo. Please set "
                          "'CALLING_ALGO' to one of " + str(valid_algos))
            sys.exit(-1)
    if not args.no_check_inputs_exist:
        check_inputs_exist(job_vars, credentials)

    # Resources dict
    zones = job_vars["ZONES"].split(',') if job_vars["ZONES"] else []
    if not zones:
        logging.error("Please supply at least one zone to run the pipeline")
    region = zones[0][:-2]
    disk = {
        "name": "local-disk",
        "type": "local-ssd",
        "sizeGb": int(job_vars["DISK_SIZE"])
    }
    vm_dict = {
        "machineType": job_vars["MACHINE_TYPE"],
        "preemptible": preemptible,
        "disks": [disk],
        "serviceAccount": {"scopes": [
            "https://www.googleapis.com/auth/cloud-platform"]},
        "cpuPlatform": job_vars["CPU_PLATFORM"]
    }
    resources_dict = {
        "zones": zones,
        "virtualMachine": vm_dict
    }

    # Environment
    env_dict = {}
    for input_var in pipeline_dict["inputParameters"]:
        env_dict[input_var["name"]] = job_vars[input_var["name"]]
        if env_dict[input_var["name"]] is None:
            env_dict[input_var["name"]] = "None"

    # Action
    if pipeline == "GERMLINE":
        _cmd = "/opt/sentieon/gc_germline.sh"
    elif pipeline == "SOMATIC":
        _cmd = "/opt/sentieon/gc_somatic.sh"
    elif pipeline == "CCDG":
        _cmd = "/opt/sentieon/gc_ccdg_germline.sh"
    else:
        logging.error("Error: Unknown pipeline " + pipeline)
        sys.exit(-1)

    run_action = {
        "containerName": "run-pipeline",
        "imageUri": job_vars["DOCKER_IMAGE"],
        "commands": ["/bin/bash", _cmd],
        "mounts": [{
            "disk": "local-disk",
            "path": "/mnt/work",
            "readOnly": False
        }],
    }

    cleanup_action = {
        "containerName": "cleanup",
        "imageUri": job_vars["DOCKER_IMAGE"],
        "commands": [
            "/bin/bash",
            "-c",
            ("gsutil cp /google/logs/action/1/stderr "
             "\"{}/worker_logs/stderr.txt\" && "
             "gsutil cp /google/logs/action/1/stdout "
             "\"{}/worker_logs/stdout.txt\"").format(
                 job_vars["OUTPUT_BUCKET"], job_vars["OUTPUT_BUCKET"])],
        "alwaysRun": True
    }

    # Run the pipeline #
    project = job_vars["PROJECT_ID"]
    service_parent = "projects/" + project + "/locations/" + region
    service = build('lifesciences', 'v2beta', credentials=credentials)
    compute_service = build("compute", "v1", credentials=credentials)
    operation = None
    counter = 0

    while non_preemptible_tries > 0 or preemptible_tries > 0:
        if operation:
            while not operation.get("done", False):
                new_op, tries = None, 0
                while tries <= 5:
                    time.sleep(polling_interval)
                    try:
                        ops = service.projects().locations().operations()
                        new_op = ops.get(name=operation['name']).execute()
                        break
                    except googleapiclient.errors.HttpError as e:
                        logging.warning(str(e))
                        tries += 1
                    except ssl.SSLError as e:
                        logging.warning(str(e))
                        tries += 1
                if not new_op:
                    logging.error("Network error while polling running "
                                  "operation.")
                    sys.exit(1)
                operation = new_op
            logging.debug(pformat(operation, indent=2))
            if "error" in operation:
                assigned_events = list(filter(
                        lambda x: "workerAssigned" in x.keys(),
                        operation["metadata"]["events"]))
                if not assigned_events:
                    logging.error("Genomics operation failed before running:")
                    logging.error(pformat(operation["error"], indent=2))
                    sys.exit(2)

                startup_event = assigned_events[-1]
                instance = startup_event["workerAssigned"]["instance"]
                zone = startup_event["workerAssigned"]["zone"]
                url = target_url_base.format(**locals())
                time.sleep(300)  # It may take some time to set the operation
                compute_ops = (
                        compute_service.zoneOperations().list(
                            project=project, zone=zone, filter=(
                                "(targetLink eq {url}) (operationType eq "
                                "compute.instances.preempted)"
                            ).format(**locals())).execute())
                if ("items" in compute_ops and
                        any([(x["operationType"] ==
                              "compute.instances.preempted")
                            for x in compute_ops["items"]])):
                    logging.warning("Run {} failed. "
                                    "Retrying...".format(counter))
                else:
                    logging.error("Run {} failed, but not due to preemption. "
                                  "Exit".format(counter))
                    operation = None
                    break
            else:
                break

        if preemptible_tries > 0:
            vm_dict["preemptible"] = True
            preemptible_tries -= 1
        else:
            vm_dict["preemptible"] = False
            non_preemptible_tries -= 1

        logging.debug("Running pipeline:")
        body = {
            "pipeline": {
                "actions": [run_action, cleanup_action],
                "resources": resources_dict,
                "environment": env_dict
            }
        }

        logging.debug(pformat(body, indent=2))
        sys.stderr.flush()
        backoff, backoff_interval = 0, 1
        while backoff < 6:
            time.sleep(backoff_interval * random.random() * (2 ** backoff - 1))
            try:
                op_pipelines = service.projects().locations().pipelines()
                request = op_pipelines.run(
                        parent=service_parent,
                        body=body)
                operation = request.execute()
                break
            except googleapiclient.errors.HttpError as e:
                logging.warning(str(e))
                backoff += 1
        if not operation:
            logging.error("Failed to launch job")
            sys.exit(3)
        else:
            logging.warning("Launched job: " + operation["name"])
        counter += 1
        logging.debug(pformat(operation, indent=2))

    if operation:
        while not operation.get("done", False):
            new_op, tries = None, 0
            while tries <= 5:
                time.sleep(polling_interval)
                try:
                    new_op = service.projects().locations().operations().get(
                        name=operation["name"]).execute()
                    break
                except googleapiclient.errors.HttpError as e:
                    logging.warning(str(e))
                    tries += 1
                except ssl.SSLError:
                    logging.warning(str(e))
                    tries += 1
            if not new_op:
                logging.error("Network error while waiting for the final "
                              "operation to finish")
                sys.exit(1)
            operation = new_op
        logging.debug(pformat(operation, indent=2))
        if "error" in operation:
            assigned_events = list(filter(
                    lambda x: "workerAssigned" in x.keys(),
                    operation["metadata"]["events"]))
            if not assigned_events:
                logging.error("Genomics operation failed before running:")
                logging.error(pformat(operation["error"], indent=2))
                sys.exit(2)

            startup_event = assigned_events[-1]
            instance = startup_event["workerAssigned"]["instance"]
            zone = startup_event["workerAssigned"]["zone"]
            url = target_url_base.format(**locals())
            compute_ops = compute_service.zoneOperations().list(
                    project=project,
                    zone=zone,
                    filter=("(targetLink eq {url}) (operationType eq "
                            "compute.instances.preempted)").format(
                        **locals())).execute()
            logging.error("Final run failed.")
        else:
            logging.warning("Operation succeeded")


if __name__ == "__main__":
    main()
