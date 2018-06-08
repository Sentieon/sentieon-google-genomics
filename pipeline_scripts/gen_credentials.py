#!/usr/bin/env python

from __future__ import print_function

import requests
import time
import argparse
import json

audience = "https://sentieon.com"
headers = {'Metadata-Flavor': 'Google'}
request_format = "full"
metadata_url = "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/identity?audience={}&format={}"
project_url = "http://metadata.google.internal/computeMetadata/v1/project/project-id?format={}"

def process_args():
    parser = argparse.ArgumentParser(description="Write fresh instance metadata credentials to a file for license authentication")
    parser.add_argument("auth_data_file", help="A file to hold the instance metadata JWT")
    parser.add_argument("sentieon_key", help="A license key string")
    return parser.parse_args()

def main(args):
    if not args:
        args = process_args()

    url = project_url.format(request_format)
    response = requests.get(url, headers=headers)
    project_id = response.text
    with open(args.auth_data_file + ".project", 'w') as f:
        print(project_id, file=f)

    url = metadata_url.format(audience, request_format)

    while True:
        out = {}
        response = None

        # Handle brief network outages
        retry_count = 0
        while retry_count < 10:
            try:
                response = requests.get(url, headers=headers)
            except requests.exceptions.ConnectionError as e:
                time.sleep(60)
                continue
            break
        if not response:
            sys.exit("Failed to connect to instance metadata server")

        out["google_session_token"] = response.text
        if args.sentieon_key:
            out["license_key"] = args.sentieon_key
        with open(args.auth_data_file, 'w') as f:
            json.dump(out, f)
        time.sleep(55 * 60) # sleep for 55 minutes before refreshing the token or until killed

if __name__ == "__main__":
    main(None)

