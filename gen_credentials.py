#!/usr/bin/env python

import requests
import time
import argparse
import json

audience = "https://sentieon.com"
headers = {'Metadata-Flavor': 'Google'}
request_format = "full"
metadata_url = "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/identity?audience={}&format={}"

def process_args():
    parser = argparse.ArgumentParser(description="Write fresh instance metadata credentials to a file for license authentication")
    parser.add_argument("auth_data_file", help="A file to hold the instance metadata JWT")
    return parser.parse_args()

def main(args):
    if not args:
        args = process_args()

    url = metadata_url.format(audience, request_format)

    while True:
        out = {}
        response = requests.get(url, headers=headers)
        out["google_session_token"] = response.text
        with open(args.auth_data_file, 'w') as f:
            json.dump(out, f)
        time.sleep(55 * 60) # sleep for 55 minutes before refreshing the token or until killed

if __name__ == "__main__":
    main(None)

