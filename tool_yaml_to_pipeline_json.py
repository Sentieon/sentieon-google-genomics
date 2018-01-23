#!/usr/bin/env python

from __future__ import print_function
import yaml
import copy
import os
import json

def add_to_yaml(a, b):
    '''
    If an object is in b but not a, add it to a.

    Fail if confilicting values are found
    '''
    if type(a) != type(b):
        raise ValueError("Incomptable yaml entries found")
    if type(a) is dict:
        for k in b.keys():
            if k in a:
                add_to_yaml(a[k], b[k])
            else:
                a[k] = b[k]
    elif type(a) is list: # In this case, order is unimportant
        try: # Only applicable with zones.
            a_set, b_set = set(a), set(b)
            intersection = list(a_set.intersection(b_set))
            for i in range(len(a) - 1, -1, -1):
                if a[i] not in intersection:
                    a.pop(i)
            if len(a) == 0:
                raise ValueError("No values left after intersection")
            return
        except TypeError:
            pass
        # item["name"] essentially acts as a key
        a_dict, b_dict = dict([(x["name"], x) for x in a]), dict([(x["name"], x) for x in b])
        for k, v in b_dict.items():
            if k in a_dict:
                if v != a_dict[k]:
                    raise ValueError("Paramters are unequal")
            else:
                a.append(v)
    else:
        if a != b:
            raise ValueError("Values are unequal")

germline = yaml.load(open(os.path.dirname(os.path.realpath(__file__)) + "/sentieon_germline.yaml"))
tn = None
try:
    tn = yaml.load(open(os.path.dirname(os.path.realpath(__file__)) + "/sentieon_tn.yaml"))
except IOError:
    pass
output = os.path.dirname(os.path.realpath(__file__)) + "/runner_default.json"

additional_input_params = {
        "ZONES": None,
        "DISK_SIZE": 300,
        "MIN_CPU": 64,
        "MIN_RAM_GB": 56,
        "PIPELINE": "DNA",
        "PROJECT_ID": None,
        "DOCKER_IMAGE": "sentieon/sentieon-google-cloud:201711.01",
        "PREEMPTIBLE_TRIES": 0,
        "N_TRIES": 1
}

out_yaml = copy.deepcopy(germline)
if tn:
    add_to_yaml(out_yaml, tn)
out_json = {}
for param in out_yaml["inputParameters"]:
    k = param["name"]
    if "defaultValue" not in param:
        v = None
    elif param["defaultValue"] == "None":
        v = None
    else:
        v = param["defaultValue"]
    out_json[k] = v
out_json.update(additional_input_params)

with open(output, 'w') as f:
    json.dump(out_json, f, indent=2)
