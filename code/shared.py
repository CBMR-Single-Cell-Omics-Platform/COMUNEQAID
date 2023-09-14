# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 15:02:40 2023

@author: vlf706
"""
import pandas as pd
from functools import reduce

def lenient_from_dict(my_dictionary):
    def convert_scalar_to_list(data):
        if isinstance(data, list):
            return data
        else:
            return [data]
        
    updated_dictionary = {key: convert_scalar_to_list(value) for \
                          key, value in my_dictionary.items()}
    return pd.DataFrame.from_dict(updated_dictionary)

def merge_tables(*argv, snake_obj):
    pd_list = [lenient_from_dict(snake_obj[key]) for key in argv]
    return reduce(pd.merge, pd_list)

def get_software_version(software):
    commands = {
        "salmon": ["salmon", "--version"],
        "alevin-fry": ["alevin-fry", "--version"],
        "bcl-convert": ["bcl-convert", "--version"]
    }

    if software not in commands:
        raise ValueError(f"Unsupported software: {software}")

    output = subprocess.check_output(commands[software], text=True).strip()
    
    # For bcl-convert, we extract the last 3 segments of the version. For others, we take the full version.
    if software == "bcl-convert":
        software_name = output.split()[0].lower()
        version_full = output.split()[2]
        version = "v" + ".".join(version_full.split(".")[-3:])
    else:
        software_name = output.split()[0].lower()
        version = "v" + output.split()[1]

    return f"{software_name}_{version}"
