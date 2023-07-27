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
