"""
This module contains utility methods for the SNP parser.
"""

import sys
import fnmatch
import os
from datetime import datetime

start_time = datetime.now()
true_values = ["TRUE", "T", "1", "YES", "Y"]

# Return the elapsed time as a tuple of minutes and seconds
def get_elapsed():
    return divmod((datetime.now() - start_time).total_seconds(), 60)

# Parse a boolean from a string value
def string_to_bool(string_val):
    if (string_val.upper() in true_values):
        value = True
    else:
        value = False
    return value

# Remove leading and trailing double quotes
def strip_quotes(val):
    if val[0:1] == "\"":
        val = val[1:]
    if val[len(val)-1:] == "\"":
        val = val[:len(val)-1]
    return val

# Increment a counter in a dictionary
def inc_dict_counter(counter, value):
    if (value in counter):
        counter[value] += 1
    else:
        counter[value] = 1
    return
