# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 15:02:05 2024

@author: trist
"""
import json

import pandas as pd

# Connor stored the values in milimetres, I want them in metres
CONVERSION_FACTOR = 1e3

beam_dict = {}
for beam_type in ('Y', "YE", "M"):
    df = pd.read_excel("./data/Bridge_Beam_Info.xlsx", sheet_name=beam_type)
    beam_dict[beam_type] = {}
    for _, row in df.iterrows():
        top_width = int(row['w']) / CONVERSION_FACTOR
        base_width = int(row['b']) / CONVERSION_FACTOR
        height = int(row['h']) / CONVERSION_FACTOR
        y0 = int(row["y0"]) / CONVERSION_FACTOR
        z0 = int(row["z0"]) / CONVERSION_FACTOR
        string_key = f"({top_width}, {base_width}, {height})"
        beam_enumeration = int(row['N'])
        beam_dict[beam_type][string_key] = {
            "N": beam_enumeration,
            "centroid": {
                "y0": y0,
                "z0": z0
            }
        }


for beam_type in ('U', "UM"):
    df = pd.read_excel("./data/Bridge_Beam_Info.xlsx", sheet_name=beam_type)
    beam_dict[beam_type] = {}
    for _, row in df.iterrows():
        opening_width  = int(row["w1"]) / CONVERSION_FACTOR
        top_width  = int(row["w2"]) / CONVERSION_FACTOR
        base_width = int(row['b']) / CONVERSION_FACTOR
        height = int(row['h']) / CONVERSION_FACTOR
        y0 = int(row["y0"]) / CONVERSION_FACTOR
        z0 = int(row["z0"]) / CONVERSION_FACTOR
        string_key = f"({opening_width}, {top_width}, {base_width}, {height})"
        beam_enumeration = int(row['N'])
        beam_dict[beam_type][string_key] = {
            "N": beam_enumeration,
            "centroid": {
                "y0": y0,
                "z0": z0
            }
        }


with open("./data/nonparametric_beam_geometries.json", 'w') as f:
    json.dump(beam_dict, f, indent=4)
    