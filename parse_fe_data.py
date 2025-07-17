# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 11:23:51 2025

@author: trist
"""
import re
import os
from collections import defaultdict
import json
from datetime import datetime, timezone
import time
from itertools import accumulate

import pymongo
from urllib.parse import quote_plus


def get_centre_span_locations(si_data):
    """
    Given some SI data, get a list of the centre span locations, this is where
    the point loads have been applied.
    """
    si_params = si_data["models"]["structuralInformation"]["parameters"]
    span_lengths = si_params["span"]["length"]["values"]
    if isinstance(span_lengths, (int, float)):
        span_lengths = [span_lengths]
    bridge_width = (si_params["beams"]["quantity"]-1) * si_params["beams"]["center_to_center"]["value"]
    
    x_temp = list(accumulate(span_lengths))
    x_temp.insert(0, 0)
    
    centre_span_locations = []
    mid_y = bridge_width / 2
    
    for i in range(len(span_lengths)):
        mid_x = (x_temp[i] + x_temp[i+1]) / 2
        centre_span_locations.append((mid_x, mid_y, 0))
    
    return centre_span_locations


def parse_loadcases(folder, filename):
    """
    Loads in the output results file csv and parses it into a more pythonic
    form.
    """
    filepath = os.path.join(folder, filename)
    loadcases = defaultdict(list)
    structure_name = filename.split('_')[0]
    structure_name = structure_name.replace("bridge-1-1", "bridge-1-1-1-0")
    loadcases["metadata"] = {"structure_name": structure_name}

    current_case = None
    recording = False
    in_header = True
    
    if "displacement" in filepath.lower():
        loadcases["metadata"]["data_type"] = "displacement"
    elif "reaction" in filepath.lower():
        loadcases["metadata"]["data_type"] = "reaction"
    else:
        raise ValueError("Unsupported data type.")
    
    with open(filepath, 'r', encoding="utf-8") as file:
        meta_data_retreived = 0
        for line in file:
            stripped = line.strip()
            if in_header:
                if "Created" in stripped:
                    dt = datetime.strptime(stripped.split('\t')[-1], "%d/%m/%y %H:%M")
                    dt = dt.replace(tzinfo=timezone.utc)
                    timestamp_seconds = dt.timestamp()
                    timestamp_nanoseconds = int(timestamp_seconds * 1e9)
                    loadcases["metadata"]["timestamps"] = {
                        "generated": timestamp_nanoseconds,
                        "stored": time.time_ns()
                    }
                    meta_data_retreived += 1
                
                elif "Model title" in stripped:
                    loadcases["metadata"]["title"] = stripped.split('\t')[-1]
                    meta_data_retreived += 1
                
                if meta_data_retreived == 2:
                    in_header = False
                continue
            
            match_reg = re.match(r'^(\d+):([A-Za-z0-9_]+)', stripped)
            
            if not stripped:
                continue
            
            if match_reg:
                if stripped.endswith('(Summary)'):
                    recording = False
                    continue
                next(file)  # Skip the blank line
                current_case = (match_reg.group(1), match_reg.group(2))
                headings = next(file).strip().split('\t')
                headings.insert(0, "Line Number")
                loadcases[current_case].append(headings)
                recording = True
                continue
                
            if recording:
                data = stripped.split('\t')
                loadcases[current_case].append(data)
    
    # Now loop over the extracted data, trim for only the relevant data, and
    # format
    for key, value in loadcases.items():
        if key == "metadata":
            loadcases["metadata"]["units"] = {}
            continue
        
        if len(value[0]) == 9:
            raise NotImplementedError("Need to implement for old case")
        elif len(value[0]) == 12:
            slice_obj = slice(2, 8, 1)
        else:
            raise NotImplementedError()
        if loadcases["metadata"]["data_type"] == "reaction":
            # To remove the sum row.
            del value[-1]
 
        for i in range(len(value)):
            temp = value[i][slice_obj]
            if i == 0:
                units = []
                clean_labels = []

                for label in temp:
                    match = re.search(r'\[(.*?)\]', label)
                    if match:
                        units.append(match.group(1))  # Get the unit inside []
                        clean_labels.append(re.sub(r'\[.*?\]', '', label))
                    else:
                        units.append(None)  # No unit present
                        clean_labels.append(label)
                    
                value[i] = clean_labels
                loadcases["metadata"]["units"][key] = units
                continue
            
            value[i] = [cast_str_to_num(val) for val in temp]
    
    return loadcases


def cast_str_to_num(val):
    """
    Handy function used for casting the results entries from strings to ints or
    floats
    """
    if '.' in val or 'e' in val or 'E' in val:
        return float(val)
    else:
        return int(val)


def reform_into_feature_schema(loadcases, save_mongo=True, save_csv=False):
    """
    """
    params_path = os.path.join(os.getcwd(), "instance", "fe_run_params.json")
    if not os.path.exists(params_path):
        raise FileNotFoundError()
    
    with open(params_path) as f:
        fe_params = json.load(f)
    
    type_root = "spatial"
    type_header = {
        "name": type_root,
        "type": {"name": loadcases["metadata"]["data_type"]}
    }
    
    si_folder = os.path.join(os.getcwd(), "si-models")
    filename = f"{loadcases['metadata']['structure_name']}.json"
    file_path = os.path.join(si_folder, filename)

    if os.path.isfile(file_path):
        with open(file_path, 'r') as f_si:
            si_data = json.load(f_si)
    else:
        raise FileNotFoundError
        
    
    for key, value in loadcases.items():
        if "DeadLoad" in key[1]:
            variant_id = "static"
            scenario_id = "1" if loadcases["metadata"]["data_type"] == "displacement" else "2"
            feature_name = f"spatial-{loadcases['metadata']['data_type']}-whole-structure-dead-load"
            environment = {
                "load": {
                    "description": "Structure is loaded under own weight in z axis in m/s^2",
                    "value": 9.81
                }
            }
        elif "Span" in key[1]:
            
            centre_span_locations = get_centre_span_locations(si_data)
            span_id = re.search(r"Span(\d+)_", key[1]).group(1)
            variant_id = "static"
            scenario_id = "3"
            feature_name = f"spatial-{loadcases['metadata']['data_type']}-whole-structure-{span_id}-midspan-point-load"
            loc = centre_span_locations[int(span_id) - 1]
            environment = {
               "loadValue": {
                   "description": f"Point load at midspan of span {span_id} in -z axis in kN",
                   "value": 400
               },
               "loadLocationX": {"description": "X coordinates of load position in meters", "value": loc[0]},
               "loadLocationY": {"description": "Y coordinates of load position in meters", "value": loc[1]},
               "loadLocationZ": {"description": "Z coordinates of load position in meters", "value": loc[2]},
            }
        else:
            continue
        parts = loadcases["metadata"]["structure_name"].split('-')
        parts.insert(-1, scenario_id)
        structure_name_filename = '-'.join(parts)
        
        # Handle the features section
        transposed = list(zip(*value))
        features = {
            "name": feature_name,
            "type": type_header,
            "variant": variant_id,
            "coordinates": { "global": { "translational": {
                axis: {
                    "indices": {
                        "start": 0,
                        "end": len(coord_row) - 1
                    },
                    "vector": list(coord_row[1:])
                }
                for coord_row, axis in zip(transposed[:3], ('x', 'y', 'z'))
            }}},
            "values": {
                axis: {
                    "indices": {
                        "start": 0,
                        "end": len(value_row) - 1
                    },
                    "vector": list(value_row[1:])
                }
                for value_row, axis in zip(transposed[3:], ('x', 'y', 'z'))
            }
        }
        features["coordinates"]["global"]["translational"]["unit"] = 'm'
        features["values"]["unit"] = loadcases["metadata"]["units"][key][-1]
        
        # Create the output JSON
        json_file = {
            "version": si_data["version"],
            "name": structure_name_filename,
            "population": si_data["population"],
            "timestamps": loadcases["metadata"]["timestamps"],
            "selection": {"name": loadcases["metadata"]["structure_name"]},
            "software": {
                "name": "IE2FE",
                "version": "0.1",
                "source": "https://github.com/TristanGowdridge/IE2FE",
                "parameters": fe_params
            },
            "environment": environment,
            "features": features
        }
        
        if save_csv:
            save_name = structure_name_filename + "-" + feature_name + ".json"
            save_folder = os.path.join(os.getcwd(), "parsed_feature_data")
            if not os.path.isdir(save_folder):
                os.mkdir(save_folder)
            save_path = os.path.join(save_folder, save_name)
            with open(save_path, "w") as f:
                json.dump(json_file, f, indent=4)
        
        if save_mongo:
            collection_obj = db_get_collection()
            collection_obj.insert_one(json_file)
        

def db_get_collection():
    """
    Returns a connection to the MongoDB collection.
    """
    credentials_path = os.path.join(os.getcwd(), "instance", "database_credentials.json")
    if os.path.isfile(credentials_path):
        with open(credentials_path, 'r') as f_cred:
            config = json.load(f_cred)

    mongodb_uri = "mongodb://"
    mongodb_uri += f"{quote_plus(config['authentication']['username'])}:"
    mongodb_uri += f"{quote_plus(config['authentication']['password'])}@"
    mongodb_uri += f"{config['hostname']}:{config['port']}/"
    mongodb_uri += f"{config['authentication']['database']}"
    return pymongo.MongoClient(mongodb_uri)[config["database"]][config["collection"]]


if __name__ == "__main__":
    t0 = time.time()
    fe_outputs_folder = r"C:\Users\trist\Desktop\University of Sheffield\ROSEHIPS\IE2FE\fe_outputs"
    for filename in os.listdir(fe_outputs_folder):
        loadcases = parse_loadcases(fe_outputs_folder, filename)
        reform_into_feature_schema(loadcases, save_csv=True)
    print(f"Time to execute: {time.time() - t0:.2f}s")
    

