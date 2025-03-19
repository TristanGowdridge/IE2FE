# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 12:48:16 2024

@author: trist
"""

import json
import time
from urllib.parse import quote_plus
import pymongo

 

def parse_loadcase_data(filename):
    """
    Parses the given file to extract only the loadcase node data of the first
    loadcase.
    
    Returns:
        A list of dictionaries, in the JSON format determined by the MongoDB
        Schema each dictionary corresponding to one node line of the loadcase.
    """
    timestamp = time.time_ns()
    loadcase_data = []
    in_headings = False
    in_loadcase = False
    summary_start_str = "1:Loadcase 1(Summary)"

    with open(f"{filename}_results.txt", 'r', encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            

            # Check if we've reached the first loadcase
            if line.startswith("1:Loadcase 1") and line != summary_start_str:
                in_model_info = False
                in_headings = True
                continue

            if in_headings:
                headings = []
                for s in [p for p in line.split('\t') if p.strip()]:
                    start = s.find('[')
                    end = s.find(']')
                    measurement_name = s[:start].strip()
                    unit = s[start+1:end].strip()
                    headings.append((measurement_name, unit))
                headings.insert(0, ("index", None))
                
                in_headings = False
                in_loadcase = True
                continue
            
            # Once we reach the summary, stop parsing loadcase data
            if line.startswith(summary_start_str):
                break
            
            if in_loadcase:
                # Attempt to parse a node data line
                readings = [p for p in line.split('\t') if p.strip()]
                
                channels = [
                    {
                        "name": measurement[0],
                        "type": "displacement",
                        "unit": measurement[1],
                        "value": float(value.replace('E', 'e'))
                    }
                    for measurement, value in list(zip(headings, readings))[2:]
                ]
                
                json_format = {
                    "version": "1.1.0",
                    "name": f"MY bridge, node {readings[0]}",
                    "population": "PEAR Bridges",
                    "timestamp": timestamp,
                    "channels": channels
                }
                loadcase_data.append(json_format)

    return loadcase_data


def get_structure_collection(collection_name):
    """
    Returns a connection to the MongoDB Database.
    """
    config = json.load(open("database.json"))
    mongodb_uri = "mongodb://"
    mongodb_uri += f"{quote_plus(config['authentication']['username'])}:"
    mongodb_uri += f"{quote_plus(config['authentication']['password'])}@"
    mongodb_uri += f"{config['hostname']}:{config['port']}/"
    mongodb_uri += f"{config['authentication']['database']}"
    
    db = pymongo.MongoClient(mongodb_uri)[config["database"]]
    db.list_collection_names()  # To raise an error if incorrect credentials.
    
    return db[collection_name]


def add_to_database(data):
    """
    
    """
    collection = get_structure_collection("IEtoFE")
    collection.insert_many(data)