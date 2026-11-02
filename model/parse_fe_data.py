# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 11:23:51 2025

@author: trist
"""

import json
import math
import os
import re
import warnings
from datetime import datetime, timezone
from typing import Any


def cast_str_to_num(val):
    """
    Handy function used for casting the results entries from strings to ints or
    floats
    """
    return float(val) if "." in val or "e" in val or "E" in val else int(val)


def create_schema_document(
    model_name: str,
    population_name: str,
    health_state: int,
    timestamp: int,
    parameters: dict,
) -> dict:
    """Create empty schema document

    Args:
        model_name (str): model name e.g bridge-1-1-1
        population_name (str): population name e.g pear-bridge-beam-and-slab-1
        health_state (int): numerical health state of features e.g 0
        timestamp (int): nanoseconds since epoch
        parameters (dict): parameters used to generate the features

    Returns:
        dict: PBSHM Schema document
    """
    return {
        "version": "1.4.1",
        "name": model_name,
        "population": population_name,
        "health": {"state": health_state},
        "source": {
            "timestamp": {"generated": timestamp, "stored": timestamp},
            "selection": [{"name": model_name}],
            "software": [
                {
                    "name": "IE2FE",
                    "version": "v1.0",
                    "source": "https://github.com/tristanGowdridge/IE2FE",
                    "parameters": parameters,
                }
            ],
            "environment": {
                "dataset": {
                    "description": "the dataset the model was taken from",
                    "value": "pear",
                }
            },
        },
        "features": [],
    }


def create_spatial_feature(
    scenario_name: str,
    type: str,
    x_pos: list[float],
    y_pos: list[float],
    z_pos: list[float],
    x_val: list[float],
    y_val: list[float],
    z_val: list[float],
) -> dict:
    """Create spatial feature

    Args:
        scenario_name (str): scenario name e.g deadload
        type (str): spatial type e.g displacement
        x_pos (list[float]): list of values for x positions
        y_pos (list[float]): list of values for y positions
        z_pos (list[float]): list of values for z positions
        x_val (list[float]): list of values for x values
        y_val (list[float]): list of values for y values
        z_val (list[float]): list of values for z values

    Raises:
        Exception: invalid lengths

    Returns:
        dict: PBSHM feature dictionary
    """
    # check equal lengths
    if not (
        len(x_pos) == len(y_pos) == len(z_pos) == len(x_val) == len(y_val) == len(z_val)
    ):
        raise Exception("X,Y,Z values do not have equal lengths")

    return {
        "name": scenario_name,
        "type": {"name": "spatial", "type": {"name": type}},
        "variant": "static",
        "coordinates": {
            "global": {
                "translational": {
                    "unit": "m",
                    "x": {"indices": {"start": 0, "end": len(x_pos)}, "vector": x_pos},
                    "y": {"indices": {"start": 0, "end": len(y_pos)}, "vector": y_pos},
                    "z": {"indices": {"start": 0, "end": len(z_pos)}, "vector": z_pos},
                }
            }
        },
        "values": {
            "unit": "m",
            "x": {"indices": {"start": 0, "end": len(x_val)}, "vector": x_val},
            "y": {"indices": {"start": 0, "end": len(y_val)}, "vector": y_val},
            "z": {"indices": {"start": 0, "end": len(z_val)}, "vector": z_val},
        },
    }


def create_eigen_mode_feature(
    scenario_name: str,
    variation: str,
    natural_frequency: float,
    x_pos: list[float],
    y_pos: list[float],
    z_pos: list[float],
    x_val: list[float],
    y_val: list[float],
    z_val: list[float],
) -> dict:
    """Create eigen mode feature

    Args:
        scenario_name (str): scenario name e.g deadload
        variation (str): eigen mode variation e.g acceleration
        natural_frequency (float): natural frequency
        x_pos (list[float]): list of values for x positions
        y_pos (list[float]): list of values for y positions
        z_pos (list[float]): list of values for z positions
        x_val (list[float]): list of values for x values
        y_val (list[float]): list of values for y values
        z_val (list[float]): list of values for z values

    Raises:
        Exception: invalid lengths

    Returns:
        dict: PBSHM feature dictionary
    """
    # check equal lengths
    if not (
        len(x_pos) == len(y_pos) == len(z_pos) == len(x_val) == len(y_val) == len(z_val)
    ):
        raise Exception("X,Y,Z values do not have equal lengths")

    return {
        "name": scenario_name,
        "type": "eigenMode",
        "variant": {"name": variation},
        "naturalFrequency": {"unit": "Hz", "value": natural_frequency},
        "modeShape": {
            "coordinates": {
                "global": {
                    "translational": {
                        "unit": "m",
                        "x": {
                            "indices": {"start": 0, "end": len(x_pos)},
                            "vector": x_pos,
                        },
                        "y": {
                            "indices": {"start": 0, "end": len(y_pos)},
                            "vector": y_pos,
                        },
                        "z": {
                            "indices": {"start": 0, "end": len(z_pos)},
                            "vector": z_pos,
                        },
                    }
                }
            },
            "values": {
                "unit": "m",
                "x": {"indices": {"start": 0, "end": len(x_val)}, "vector": x_val},
                "y": {"indices": {"start": 0, "end": len(y_val)}, "vector": y_val},
                "z": {"indices": {"start": 0, "end": len(z_val)}, "vector": z_val},
            },
        },
    }


def process_results(
    filename: str,
) -> tuple[str, int, list[int, int, str, list[str], list[list[Any]]]]:
    """process results into tuple

    Args:
        filename (str): results filename

    Returns:
        tuple[str, int, list[int, int, str, list[str], list[list[Any]]]]: results object
    """
    # headers
    header_regex = [
        (re.compile(r"^[0-9]+\tModel title\t([a-zA-Z0-9-]+)$"), 1),
        (
            re.compile(r"^[0-9]+\tDate\t([0-9]{1,2})/([0-9]{1,2})/([0-9]{1,2})$"),
            (1, 2, 3),
        ),
        (re.compile(r"^[0-9]+\tTime\t([0-9]{1,2}):([0-9]{1,2})$"), (1, 2)),
    ]
    header_values = [None for _ in header_regex]

    # read file line by line
    sections = []
    section_re = re.compile(r"^([0-9]+):([)(_ .=a-zA-Z0-9+-]+)$")
    number_re = re.compile(r"^-?[0-9]+(\.[0-9]+(E[+-]?[0-9]+)?)?$")
    with open(filename, "r") as fs:
        index = 0
        for line in fs:

            # ignore empty lines
            contents = line.strip()
            if len(contents) == 0:
                continue

            # process file header info
            if len(sections) == 0:
                for h_index, regex in enumerate(header_regex):
                    if header_values[h_index] == None:
                        re_match = regex[0].fullmatch(contents)
                        if re_match:
                            header_values[h_index] = (
                                re_match.group(*regex[1])
                                if isinstance(regex[1], tuple)
                                else re_match.group(regex[1])
                            )

            # process start of sections
            section_match = section_re.fullmatch(contents)
            if section_match:
                if len(sections) > 0:
                    sections[-1]["end"] = index - 1
                sections.append(
                    {
                        "start": index,
                        "end": None,
                        "name": section_match.group(2),
                        "headers": [],
                        "values": [],
                    }
                )

            # process section headers
            elif len(sections) > 0 and index - 1 == sections[-1]["start"]:
                sections[-1]["headers"] = contents.split("\t")

            # process section values
            elif len(sections) > 0:
                values = contents.split("\t")
                sections[-1]["values"].append(
                    [
                        (
                            cast_str_to_num(value)
                            if number_re.fullmatch(value)
                            else value
                        )
                        for value in values[1:]
                    ]
                )

            # only increment index on valid lines
            index += 1

        # close of final section
        if len(sections) > 0:
            sections[-1]["end"] = index

    # convert date time object
    timestamp = 0
    if len(header_values[1]) == 3 and len(header_values[2]) == 2:
        dt = datetime(
            int(header_values[1][2]) + 2000,  # year, month, day
            int(header_values[1][1]),
            int(header_values[1][0]),
            int(header_values[2][0]),  # hours, minites
            int(header_values[2][1]),
        )
        dt = dt.replace(tzinfo=timezone.utc)
        timestamp = int(dt.timestamp() * 1e9)

    # return data
    return (header_values[0], timestamp, sections)


class FeatureWriter:

    def __init__(
        self,
        model_name: str,
        population_name: str,
        health_state: int,
        parameters: dict,
        scenarios: list[str] = [],
    ) -> "FeatureWriter":
        """create new feature writer to convert fe results into PBSHM features

        Args:
            model_name (str): model name
            population_name (str): population name
            health_state (int): document health state
            parameters (dict): parameters dictionary
            scenarios (list[str], optional): scenario filter. Defaults to [].

        Returns:
            FeatureWriter: _description_
        """
        self.model_name = model_name
        self.population_name = population_name
        self.health_state = health_state
        self.parameters = parameters
        self.scenarios = scenarios
        self.features = []
        self.timestamps = []

    def add_displacement(self, filename: str) -> None:
        """add displacement data to features to be written

        Args:
            filename (str): filename for data to be read from
        """
        results = process_results(filename)
        if results[0] != self.model_name:
            warnings.warn(
                f"Displacement file model conflict {self.model_name}:{results[0]}"
            )
        for segment in results[2]:
            # exclude summary lines
            if segment["name"].endswith("(Summary)"):
                continue
            # exlucde scenarios not in filter
            elif len(self.scenarios) > 0 and segment["name"] not in self.scenarios:
                continue
            x_pos, y_pos, z_pos = [], [], []
            x_val, y_val, z_val = [], [], []
            for line in segment["values"]:
                x_pos.append(line[1])
                y_pos.append(line[2])
                z_pos.append(line[3])
                x_val.append(line[4])
                y_val.append(line[5])
                z_val.append(line[6])
            self.timestamps.append(results[1])
            self.features.append(
                create_spatial_feature(
                    segment["name"].lower().replace("_", "-"),
                    "displacement",
                    x_pos,
                    y_pos,
                    z_pos,
                    x_val,
                    y_val,
                    z_val,
                )
            )

    def add_reaction(self, filename: str) -> None:
        """add reaction data to features to be written

        Args:
            filename (str): filename for data to be read from
        """
        results = process_results(filename)
        if results[0] != self.model_name:
            warnings.warn(
                f"Reaction file model conflict {self.model_name}:{results[0]}"
            )
        for segment in results[2]:
            # exclude summary lines
            if segment["name"].endswith("(Summary)"):
                continue
            # exlucde scenarios not in filter
            elif len(self.scenarios) > 0 and segment["name"] not in self.scenarios:
                continue
            x_pos, y_pos, z_pos = [], [], []
            x_val, y_val, z_val = [], [], []
            for line in segment["values"]:
                # ignore total line
                if line[0] == isinstance(line[0], str) and line[0].lower() == "total":
                    continue
                # append data
                x_pos.append(line[1])
                y_pos.append(line[2])
                z_pos.append(line[3])
                x_val.append(line[4])
                y_val.append(line[5])
                z_val.append(line[6])
            self.timestamps.append(results[1])
            self.features.append(
                create_spatial_feature(
                    segment["name"].lower().replace("_", "-"),
                    "reaction",
                    x_pos,
                    y_pos,
                    z_pos,
                    x_val,
                    y_val,
                    z_val,
                )
            )

    def add_eigen_mode(
        self, frequencies_filename: str, modeshape_filename: str
    ) -> None:
        # get frequency details
        results = process_results(frequencies_filename)
        if results[0] != self.model_name:
            warnings.warn(
                f"Frequency file model conflict {self.model_name}:{results[0]}"
            )
        frequencies = {}
        for segment in results[2]:
            if segment["name"].lower().endswith("(eigenvalues)"):
                for line in segment["values"]:
                    if line[0] in frequencies:
                        frequencies[line[0]][0] = line[2]
                    else:
                        frequencies[line[0]] = [line[2], None]
            elif segment["name"].lower().endswith("(mass participation factors)"):
                for line in segment["values"]:
                    max_value = max(
                        [line[1], line[2], line[3], line[4], line[5], line[6]]
                    )
                    if line[0] in frequencies:
                        frequencies[line[0]][1] = max_value
                    else:
                        frequencies[line[0]] = [None, max_value]
        # process mode shapes
        results = process_results(modeshape_filename)
        if results[0] != self.model_name:
            warnings.warn(
                f"Modeshape file model conflict {self.model_name}:{results[0]}"
            )
        scenario_re = re.compile(
            r"^[Mm]ode ([0-9]+) [Ff]requency = ([0-9]+(\.[0-9]+(E[+-]?[0-9]+)?)?)$"
        )
        for segment in results[2]:
            # exclude any scenario names that don't match
            scenario_match = scenario_re.fullmatch(segment["name"])
            if not scenario_match:
                continue
            # ensure a mode is linked in the frequencies
            mode_number = int(scenario_match.group(1))
            if mode_number not in frequencies:
                continue
            # ensure the frequencies are close enough
            if not math.isclose(
                frequencies[mode_number][0],
                float(scenario_match.group(2)),
                rel_tol=1e-8,
            ):
                continue
            # exclude any below frequency threshold
            if (
                "freqThreshold" in self.parameters
                and frequencies[mode_number][0] < self.parameters["freqThreshold"]
            ):
                continue
            # enxclude any below mass threshold
            if (
                "massThreshold" in self.parameters
                and frequencies[mode_number][1] < self.parameters["massThreshold"]
            ):
                continue
            # process data
            x_pos, y_pos, z_pos = [], [], []
            x_val, y_val, z_val = [], [], []
            for line in segment["values"]:
                x_pos.append(line[1])
                y_pos.append(line[2])
                z_pos.append(line[3])
                x_val.append(line[4])
                y_val.append(line[5])
                z_val.append(line[6])
            # create feature
            self.timestamps.append(results[1])
            self.features.append(
                create_eigen_mode_feature(
                    f"mode {mode_number}",
                    "acceleration",
                    frequencies[mode_number][0],
                    x_pos,
                    y_pos,
                    z_pos,
                    x_val,
                    y_val,
                    z_val,
                )
            )

    def write(self, features_directory: str) -> int:
        """write features to directory

        Args:
            features_directory (str): directiory path

        Returns:
            int: number of files written
        """
        # ensure directory
        if not os.path.exists(features_directory):
            os.mkdir(features_directory)

        # create feature groups
        feature_groups = {}
        manual_groupings = self.parameters["featureGroups"]
        if len(manual_groupings) > 0 and len(manual_groupings) % 2 == 0:
            for i in range(0, int(len(manual_groupings) / 2)):
                feature_groups[manual_groupings[i * 2]] = manual_groupings[(i * 2) + 1]
        else:
            for feature in self.features:
                key = feature["name"].lower().replace(" ", "-")
                if key not in feature_groups:
                    feature_groups[key] = f"^{key}$"

        # output features into directories
        for group in feature_groups:
            group_dir = os.path.join(features_directory, group.lower().replace(" ", ""))
            if not os.path.exists(group_dir):
                os.mkdir(group_dir)
            group_re = re.compile(feature_groups[group])
            group_feature_keys = {}
            for index, feature in enumerate(self.features):
                key = feature["name"].lower().replace(" ", "-")
                if not group_re.fullmatch(key):
                    continue
                if key not in group_feature_keys:
                    group_feature_keys[key] = [index]
                else:
                    group_feature_keys[key].append(index)
            for key in group_feature_keys:
                document = create_schema_document(
                    self.model_name,
                    self.population_name,
                    self.health_state,
                    self.timestamps[group_feature_keys[key][0]],
                    self.parameters,
                )
                for index in group_feature_keys[key]:
                    document["features"].append(self.features[index])
                with open(
                    os.path.join(group_dir, f"{self.model_name}-{key}.json"), "w"
                ) as fs:
                    json.dump(document, fs)

        return len(feature_groups.items())
