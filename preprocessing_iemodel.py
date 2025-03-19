# -*- coding: utf-8 -*-
"""
Created on Mon Nov 4 15:14:07 2024

@author: Tristan Gowdridge

The functions contained within this script serve as a preprocessing layer to
extract the IE model information in a Python-friendly format which can then be
easily passed on to any Finite Element API. Currently, the only FE API
implemented is for LUSAS, therefore references to LUSAS are heavy, despite this
code being nonspecific to LUSAS.
"""
import json
import os
from os.path import exists, join
import warnings

from ie_to_fe_funcs import distance

GEOMETRY_TYPES = {
    "beam": "Line",
    "column": "Line",
    "slab": "Surface"
}


class ElementInfo():
    """
    A class used to represent all the information per element within the IE
    model.
    """
    
    def __init__(self, name, element_type):
        self.name = name
        self.element_type = element_type
        self.coords = set()
        self.geometry = ""
        self.ground_coord = set()

    def __str__(self):
        my_string = f"{self.name}: {self.element_type}"
        for attribute, value in self.__dict__.items():
            my_string += f"\n{attribute} = {value}"
        
        return my_string
    
    def __repr__(self):
        return self.__str__()
            

def get_material_properties(element):
    """
    Extract material properties from the element.
    
    LUSAS requires a minimum of a Young's modulus, density, Poisson's ratio,
    and linear thermal expansion to be defined. These values are extracted from
    the IE model. However, if the IE model does not have this information, the
    missing values are supplied via some default values defined in
    ./data/materials.json.
    
    Needs scaling
    """
    if element["type"] == "ground":
        # If it's a ground elenent, then it has no material properties.
        return {}
    
    element_material = element.get("material")
    try:
        material_name = element_material["matrix"]["base"]["type"]["name"]
    except KeyError:
        material_name = element_material["type"]["type"]["type"]["name"]
    
    lusas_requirements = (
        "youngsModulus",
        "density",
        "poissonsRatio",
        "linearThermalExpansionCoefficient"
    )
    material_properties = {}
    if "properties" in element_material:
        for required_property in lusas_requirements:
            for available_property in element_material["properties"]:
                if required_property == available_property["type"]:
                    temp_value = available_property["value"]
                    # temp_unit = available_property["unit"]
                    # material_properties[required_property] = unit_scaler(temp_value, temp_unit)
                    material_properties[required_property] = temp_value
                    break
    
    # If there is not enough information regarding the minimum LUSAS material
    # requirements, then use a look-up table fill in the blanks.
    if any(prop not in material_properties.keys() for prop in lusas_requirements):
        # Find the absolute path to the materials.json file used to supply the
        # default values.
        materials_json_path = join(os.getcwd(), "data", "materials.json")

        # Make sure the path is correct and the file exists
        if not exists(materials_json_path):
            raise FileNotFoundError(
                f"Material data file not found: {materials_json_path}."
            )

        with open(materials_json_path, 'r') as f:
            material_defaults = json.load(f)
        
        # Check if the material is present in the default values found in 
        # materials.json
        if material_name not in material_defaults:
            raise KeyError(
                "Not enough material information has been supplied in the IE "
                "model. The requirements are at least "
                ", ".join(lusas_requirements) + ", and default values are not"
                f" provided for `{material_name}'."
            )
        
        for attribute, value in material_defaults[material_name].items():
            if attribute not in material_properties:
                #
                # Shouldn't need scaling but if LUSAS sesh loaded with
                # nondefault values, the material defaults might be wrong.
                #
                material_properties[attribute] = value
    
    return material_name, material_properties


def get_geometric_properties(element):
    """
    Extract geometric properties from the element.
    
    Needs scaling
    """
    dimensions = {
        k: v for k, v in element["geometry"]["dimensions"].items()
        if k != "length"
    }
    dimensions["profile"] = (
        element["geometry"]["type"]["type"]["type"]["name"],
        element["geometry"]["type"]["type"]["name"]
    )
    
    return dimensions

    
def extract_coordinates(coord_data):
    """
    Need to add scaling here
    """
    return tuple(coord_data[axis]["value"] for axis in ('x', 'y', 'z'))


def extract_profile_info(element):
    """
    Extract profile information from the element.
    
    Needs scaling
    """
    profile_info = []
    for key, measurement in element.profile.items():
        if key == "profile":
            continue
        value = measurement["value"]
        #
        # Need to apply scaling
        #
        profile_info.append((key, value))
    profile_info = sorted(profile_info, key=lambda x: x[0])
    profile_info.insert(0, element.profile["profile"])
    return tuple(profile_info)


def get_line_end_coords(coord_center, connected_axis, length_value):
    """
    If there is a Line element connected to a Surface, the midpoint of their
    connection is used as the relationship within the IE model. Therefore, need
    to determine the end coordinates along the length of the surface.
    
    Need to add scaling here
    """
    coord1 = list(extract_coordinates(coord_center))
    coord2 = list(extract_coordinates(coord_center))
    axis_index = {'x': 0, 'y': 1, 'z': 2}[connected_axis]
    half_length = length_value / 2
    coord1[axis_index] -= half_length
    coord2[axis_index] += half_length
    return tuple(coord1), tuple(coord2)


def join_close_points(ie_table, threshold=1e-4):
    """
    When creating the IE models, some numerical precision issues are being
    introduced, therefore points are being created very close to one another
    when in fact they are the same point. This function combines the points
    that are within some threshold value and keeps the `simplest' point (the
    one with the fewest significant figures, therefore the one less likely to
    be affected by numerical precision errors).
    """
    # Post processing step to combine very close points ** BODGY **
    point_set = set()
    for element in ie_table.values():
        point_set.update(element.coords)
    
    point_map = dict()
    point_list = list(point_set)
    for point1 in point_set:
        close_points = []
        for point2 in point_list:
            if distance(point1, point2) < threshold:
                close_points.append(point2)
        point_map[point1] = min(close_points, key=lambda x: len(str(x)))
    
    for old_point, new_point in point_map.items():
        if old_point != new_point:  # Only different points are processed
            for element in ie_table.values():
                if old_point in element.coords:
                    element.coords.remove(old_point)
                    element.coords.add(new_point)


def extract_ie_information(filename):
    """
    Process the IE model JSON and extract relevant information into a
    dictionary keyed by the name of the element with the value of the
    ElementInfo class.
    """
    # Load the IE model JSON file
    if not exists(filename):
        raise FileNotFoundError(f"Cannot find the IE model {filename}")
    
    with open(filename) as _f:
        ie_model = json.load(_f)["models"]["irreducibleElement"]
    
    # Use a dictionary for speedy lookups and easier name-based indexing at a
    # later stage.
    elements = {
        element["name"]: dict(element)
        for element in ie_model["elements"]
    }
    relationships = {
        relationship["name"]: dict(relationship)
        for relationship in ie_model["relationships"]
    }
    
    # From the relationships, determine all the connections and the geometry of
    # the structure.
    ie_table = {}
    for relationship in relationships.values():
        relationship_element1, relationship_element2 = relationship["elements"]
        
        # Get the names involved in the relationships, so more information can
        # be obtained from the elements dictionary.
        element1_name = relationship_element1["name"]
        element2_name = relationship_element2["name"]
        
        # Use the element names from the relationship to obtain the element
        # objects.
        element1 = elements[element1_name]
        element2 = elements[element2_name]
        
        # If the element isn't already being tracked, and is not a ground
        # element, add to the ie_table.
        for element in [element1, element2]:
            ele_name, ele_type = element["name"], element["type"]
            if (ele_name not in ie_table) and (ele_type != "ground"):
                ie_table[ele_name] = ElementInfo(ele_name, ele_type)
        
        # Determine the geometry of the relationship element (Point, Line,
        # Surface, Volume). Coordinates are appended to a set, these are later
        # used to extract the geometry of each element. The relationships
        # between different geometries have different processes of extracting
        # the coordinates defining the element geometry.
        geometries = []
        for element in [element1, element2]:
            if element["type"] == "ground":
                if element["name"] == element1_name:
                    grounded_element_name = element2_name
                else:
                    grounded_element_name = element1_name
                ground_coord = extract_coordinates(
                    relationship["coordinates"]["global"]["translational"]
                )
                ie_table[grounded_element_name].ground_coord.add(ground_coord)
                continue
            element_geometry = GEOMETRY_TYPES[element["contextual"]["type"]]
            ie_table[element["name"]].geometry = element_geometry
            geometries.append(element_geometry)
        
        if all(geometry == "Line" for geometry in geometries):
            if "coordinates" in relationship:
                for_coord = relationship["coordinates"]["global"]["translational"]
            else:
                for_coord = relationship_element1["coordinates"]["global"]["translational"]
            coord = extract_coordinates(for_coord)
            
            for element in (element1, element2):
                if element["type"] != "ground":
                    ie_table[element["name"]].coords.add(coord)
        
        elif geometries.count("Line") == 1:
            # In the case of a line intersecting a higher-dimensional geometry,
            # we need to dtermine the end point of the line for its
            # coordinates.
            line_index = geometries.index("Line")
            line_element = (element1, element2)[line_index]
            length_axis = line_element["geometry"]["dimensions"]["length"]["axis"]
            length_value = line_element["geometry"]["dimensions"]["length"]["value"]
            
            if "coordinates" in relationship:
                rel_coord = relationship["coordinates"]["global"]["translational"]
            else:
                rel_coord = relationship["elements"][line_index]["coordinates"]["global"]["translational"]
            
            p1, p2 = get_line_end_coords(rel_coord, length_axis, length_value)

            for element in (element1, element2):
                if element["type"] != "ground":
                    ie_table[element["name"]].coords.add(p1)
                    ie_table[element["name"]].coords.add(p2)
        
        elif geometries.count("Surface") == 2:
            # Connection between two surfaces has to be extracted more manually
            # by taking out from the element geometries, this process is more
            # involved, error prone and currently not required since all the
            # line geometries in the IE models being created for PEAR have an
            # attributed geometry. However, this is not neccessarily the case,
            # and this will almost certainly cause issues. (tl;dr a later
            # problem).
            warnings.warn(
                "There is a surface to surface connection here. If the surface"
                " to surface connection does not have a line geometry, this"
                " will be a problem. There is a solution but it is not needed"
                " just yet (i.e. a later problem)"
            )
            continue
        else:
            raise NotImplementedError("Unknown geometry combination.")
    
    # Combines the points that are within a some threshold distance from one
    # another.
    join_close_points(ie_table, threshold=1e-4)
    
    # Extract the material, and geometric properties foe each element.
    for name, element in ie_table.items():
        element.material = get_material_properties(elements[name])
        element.profile = get_geometric_properties(elements[name])
    
    return ie_table
