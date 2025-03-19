# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 14:21:19 2024

@author: trist
"""
import os
from collections import defaultdict

from preprocessing_iemodel import extract_ie_information, extract_profile_info
from LUSAS_classes import LUSASSession
from ie_to_fe_funcs import order_coords
     
      
def run_LUSAS(filename: str):
    """
    This is the function that is runs the whole analysis, this is passed into
    the batch script.
    """

    ie_path = os.path.join(".", "ie_models", filename)
    if not os.path.exists(ie_path):
        raise FileNotFoundError(f"The file {ie_path} does not exist.")
        
    # all_in_one extracts all the information from a JSON IE model, and
    # converts it into a more handy table structure, which is used internally
    # in Python.
    ie_table = extract_ie_information(ie_path)
    
    # Instantiate the LUSAS project.
    lusas_sesh = LUSASSession(filename)
    
    geometry_collection = defaultdict(dict)
    all_materials, all_profiles = set(), set()
    for info in ie_table.values():
        geometry_collection[info.geometry][order_coords(info.coords)] = info
        all_materials.add(info.material[0])
        
        profile_info = extract_profile_info(info)
        all_profiles.add(profile_info)
        
        if info.ground_coord:
            lusas_sesh.ground_elements.extend(info.ground_coord)
    
    # Create the base of the structure, to which all subsequent information is
    # added to.
    lusas_sesh.create_geometry(geometry_collection)
    
    # From the wireframe, we now add specific information.
    lusas_sesh.assign_materials(ie_table, all_materials)
    lusas_sesh.assign_profiles(ie_table, all_profiles)
    lusas_sesh.assign_ground_elements()
    
    # Bridge model fine tuning
    lusas_sesh.apply_surface_eccentricity()
    lusas_sesh.correct_perimeter_line_direction_for_geometry()
 
    if "ladder" in bridge_identifier:
        obj = lusas_sesh.selection.add("Surface", "All").getObjects("Surface", "All")[0]
        ladder_beam_profiles = set()
        main_span_profiles = set()
        for lof in obj.getLOFs():
            if lof.getTypeCode() != 2:
                continue
            facet_coords = lof.getFacetCoordinates()
            p1, p2 = facet_coords[:3], facet_coords[3:]
            dx, dy = p2[0] - p1[0], p2[1] - p1[1]
            
            if abs(dy) < 1e-5:
                profile_name = lof.getAssignments("Geometric")[0].getAttribute().getName()
                main_span_profiles.add(profile_name)
            
            if abs(dx) < 1e-5:
                profile_name = lof.getAssignments("Geometric")[0].getAttribute().getName()
                ladder_beam_profiles.add(profile_name)
            
        if not set.intersection(ladder_beam_profiles, main_span_profiles):
            if len(main_span_profiles) != 1:
                raise ValueError("too many main spans")
            main_prof = tuple(main_span_profiles)[0]
            main_span_height = lusas_sesh.database.getAttribute("Geometric", main_prof).getValue("D")
            for prof in tuple(ladder_beam_profiles):
                attr = lusas_sesh.database.createGeometricLine(prof)
                
                attr.setValue("elementType", "3D Thick Beam")
                attr.setFromLibrary("Utilities", "", prof[5:], 0, 0, 0)
                attr.setEccentricityOrigin("Centroid", "Fibre", "", "I3")
                ladder_height = attr.getValue("D")
                attr.setValue("ez0", -(main_span_height-ladder_height), 0)
                attr.setAnalysisCategory("3D")
    
    # # Make basic meshes.
    lusas_sesh.create_basic_line_mesh(2)
    # lusas_sesh.create_equal_line_mesh(1)
    lusas_sesh.create_basic_surface_mesh(0, 0)
    lusas_sesh.create_basic_volume_mesh(2, 2, 2)
    
    # Assign the loadcases.
    lusas_sesh.assign_loadcase_from_lvb(bridge_identifier)
    
    # Save the results to txt.
    lusas_sesh.save_results(bridge_identifier)
    
    # # These need modifying...................................................
    # # parsed_data = parse_loadcase_data(project_name)
    # # add_to_database(parsed_data)
    
    lusas_sesh.save_model()
    # lusas_sesh.modeller.quit()
    
    
    
