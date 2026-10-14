# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 12:09:43 2025

@author: Tristan Gowdridge

This file contains all the logic for the translate and scale IE model profiles
when automatically creating the FE models from IE models in LUSAS.
"""


def tapered_cylindrical_shell(lusas_sesh, **measurements):
    """
    Constructs a tapered cylindrical shell profile by linearly interpolating
    between the start and end circular faces using the given measurements.

    :param lusas_sesh: The LUSAS session object.
    :type lusas_sesh: object
    :param measurements: A dictionary of geometry parameters.
    :type measurements: dict

    The expected keys in ``measurements`` are:
        - ``r1``: Radius of the starting circle.
        - ``t1``: Thickness of the starting circle.
        - ``r2``: Radius of the ending circle.
        - ``t2``: Thickness of the ending circle.

    :return: Identifier string for the created tapered cylinder.
    :rtype: str
    """
    # dim_values is used to create a unique profile name for reference in LUSAS
    dim_values = ', '.join(map(str, map(float, measurements.values())))
    geom_identifier = f"({dim_values})"
    cylinder_identifier = f"Tapered cylinder {geom_identifier}"
    
    # Section 1: define the start face by specifying its coordinates.
    # Reusing a non-unique name will overwrite it in subsequent profiles.
    names = ["D", "t"]
    values = [2*measurements["r1"], measurements["t1"]]
    dim_values = ', '.join(map(str, zip(names, values)))
    geom_identifier1 = f"Sct ({dim_values})"
    attr = lusas_sesh.database.createParametricSection(geom_identifier1)
    attr.setType("Circular Hollow")
    attr.setDimensions(names, values)
    
    # Section 2: define the end face.
    names = ["D", "t"]
    values = [2*measurements["r2"], measurements["t2"]]
    dim_values = ', '.join(map(str, zip(names, values)))
    geom_identifier2 = f"Sct ({dim_values})"
    attr = lusas_sesh.database.createParametricSection(geom_identifier2)
    attr.setType("Circular Hollow")
    attr.setDimensions(names, values)
    
    # Create a geometric line interpolating between the two sections
    attr = lusas_sesh.database.createGeometricLine(cylinder_identifier)
    attr.setValue("elementType", "3D Thick Beam")
    attr.setNumberOfSections(2)
    attr.setValue("interpMethod", "Use Section Calculator")
    attr.setFromLibrary("Utilities", "", geom_identifier1, 0, 0, 0)
    attr.setFromLibrary("Utilities", "", geom_identifier2, 0, 0, 1)
    attr.setVerticalAlignment("CenterToCenter")
    attr.setHorizontalAlignment("CenterToCenter")
    attr.setAlignmentSection(0)
    attr.setAnalysisCategory("3D")
    
    return cylinder_identifier
