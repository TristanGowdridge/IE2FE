# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 12:09:43 2025

@author: trist
"""


def tapered_cylindrical_shell(lusas_sesh, **measurements):
    """
    """
    dim_values = ', '.join(map(str, map(float, measurements.values())))
    geom_identifier = f"({dim_values})"
    cylinder_identifier = f"Tapered cylinder {geom_identifier}"
    
    # Section 1
    names = ["D", "t"]
    values = [2*measurements["r1"], measurements["t1"]]
    dim_values = ', '.join(map(str, zip(names, values)))
    geom_identifier1 = f"Sct ({dim_values})"
    attr = lusas_sesh.database.createParametricSection(geom_identifier1)
    attr.setType("Circular Hollow")
    attr.setDimensions(names, values)
    
    # Section 2
    names = ["D", "t"]
    values = [2*measurements["r2"], measurements["t2"]]
    dim_values = ', '.join(map(str, zip(names, values)))
    geom_identifier2 = f"Sct ({dim_values})"
    attr = lusas_sesh.database.createParametricSection(geom_identifier2)
    attr.setType("Circular Hollow")
    attr.setDimensions(names, values)
    
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
