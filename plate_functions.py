# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:29:04 2024

@author: Tristan Gowdridge

The code to create all of the Surface geometries used in the IE to FE models.
"""


def rectangular(lusas_sesh, **measurements):
    """
    Create a rectangular solid plate with the given measurements.
    """
    # Create the length of the beam
    thickness = measurements["thickness"]
    geom_identifier = f"LGeo Rectangular Solid Plate {thickness}"
    attr = lusas_sesh.database.createGeometricSurface(geom_identifier)
    # Thickness and eccentricity are the parameters passed into the set surface
    # for now. The eccentricity is changed further down the line, but this
    # cannot be done at this stage as information regarding the beam geometires
    # is required.
    attr.setSurface(f"{measurements['thickness']}", 0.0)

    return geom_identifier