# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:29:04 2024

@author: Tristan Gowdridge

Defines all surface geometries used in the IE to FE model conversion.
"""


def rectangular(lusas_sesh, **measurements):
    """
    Create a rectangular solid plate geometry using the specified thickness.
    This function is primarily used for the decks of bridges.

    :param lusas_sesh: The active LUSAS session object.
    :type lusas_sesh: object
    :param measurements: Dictionary of geometric parameters.
    :type measurements: dict

    The expected key in ``measurements`` is:
        - ``thickness``: The thickness of the plate.

    :return: Identifier string for the created rectangular plate geometry.
    :rtype: str
    """
    # Create the plate geometry using the specified thickness
    thickness = measurements["thickness"]
    geom_identifier = f"LGeo Rectangular Solid Plate {thickness}"
    attr = lusas_sesh.database.createGeometricSurface(geom_identifier)
    # Set the surface using thickness and an initial eccentricity of 0.0.
    # Eccentricity is later updated once beam geometries are defined.
    attr.setSurface(f"{measurements['thickness']}", 0.0)

    return geom_identifier
