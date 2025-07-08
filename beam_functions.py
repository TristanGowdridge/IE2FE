# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 14:59:42 2024

@author: Tristan Gowdridge

Defines the code to create all Line geometries used in IE to FE model
conversion.
"""
import json
import os
from os.path import exists, join


def _itc_beams_base(lusas_sesh, beam_geom, **measurements):
    """
    A private function used internally to handle the bulk of the logic to
    generate a geometric line for I, T, or C beams. The respective beams call
    this function with the relevant beam_geom.

    These beams share the same structure with minimal variation.

    :param lusas_sesh: Active LUSAS session.
    :param beam_geom: Beam type ("I", "T", or "C").
    :param measurements: Required keys:
        - width
        - height
        - webThickness
        - flangeThickness
    :return: Unique geometry identifier string.
    """
    beam_dims = {}
    ie_names = ("width", "height", "webThickness", "flangeThickness")
    lusas_names = ('B', 'D', 'tw', 'tf')
    for ie_name, lusas_name in zip(ie_names, lusas_names):
        beam_dims[lusas_name] = measurements[ie_name]
    beam_dims['r'] = 0  # Edge radius between web and flange

    # Create the beam geometry attribute
    dim_values = ', '.join(map(str, map(float, beam_dims.values())))
    geom_identifier = f"({dim_values})"
    ibeam_identifier = f"{beam_geom} Beam {geom_identifier}"
    
    beam_face = lusas_sesh.database.createParametricSection(ibeam_identifier)
    beam_face.setType(f'{beam_geom}')
    beam_face.setDimensions(tuple(beam_dims.keys()), tuple(beam_dims.values()))
    
    # Create the length of the beam
    geometry_identifier = f"LGeo {ibeam_identifier}"
    attr = lusas_sesh.database.createGeometricLine(geometry_identifier)
    attr.setValue("elementType", "3D Thick Beam")
    attr.setFromLibrary("Utilities", "", ibeam_identifier, 0, 0, 0)
    attr.setEccentricityOrigin("Centroid", "Fibre", "", "I3")
    attr.setAnalysisCategory("3D")

    return geometry_identifier


def i_beam(lusas_sesh, **measurements):
    """
    Create a geometric line for an I-beam.

    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 461

    :param measurements: Dictionary with:
        - width -> B
        - height -> D
        - webThickness -> tw
        - flangeThickness -> tf
    :return: Geometry identifier string.
    """
    return _itc_beams_base(lusas_sesh, 'I', **measurements)
        

def t_beam(lusas_sesh, **measurements):
    """
    Create a geometric line for a T-beam.
    
    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 462
    
    Required conditions:
        width > flangeThickness + 2r (r assumed 0)
        height > webThickness + r (r assumed 0)
    
    :param measurements: Dictionary with required keys.
    :return: Geometry identifier string.
    """
    return _itc_beams_base(lusas_sesh, 'T', **measurements)

    
def c_beam(lusas_sesh, **measurements):
    """
    Create a geometric line for a C-beam.
    
    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 46363
    
    :param measurements: Dictionary with:
        - width
        - height
        - webThickness
        - flangeThickness
    :return: Geometry identifier string.
    """
    return _itc_beams_base(lusas_sesh, 'C', **measurements)
       

def l_beam(lusas_sesh, **measurements):
    """
    Create a geometric line for an L-beam.
    
    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 507
    
    :param measurements: Dictionary with:
        - width
        - height
        - thickness
        - angle (not currently handled)
    :return: Geometry identifier string.
    """
    beam_dims = {}
    
    ie_names = ("width", "height", "thickness", "thickness")
    lusas_names = ('B', 'D', 'tw', 'tf')
    for ie_name, lusas_name in zip(ie_names, lusas_names):
        beam_dims[lusas_name] = measurements[ie_name]
    beam_dims['r1'] = 0
    beam_dims['r2'] = 0

    # Create the beam geometry attribute
    ibeam_identifier = f"L Beam {str(measurements)}"
    beam_face = lusas_sesh.database.createParametricSection(ibeam_identifier)
    beam_face.setType('L')
    beam_face.setDimensions(tuple(beam_dims.keys()), tuple(beam_dims.values()))
    
    # Create the length of the beam
    geometry_identifier = f"LGeo L Beam {str(measurements)}"
    attr = lusas_sesh.database.createGeometricLine(geometry_identifier)
    attr.setValue("elementType", "3D Thick Beam")
    attr.setFromLibrary("Utilities", "", ibeam_identifier, 0, 0, 0)
    attr.setEccentricityOrigin("Centroid", "Fibre", "", "I3")
    attr.setAnalysisCategory("3D")

    return geometry_identifier


def _yyem_base(lusas_sesh, beam_geom, **measurements):
    """
    Base function for Y, YE, and M beams using external JSON profiles.

    :param beam_geom: Beam type string, e.g., "Y", "YE", or "M".
    :param measurements: Dict with:
        - topWidth
        - baseWidth
        - height
    :return: Geometry identifier string.
    """
    beam_dims = []
    # Order of the key in the beam data json, vital for keying later
    for ie_name in ("topWidth", "baseWidth", "height"):
        beam_dims.append(measurements[ie_name])
    
    beam_path = join(os.getcwd(), "data", "nonparametric_beam_geometries.json")
    if not exists(beam_path):
        raise FileNotFoundError(
            "Cannot locate the file 'nonparametric_beam_geometries.json'"
        )
    
    with open(beam_path, 'r') as f:
        beam_data = json.load(f)
    dimension_identifier = f"({', '.join(map(str, map(float, beam_dims)))})"
    number = beam_data[beam_geom].get(dimension_identifier, None)['N']
    
    if number:
        beam_profile = f"{beam_geom}{number}"
    else:
        raise KeyError("This beam section is not supported.")
    
    # Create the length of the beam
    geometry_identifier = f"LGeo {beam_profile} Beam {str(measurements)}"
    attr = lusas_sesh.database.createGeometricLine(geometry_identifier)
    attr.setValue("elementType", "3D Thick Beam")
    attr.setFromLibrary(
        "UK Sections", f"Precast {beam_geom} Beams", beam_profile, 0, 0, 0
    )
    attr.setEccentricityOrigin("Centroid", "Fibre", "", "A3")
    # Bug with M-beam fibre locations not loading, need to manually set their
    # eccentricities.
    if beam_geom == 'M':
        z_offset = beam_data['M'][dimension_identifier]["centroid"]["z0"]
        attr.setValue("ez0", -z_offset, 0)
    attr.setAnalysisCategory("3D")
        
    return geometry_identifier
    
      
def y_beam(lusas_sesh, **measurements):
    """
    Create a geometric line for a Y-beam.

    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 511

    :param lusas_sesh: LUSAS session object.
    :param measurements: Dictionary with:
        - height
        - topWidth
        - baseWidth
    :return: Geometry identifier string.
    """
    return _yyem_base(lusas_sesh, 'Y', **measurements)
    

def ye_beam(lusas_sesh, **measurements):
    """
    Create a geometric line for a YE-beam.

    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 552

    :param lusas_sesh: LUSAS session object.
    :param measurements: Dictionary with:
        - topWidth
        - baseWidth
        - height
    :return: Geometry identifier string.
    """
    return _yyem_base(lusas_sesh, "YE", **measurements)


def m_beam(lusas_sesh, **measurements):
    """
    Create a geometric line for a YE-beam.

    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 553

    :param lusas_sesh: LUSAS session object.
    :param measurements: Dictionary with:
        - topWidth
        - baseWidth
        - height
    :return: Geometry identifier string.
    """
    return _yyem_base(lusas_sesh, 'M', **measurements)


def _uume_base(lusas_sesh, beam_geom, **measurements):
    """
    Internal base function for generating U and UM beam geometry.

    :param lusas_sesh: LUSAS session object.
    :param beam_geom: Beam type string, e.g., "U" or "UM".
    :param measurements: Dictionary with:
        - openingWidth
        - topWidth
        - baseWidth
        - height
    :return: Geometry identifier string.
    """
    beam_dims = []
    # Order of the key in the beam data json, vital for keying later
    for ie_name in ("openingWidth", "topWidth", "baseWidth", "height"):
        beam_dims.append(measurements[ie_name])
   
    beam_path = join(os.getcwd(), "data", "nonparametric_beam_geometries.json")
    if not exists(beam_path):
        raise FileNotFoundError(
            "Cannot locate the file 'nonparametric_beam_geometries.json'"
        )
    
    with open(beam_path, 'r') as f:
        beam_data = json.load(f)
    dimension_identifier = f"({', '.join(map(str, map(float, beam_dims)))})"
    number = beam_data[beam_geom].get(dimension_identifier, None)['N']
    
    if number:
        beam_profile = f"{beam_geom}{number}"
    else:
        raise KeyError("This beam section is not supported.")
    
    # Create the length of the beam
    geometry_identifier = f"LGeo {beam_profile} Beam {str(measurements)}"
    attr = lusas_sesh.database.createGeometricLine(geometry_identifier)
    attr.setValue("elementType", "3D Thick Beam")
    attr.setFromLibrary(
        "UK Sections", f"Precast {beam_geom} Beams", beam_profile, 0, 0, 0
    )
    attr.setEccentricityOrigin("Centroid", "Fibre", "", "A3")
    attr.setAnalysisCategory("3D")

    return geometry_identifier


def u_beam(lusas_sesh, **measurements):
    """
    Create a geometric line for a U-beam.

    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 593

    :param lusas_sesh: LUSAS session object.
    :param measurements: Dictionary with:
        - height
        - baseWidth
        - topWidth
        - openingWidth
    :return: Geometry identifier string.
    """
    return _uume_base(lusas_sesh, 'U', **measurements)


def um_beam(lusas_sesh, **measurements):
    """
    Create a geometric line for a UM-beam.

    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 594

    :param lusas_sesh: LUSAS session object.
    :param measurements: Dictionary with:
        - height
        - baseWidth
        - topWidth
        - openingWidth
    :return: Geometry identifier string.
    """
    return _uume_base(lusas_sesh, 'UM', **measurements)


def rectangular(lusas_sesh, **measurements):
    """
    Create a geometric line for a rectangular solid beam.

    :param lusas_sesh: LUSAS session object.
    :param measurements: Dictionary with:
        - width
        - height
    :return: Geometry identifier string.
    """
    beam_dims = {}

    ie_names = ("width", "height")
    lusas_names = ('B', 'D')
    for ie_name, lusas_name in zip(ie_names, lusas_names):
        beam_dims[lusas_name] = measurements[ie_name]

    # Create the beam geometry attribute
    ibeam_identifier = f"Rectangular Solid Beam {str(measurements)}"
    beam_face = lusas_sesh.database.createParametricSection(ibeam_identifier)
    beam_face.setType("Rectangular Solid")
    beam_face.setDimensions(tuple(beam_dims.keys()), tuple(beam_dims.values()))
    
    # Create the length of the beam
    geometry_identifier = f"LGeo Rectangular Solid Beam {ibeam_identifier}"
    attr = lusas_sesh.database.createGeometricLine(geometry_identifier)
    attr.setValue("elementType", "3D Thick Beam")
    attr.setEccentricityOrigin("Centroid", "Fibre", "", "S3")
    attr.setFromLibrary("Utilities", "", ibeam_identifier, 0, 0, 0)
    attr.setAnalysisCategory("3D")

    return geometry_identifier


def circular(lusas_sesh, **measurements):
    """
    Create a geometric line for a circular solid beam.

    https://github.com/dynamics-research-group/pbshm-schema/blob/main/
        model-irreducible-element-geometry-data.json
    line 920

    :param lusas_sesh: LUSAS session object.
    :param measurements: Dictionary with:
        - radius
    :return: Geometry identifier string.
    """
    beam_dims = {}

    ie_names = ("radius",)
    lusas_names = ('D',)
    for ie_name, lusas_name in zip(ie_names, lusas_names):
        beam_dims[lusas_name] = measurements[ie_name]

    # Create the beam geometry attribute
    ibeam_identifier = f"Circular Solid Beam {str(measurements)}"
    beam_face = lusas_sesh.database.createParametricSection(ibeam_identifier)
    beam_face.setType("Circular Solid")
    beam_dims['D'] *= 2  # LUSAS wants the diameter
    beam_face.setDimensions(tuple(beam_dims.keys()), tuple(beam_dims.values()))
    
    # Create the length of the beam
    geometry_identifier = f"LGeo Circular Solid Beam {str(measurements)}"
    attr = lusas_sesh.database.createGeometricLine(geometry_identifier)
    attr.setValue("elementType", "3D Thick Beam")
    attr.setFromLibrary("Utilities", "", ibeam_identifier, 0, 0, 0)
    attr.setAnalysisCategory("3D")

    return geometry_identifier
