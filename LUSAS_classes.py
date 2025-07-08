# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 09:42:57 2025

@author: Tristan Gowdridge

This module defines structure-type-specific logic for constructing finite
element (FE) models in LUSAS. It provides a framework for interpreting
structured filenames and instantiating appropriate subclasses representing
different civil infrastructure types, such as bridges, masts, and wind
turbines.

Each structure type can have multiple subclasses, each with its own geometric
and topological rules. The correct subclass is selected based on a structured
identifier encoded in the filename.

The core FE logic is handled via a `LUSASSession` instance, which manages:
    - Creating and configuring a new LUSAS model.
    - Structure wireframe generation
    - Parsing a structured identifier from the filename.
    - Loading and converting Irreducible Element JSON models.
    - Assigning materials and profiles to the model.
    - Meshing geometry elements.
    - Applying subclass-specific modelling logic and loadcases.

The `LUSASSession` class acts as a backend interface and is passed to
structure-specific subclasses that implement their own logic for geometry and
model construction.
"""
import re
import os
from os.path import exists, join
import ast
from collections import defaultdict
from abc import ABC, abstractmethod
import json

from win32com.client.dynamic import Dispatch

from ie_to_fe_funcs import order_coords, greedy_cyclic_sort
from preprocessing_iemodel import extract_profile_info, extract_ie_information
import beam_functions as bf
import plate_functions as pf
import translateAndScale_functions as tasf


def extract_structure_identifiers(filename):
    """
    Extract structure identifiers from a filename. The structure identifiers
    are later used for the LUSAS model name, and all the respective output
    data.
    
    Filenames must contain a structured identifier in the form `a-b-c-d-e-f`,
    where each component conveys specific metadata about the structure:

        a: Structure type
           (e.g., 1 = bridge, 2 = mast, 3 = wind turbine)
        b: Structure subclass
           (type-specific; e.g., 1-1 = beam and slab bridge, 1-2 = ladder deck)
        c: Dataset ID
           (used to group files by test campaign, location, etc.)
        d: State ID
           (represents a structural state, e.g., damaged vs undamaged)
        e: Structure instance ID
           (a unique number for the instance of this subclass)
        f: Variant ID (optional)
           (used to track slight variations in configuration, materials, etc.)

    The first five identifiers (a-e) are mandatory. If a sixth (variant ID) is
    provided, it will be included in the output tuple; otherwise, a default
    value should be handled elsewhere if needed.

    :param filename: The filename containing the structured identifier.
    :type filename: str
    :return: Tuple containing:
        - structure_type (str)
        - structure_subclass (str)
        - dataset_id (str)
        - state_id (str)
        - structure_id (str)
        - variant_id (str)
    :rtype: tuple
    :raises ValueError: If the filename does not contain all required
    identifiers.
    """
    match = re.search(r"(\d+)-(\d+)-(\d+)-(\d+)-(\d+)-?(\d+)?", filename)
    
    if match:
        structure_type = match.group(1)
        structure_subclass = match.group(2)
        dataset_id = match.group(3)
        state_id = match.group(4)
        structure_id = match.group(5)
        variant_id = match.group(6)
    
    if (structure_type and structure_subclass and dataset_id and structure_id
            and state_id):
        return (structure_type, structure_subclass, dataset_id, state_id,
                structure_id, variant_id)
    else:
        raise ValueError(
            "Filename must contain a structure identifier of the form a-b-c-d-e-f."
        )
    
    
def allocate_subclass(filename):
    """
    Allocate a subclass object based on the structure identifiers in the filename.

    :param filename: Filename containing the structured identifier.
    :type filename: str
    :return: Instance of a subclass representing the specific structure type.
    :rtype: Structure subclass instance
    :raises ValueError: If the structure type or subclass is invalid.
    """
    structure_type, structure_subclass, _, _, _, _ = extract_structure_identifiers(filename)

    # Define a dictionary mapping structure types to their corresponding
    # subclasses
    subclasses = {
        "1": {  # Bridges
            "1": BeamAndSlab,
            "2": LadderDeck,
        },
        "2": {  # Masts
            "1": LatticeTower,
            "2": Monopole
        },
        "3": {  # Wind Turbines
            "1": WindTurbine
        }
    }

    # Check if the struct_type and struct_subclass exist in the dictionary
    try:
        subclass_class = subclasses[structure_type][structure_subclass]
    except KeyError:
        raise ValueError(
            f"Invalid structure type '{structure_type}' or subclass "
            f"'{structure_subclass}' for filename '{filename}'"
        )
    # Instantiate the appropriate subclass with the filename
    return subclass_class(filename)


def is_clockwise_about(p1, p2, centre):
    """
    Determine if the vector from `p1` to `p2` rotates clockwise around a
    `centre`.

    This is done using the 2D cross product in the XY plane.

    :param p1: First point (x1, y1, z1)
    :type p1: tuple
    :param p2: Second point (x2, y2, z2)
    :type p2: tuple
    :param centre: Centre point for rotation (cx, cy, cz)
    :type centre: tuple
    :return: True if the vector from p1 to p2 is clockwise about the centre.
    :rtype: bool
    """
    x1, y1, _ = p1
    x2, y2, _ = p2
    cx, cy, _ = centre

    # Compute 2D cross product to determine rotational direction
    cross = (x2 - x1) * (y1 - cy) - (y2 - y1) * (x1 - cx)
    
    return cross < 0


class LUSASSession(ABC):
    """
    Abstract base class for constructing LUSAS FE models from structured data.

    This class encapsulates the shared logic needed to:

    - Parse structure identifiers from the filename.
    - Load and interpret instance-exchange (IE) JSON models.
    - Create a hierarchical wireframe model in LUSAS (volumes > surfaces >
        lines > points).
    - Assign materials, profiles, and ground elements.
    - Apply surface vector orientation and geometric corrections.
    - Create line, surface, and volume meshes.
    - Load VBA scripts and assign loadcases.
    - Save model outputs for downstream analysis.

    Structure-specific subclasses must implement the abstract method
    :meth:`subclass_specific_logic` to define any behaviour unique to the
    given structure type or subclass (e.g., beam-and-slab bridge, monopile).

    Geometry and other modelling parameters are read from
    ``fe_run_params.json`` if it exists; otherwise, a set of default parameters
    is used.
    """
    
    VERSION = "21.1"
    FE_FOLDER = "fe_models"  # Directory where the FE models will be saved.
    IE_FOLDER = "ie_models"  # Directory from where the IE models are loaded.
    LOADCASE_FOLDER = "loadcase_vba_scripts"
    
    def __init__(self, filename):
        """
        Initialise a LUSAS session with the given IE path and project name.
        This class is used as a base class for all the structure subclasses,
        and contains all of the logic for creating FE models. The children of
        this class apply the specific logic to creating that subclass.
        
        The logic within this class leverages the fact that points can be
        inherited from lines, which can be inherited from surfaces, which can
        be inherited from volumes.
        """
        
        fe_params_path = join(os.getcwd(), "instance", "fe_run_params.json")
        if exists(fe_params_path):
            print("Loading from JSON")
            with open(fe_params_path, "r") as f:
                self.fe_params = json.load(f)
        else:
            self.fe_params = {
                "verticalDirection": "Z",
                "analysisCategory": "3D",
                "modelUnits": "kN,m,t,s,C",
                "useLineDistance": True,
                "lineMeshSpacing": [2],
                "lineMeshDistance": [0.5],
                "surfaceMeshSpacing": [0, 0],
                "volumeMeshSpacing": [2, 2, 2],
                "performDynamicAnalysis": True,
                "performStaticAnalysis": True,
                "nModalFeatures": 10,
                "saveData": {
                    "displacement": True,
                    "reaction": True,
                    "forceMomentBeam": True,
                    "forceMomentShell": True,
                    "loading": True,
                    "modal": True
                }
            }
        
        self.structure_identifiers = extract_structure_identifiers(filename)
        self.bridge_identifier = "-".join(self.structure_identifiers)
        
        # Create a new LUSAS modeller instance, this boots up a LUSAS session.
        self.modeller = Dispatch("Lusas.Modeller." + LUSASSession.VERSION)
        
        # Folder to save the FE models
        if not exists(join(os.getcwd(), LUSASSession.FE_FOLDER)):
            os.mkdir(join(os.getcwd(), LUSASSession.FE_FOLDER))

        # Create a new project inside the LUSAS session, defining the analysis
        # type and where to save the FE mode.
        self.database = self.modeller.newDatabase(
            "structural",
            join(os.getcwd(), LUSASSession.FE_FOLDER, self.bridge_identifier)
        )
        
        # Defining some traits of the model behaviour.
        self.database.setModelTitle(self.bridge_identifier)
        self.database.setVerticalDir(self.fe_params["verticalDirection"])
        self.database.setAnalysisCategory(self.fe_params["analysisCategory"])
        self.database.setModelUnits(self.fe_params["modelUnits"])
        
        # The following are used repeatedly to store different geometries and
        # assignments.
        self.selection = self.modeller.getSelection()
        self.assignment = self.modeller.getSessionFileAssignment()
        self.assignment.setAllDefaults()
        self.geom = self.modeller.getSessionFileGeometryData()
        self.geom.setAllDefaults()
        
        # Used to store the coordinates of the ground elements.
        self.ground_elements = []
        
        # Used as a handy Python-side lookup table to reduce the number of time
        # consuming LUSAS API calls.
        self.coords_name_map = {}
        
        # Used as a flag to run the the simulation if a relevant loading script
        # is found.
        self.loading_script_found = False
        
        ie_path = join(os.getcwd(), LUSASSession.IE_FOLDER, filename)
        if not exists(ie_path):
            raise FileNotFoundError(f"The file {ie_path} does not exist.")
            
        # Extracts all the information from a JSON IE model, and converts it
        # into a more handy table structure, which is used internally in
        # Python.
        self.ie_table = extract_ie_information(ie_path)
        
        # Orders by the geometry types to make use of geometric inheritance.
        # For instance, points are inherited from lines, lines from surfaces,
        # etc...
        self.geometry_collection = defaultdict(dict)
        self.all_materials, self.all_profiles = set(), set()
        for info in self.ie_table.values():
            geom = info.geometry
            coords = order_coords(info.coords)
            self.geometry_collection[geom][coords] = info
            self.all_materials.add(info.material[0])
            
            profile_info = extract_profile_info(info)
            self.all_profiles.add(profile_info)
            
            if info.ground_coord:
                self.ground_elements.extend(info.ground_coord)
                
        self.create_wireframe()
        
        # From the wireframe, we now add specific information.
        self.assign_materials()
        self.assign_profiles()
        self.assign_ground_elements()
        
        # Bridge model fine tuning
        self.apply_surface_eccentricity()
        self.correct_perimeter_line_direction_for_geometry()
        
        # Apply the model specific logic
        self.subclass_specific_logic()
        
        # Make basic meshes.
        if self.fe_params["useLineDistance"]:
            self.create_equal_line_mesh(*self.fe_params["lineMeshDistance"])
        else:
            self.create_basic_line_mesh(*self.fe_params["lineMeshSpacing"])
        
        self.create_basic_surface_mesh(*self.fe_params["surfaceMeshSpacing"])
        self.create_basic_volume_mesh(*self.fe_params["volumeMeshSpacing"])
        
        # Assign the loadcases.
        self.assign_loadcase_from_lvb()
        
        # Save the results to txt.
        self.save_results()
        
        # These need modifying.................................................
        # parsed_data = parse_loadcase_data(project_name)
        # add_to_database(parsed_data)
        
        self.save_model()
        # self.modeller.quit()
                
    @abstractmethod
    def subclass_specific_logic(self):
        """
        Abstract method requried to be overwritten to provide structure
        specific behaviour for FE model construction.
        """
        pass
        
    def create_wireframe(self):
        """
        Create the wireframe model based on the geometry collection.

        This method processes the geometry elements starting from the highest
        dimension (Volume) and working down to the lowest dimension (Point).
        Lower-dimensional elements are likely subsets of higher-dimensional
        ones, so this processing order is eliminates inefficiency.
        """
        for coords, element in self.geometry_collection["Volume"].items():
            raise NotImplementedError
        
        for coords, element in self.geometry_collection["Surface"].items():
            self.create_surface(element)
        
        for coords, element in self.geometry_collection["Line"].items():
            self.create_line(element)
              
        for coords, element in self.geometry_collection["Point"].items():
            raise NotImplementedError
    
    def create_surface(self, element):
        """
        This method processes the given surface element by sorting its
        coordinates, adding them to the model, and creating the surface in the
        LUSAS database. It also applies surface vector orientation, assigns a
        unique LUSAS identifier, and processes associated lines.
        """
        surface_key = order_coords(element.coords)
        if surface_key in self.coords_name_map:
            print("Surface already exists could optimise for lookup")
        
        sorted_points = greedy_cyclic_sort(list(element.coords))
        self.geom.setAllDefaults()
        self.geom.setCreateMethod("coons")
        for sorted_point in sorted_points:
            self.geom.addCoords(*sorted_point)
        self.geom.setLowerOrderGeometryType("coordinates")
        pointer = self.database.createSurface(self.geom)
        
        for obj in pointer.getObjects("All"):
            name = obj.getName()  # LUSAS index
            if obj.getTypeCode() == 4:  # Surface type code
                self.orient_surface_vectors(obj, name)
                element.lusas_identifier = (element.geometry, name)
                self.coords_name_map[surface_key] = ("Surface", name)
                self.get_lines_of_surface(obj)
    
    def orient_surface_vectors(self, obj, name):
        """
        Reverse the orientation of a surface if its normal vector points
        downwards.
        """
        normal_vector = obj.getNormal()
        if normal_vector[2] < 0:
            self.selection.remove("All")
            self.selection.add("Surface", name)
            self.geom.setAllDefaults()
            self.geom.cycleReverse(True)
            self.geom.setReverseKeepX(True)
            self.geom.cycleReset(False)
            tempSetA = self.modeller.newObjectSet()
            tempSetA.add(self.selection, "Surface")
            tempSetA.reverse(self.geom)
        self.selection.remove("All")
        
    def get_lines_of_surface(self, pointer):
        """
        Get lines of a surface from the pointer. This is used to store the
        lines created from the faces of higher-dimensional geometries. Also
        stores the points at the end of the lines.
        """
        for line in pointer.getLOFs():
            start_point = line.getStartPoint()
            end_point = line.getEndPoint()
            p1, p2 = start_point.getXYZ(), end_point.getXYZ()
            line_key = order_coords((p1, p2))
            if line_key not in self.coords_name_map:
                self.coords_name_map[line_key] = ("Line", line.getName())
                self.coords_name_map[p1] = ("Point", start_point.getName())
                self.coords_name_map[p2] = ("Point", end_point.getName())

    def create_line(self, element):
        """
        Create a line element and store its endpoints as points in the model.
        """
        line_key = order_coords(element.coords)
        if line_key in self.coords_name_map:
            # If the line already exists, ie it is the face of a surface, then
            # skip creating.
            element.lusas_identifier = self.coords_name_map[line_key]
            return
        
        self.geom.setAllDefaults()
        self.geom.setCreateMethod("straight")
        
        for point in element.coords:
            self.geom.addCoords(*point)
        self.geom.setLowerOrderGeometryType("coordinates")
        
        objects = self.database.createLine(self.geom)

        for obj in objects.getObjects("All"):
            if obj.getTypeCode() == 2:  # Type code for a line
                line = obj
                break
        else:
            raise TypeError("No lines in object set?")
        
        self.coords_name_map[line_key] = ("Line", line.getName())
        element.lusas_identifier = (element.geometry, line.getName())
        
        p1 = line.getStartPoint()
        p2 = line.getEndPoint()
        self.coords_name_map[p1.getXYZ()] = ("Point", p1.getName())
        self.coords_name_map[p2.getXYZ()] = ("Point", p2.getName())
        
    def create_basic_volume_mesh(self, sx, sy, sz):
        """
        Create and assign a basic volume mesh, with a specified spacing in each
        axis.
        """
        self.selection.remove("All")
        # Create and assign a mesh
        self.selection.add("Volume", "All")
    
        attr = self.database.createMeshVolume("VMsh1")
        attr.setRegular("HX8M", sx, sy, sz, False)
        attr.setValue("allowIrregular", True)
        attr.setAnalysisCategory("3D")
    
        self.assignment.setSelectionNone()
        self.assignment.addToSelection("Volume")
        self.assignment.setLoadset("Loadcase 1")
        self.database.getAttribute("Volume Mesh", "VMsh1").assignTo(
            self.selection,
            self.assignment
        )
        self.database.updateMesh()
        self.selection.remove("All")
    
    def create_basic_surface_mesh(self, sx, sy):
        """
        Create and assign a basic surface mesh, with a specified spacing in
        each axis.
        """
        self.selection.remove("All")
        # Create and assign a mesh
        self.selection.add("Surface", "All")
    
        attr = self.database.createMeshSurface("SMsh1")
        attr.setRegular("QTS4", sx, sy, True)
        attr.setValue("allowIrregular", True)
        attr.setAnalysisCategory("3D")
    
        self.assignment.setAllDefaults()
        self.assignment.setLoadset("Loadcase 1")
        self.database.getAttribute("Surface Mesh", "SMsh1").assignTo(
            self.selection, self.assignment
        )
        self.database.updateMesh()
        self.selection.remove("All")
    
    def create_basic_line_mesh(self, n_divisions):
        """
        Create and assign a basic line mesh with divisions specified by
        n_divisions.
        """
        self.selection.remove("All")
        attr = self.database.createMeshLine("LMsh1")
        attr.setSpacing("BMI21")
        attr.addSpacing(n_divisions, 1.0)
        attr.setEndReleasesSameAsStart(True)
        attr.setValue("uiSpacing", "uniform")
        attr.setAnalysisCategory("3D")
        
        self.selection.add("Line", "All")
        
        self.assignment.setSelectionNone()
        self.assignment.addToSelection("Line")
        self.assignment.setLoadset("Loadcase 1")
        self.database.getAttribute("Line Mesh", "LMsh1").assignTo(
            self.selection,
            self.assignment
        )
        self.database.updateMesh()
        self.selection.remove("All")
    
    def create_equal_line_mesh(self, mesh_length):
        """
        Creates a line mesh with nodes spaced at mesh_length.
        """
        self.selection.remove("All")
        attr = self.database.createMeshLine("LMsh1")
        attr.setSize("BMI21", mesh_length)
        attr.setEndReleasesSameAsStart(True)
        attr.setValue("uiSpacing", "irregular")
        attr.setAnalysisCategory("3D")
        
        self.selection.add("Line", "All")
        
        self.assignment.setSelectionNone()
        self.assignment.addToSelection("Line")
        self.assignment.setLoadset("Loadcase 1")
        self.database.getAttribute("Line Mesh", "LMsh1").assignTo(
            self.selection,
            self.assignment
        )
        self.database.updateMesh()
        self.selection.remove("All")
        
    def assign_line_geometries(self, lines, geometry_identifier):
        """
        Add a cross-section to the specified lines.
        """
        self.selection.remove("All")
        self.assignment.setAllDefaults()
        self.assignment.setSelectionNone()
        self.selection.add("Line", lines)
        self.assignment.addToSelection("Line")
        self.assignment.setLoadset("Loadcase 1")
        self.assignment.setAnalysis("Analysis 1")
        self.database.getAttribute("Line Geometric", geometry_identifier).assignTo(
            self.selection, self.assignment
        )
        self.selection.remove("All")
    
    def assign_surface_geometries(self, surfaces, geometry_identifier):
        """
        Add a cross-section to the specified surfaces.
        """
        self.selection.remove("All")
        self.assignment.setAllDefaults()
        self.assignment.setSelectionNone()
        self.selection.add("Surface", surfaces)
        self.assignment.addToSelection("Surface")
        self.assignment.setLoadset("Loadcase 1")
        self.assignment.setAnalysis("Analysis 1")
        self.database.getAttribute("Surface Geometric", geometry_identifier).assignTo(
            self.selection, self.assignment
        )
        self.selection.remove("All")
    
    def create_material(self, material_name, material_properties):
        """
        Create a material with the given properties.
        """
        attr = self.database.createIsotropicMaterial(
            material_name,
            material_properties["youngsModulus"],
            material_properties["poissonsRatio"],
            material_properties["density"]
        ).setValue(
            "alpha",  # Parameter being set
            material_properties["linearThermalExpansionCoefficient"],  # value
            0
        )
        attr.setDefinitionMenuID(1, None, True)
        attr.setDescription(f"Description for {material_name.title()}")
    
    def assign_materials(self):
        """
        Assign materials to their corresponding FE elements.
        """
        for material in self.all_materials:
            material_properties = None
            elements_to_select = defaultdict(list)
            for element in self.ie_table.values():
                if element.material[0] == material:
                    if not material_properties:
                        # This does assume that every material by the same name
                        # has the same property values.
                        material_properties = element.material[1]
                    elements_to_select[element.lusas_identifier[0]].append(element.lusas_identifier[1])
            
            self.selection.remove("All")
            # Select all the elements that have this material.
            for geometry, ids in elements_to_select.items():
                self.selection.add(geometry, ids)
            self.create_material(material, material_properties)
            self.assignment.setAllDefaults()
            self.assignment.setLoadset("Loadcase 1")
            self.database.getAttribute("Isotropic Material", material).assignTo(self.selection, self.assignment)
            self.selection.remove("All")

    def assign_profiles(self):
        """
        Assign the geometric profiles to all the corresponding FE elements.
        """
        for profile in self.all_profiles:
            profile_properties = None
            add_profile = []
            for element in self.ie_table.values():
                if extract_profile_info(element) == profile:
                    if not profile_properties:
                        profile_properties = profile
                    add_profile.append(element.lusas_identifier[1])
              
            if profile_properties[0][1] == "beam":
                beam_type = profile_properties[0][0]
                measurements = dict(profile_properties[1:])
                beam_type = beam_type.replace('-', '_')
                beam_function = getattr(bf, beam_type, None)
                
                if beam_function:
                    beam_identifier = beam_function(self, **measurements)
                else:
                    raise IndexError("Beam function does not exist")
                
                self.assign_line_geometries(add_profile, beam_identifier)
                self.selection.remove("All")
            
            elif profile_properties[0][1] == "plate":
                plate_type = profile_properties[0][0]
                measurements = dict(profile_properties[1:])
                plate_type = plate_type.replace('-', '_')
                plate_function = getattr(pf, plate_type, None)
                
                if plate_function:
                    plate_identifier = plate_function(self, **measurements)
                else:
                    raise IndexError("Plate function does not exist")
                
                self.assign_surface_geometries(add_profile, plate_identifier)
                self.selection.remove("All")
            
            elif profile_properties[0][1] == "translateAndScale":
                tas_type = profile_properties[0][0]
                measurements = dict(profile_properties[1:])
                tas_type = tas_type.replace('-', '_')
                tas_function = getattr(tasf, tas_type, None)
                
                if tas_function:
                    tas_identifier = tas_function(self, **measurements)
                else:
                    raise IndexError("translateAndScale function does not exist")
                
                self.assign_line_geometries(add_profile, tas_identifier)
                self.selection.remove("All")

    def assign_ground_elements(self):
        """
        Assigns a ground support to all of the listed ground elements.
        """
        # Create the ground attribute
        attr = self.database.createSupportStructural("Ground")
        attr.setStructural("R", "R", "R", "F", "F", "F", "F", "F", "C", "F")
        attr.setAnalysisCategory("3D")
        ground_points = []
        for ground_coord in self.ground_elements:
            ground_points.append(self.coords_name_map[ground_coord][1])
        self.selection.add("Point", ground_points)
        self.assignment.setAllDefaults()
        self.assignment.setLoadset("Loadcase 1")
        self.database.getAttribute(
            "Structural Support", "Ground"
        ).assignTo(self.selection, self.assignment)

    def save_results(self):
        """
        Saves model results based on the analysis settings.
        """
        filename = self.bridge_identifier
        
        if self.fe_params["saveData"]["modal"] and self.fe_params["performDynamicAnalysis"]:
            self.setup_modal_analysis(filename)
            self.save_modal_data(filename)
        
        if not self.loading_script_found:
            print(f"Cannot save results for {filename} as no loading script was found.")
            return
        elif not self.fe_params["performStaticAnalysis"]:
            return
        
        if self.fe_params["saveData"]["displacement"]:
            self.save_displacement_data(filename)
        if self.fe_params["saveData"]["reaction"]:
            self.save_reaction_data(filename)
        if self.fe_params["saveData"]["forceMomentBeam"]:
            self.save_force_moment_beam_data(filename)
        if self.fe_params["saveData"]["forceMomentShell"]:
            self.save_force_moment_shell_data(filename)
        if self.fe_params["saveData"]["loading"]:
            self.save_loading_data(filename)

    def common_save_params(self, attr, set_results_type, set_analysis_results_type, loadcase_option):
        """
        Sets common save parameters to avoid repetition in result configuration.
        """
        attr.setUnits(None)
        attr.setResultsType(set_results_type)
        attr.setResultsOrder("Mesh")
        attr.setResultsContent("Tabular and Summary")
        attr.setExtent("Elements showing results", "")
        attr.setResultsLocation("Nodal")
        attr.setLoadcasesOption(loadcase_option)
        primaryComponents = ["All"]
        primaryEntities = ["Displacement"]
        attr.setPrimaryResultsData(primaryComponents, primaryEntities)
        attr.setAnalysisResultTypes(set_analysis_results_type)
        attr.setResultsTransformNone()
        attr.showCoordinates(True)
        attr.showExtremeResults(False)
        attr.setSlice(False)
        attr.setAllowDerived(False)
        attr.setDisplayNow(True)
        attr.showActiveNodesOnly(True)
        attr.setSigFig(6, False)
        attr.setThreshold(None)
        attr.showResults(False)
    
    def create_save_folder_and_save(self, filename, folder_name):
        """
        Creates a folder and saves the model data to a text file.
        """
        folder_path = join(
            os.getcwd(), LUSASSession.FE_FOLDER, folder_name
        )
        if not exists(folder_path):
            os.mkdir(folder_path)
        
        self.modeller.getGridWindowByID(1).saveAllAs(
            join(folder_path, f"{filename}.txt"), "Text"
        )
        self.modeller.getGridWindowByID(1).close()
        
    def save_displacement_data(self, filename):
        """
        Save nodal displacement results to file.
        """
        attr = self.database.createPrintResultsWizard("PRW Displacement")
        attr.setResultsEntity("Displacement")
        attr.setComponents(["DX", "DY", "DZ", "THX", "THY", "THZ", "RSLT"])
        self.common_save_params(attr, "Components", None, "All")
        self.create_save_folder_and_save(filename, "displacement_results")
    
    def save_reaction_data(self, filename):
        """
        Save reaction force results to file.
        """
        attr = self.database.createPrintResultsWizard("PRW Reaction")
        attr.setResultsEntity("Reaction")
        attr.setComponents(["FX", "FY", "FZ", "MX", "MY", "MZ", "RSLT"])
        self.common_save_params(attr, "Components", None, "All")
        self.create_save_folder_and_save(filename, "reaction_results")

    def save_force_moment_beam_data(self, filename):
        """
        Save beam force moment results to file.
        """
        attr = self.database.createPrintResultsWizard("Force-Moment_Beam")
        attr.setResultsEntity("Force/Moment - Thick 3D Beam")
        attr.setComponents(["Fx", "Fy", "Fz", "Mx", "My", "Mz"])
        self.common_save_params(attr, "Components", None, "All")
        self.create_save_folder_and_save(filename, "force_moment_beam_results")
    
    def save_force_moment_shell_data(self, filename):
        """
        Save shell force/moment results to file.
        """
        attr = self.database.createPrintResultsWizard("Force-Moment_Shell")
        attr.setResultsEntity("Force/Moment - Thick Shell")
        attr.setComponents(["Nx", "Ny", "Nxy", "Mx", "My", "Mxy", "Sx", "Sy"])
        self.common_save_params(attr, "Components", None, "All")
        self.create_save_folder_and_save(filename, "force_moment_shell_results")
    
    def save_loading_data(self, filename):
        """
        Save applied loading results to file.
        """
        attr = self.database.createPrintResultsWizard("Loading")
        attr.setResultsEntity("Loading")
        attr.setComponents(["FX", "FY", "FZ", "MX", "MY", "MZ", "RSLT"])
        self.common_save_params(attr, "Components", None, "All")
        self.create_save_folder_and_save(filename, "loading_results")
    
    def setup_modal_analysis(self, filename):
        """
        Configure and solve the modal analysis.
        """
        analysis = self.database.createAnalysisStructural("Modal_analysis", True, "Nonlinear and transient")
        analysis.setInheritFromBase("None")
        analysis.setInheritFromBase("All")
        analysis.setUndeformedMeshStart()
        analysis.setCoupled(False)
        analysis.setSelectedResultsGroup("all")
        analysis.setSelectedElementOutputGroup("all")
        analysis.setSelectedNodeOutputGroup("all")
        
        # Modify loadcase/control
        self.database.options().setBoolean("Option 311", False, False, "Modal_analysis")
        self.database.options().setBoolean("Option 432", False, False, "Modal_analysis")
        loadcase = self.database.getLoadset("Loadcase 1", 0)
        loadcase.setEigenvalueMaxMinControl("Frequency", "Minimum", 10)
        loadcase.getEigenvalueControl().setValue("nivc", 0).setValue("nitem", 30).setValue("namast", 0).setValue("shift", 0.0)
        loadcase.getEigenvalueControl().setValue("rtol", 0.1E-3).setValue("sturm", True).setValue("Guyan", False).setValue("buckl", False)
        loadcase.getEigenvalueControl().setValue("NormalisationProcedure", "GlobalMass").setValue("Eigensolver", "Default")
        
        # Solve
        self.database.closeAllResults()
        self.database.updateMesh()
        self.database.save()
        self.database.getAnalysis("Modal_analysis").solve(True)
        self.database.openAllResults(True)
        
    def save_modal_data(self, filename):
        """
        Save modal frequencies and mode shapes.
        """
        # Handle Natural Frequencies
        attr = self.database.createPrintResultsWizard("Frequencies")
        attr.setResultsEntity("Loading")
        attr.setComponents(["FX", "FY", "FZ", "MX", "MY", "MZ", "RSLT"])
        self.common_save_params(attr, "Eigenvalues", ["Eigenvalues"], "All")
        self.create_save_folder_and_save(filename, "frequencies")
        
        # Handle Mode Shapes
        attr = self.database.createPrintResultsWizard("Mode_shapes")
        
        for loadcase in self.database.getLoadsets("All"):
            if "Mode" in loadcase.getName():
                loadcase_id = loadcase.getID()
                break
        n_modes = self.fe_params["nModalFeatures"]
        lcIDs = [loadcase_id] * n_modes
        lcResFileIDs = [2] * n_modes  # self.database.getAnalyses()[1].getAnalysisResultsFileIDs() to get 2
        lcEigenvalueIDs = list(range(1, n_modes+1))
        lcHarmonicIDs = [-1] * n_modes  # Not entirely sure what the -1 is

        attr.setResultsEntity("Displacement")
        attr.setLoadcases(lcIDs, lcResFileIDs, lcEigenvalueIDs, lcHarmonicIDs)
        attr.setComponents(["DX", "DY", "DZ", "THX", "THY", "THZ", "RSLT"])
        self.common_save_params(attr, "Components", ["Eigenvalues"], "Selected")
        self.create_save_folder_and_save(filename, "modalshapes")

    def assign_gravity(self):
        """
        Assigns gravity to the model. Originally used for debugging, but
        retained for convenience in testing or standalone runs where load
        scripts are not used.
        """
        self.loadcase = self.database.getLoadset("Loadcase 1", 0)
        self.loadcase.addGravity(True)
        self.loadcase.setGravityFactor(1)
    
    def run_analysis(self, version=1):
        """
        Runs the analysis for the model. Originally used for debugging, but
        retained as a convenient way to manually trigger an analysis outside of
        automated scripts.
        """
        self.database.closeAllResults()
        self.database.updateMesh()
        self.database.save()
        self.database.getAnalysis(f"Analysis {version}").solve(True)
        self.database.openAllResults(False)

    def assign_loadcase_from_lvb(self):
        """
        Searches for the appropriate loading script and applies it to the model
        if a matching identifier is found.
        """
        if not self.fe_params["performStaticAnalysis"]:
            return
        
        loadcase_path = join(os.getcwd(), LUSASSession.LOADCASE_FOLDER)
        if not exists(loadcase_path):
            return
        
        for filename in os.listdir(loadcase_path):
            print("This is now different, check with Connor")
            match = re.search(r"(\d+)-(\d+)-(\d+)", filename)
            if not match:
                continue
            if match.group() == self.bridge_identifier:
                self.modeller.fileOpen(join(loadcase_path, filename))
                self.loading_script_found = True
                break
         
    def save_model(self):
        """
        Saves the current model to the FE folder using the bridge identifier.
        """
        self.database.saveAs(
            join(os.getcwd(), LUSASSession.FE_FOLDER, self.bridge_identifier)
        )
        
    def apply_surface_eccentricity(self):
        """
        Applies a vertical eccentricity to deck surfaces so that they align
        flush with adjacent supporting side beams. The shift equals the
        full beam height plus half the deck thickness, moving the deck
        upward relative to its original position.
        """
        # Eccentric shift for surfaces respective to their beam geoms.
        surfaces_fixed = set()
        for obj in self.selection.add("Surface", "All").getObjects("Surface", "All"):
            if obj.getID() in surfaces_fixed:
                continue
            
            surf_geom_ass = obj.getAssignments("Geometric")[0]
            surf_geom = surf_geom_ass.getAttribute()
            
            for in_selection in self.selection.add(surf_geom).getObjects("All"):
                if in_selection.getTypeCode() == 4:
                    surfaces_fixed.add(in_selection.getID())
            
            surface_name = surf_geom.getName()
            thickness = surf_geom.getValue("t")
            
            height = 0
            for lof in obj.getLOFs():
                if lof.getTypeCode() != 2:
                    continue

                attr = lof.getAssignments("Geometric")[0].getAttribute()
                if 'D' in attr.getValueNames():
                    height = max(height, attr.getValue("D"))
                
                else: ## This is pretty bodgy, cause getting dim from the name of the attribute.
                    attr_name = attr.getName()
                    re_match = re.search(r"\{.*\}", attr_name)
                    if re_match:
                        measurements = ast.literal_eval(re_match[0])
                        height = max(height, measurements["height"])

            attr = self.database.createGeometricSurface(surface_name)
            attr.setSurface(thickness, -(height + thickness/2))

    def correct_perimeter_line_direction_for_geometry(self):
        """
        Corrects the direction of perimeter lines for asymmetric beam profiles,
        specifically YE and UM sections, to ensure proper orientation. For
        these profiles, the exposed (flat) face must align consistently --
        clockwise for YE and anti-clockwise for UM -- so that they render
        correctly in the FE model.
        """
        for attr in self.database.getAttributes("Line Geometric"):
            attr_id_and_name = attr.getIDandName()
            attr_id, attr_name = attr_id_and_name.split(':')[:2]
            beam_geometry = attr_name.split()[1]
            
            match = re.match(r"^([A-Z]{1,2})(\d{1,2})$", beam_geometry)
            if match:
                letters = match.group(1)
                if letters in ("YE",):  # Clockwise (negative dx at y=0)
                    reverse_condition = -1
                elif letters in ("UM",):  # Antilockwise (positive dx at y=0)
                    reverse_condition = 1
                else:
                    continue
                
                lines_to_be_reversed = []
                for beam in self.selection.add(attr).getObjects("Line", "All"):
                    facet_coords = beam.getFacetCoordinates()
                    p1, p2 = facet_coords[:3], facet_coords[3:]
                    dx = p2[0] - p1[0]
            
                    if (
                            ((p1[1] == 0) and (dx * reverse_condition < 0)) or
                            ((p1[1] != 0) and (dx * reverse_condition > 0 ))
                        ):
                        lines_to_be_reversed.append(beam.getID())
                    
                # Perform reversal
                self.geom.setAllDefaults()
                self.geom.cycleReverse(True)
                self.geom.cycleReset(False)
                self.selection.remove("All")
                
                self.selection.add("Line", lines_to_be_reversed)
                self.selection.reverse(self.geom)
                self.selection.remove("All")


class BeamAndSlab(LUSASSession):
    def __init__(self, filename):
        super(BeamAndSlab, self).__init__(filename)
    
    def subclass_specific_logic(self):
        """
        There is no additional specific logic for the beam and slab bridge, but
        need to overwrite the abstract method in the inherited class.
        """
        pass


class LadderDeck(LUSASSession):
    def __init__(self, filename):
        super(LadderDeck, self).__init__(filename)
    
    def subclass_specific_logic(self):
        """
        Ladder deck specific implementation, in this case it is required that
        the topside of the cross beams are flush with the surface of the
        bridge.
        """
        obj = self.selection.add("Surface", "All").getObjects("Surface", "All")[0]
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
            main_span_height = self.database.getAttribute("Geometric", main_prof).getValue("D")
            for prof in tuple(ladder_beam_profiles):
                attr = self.database.createGeometricLine(prof)
                
                attr.setValue("elementType", "3D Thick Beam")
                attr.setFromLibrary("Utilities", "", prof[5:], 0, 0, 0)
                attr.setEccentricityOrigin("Centroid", "Fibre", "", "I3")
                ladder_height = attr.getValue("D")
                attr.setValue("ez0", -(main_span_height-ladder_height), 0)
                attr.setAnalysisCategory("3D")


class LatticeTower(LUSASSession):
    def __init__(self, filename):
        super(LatticeTower, self).__init__(filename)

    def subclass_specific_logic(self):
        """
        Need to:
            1) Orient all the lines to point upwards.
            2) For the horizontal bracing, we need to make the line go
                clockwise (in plan view)
            3) Create 4 local coordinate items for each leg rotation.
    
        """
        # Used to determine the centre of the structure, upon which the
        # rotation of the horizontal elements is tested against.
        x_min, y_min = float("inf"), float("inf")
        x_max, y_max = float("-inf"), float("-inf")
        for element in self.ie_table.values():
            for coord in element.coords:
                x_max, y_max = max(x_max, coord[0]), max(y_max, coord[1])
                x_min, y_min = min(x_min, coord[0]), min(y_min, coord[1])
        centre = ((x_max-x_min)/2, (y_max-y_min)/2, 0)
        
        lines_to_be_reversed = []
        for line in self.selection.add("Line", "All").getObjects("Line", "All"):
            start_point = line.getStartPoint()
            end_point = line.getEndPoint()
            p1, p2 = start_point.getXYZ(), end_point.getXYZ()
            if p2[2] - p1[2] < 0:  # 1) Orient Lines Upwards
                lines_to_be_reversed.append(line.getID())
            
            elif abs(p2[2] - p1[2]) < 1e-6:  # 2) Make horizontal clockwise
                if is_clockwise_about(p1, p2, centre):
                    lines_to_be_reversed.append(line.getID())
        
        # Perform reversal
        self.geom.setAllDefaults()
        self.geom.cycleReverse(True)
        self.geom.cycleReset(False)
        self.selection.remove("All")
        
        self.selection.add("Line", lines_to_be_reversed)
        self.selection.reverse(self.geom)
        self.selection.remove("All")
        
        rotation_map = {
            1: 180,
            2: 90,
            3: 0,
            4: 270
        }
        
        legs = defaultdict(list)
        for key in self.ie_table:
            if key.startswith("leg-"):
                parts = key.split("-")
                if len(parts) == 3:
                    legs[parts[2]].append(key)
        
        origin = (0, 0, 0)
        for i in range(1, 5):
            lines_to_rotate = []
            rotat = rotation_map[i]
            for leg_name in legs[str(i)]:
                lines_to_rotate.append(self.ie_table[leg_name].lusas_identifier[1])
            attr = self.database.createLocalCartesianXYAttr(f"Rotation_{rotat}", rotat, origin).setAxesType("Cartesian")
            self.selection.add("Line", lines_to_rotate)
            self.assignment.setAllDefaults()
            self.assignment.setLoadset("Loadcase 1")
            self.database.getAttribute("Local Coordinates", f"Rotation_{rotat}").assignTo(self.selection, self.assignment)
            self.selection.remove("All")
    
        # Rotating each leg so corner faces outwards
        rotation_map = {
            1: 180,
            2: 90,
            3: 0,
            4: 270
        }
        
        legs = defaultdict(list)
        for key in self.ie_table:
            if key.startswith("leg-"):
                parts = key.split("-")
                if len(parts) == 3:
                    legs[parts[2]].append(key)
        
        origin = (0, 0, 0)
        for i in range(1, 5):
            lines_to_rotate = []
            rotat = rotation_map[i]
            for leg_name in legs[str(i)]:
                lines_to_rotate.append(self.ie_table[leg_name].lusas_identifier[1])
            attr = self.database.createLocalCartesianXYAttr(f"Rotation_{rotat}", rotat, origin).setAxesType("Cartesian")
            self.selection.add("Line", lines_to_rotate)
            self.assignment.setAllDefaults()
            self.assignment.setLoadset("Loadcase 1")
            self.database.getAttribute("Local Coordinates", f"Rotation_{rotat}").assignTo(self.selection, self.assignment)
            self.selection.remove("All")
        
        # Setting the cross-section centroid as L2 (the corner of the L)
        for line_geom in self.database.getAttributes("Line Geometric"):
            line_geom.setEccentricityOrigin("Fibre", "Fibre", "L2", "L2")

    def assign_loadcase_from_lvb(self):
        """
        All Mast loadcases are the same, so only creating one file, and apply
        this to all.
        """
        if not self.fe_params["performStaticAnalysis"]:
            return
        
        loadcase_path = join(
            os.getcwd(), LUSASSession.LOADCASE_FOLDER, "mast_loading.lvb"
        )
        if not exists(loadcase_path):
            return

        self.modeller.fileOpen(loadcase_path)
        self.loading_script_found = True


class Monopole(LUSASSession):
    def __init__(self, filename):
        super(Monopole, self).__init__(filename)
        
    def subclass_specific_logic(self):
        """
    
        """
        lines_to_be_reversed = []
        for line in self.selection.add("Line", "All").getObjects("Line", "All"):
            start_point = line.getStartPoint()
            end_point = line.getEndPoint()
            p1, p2 = start_point.getXYZ(), end_point.getXYZ()
            if p2[2] - p1[2] < 0:  # 1) Orient Lines Upwards
                lines_to_be_reversed.append(line.getID())
        
        # Perform reversal
        self.geom.setAllDefaults()
        self.geom.cycleReverse(True)
        self.geom.cycleReset(False)
        self.selection.remove("All")
        
        self.selection.add("Line", lines_to_be_reversed)
        self.selection.reverse(self.geom)
        self.selection.remove("All")

    def assign_loadcase_from_lvb(self):
        """
        All Mast loadcases are the same, so only creating one file, and apply
        this to all.
        """
        if not self.fe_params["performStaticAnalysis"]:
            return
        
        loadcase_path = join(
            os.getcwd(), LUSASSession.LOADCASE_FOLDER, "mast_loading.lvb"
        )
        if not exists(loadcase_path):
            return

        self.modeller.fileOpen(loadcase_path)
        self.loading_script_found = True


class WindTurbine(LUSASSession):
    def __init__(self):
        raise NotImplementedError("Wind Turbine not implemented.")
