# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 11:03:47 2024

@author: trist
"""
import os
import ast
import re
import math
import json
import warnings
from decimal import Decimal
from collections import defaultdict


def order_coords(coords):
    """
    Order coordinates for consistent processing. The ordered coordinates are
    used as a dictionary key to reference their geometry.
    """
    coords = list(coords)
    ordered = sorted(coords, key=lambda coord: (coord[0], coord[1], coord[2]))
    return tuple(ordered)


def distance(p1, p2):
    """
    Calculate the distance between two 3D points in their tuple form.
    """
    return math.sqrt(
        (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
    )


def greedy_cyclic_sort(points):
    """
    Sort the points in a greedy cyclic manner. The reasoning for sorting by the
    closest point is because LUSAS requires the points in order to define a
    surface. Applying this greedy methodology assumes that the points are
    enclosing a convex hull. An issue arises later on where the normal vectors
    of surfaces are defined by the perimeter points being clockwise or
    anticlockwise.
    """
    sorted_points = [points[0]]
    remaining_points = points[1:]
    
    # Continue until all points are added
    while remaining_points:
        last_point = sorted_points[-1]
        # Find the nearest neighbor in the remaining points
        nearest_point = min(remaining_points, key=lambda point: distance(last_point, point))
        sorted_points.append(nearest_point)
        remaining_points.remove(nearest_point)
    
    return sorted_points


def unit_scaler(value, unit):
    """
    Scale units to the set standard units.
    """
    # Convert value to Decimal for higher precision in calculations
    value = Decimal(value)
    
    match unit:
        # Force
        case "kgf":  # kilograms of force
            scaled_value = Decimal("0.001") * value
        
        # Length
        case "mm":  # millimetre
            scaled_value = Decimal("0.001") * value
        case "cm":  # centimetre
            scaled_value = Decimal("0.01") * value
        case "m":  # metre
            scaled_value = 1 * value  # Already in meters
        case "km":  # kilometre
            scaled_value = Decimal("1000") * value
        
        # Mass
        case "kg":  # kilogram
            warnings.warn("Needs converting to tonnes? need to double check")
            scaled_value = Decimal("0.001") * value  # Convert to metric tonnes if needed
        
        # Time
        case "s":  # second
            scaled_value = 1 * value  # Already in seconds
        
        # Temperature
        case "C":  # celsius
            scaled_value = 1 * value  # Already in Celsius
        case "F":  # fahrenheit
            scaled_value = (value - Decimal("32")) * Decimal("5") / Decimal("9")
        case "K":  # kelvin
            scaled_value = value - Decimal("273.15")
        
        # Pressure
        case "Pa":  # Pascals
            scaled_value = 1 * value
        case "kPa":  # Kilopascals
            scaled_value = 1 * value
        case "MPa":  # Megapascals
            scaled_value = 1 * value
        case "GPa":  # Gigapascals
            scaled_value = 1 * value
        case "psi":  # Pound per square inch
            scaled_value = 1 * value
        case "ksi":  # Kilopound per square inch
            scaled_value = 1 * value
        case "Mpsi":  # Megapound per square inch
            scaled_value = 1 * value
        
        # Density
        case "kg/m^3":  # Kilogram per metre cubed
            scaled_value = 1 * value
        case "g/cm^3":  # Gram per centimetre cubed
            scaled_value = 1 * value
        case "kg/L":  # Kilogram per litre
            scaled_value = 1 * value
        case "g/mL":  # Gram per millilitre
            scaled_value = 1 * value
        case "t/m^3":  # Tonnes per metre cubed
            scaled_value = 1 * value
        case "kg/dm^3":  # Kilogram per decimetre cubed
            scaled_value = 1 * value
        case "oz/cu in":  #
            print("Possibly ounces per cubic inch, checking with Dan")
            scaled_value = 1 * value
        
        # Miscellaneous
        case "degrees":
            warnings.warn("Need to check if degrees is standard")
            scaled_value = 1 * value
        case "radians":
            warnings.warn("Need to check if degrees is standard")
            scaled_value = value * Decimal("57.2958")  # Convert to degrees
        
        # If unit is unknown
        case _:
            warnings.warn(f"The input unit '{unit}' was not detected, and the value has not been scaled.")
            scaled_value = value

    return float(scaled_value)  # Convert back to float if preferred