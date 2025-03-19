# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 08:08:37 2024

@author: Tristan Gowdridge


To run this script, make sure there is an IE model JSON in the ie_models/
folder and an appropriates names loadcase script.
"""
from LUSAS_classes import allocate_subclass

ie_model = f"bridge-1-2-10_v10_dev.json"
allocate_subclass(ie_model)