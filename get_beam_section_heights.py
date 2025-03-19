# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 08:23:05 2025

@author: trist
"""
import json

from ie_to_fe_funcs import LUSASSession


with open("./data/nonparametric_beam_geometries.json", "r") as f:
    previous_info = json.load(f)


bridge_identifier = f"bridge-ladder-1_v10_Test_dev"
ie_path = f"./ie_models/{bridge_identifier}.json"
lusas_sesh = LUSASSession(ie_path, bridge_identifier)

for profile_letter, v in previous_info.items():
    for dimensions, info in v.items():
        beam = f"{profile_letter}{info['N']}"
        if "M" == profile_letter:  # M is broken as explained below.
            continue
        attr = lusas_sesh.database.createGeometricLine("LGeo1")
        attr.setValue("elementType", "3D Thick Beam")
        attr.setFromLibrary("UK Sections", f"Precast {profile_letter} Beams", beam, 0, 0, 0)
        attr.setAnalysisCategory("3D")
        
        top = attr.getFibrePosition("A2")[0]  # coordinates are listed (z, y)
        bottom = attr.getFibrePosition("A3")[0]

        info["height"] = round(top - bottom, 3)


# For some reason, the fibre locations for M cannot be returned, they exists in
# the GUI. I'm assuming they're called something other than A2 and A3 under the
# hood. I'll manually enter these.
m_dims = {
    2:  0.454673 + 0.265327,
    3:  0.490269 + 0.309731,
    4:  0.526607 + 0.353393,
    5:  0.603151 + 0.356849,
    6:  0.630852 + 0.409148,
    7:  0.660437 + 0.459563,
    8:  0.746053 + 0.453947,
    9:  0.767710 + 0.512290,
    10: 0.792525 + 0.567474
}

for dimension, info in previous_info["M"].items():
    info["height"] = round(m_dims[info["N"]], 3)
    
with open("./data/nonparametric_beam_geometries.json", "w") as f:
    json.dump(previous_info, f, indent=4)
