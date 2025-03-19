import sys
import os
import pytest
from unittest.mock import patch, mock_open
from FE_to_IE.preprocessing_iemodel import (
    distance, ElementInfo, get_material_properties, get_geometric_properties,
    extract_coordinates, create_line_end_coords
)

# Add the parent directory (root of the project) to the system path.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))


class TestDistance:
    def test_distance1(self):
        result = distance(
            (0, 0, 0),
            (0, 0, 0)
        )
        assert result == 0

    def test_distance2(self):
        result = distance(
            (0, 0, 1),
            (0, 0, 0)
        )
        assert result == 1

    def test_distance3(self):
        point1 = (1, 2, 3)
        point2 = (4, 5, 6)
        result = distance(point1, point2)
        expected_result = 5.196152422706632
        assert result == pytest.approx(expected_result, rel=1e-9)


class TestElementInfo:
    def test_element_info_str(self):
        element = ElementInfo("Element1", "beam")
        element.coords.add((1, 2, 3))
        element.coords.add((4, 5, 6))
        
        expected_str = "Element1: beam\n(1, 2, 3)\n(4, 5, 6)"
        assert str(element) == expected_str


class TestGetMaterialProperties:
    @patch("builtins.open", mock_open(read_data='{"materialName": {"youngsModulus": 210000, "density": 7800}}'))
    def test_get_material_properties(self):
        element = {
            "type": "beam",
            "material": {
                "matrix": {
                    "base": {
                        "type": {
                            "name": "materialName"
                        }
                    }
                },
                "properties": [
                    {"type": "youngsModulus", "value": 210000},
                    {"type": "density", "value": 7800}
                ]
            }
        }
        material_name, material_properties = get_material_properties(element)
        
        # Expected values from mock material.json data
        assert material_name == "materialName"
        assert material_properties["youngsModulus"] == 210000
        assert material_properties["density"] == 7800

    def test_get_material_properties_ground(self):
        element = {"type": "ground"}
        
        result = get_material_properties(element)
        
        assert result == {}

    def test_get_material_properties_complete(self):
        element = {
            "type": "beam",
            "material": {
                "matrix": {
                    "base": {
                        "type": {
                            "name": "steel"
                        }
                    }
                },
                "properties": [
                    {"type": "youngsModulus", "value": 210e9},
                    {"type": "density", "value": 7800},
                ]
            }
        }
        
        material_name, material_properties = get_material_properties(element)
        
        # Check if the material name is correctly extracted
        assert material_name == "steel"
        # Check if the required properties are in the material properties
        assert material_properties["youngsModulus"] == 210e9
        assert material_properties["density"] == 7800
        assert material_properties["poissonsRatio"] == 0.3
        assert material_properties["linearThermalExpansionCoefficient"] == 12e-6

    def test_get_material_properties_with_defaults(self):
        element = {
            "type": "beam",
            "material": {
                "matrix": {
                    "base": {
                        "type": {
                            "name": "steel"
                        }
                    }
                },
                "properties": [
                    {"type": "youngsModulus", "value": 210000000},
                    # Missing density, poissonsRatio, and linearThermalExpansionCoefficient
                ]
            }
        }
        
        material_name, material_properties = get_material_properties(element)
        
        # Check if the material name is correctly extracted
        assert material_name == "steel"
        # Check if the missing properties were filled with default values
        assert material_properties["youngsModulus"] == 210000000
        assert material_properties["density"] == 7850
        assert material_properties["poissonsRatio"] == 0.3
        assert material_properties["linearThermalExpansionCoefficient"] == 1.2e-5

    def test_get_material_properties_material_not_found(self):
        element = {
            "type": "beam",
            "material": {
                "matrix": {
                    "base": {
                        "type": {
                            "name": "unknownMaterial"
                        }
                    }
                },
                "properties": [
                    {"type": "youngsModulus", "value": 210000},
                    {"type": "density", "value": 7800}
                ]
            }
        }
        
        with pytest.raises(KeyError):
            get_material_properties(element)

    def test_get_material_properties_missing_properties(self):
        element = {
            "type": "beam",
            "material": {
                "matrix": {
                    "base": {
                        "type": {
                            "name": "concrete"
                        }
                    }
                },
                # No properties, just use defaults
            }
        }
        
        material_name, material_properties = get_material_properties(element)
        
        # Check if the material name is correctly extracted
        assert material_name == "concrete"
        # Check if the properties are filled with the defaults
        assert material_properties["youngsModulus"] == 35e9
        assert material_properties["density"] == 2400
        assert material_properties["poissonsRatio"] == 0.2
        assert material_properties["linearThermalExpansionCoefficient"] == 1e-5


class TestGetGeometricProperties():
    def test_need_to_write(self):
        assert True


class TestExtractCoordinates():
    def test_need_to_write(self):
        assert True


class TestCreateLineEndCoordinates():
    def test_need_to_write(self):
        assert True


class TestExtractIEInformation():
    def test_need_to_write(self):
        assert True


if __name__ == "__main__":
    pytest.main()
