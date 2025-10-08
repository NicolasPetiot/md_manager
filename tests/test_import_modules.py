import pytest
import importlib

modules_packages = [
    ("md_manager", None), 
    ("md_manager", "CUTABI"), 
    ("md_manager", "NMA"), 
    ("md_manager", "Graph"), 
    ("md_manager", "PCA")
]

@pytest.mark.parametrize("module, package", modules_packages)
def test_import(module:str, package:str) -> bool:
    """
    
    """
    assert can_import(module, package)
    
def can_import(module:str, package:str = None) -> bool:
    """
    
    """
    try:
        importlib.import_module(module, package)
        return True

    except ImportError:
        return False