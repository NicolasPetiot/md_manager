from ..utils._ANM_utils import _jit_ANM_hessian, _jit_pfANM_hessian, _jit_local_MSF
from ..utils._MD_utils import pdb2df

import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.spatial import distance_matrix

__all__ = [
    "ANM_hessian",
    "pfANM_hessian",
    "collective_modes",
    "predicted_Bfactors"
]

def ANM_hessian(nodes_position = None, nodes_mass = None, df = None, distance_inter_nodes = None, cutoff_radius = 5.0, spring_constant = 1.0) -> np.ndarray:
    """
    
    """
    NoneType = type(None)
    if type(df) != NoneType:
        nodes_position = df[["x", "y", "z"]].to_numpy(dtype = float)
        nodes_mass = df.m.to_numpy(dtype = float)

    elif type(nodes_position) == NoneType or type(nodes_mass) == NoneType:
        raise ValueError("`nodes_position` and/or `nodes_mass` are NoneType.")
    
    if type(nodes_position) != np.ndarray:
        nodes_position = np.array(nodes_position, dtype=float)

    if type(nodes_mass) != np.ndarray:
        nodes_mass = np.array(nodes_mass, dtype=float)
    
    if type(distance_inter_nodes) == NoneType:
        distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

    return _jit_ANM_hessian(nodes_position, nodes_mass, distance_inter_nodes, spring_constant, cutoff_radius)

def pfANM_hessian(df:pd.DataFrame = None, nodes_position:np.ndarray = None, nodes_mass:np.ndarray = None, distance_inter_nodes:np.ndarray = None, spring_constant = 1.0) -> np.ndarray:
    """
    
    """
    if df is not None:
        nodes_position = df[["x", "y", "z"]].to_numpy(dtype = float)
        nodes_mass = df.m.to_numpy(dtype = float)

    elif (nodes_position is None) or (nodes_mass is None):
        raise ValueError("Missing input. Please specify either `df` or `nodes_position` and `nodes_mass`")

    if distance_inter_nodes is None:
        distance_inter_nodes = distance_matrix(nodes_position, nodes_position)
        
    return _jit_pfANM_hessian(nodes_position, nodes_mass, distance_inter_nodes, spring_constant)

def collective_modes(hessian):
    """
    Returns a tuple containing the eigenvalues and eigenvectors of the inputed mass-weighted hessian. 

    The mode of frequency `eigenvals[k]` is given by `mode = eigenvecs[:, k]`.
    
    Caution : 6 firsts null eigenmodes have been removed.
    """
    eigenvals, eigenvecs = eigh(hessian)
    eigenvals = eigenvals[6:]
    eigenvecs = eigenvecs[:, 6:]
    
    return eigenvals, eigenvecs

def predicted_Bfactors(df:pd.DataFrame = None, eigenvalues:np.ndarray = None, eigenvectors:np.ndarray = None, convert2bfactors = True):
    if df is None:
        raise ValueError("Missing input. Please specify `df`.")

    if (eigenvalues is None) or (eigenvectors is None):        
        H = pfANM_hessian(df=df)
        eigenvalues, eigenvectors = collective_modes(H)
        
    
    # B-factors computation:
    nodes_mass = df.m.to_numpy()

    return _jit_local_MSF(eigenvalues, eigenvectors, nodes_mass, convert2bfactors=True)