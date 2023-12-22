from ..utils._ANM_utils import _jit_ANM_hessian, _jit_pfANM_hessian, _jit_local_MSF
from .._PDBfile import pdb

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
    if df != None:
        nodes_position = df[["x", "y", "z"]].to_numpy(dtype = float)
        nodes_mass = df.m.to_numpy(dtype = float)

    elif nodes_position == None or nodes_mass == None:
        raise ValueError("`nodes_position` and/or `nodes_mass` are NoneType.")
    
    if distance_inter_nodes == None:
        distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

    return _jit_ANM_hessian(nodes_position, nodes_mass, distance_inter_nodes, spring_constant, cutoff_radius)

def pfANM_hessian(nodes_position = None, nodes_mass = None, df = None, distance_inter_nodes = None, spring_constant = 1.0) -> np.ndarray:
    """
    
    """
    NoneType = type(None)
    if type(df) != NoneType:
        nodes_position = df[["x", "y", "z"]].to_numpy(dtype = float)
        nodes_mass = df.m.to_numpy(dtype = float)

    elif type(nodes_position) == NoneType or type(nodes_mass) == NoneType:
        raise ValueError("`nodes_position` and/or `nodes_mass` are NoneType.")
    
    if type(distance_inter_nodes) == NoneType:
        distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

    return _jit_pfANM_hessian(nodes_position, nodes_mass, distance_inter_nodes, spring_constant)

def collective_modes(hessian):
    """
    Returns a tuple containing the eigenvalues and eigenfrequencies of the inputed mass-weighted hessian. 

    The mode of frequency `eigenvals[k]` is given by `mode = eigenvecs[:, k]`.
    
    Caution : 6 firsts null eigenmodes have been removed.
    """
    eigenvals, eigenvecs = eigh(hessian)
    eigenvals = eigenvals[6:]
    eigenvecs = eigenvecs[:, 6:]
    
    return eigenvals, eigenvecs

def predicted_Bfactors(
        eigenvals:np.ndarray=None, eigenvecs:np.ndarray=None, nodes_mass:np.ndarray = None,
        hessian:np.ndarray=None, PDB:pdb = None, df:pd.DataFrame=None,
        spring_constant=1.0, cutoff_radius=None, 
        convert2bfactors=True
    ):
    """
    Compute the thermal B-factors from the input informations.

    The minimal use case is :
    ```
    import md_manager as md

    with md.pdb(filename) as pdb:
        b = md.NMA.predicted_Bfactors(pdb)
    ```

    To have a deeper overview of the use of this function. See "https://github.com/NicolasPetiot/md_manager/"
    """
    # Compute collective modes if needed.
    NoneType = type(None)
    if type(eigenvals) == NoneType or type(eigenvecs) == NoneType:
        if type(hessian) == NoneType:
            if type(PDB)==NoneType:
                raise ValueError("Not enougth information. Unable to compute collective modes.")
            
            df = PDB.read2df()
            if type(cutoff_radius) == NoneType:
                hessian = pfANM_hessian(df=df, spring_constant=spring_constant)
            else :
                hessian = ANM_hessian(df=df, spring_constant=spring_constant, cutoff_radius=cutoff_radius)
        eigenvals, eigenvecs = collective_modes(hessian)

    # compute nodes mass if needed :
    if type(nodes_mass) == NoneType:
        if type(df) == NoneType :
            if type(PDB) == NoneType:
                raise ValueError("Not enougth information. Unable to compute nodes mass.")
            df = PDB.read2df()
        nodes_mass = df.m.to_numpy(dtype = float)

    return _jit_local_MSF(eigenvals, eigenvecs, nodes_mass, convert2bfactors)
    