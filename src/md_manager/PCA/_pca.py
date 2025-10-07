from .._Traj import Traj
from ..df_utils import dihedral_angles, angles

import numpy as np
import pandas as pd
from scipy.linalg import eigh, ishermitian

def cartesian_pca(traj:Traj) -> tuple[np.ndarray]:
    pos = load_atom_coordinates(traj)
    cov = get_cartesian_covariance_matrix(pos)
    
    if not ishermitian(cov):
        raise ValueError("Covariance matrix should be hermisian")
    
    eigenvalues, eigenvectors = eigh(cov)
    return eigenvalues, eigenvectors

def dihedral_pca(traj:Traj) -> tuple[np.ndarray]:
    angles = load_local_conformation(traj)


def get_cartesian_covariance_matrix(pos:np.ndarray) -> np.ndarray:
    """
    
    """
    Nframe, Natom, Ndim = pos.shape
    if Ndim != 3:
        raise ValueError(f"Positions should be three dimentional")
    
    pos = np.reshape(pos, (Nframe, 3*Natom))
    return np.cov(pos)

def load_atom_coordinates(traj:Traj) -> np.ndarray:
    """
    Returns a numpy array of shape (Nframe, Natom, 3) that contains xyz coordinates of atoms
    """
    Nframe = len(traj)
    Natom = len(traj.top)

    data = np.zeros((Nframe, Natom, 3))
    xyz = ["x", "y", "z"]
    for i, df in enumerate(traj):
        data[i, :, :] = df[xyz].values
    return data

def load_local_conformation(traj:np.ndarray) -> np.ndarray:
    """
    
    """
    Nframe = len(traj)

    Natom  = len(traj.top)
    Nca    = len(traj.top.query("name == 'CA'"))
    if Natom != Nca:
        use_slice = True
        idx = traj.top.query("name == 'CA'").index
    else:
        use_slice = False

    Nchain = 1
    if "chain" in traj.top.columns:
        Nchain = len(traj.top.chain.unique())

    

    

    

    
