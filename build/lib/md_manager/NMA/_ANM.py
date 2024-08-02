from .._params import BOLTZMANN, TEMPERATURE

import numpy as np
import pandas as pd
from numba import njit
from scipy.linalg import eigh
from scipy.spatial import distance_matrix

__all__ = [
    "ANM_hessian",
    "pfANM_hessian",
    "local_MSF",
    "predicted_Bfactors"
]

@njit(cache=True)
def ANM_hessian(node_position:np.ndarray, node_mass:np.ndarray, distance_inter_node:np.ndarray, cutoff_radius:float, spring_constant = 1.0) -> np.ndarray:
    """Returns a mass-weighted Hessian based on ANM model"""
    Nnodes = len(node_position)
    hessian = np.zeros((3*Nnodes, 3*Nnodes))

    # loop over nodes :
    for i in range(Nnodes):
        Ri = node_position[i, :]
        Mi = node_mass[i]
        for j in range(i+1, Nnodes):
            if distance_inter_node[i, j] < cutoff_radius:
                Rj = node_position[j, :]
                Mj = node_mass[j]

                Hij = jit_Hij(coordI=Ri, coordJ=Rj)
                Hij = -Hij / distance_inter_node[i, j]**2 

                hessian[3*i:3*(i+1), 3*j:3*(j+1)] = Hij / np.sqrt(Mi*Mj)
                hessian[3*j:3*(j+1), 3*i:3*(i+1)] = Hij / np.sqrt(Mi*Mj)
                hessian[3*i:3*(i+1), 3*i:3*(i+1)] -= Hij / Mi
                hessian[3*j:3*(j+1), 3*j:3*(j+1)] -= Hij / Mj
    
    return spring_constant * hessian

@njit(cache=True)
def pfANM_hessian(node_position:np.ndarray, node_mass:np.ndarray, distance_inter_node:np.ndarray, spring_constant = 1.0) -> np.ndarray:
    """Returns a mass-weighted Hessian based on pfANM model"""
    Nnodes = len(node_position)
    hessian = np.zeros((3*Nnodes, 3*Nnodes))

    # loop over nodes :
    for i in range(Nnodes):
        Ri = node_position[i, :]
        Mi = node_mass[i]
        for j in range(i+1, Nnodes):
            Rj = node_position[j, :]
            Mj = node_mass[j]

            Hij = jit_Hij(coordI=Ri, coordJ=Rj)
            Hij = -Hij / distance_inter_node[i, j]**4 

            hessian[3*i:3*(i+1), 3*j:3*(j+1)] = Hij / np.sqrt(Mi*Mj)
            hessian[3*j:3*(j+1), 3*i:3*(i+1)] = Hij / np.sqrt(Mi*Mj)
            hessian[3*i:3*(i+1), 3*i:3*(i+1)] -= Hij / Mi
            hessian[3*j:3*(j+1), 3*j:3*(j+1)] -= Hij / Mj
    
    return spring_constant * hessian

@njit(cache=True)
def local_MSF(eigenvals:np.ndarray, eigenvecs:np.ndarray, node_mass:np.ndarray, convert2bfactors = False) -> np.ndarray:
    """Computes the Mean Squared Fluctuations of the atoms arround their equilibrium position based on normal modes"""
    Nnode = len(node_mass)
    b = np.zeros(Nnode)

    # scaling factor :
    KT = BOLTZMANN * TEMPERATURE
    if convert2bfactors :
        KT = KT * 8*np.pi**2/3
    
    for i in range(Nnode):
        vi = eigenvecs[3*i:3*(i+1), :]
        mi = node_mass[i]

        bi = np.sum(np.sum(vi**2, axis=0) / eigenvals)
        b[i] = bi / mi

    return KT * b

def predicted_Bfactors(df:pd.DataFrame, spring_constant = 1.0) -> pd.Series:
    """Use the pfANM model to compute Normal Modes and deduce the predicted B-factors of the input structure"""
    # Hessian:
    xyz = ["x", "y", "z"]
    node_position = df[xyz].to_numpy()
    node_mass = df.m.to_numpy()
    hessian = pfANM_hessian(node_position, node_mass, distance_matrix(node_position, node_position), spring_constant)

    # Normal Modes:
    eigenvals, eigenvecs = eigh(hessian)
    eigenvals = eigenvals[6:]
    eigenvecs = eigenvecs[:, 6:]

    # Thermal B-factors:
    return pd.Series(
        local_MSF(eigenvals, eigenvecs, node_mass, convert2bfactors=True), 
        index=df.index, name="b"
    )

@njit(cache=True)
def jit_Hij(coordI:np.ndarray, coordJ:np.ndarray) -> np.ndarray:
    """
    Is the compiled version of the dot product Rij[:, None] @ Rij[None, :].
    """
    Rij = coordJ - coordI
    Hij = np.zeros((3, 3))
    for xi in range(3):
        x = Rij[xi]
        for yi in range(xi, 3):
            y = Rij[yi]

            Hij[xi, yi] = x*y
            Hij[yi, xi] = x*y
    return Hij