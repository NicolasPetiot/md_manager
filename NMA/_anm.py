from .._parameters import *

import numpy as np
from numba import jit
from scipy.linalg import eigh
from scipy.spatial import distance_matrix

__all__ = ["ANM_hessian", "pfANM_hessian", "collective_modes", "local_MSF", "non_local_MSF"]

### Hessian computations ###
def ANM_hessian(nodes_position = None, nodes_mass = None, df = None, distance_inter_nodes = None, cutoff_radius = 5.0, spring_constant = 1.0):
    """
    Returns a np.ndarray that corresponds to the mass-weighted Hessian for the set of input nodes.

    Input can be 'nodes_position' and 'nodes_mass' or 'df'. By default the spring constant is set to 1.0 kcal/mol/AA² and the cutoff radius is set to 5.0AA.
    """
    if df != None :
        nodes_position = df[["x", "y", "z"]].to_numpy(dtype = float)
        nodes_mass = df.m.to_numpy(dtype = float)

    elif nodes_position == None or nodes_mass == None:
        raise ValueError("Please use nodes_position and nodes_mass or df as argument(s).")
    
    if distance_inter_nodes == None :
        distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

    return compiled_ANM_hessian_builder(nodes_position, nodes_mass, distance_inter_nodes, cutoff_radius, spring_constant)
    

def pfANM_hessian(nodes_position = None, nodes_mass = None, df = None, distance_inter_nodes = None, spring_constant = 1.0):
    """
    Returns a np.ndarray that corresponds to the parameter-free mass-weighted Hessian for the set of input nodes.

    Input can be 'nodes_position' and 'nodes_mass' or 'df'. By default the spring constant is set to 1.0 kcal/mol/AA².
    """
    if df != None :
        nodes_position = df[["x", "y", "z"]].to_numpy(dtype = float)
        nodes_mass = df.m.to_numpy(dtype = float)

    elif nodes_position == None or nodes_mass == None:
        raise ValueError("Please use nodes_position and nodes_mass or df as argument(s).")

    if distance_inter_nodes == None :
        distance_inter_nodes = distance_matrix(nodes_position, nodes_position)
    
    return compiled_pfANM_hessian_builder(nodes_position, nodes_mass, distance_inter_nodes, spring_constant)

def collective_modes(hessian):
    """
    Returns a tuple containing the eigenvalues and eigenfrequencies of the inputed mass-weighted hessian. Note that the 6 firsts null eigenmodes have been removed.
    """
    eigenvals, eigenvecs = eigh(hessian)
    eigenvals = eigenvals[6:]
    eigenvecs = eigenvecs[:, 6:]
    
    return eigenvals, eigenvecs

### locan & non-local Mean Squared Fluctuations ###
@jit(nopython = True, cache=True)
def local_MSF(eigenvals:np.ndarray, eigenvecs:np.ndarray, nodes_mass:np.ndarray, convert2bfactors = False):
    """
    Returns a numpy.ndarray containing the predicted local mean squared fluctuations associated to the inputed collective modes. 
    If convert2bfactors is set to True, the local MSF is multiplied by 8*pi^2/3 which is the convertion to Bfactors coefficients.
    The node_mass array must have a length of Nnodes, the eigenvals array must have a length of 3*Nnodes - 6 and the eigenvecs array must have a shape of (3*Nnodes, 3*Nnodes - 6).
    The output array will have a length of Nnodes.
    """
    Nnode = len(nodes_mass)
    b = np.zeros(Nnode)

    # scaling factor :
    KT = BOLTZMANN * TEMPERATURE
    if convert2bfactors :
        KT = KT * 8*np.pi**2/3
    
    for i in range(Nnode):
        vi = eigenvecs[3*i:3*(i+1), :]
        mi = nodes_mass[i]

        bi = np.sum(np.sum(vi**2, axis=0) / eigenvals)
        b[i] = bi / mi

    return KT * b

@jit(nopython = True, cache=True)
def non_local_MSF(eigenvals:np.ndarray, eigenvecs:np.ndarray, nodes_mass:np.ndarray):
    """

    """
    Nnode = len(nodes_mass)
    d = np.zeros((Nnode, Nnode))

    # scaline factor :
    KT = BOLTZMANN * TEMPERATURE

    for i in range(Nnode):
        mi = nodes_mass[i]
        vi = eigenvecs[3*i:3*(i+1), :] / np.sqrt(mi)
        for j in range(i+1, Nnode):
            mj = nodes_mass[j]
            vj = eigenvecs[3*j:3*(j+1), :] / np.sqrt(mj)

            dij = np.sum(np.sum((vi - vj)**2, axis = 0) / eigenvals)
            d[i, j] = dij
            d[j, i] = dij

    return d

@jit(nopython = True, cache = True)
def compiled_pfANM_hessian_builder(nodes_position:np.ndarray, nodes_mass:np.ndarray, distance_inter_nodes:np.ndarray, spring_constant:float):
    """
    
    """
    Nnodes = len(nodes_position)
    hessian = np.zeros((3*Nnodes, 3*Nnodes))
    squared_mass = np.zeros((3*Nnodes, 3*Nnodes))

    # loop over nodes :
    for i in range(Nnodes):
        Ri = nodes_position[i, :]
        mi = nodes_mass[i]
        for j in range(i, Nnodes):
            mj = nodes_mass[j]

            # update squared masses :
            squared_mass[3*i:3*(i+1), 3*j:3*(j+1)] = mi*mj
            squared_mass[3*j:3*(j+1), 3*i:3*(i+1)] = mi*mj

            # update hessian :
            if i != j :
                Rj = nodes_position[j, :]
                Rij = Rj - Ri
                Rij4 = distance_inter_nodes[i, j]**4

                Hij = np.zeros((3, 3))
                for xi in range(3):
                    x = Rij[xi]
                    for yi in range(xi, 3):
                        y = Rij[yi]

                        Hij[xi, yi] = x*y
                        Hij[yi, xi] = x*y

                Hij = -spring_constant * Hij / Rij4
                
                hessian[3*i:3*(i+1), 3*j:3*(j+1)] = Hij
                hessian[3*j:3*(j+1), 3*i:3*(i+1)] = Hij
                hessian[3*i:3*(i+1), 3*i:3*(i+1)] -= Hij
                hessian[3*j:3*(j+1), 3*j:3*(j+1)] -= Hij

    return hessian / np.sqrt(squared_mass)

@jit(nopython = True, cache = True)
def compiled_ANM_hessian_builder(nodes_position, nodes_mass, distance_inter_nodes, cutoff_radius = 5.0, spring_constant = 1.0):
    Nnodes = len(nodes_position)
    hessian = np.zeros((3*Nnodes, 3*Nnodes))
    squared_mass = np.zeros((3*Nnodes, 3*Nnodes))

    # loop over nodes :
    for i in range(Nnodes):
        Ri = nodes_position[i, :]
        mi = nodes_mass[i]
        for j in range(i, Nnodes):
            mj = nodes_mass[j]

            # update squared masses :
            squared_mass[3*i:3*(i+1), 3*j:3*(j+1)] = mi*mj
            squared_mass[3*j:3*(j+1), 3*i:3*(i+1)] = mi*mj

            # update hessian :
            if i != j and distance_inter_nodes[i, j] <= cutoff_radius:
                Rj = nodes_position[j, :]
                Rij = Rj - Ri
                Rij2 = distance_inter_nodes[i, j]**2

                Hij = np.zeros((3, 3))
                for xi in range(3):
                    x = Rij[xi]
                    for yi in range(xi, 3):
                        y = Rij[yi]

                        Hij[xi, yi] = x*y
                        Hij[yi, xi] = x*y

                Hij = -spring_constant * Hij / Rij2
                
                hessian[3*i:3*(i+1), 3*j:3*(j+1)] = Hij
                hessian[3*j:3*(j+1), 3*i:3*(i+1)] = Hij
                hessian[3*i:3*(i+1), 3*i:3*(i+1)] -= Hij
                hessian[3*j:3*(j+1), 3*j:3*(j+1)] -= Hij

    return hessian / np.sqrt(squared_mass)