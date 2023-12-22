from ._params import BOLTZMANN, TEMPERATURE

import numpy as np
from numba import jit

__all__ = [
    "_jit_ANM_hessian",
    "_jit_pfANM_hessian",
    "_jit_local_MSF"
]

@jit(nopython = True, cache = True)
def _jit_ANM_hessian(nodes_position:np.ndarray, nodes_mass:np.ndarray, distance_inter_nodes:np.ndarray, spring_constant:float, cutoff_radius:float):
    """
    
    """
    Nnodes = len(nodes_position)
    hessian = np.zeros((3*Nnodes, 3*Nnodes))

    # loop over nodes :
    for i in range(Nnodes):
        Ri = nodes_position[i, :]
        Mi = nodes_mass[i]
        for j in range(1, Nnodes):
            if distance_inter_nodes[i, j] < cutoff_radius:
                Rj = nodes_position[j, :]
                Mj = nodes_mass[j]

                Hij = _jit_Hij(coordI=Ri, coordJ=Rj)
                Hij = -Hij / distance_inter_nodes[i, j]**2 / np.sqrt(Mi*Mj)

                hessian[3*i:3*(i+1), 3*j:3*(j+1)] = Hij
                hessian[3*j:3*(j+1), 3*i:3*(i+1)] = Hij
                hessian[3*i:3*(i+1), 3*i:3*(i+1)] -= Hij
                hessian[3*j:3*(j+1), 3*j:3*(j+1)] -= Hij
    
    return spring_constant * hessian

@jit(nopython = True, cache = True)
def _jit_pfANM_hessian(nodes_position:np.ndarray, nodes_mass:np.ndarray, distance_inter_nodes:np.ndarray, spring_constant:float):
    """
    
    """
    Nnodes = len(nodes_position)
    hessian = np.zeros((3*Nnodes, 3*Nnodes))

    # loop over nodes :
    for i in range(Nnodes):
        Ri = nodes_position[i, :]
        Mi = nodes_mass[i]
        for j in range(i+1, Nnodes):
            Rj = nodes_position[j, :]
            Mj = nodes_mass[j]

            Hij = _jit_Hij(coordI=Ri, coordJ=Rj)
            Hij = -Hij / distance_inter_nodes[i, j]**4 / np.sqrt(Mi*Mj)

            hessian[3*i:3*(i+1), 3*j:3*(j+1)] = Hij
            hessian[3*j:3*(j+1), 3*i:3*(i+1)] = Hij
            hessian[3*i:3*(i+1), 3*i:3*(i+1)] -= Hij
            hessian[3*j:3*(j+1), 3*j:3*(j+1)] -= Hij
    
    return spring_constant * hessian


# not in __all__
@jit(nopython = True, cache = True)
def _jit_Hij(coordI, coordJ):
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

@jit(nopython = True, cache=True)
def _jit_local_MSF(eigenvals:np.ndarray, eigenvecs:np.ndarray, nodes_mass:np.ndarray, convert2bfactors = False):
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
def _jit_non_local_MSF(eigenvals:np.ndarray, eigenvecs:np.ndarray, nodes_mass:np.ndarray):
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

    return KT * d
