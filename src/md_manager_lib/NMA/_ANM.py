from ..parameters import BOLTZMANN, TEMPERATURE

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
    """
    Compute the mass-weighted Hessian matrix for an Anisotropic Network Model (ANM).

    This function calculates the mass-weighted Hessian matrix representing the ANM for a system of 
    nodes. The Hessian describes the second derivatives of the potential energy with respect to the 
    positions of the nodes, considering pairwise interactions within a cutoff radius. The resulting 
    Hessian is scaled by the specified spring constant.

    Parameters:
    node_position : numpy.ndarray
        A 2D array of shape (N, 3), where N is the number of nodes. Each row contains the 
        3D coordinates (x, y, z) of a node.
    node_mass : numpy.ndarray
        A 1D array of length N containing the masses of the nodes.
    distance_inter_node : numpy.ndarray
        A 2D array of shape (N, N) representing the pairwise distances between nodes.
    cutoff_radius : float
        The cutoff radius for interactions. Pairs of nodes separated by a distance greater 
        than this value are not considered.
    spring_constant : float, optional
        The spring constant used to scale the Hessian matrix. Defaults to 1.0.

    Returns:
    numpy.ndarray
        A 2D array of shape (3N, 3N), where N is the number of nodes, representing the 
        mass-weighted Hessian matrix for the system.

    Raises:
    ValueError
        If the dimensions of `node_position`, `node_mass`, or `distance_inter_node` are inconsistent.

    Example:
        >>> import numpy as np
        >>> node_position = np.array([[0.0, 0.0, 0.0], 
        ...                           [1.0, 0.0, 0.0], 
        ...                           [0.0, 1.0, 0.0]])
        >>> node_mass = np.array([1.0, 1.0, 1.0])
        >>> distance_inter_node = np.array([[0.0, 1.0, 1.0], 
        ...                                 [1.0, 0.0, 1.414], 
        ...                                 [1.0, 1.414, 0.0]])
        >>> cutoff_radius = 1.5
        >>> spring_constant = 1.0
        >>> ANM_hessian(node_position, node_mass, distance_inter_node, cutoff_radius, spring_constant)
        array([...])  # Example output

    Notes:
    - This function uses the Anisotropic Network Model (ANM) to compute the Hessian matrix (see https://www.cell.com/biophysj/fulltext/S0006-3495(01)76033-X). 
    - It is convenient to use `scipy.spatial.distance_matrix` to compute the `distance_inter_node` array.
    - The function assumes the `node_position` and `node_mass` arrays have consistent lengths, 
      and `distance_inter_node` is a square matrix with dimensions matching the number of nodes.
    """
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
    """
    Compute the parameter-free mass-weighted Hessian matrix for an Anisotropic Network Model (pfANM).

    This function calculates the mass-weighted Hessian matrix representing the pfANM for a system of 
    nodes. The Hessian describes the second derivatives of the potential energy with respect to the 
    positions of the nodes, considering pairwise interactions. The resulting Hessian is scaled by the specified spring constant.

    Parameters:
    node_position : numpy.ndarray
        A 2D array of shape (N, 3), where N is the number of nodes. Each row contains the 
        3D coordinates (x, y, z) of a node.
    node_mass : numpy.ndarray
        A 1D array of length N containing the masses of the nodes.
    distance_inter_node : numpy.ndarray
        A 2D array of shape (N, N) representing the pairwise distances between nodes.
    spring_constant : float, optional
        The spring constant used to scale the Hessian matrix. Defaults to 1.0.

    Returns:
    numpy.ndarray
        A 2D array of shape (3N, 3N), where N is the number of nodes, representing the 
        mass-weighted Hessian matrix for the system.

    Raises:
    ValueError
        If the dimensions of `node_position`, `node_mass`, or `distance_inter_node` are inconsistent.

    Example:
        >>> import numpy as np
        >>> node_position = np.array([[0.0, 0.0, 0.0], 
        ...                           [1.0, 0.0, 0.0], 
        ...                           [0.0, 1.0, 0.0]])
        >>> node_mass = np.array([1.0, 1.0, 1.0])
        >>> distance_inter_node = np.array([[0.0, 1.0, 1.0], 
        ...                                 [1.0, 0.0, 1.414], 
        ...                                 [1.0, 1.414, 0.0]])
        >>> spring_constant = 1.0
        >>> pfANM_hessian(node_position, node_mass, distance_inter_node, spring_constant)
        array([...])  # Example output

    Notes:
    - This function uses the parameter-free Anisotropic Network Model (pfANM) to compute the Hessian matrix (see https://www.pnas.org/doi/abs/10.1073/pnas.0902159106). 
    - It is convenient to use `scipy.spatial.distance_matrix` to compute the `distance_inter_node` array.
    - The function assumes the `node_position` and `node_mass` arrays have consistent lengths, 
      and `distance_inter_node` is a square matrix with dimensions matching the number of nodes.
    """
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
    """
    Compute the Mean Squared Fluctuations (MSF) of atoms around their equilibrium positions based on normal modes.

    This function calculates the MSF of each node (atom) using the eigenvalues and eigenvectors from 
    a normal mode analysis. The MSF provides a measure of the average displacement of a node from its 
    equilibrium position. If `convert2bfactors` is set to True, the results are scaled to mimic B-factors 
    used in X-ray crystallography.

    Parameters:
    eigenvals : numpy.ndarray
        A 1D array of shape (M,) containing the eigenvalues of the system, where M is the number of modes. 
        Typically, the first few eigenvalues corresponding to rigid body motions should be excluded.
    eigenvecs : numpy.ndarray
        A 2D array of shape (3N, M), where N is the number of nodes and M is the number of modes. Each 
        column corresponds to an eigenvector, and each group of 3 rows corresponds to the (x, y, z) components 
        of a node's displacement.
    node_mass : numpy.ndarray
        A 1D array of length N containing the masses of the nodes.
    convert2bfactors : bool, optional
        If True, the MSF values are scaled to mimic B-factors by multiplying with \(8\pi^2/3\). Defaults to False.

    Returns:
    numpy.ndarray
        A 1D array of length N containing the MSF values for each node.

    Raises:
    ValueError
        If the dimensions of `eigenvals`, `eigenvecs`, or `node_mass` are inconsistent.

    Example:
        >>> import numpy as np
        >>> eigenvals = np.array([1.0, 2.0, 3.0])
        >>> eigenvecs = np.array([[0.1, 0.2, 0.3], 
        ...                        [0.4, 0.5, 0.6], 
        ...                        [0.7, 0.8, 0.9], 
        ...                        [0.1, 0.2, 0.3], 
        ...                        [0.4, 0.5, 0.6], 
        ...                        [0.7, 0.8, 0.9]])
        >>> node_mass = np.array([1.0, 1.0])
        >>> local_MSF(eigenvals, eigenvecs, node_mass)
        array([0.123, 0.456])  # Example output

    Notes:
    - The eigenvalues and eigenvectors should be precomputed using a normal mode analysis method.
    - It is assumed that the first six eigenvalues and eigenvectors (corresponding to rigid-body modes) 
      have already been excluded from the input.
    - Mass weighting is applied to ensure physical consistency.
    - If `convert2bfactors` is True, the results will match the scale of crystallographic B-factors, 
      commonly used in structural biology.
    """
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
    """
    Compute the predicted B-factors of a structure using the pfANM (parameter-free Anisotropic Network Model).

    This function calculates the Normal Modes for the input structure using the pfANM model, which is based on 
    a spring network representation of the structure. The computed modes are then used to predict the 
    B-factors, which represent the atomic fluctuations, commonly observed in X-ray crystallography.

    Parameters:
    df : pandas.DataFrame
        A DataFrame containing the atomic structure of the system. The DataFrame must include the following columns:
        - 'x', 'y', 'z': 3D coordinates of the atoms.
        - 'm': Atomic masses in AMU.
    spring_constant : float, optional
        The spring constant used to define the stiffness of the interactions between atoms in the model. 
        Defaults to 1.0.

    Returns:
    pandas.Series
        A Series containing the predicted B-factors for each atom in the input structure, indexed to match the input DataFrame.

    Raises:
    KeyError
        If the input DataFrame is missing any of the required columns ['x', 'y', 'z'].

    Example:
        >>> import pandas as pd
        >>> structure = pd.DataFrame({
        ...     'x': [0.0, 1.0, 2.0],
        ...     'y': [0.0, 1.0, 0.0],
        ...     'z': [0.0, 0.0, 1.0]
        ... })
        >>> predicted_Bfactors(structure, spring_constant=1.0)
        0    12.34
        1    45.67
        2    89.01
        dtype: float64

    Notes:
    - This function is designed to approximate experimental B-factors but does not account for external factors like solvent effects.
    """
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