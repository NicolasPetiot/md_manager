import MDAnalysis as mda
import numpy as np

from numpy.typing import NDArray

from .constants import BOLTZMANN

def pfANM(ag:mda.AtomGroup, gamma = 1.0) -> NDArray:
    """
    Builds a mass-weighted Hessian matrix based on the parameter-free Anisotropic Network Model.

    Coordinates are extracted from the input AtomGroup and the matrix is scaled by a factor gamma with unit [kJ/mol].
    """
    # Read AtomGroup Data
    Natom = len(ag)
    pos  = ag.atoms.positions # Angstrom
    mass = ag.atoms.masses    # g/mol
    # [gamma] = kJ/mol

    # Distance matrix:
    diff  = pos[:, None, :] - pos[None, :, :] # Shape (N, N, 3)
    dists = np.linalg.norm(diff, axis=2)      # Shape (N, N)

    hessian = diff[:, :, :, None] * diff[:, :, None, :] # Shape (N, N, 3, 3)

    with np.errstate(invalid='ignore'):
        hessian = -gamma * hessian / (dists[:, :, None, None])**4

    # Self-term:
    diag = -np.nan_to_num(hessian).sum(axis = 0)
    hessian[range(Natom), range(Natom)] = diag

    # Mass Weighting
    MiMj = mass[:, None] * mass[None, :]
    hessian = hessian / np.sqrt(MiMj[:, :, None, None])

    # Reshape:
    # Expected shape (3N, 3N) where hessian[0, :] = [x0, y0, z0, x1, y1, z1, (...)]
    hessian = np.permute_dims(hessian, [0, 2, 1, 3]).reshape((3*Natom, 3*Natom))

    return hessian

def normal_modes(hessian:NDArray) -> tuple[NDArray, NDArray]:
    """
    Returns the eigenfrequencies and eigennmodes associated with input hessian matrix.

     - eigenfreqs: shape (Nmode)
     - eigenmodes: shape (Nmode, Natom, 3)

    Note: the six null eigenfrequencies are removed automatically without checking the values.
    Please make sure that the input hessian checks the translation invariance relation.
    """
    # Compute eigenmodes and remove six firsts zeros eigenfrequencies
    eigenfreqs, eigenmodes = np.linalg.eigh(hessian)
    eigenfreqs = eigenfreqs[6:]
    eigenmodes = eigenmodes[:, 6:]
    Nmode = len(eigenfreqs)

    # Reshape eigenmodes:
    eigenmodes = eigenmodes.reshape(-1, 3, Nmode) # Shape (Natom, 3, Nmode)
    eigenmodes = np.permute_dims(eigenmodes, axes=[2, 0, 1]) # Shape (Nmode, Natom, 3)

    return eigenfreqs, eigenmodes

def thermal_factors(eigenfreqs:NDArray, eigenmodes:NDArray, mass:NDArray, temp = 2.5) -> NDArray:
    """
    Implements the relation for the predicted debey thermal factors from normal modes
    """
    debey = (eigenmodes ** 2).sum(axis = 2) # Shape (Nmode, Natom)
    debey /= mass[None, :]
    debey /= (eigenfreqs**2)[:, None]
    debey = np.sum(debey, axis=0) # Sum over modes

    return debey * 8*np.pi**2/3 * BOLTZMANN * temp

def predict_thermal_factors(ag:mda.AtomGroup, model = pfANM, temp = 2.5, **kwargs) -> NDArray:
    """
    Implementation of the thermal factor prediction pipeline:
        - The hessian matrix is built from the input AtomGroup and model (default: pfANM)
        - The modes are computed from diagonalization of the hessian
        - Thermal factors are computed from normal modes
    """

    hessian = model(ag, **kwargs)
    eigenfreqs, eigenmodes = normal_modes(hessian)
    return thermal_factors(eigenfreqs, eigenmodes, ag.atoms.masses, temp)
