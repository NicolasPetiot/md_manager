import MDAnalysis as mda
import numpy as np
import pandas as pd
from numpy.typing import NDArray

from .core import universe_to_top


def backbone_theta_gamma(u:mda.Universe, return_theta = True, return_gamma = True) -> pd.DataFrame:
    if not return_theta and not return_gamma:
        raise ValueError("return_theta or return_gamma should be set to True")

    expected_cols = {"record_name", "alt", "resn", "chain", "resi", "segi"}

    ca = u.select_atoms("name CA")
    top = universe_to_top(ca)

    cols = [col for col in top.columns if col in expected_cols]
    top[["x", "y", "z"]] = ca.atoms.positions

    if "chain" not in cols:
        top["chain"] = 'A'

    if return_theta:
        top["Theta"] = theta_angles(top)
        cols.append("Theta")
    if return_gamma:
        top["Gamma"] = gamma_angles(top)
        cols.append("Gamma")

    df = pd.DataFrame(top[cols])
    return df

def theta_angles(df:pd.DataFrame) -> pd.Series:
    for _, chain in df.groupby("chain")[["x", "y", "z"]]:
        idx = chain.index
        pos = np.array(chain)
        df.loc[idx[1:-1], "Theta"] = bend_angles(pos)

    return df.Theta

def gamma_angles(df:pd.DataFrame) -> pd.Series:
    for _, chain in df.groupby("chain")[["x", "y", "z"]]:
        idx = chain.index
        pos = np.array(chain)
        df.loc[idx[1:-2], "Gamma"] = dihedral_angles(pos)

    return df.Gamma

def bend_angles(atom_position:NDArray) -> np.ndarray:
    """
    Calculate the bond angles for an ensemble of atomic positions.

    Given a set of atomic positions in a 3D space, this function computes the
    bond angles between successive triplets of atoms, based on the vectors formed
    by consecutive atoms. The angles are returned in degrees.

    The length of the output array is `len(atom_position) - 2`, since the bond angles
    cannot be defined for the first and last atoms of a chain.
    """
    B1 = atom_position[1:-1, :] - atom_position[:-2, :]
    B2 = atom_position[2:, :] - atom_position[1:-1, :]

    B1 = B1 / np.sqrt(np.sum(B1**2, axis = 1, keepdims=True))
    B2 = B2 / np.sqrt(np.sum(B2**2, axis = 1, keepdims=True))

    return np.rad2deg(np.arccos(np.sum(-B1 * B2, axis=1)))

def dihedral_angles(atom_position:np.ndarray) -> np.ndarray:
    """
    Calculate the dihedral angles for an ensemble of atomic positions.

    Given a set of atomic positions in a 3D space, this function computes the
    dihedral angles between consecutive sets of four atoms, based on the vectors
    formed by these atoms. The angles are returned in degrees.

    The length of the output array is `len(atom_position) - 3`, since the dihedral
    angles cannot be defined for the first, last, and second-last atoms of a chain.
    """
    B1 = atom_position[1:-2, :] - atom_position[:-3, :]
    B2 = atom_position[2:-1, :] - atom_position[1:-2, :]
    B3 = atom_position[3:, :] - atom_position[2:-1, :]

    B1 = B1 / np.sqrt(np.sum(B1**2, axis = 1, keepdims=True))
    B2 = B2 / np.sqrt(np.sum(B2**2, axis = 1, keepdims=True))
    B3 = B3 / np.sqrt(np.sum(B3**2, axis = 1, keepdims=True))

    N1 = np.cross(B1, B2)
    N2 = np.cross(B2, B3)

    cos = np.sum(N1 * N2, axis = 1)
    sin = np.sum(np.cross(N1, N2)*B2, axis = 1)

    return np.rad2deg(np.arctan2(sin, cos))
