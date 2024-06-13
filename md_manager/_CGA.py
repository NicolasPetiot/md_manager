import numpy as np
import pandas as pd

from ._params import ATOM_NAME_SELECTION_CHI

__all__ = [
    "angles",
    "dihedral_angles",
    "chain_theta_angles",
    "chain_gamma_angles",
    "residue_chi_angles"
]

def chain_theta_angles(chain:pd.DataFrame, CA_only = True) -> pd.Series:
    """
    Returns a pandas.Series that contains the computed Theta angles of the input chain. The series index is the `resi` number 
    and the associated value is the computed theta angles formed by the surounding "CA" atoms.

    Can be used as `frame.groupby(...)[["name", "resi", "x", "y", "z"]].apply(chain_theta_angles)` for multimers.

    By default, the method considers that the input chain contains only 'CA' atoms. To apply query methods, please set `CA_only=False`.

    CAUTION : This methods will NOT returns an error if two residue index are not consecutive. Please make sure to avoid missing residues in the `chain` DataFrame.
    """
    xyz = ["x", "y", "z"]

    if not CA_only:
        atm_selection = ["CA"]
        chain = chain.query("name in @atm_selection")

    thetas = angles(chain[xyz].to_numpy())

    s = pd.Series(index = chain.resi, name = "Theta", dtype = float)
    s[s.index[1:-1]] = thetas
    return s

def chain_gamma_angles(chain:pd.DataFrame, CA_only = True) -> pd.Series:
    """
    Returns a pandas.Series that contains the computed Gamma angles of the input chain. The series index is the `resi` number 
    and the associated value is the computed gamma angles formed by the surounding "CA" atoms.

    Can be used as `frame.groupby(...)[["name", "resi", "x", "y", "z"]].apply(chain_gamma_angles)` for multimers.

    By default, the method considers that the input chain contains only 'CA' atoms. To apply query methods, please set `CA_only=False`.

    CAUTION : This methods will NOT returns an error if two residue index are not consecutive. Please make sure to avoid missing residues in the `chain` DataFrame.
    """
    xyz = ["x", "y", "z"]

    if not CA_only:
        atm_selection = ["CA"]
        chain = chain.query("name in @atm_selection")
        
    gammas = dihedral_angles(chain[xyz].to_numpy())

    s = pd.Series(index = chain.resi, name = "Gamma", dtype = float)
    s[s.index[1:-2]] = gammas
    return s

def residue_chi_angles(res:pd.DataFrame) -> pd.Series:
    """
    Returns a pandas.Series that contains the computed Chi angles of the input residue. The series index in the 'ChiN' identifier 
    and the associated value is the computed chi angle formed by the surounding atoms.

    Can be used as `frame.groupby(...)[["name", "resn", "x", "y", "z"]].apply(residue_chi_angles)`.
    """
    xyz = ["x", "y", "z"]

    atm_selection = ATOM_NAME_SELECTION_CHI[res.resn.unique()[0]]
    res = res.query("name in @atm_selection")
    chis = dihedral_angles(res[xyz].to_numpy())
    chis = [chi for chi in chis] + [np.nan for _ in range(5 - len(chis))]
    
    return pd.Series(chis, index=["Chi%d"%i for i in range(1, 6)])

def angles(atom_position:np.ndarray) -> np.ndarray:
    """
    Returns the angles associated to the ensemble of positions in 'atom_position'.
    The length of the output is len(atom_position) - 2 because angles are not defined for the first and last atom of a chain.
    """
    B1 = atom_position[1:-1, :] - atom_position[:-2, :] 
    B2 = atom_position[2:, :] - atom_position[1:-1, :] 

    B1 = B1 / np.sqrt(np.sum(B1**2, axis = 1, keepdims=True))
    B2 = B2 / np.sqrt(np.sum(B2**2, axis = 1, keepdims=True))

    return np.rad2deg(np.arccos(np.sum(-B1 * B2, axis=1)))

def dihedral_angles(atom_position:np.ndarray) -> np.ndarray:
    """
    Returns the dihedral angles associated to the ensemble of positions in 'atom_position'.
    The length of the output is len(atom_position) - 3 because dihedral angles are not defined for the first, last and second last atom of a chain.
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