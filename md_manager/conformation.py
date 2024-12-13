from .parameters import ATOM_NAME_CHI

import numpy as np
import pandas as pd

__all__ = [
    "angles",
    "dihedral_angles",
    "backbone_conformation",
    "side_chain_conformation",
    #"chain_phi_angles",
    #"chain_psi_angles",
    "chain_theta_angles",
    "chain_gamma_angles",
    "residue_chi_angles"
]

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

def backbone_conformation(df:pd.DataFrame) -> pd.DataFrame:
    """
    Returns a DataFrame containing the conformational angles Theta and Gamma. 

    Will raise an KeyError if chain is missing one of the ['x', 'y', 'z'] columns. 
    """
    CA = df.query("name == 'CA'").copy()
    CA.index = CA.resi
    is_multimer = "chain" in df.columns and len(df["chain"].unique()) > 1

    if is_multimer:
        groups = CA.groupby("chain")[["x", "y", "z"]]
        return groups.apply(lambda df:
            pd.DataFrame(dict(
                theta = chain_theta_angles(df),
                gamma = chain_gamma_angles(df)
            ))                              
        )

    else:
        return pd.DataFrame(dict(
                theta = chain_theta_angles(CA),
                gamma = chain_gamma_angles(CA)
            ))

def side_chain_conformation(df:pd.DataFrame) -> pd.DataFrame:
    """
    
    """
    # Select atoms and order them.
    # Especially important because of prolines.
    sele = df.query("name in @ATOM_NAME_CHI").copy()
    sele["name"] = pd.Categorical(sele["name"], categories=ATOM_NAME_CHI, ordered=True)

    is_multimer = "chain" in df.columns and len(df["chain"].unique()) > 1
    if is_multimer:
        sele = sele.sort_values(by=["chain", "resi", "name"])
        return sele.groupby(["chain", "resi"])[["x", "y", "z"]].apply(residue_chi_angles)

    else:
        sele = sele.sort_values(by=["resi", "name"])
        return sele.groupby(["resi"])[["x", "y", "z"]].apply(residue_chi_angles)

def chain_phi_angles(chain:pd.DataFrame) -> pd.Series:
    raise NotImplemented

def chain_psi_angles(chain:pd.DataFrame) -> pd.Series:
    raise NotImplemented

def chain_theta_angles(chain:pd.DataFrame) -> pd.Series:
    """
    Returns a Series that contains the computed Theta angles of the input chain.

    Will raise an KeyError if chain is missing one of the ['x', 'y', 'z'] columns. 

    Caution: Will not check for missing/extra atoms.
    """

    theta = pd.Series(index=chain.index, dtype=float, name="Theta")
    theta.iloc[1:-1] = angles(chain[["x", "y", "z"]].values)

    return theta

def chain_gamma_angles(chain:pd.DataFrame) -> pd.Series:
    """
    Returns a Series that contains the computed Gamma angles of the input chain. 

    Will raise an KeyError if chain is missing one of the ['x', 'y', 'z'] columns.

    Caution: Will not check for missing/extra atoms.
    """

    gamma = pd.Series(index=chain.index, dtype=float, name="Gamma")
    gamma.iloc[1:-2] = dihedral_angles(chain[["x", "y", "z"]].values)

    return gamma

def residue_chi_angles(res:pd.DataFrame) -> pd.Series:
    """
    Returns a Series that contains the computed Chi angles of the input residue. 

    Will raise an KeyError if chain is missing one of the ['x', 'y', 'z'] columns.

    Caution: Will not check for missing/extra atoms.
    """
    chi = np.array([np.nan for _ in range(5)])
    tmp = dihedral_angles(res[["x", "y", "z"]].values)
    chi[:len(tmp)] = tmp

    return pd.Series(chi, index=["Chi%d"%(i+1) for i in range(5)])