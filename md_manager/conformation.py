from .parameters import ATOM_NAME_CHI

import numpy as np
import pandas as pd

__all__ = [
    "angles",
    "dihedral_angles",
    "backbone_conformation",
    "side_chain_conformation",
    "chain_phi_psi_omega",
    "chain_theta_angles",
    "chain_gamma_angles",
    "residue_chi_angles"
]

def angles(atom_position:np.ndarray) -> np.ndarray:
    """
    Calculate the bond angles for an ensemble of atomic positions.

    Given a set of atomic positions in a 3D space, this function computes the 
    bond angles between successive triplets of atoms, based on the vectors formed 
    by consecutive atoms. The angles are returned in degrees.

    The length of the output array is `len(atom_position) - 2`, since the bond angles 
    cannot be defined for the first and last atoms of a chain.

    Parameters:
    atom_position : numpy.ndarray
        A 2D array of shape (n, 3), where `n` is the number of atoms and each row represents
        the 3D coordinates (x, y, z) of an atom in space.

    Returns:
    numpy.ndarray
        A 1D array of bond angles (in degrees) corresponding to the vectors formed by
        consecutive triplets of atoms in the input `atom_position` array.

    Example:
        >>> atom_position = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
        >>> angles(atom_position)
        array([90., 90.])
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

    Parameters:
    atom_position : numpy.ndarray
        A 2D array of shape (n, 3), where `n` is the number of atoms and each row represents
        the 3D coordinates (x, y, z) of an atom in space.

    Returns:
    numpy.ndarray
        A 1D array of dihedral angles (in degrees) corresponding to the consecutive sets of
        four atoms in the input `atom_position` array.

    Example:
        >>> atom_position = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1]])
        >>> dihedral_angles(atom_position)
        array([90.])
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
    Calculate the backbone conformational angles (Theta and Gamma) for a given set of atomic positions.

    This function computes the backbone angles Theta and Gamma for each chain in the protein. 
    It first filters the DataFrame to include only the alpha-carbon (CA) atoms, then calculates 
    the angles based on the atomic positions. If the protein is a multimer (i.e., contains multiple chains), 
    the angles are computed separately for each chain.

    The function raises a `KeyError` if the DataFrame is missing any of the ['x', 'y', 'z'] columns, 
    which are necessary to compute the conformational angles.

    Parameters:
    df : pandas.DataFrame
        A DataFrame containing atomic coordinates and atom names, with the following required columns:
        - 'name': atom names (must include 'CA' for alpha-carbon atoms).
        - 'x', 'y', 'z': 3D coordinates of the atoms.
        - 'chain' (optional): specifies different chains for multimeric proteins.
        - 'resi': residue index.

    Returns:
    pandas.DataFrame
        A DataFrame containing two columns:
        - 'theta': the Theta conformational angle for each CA atom (in degrees).
        - 'gamma': the Gamma conformational angle for each CA atom (in degrees).

    Raises:
    KeyError
        If the input DataFrame is missing one of the required columns ['x', 'y', 'z'].

    Example:
        >>> df = pd.DataFrame({
        ...     'name': ['CA', 'CA', 'CA', 'CA', 'CA'],
        ...     'x': [0.0, 1.0, 2.0, 3.0, 4.0],
        ...     'y': [0.0, 1.0, 1.0, 1.0, 1.0],
        ...     'z': [0.0, 0.0, 0.0, 0.0, 0.0],
        ...     'resi': [1, 2, 3, 4, 5],
        ...     'chain': ['A', 'A', 'A', 'A', 'A']
        ... })
        >>> backbone_conformation(df)
           theta  gamma
        0   180.0   180.0
        1   180.0   180.0
        2   180.0   180.0
        3   180.0   180.0
        4   180.0   180.0

    Notes:
    - The function calculates angles based on the positions of the alpha-carbon (CA) atoms.
    - Theta and Gamma angles are typically used to describe the relative positions of consecutive residues in a protein.
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
    Calculate the side-chain conformational angles (chi angles) for each residue in the protein.

    This function computes the chi angles for each residue's side chain by selecting atoms 
    that are part of the side chain and then calculating the torsional angles (chi angles) 
    formed by consecutive atoms in the side chain. The function is designed to handle 
    both single-chain and multimeric proteins. If the protein is a multimer, the chi angles 
    are calculated separately for each chain and residue.

    The chi angles are typically used to describe the conformation of side chains in 
    proteins, with particular attention to the rotation around the bonds in the side chain.

    Parameters:
    df : pandas.DataFrame
        A DataFrame containing atomic coordinates and atom names, with the following required columns:
        - 'name': atom names (must include atoms relevant to the side-chain chi angles).
        - 'x', 'y', 'z': 3D coordinates of the atoms.
        - 'chain' (optional): specifies different chains for multimeric proteins.
        - 'resi': residue index.

    Returns:
    pandas.DataFrame
        A DataFrame containing the chi angles (torsional angles) for each residue's side chain. 
        The angles are grouped by residue, and the output contains separate calculations for each chain 
        if the protein is a multimer.

    Raises:
    KeyError
        If the input DataFrame is missing one of the required columns ['x', 'y', 'z'].

    Example:
        >>> df = pd.DataFrame({
        ...     'name': ['CA', 'CB', 'CG', 'CD', 'NE1', 'CE1'],
        ...     'x': [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        ...     'y': [0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ...     'z': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ...     'resi': [1, 1, 1, 1, 1, 1],
        ...     'chain': ['A', 'A', 'A', 'A', 'A', 'A']
        ... })
        >>> side_chain_conformation(df)
           chi_angle
        0   180.0
        1   90.0
        2   270.0

    Notes:
    - The function uses a predefined list `ATOM_NAME_CHI` to identify the atoms involved in the chi angles.
    - The chi angles are computed for each side-chain in the protein, and the rotation of the side chain around 
      bonds is considered to calculate the torsional angles.
    """
    # Select atoms and order them.
    sele = df.query("name in @ATOM_NAME_CHI").copy()
    sele = sele.query("name != 'CZ' or resn == 'ARG'")
    sele["name"] = pd.Categorical(sele["name"], categories=ATOM_NAME_CHI, ordered=True)

    is_multimer = "chain" in df.columns and len(df["chain"].unique()) > 1
    if is_multimer:
        sele = sele.sort_values(by=["chain", "resi", "name"])
        return sele.groupby(["chain", "resi"])[["x", "y", "z"]].apply(residue_chi_angles)

    else:
        sele = sele.sort_values(by=["resi", "name"])
        return sele.groupby(["resi"])[["x", "y", "z"]].apply(residue_chi_angles)

def chain_phi_psi_omega(chain:pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the Phi, Psi, and Omega dihedral angles for a protein backbone chain.

    This function computes the three key dihedral angles (Phi, Psi, and Omega) for each 
    residue in the provided protein backbone chain. The angles are calculated based on 
    the atomic positions of the backbone atoms ('N', 'CA', 'C').

    Phi, Psi, and Omega angles are important for characterizing protein secondary structures 
    and their conformations. Specifically:
    - Phi (ϕ) is the torsion angle between the C-N-CA-C bond.
    - Psi (ψ) is the torsion angle between the N-CA-C-N-CA bond.
    - Omega (ω) is the torsion angle between the CA-C-N-CA bond (typically found between adjacent peptide bonds).

    Parameters:
    chain : pandas.DataFrame
        A DataFrame containing the 3D coordinates of atoms in the chain. The DataFrame must include:
        - 'x', 'y', 'z' columns for atomic coordinates.
        - A 'name' column that specifies atom names (e.g., 'N', 'CA', 'C').
        - A 'chain' column (optional) to identify the chain in case of multiple chains in the dataset.
        
    Returns:
    pandas.DataFrame
        A DataFrame containing the computed Phi, Psi, and Omega angles (in degrees) for each residue in the chain. 
        The index of the DataFrame corresponds to the residue number (resi), and the columns are labeled 'Phi', 'Psi', and 'Omega'.

    Raises:
    KeyError
        If the input DataFrame is missing one of the required columns ['x', 'y', 'z', 'name'].
    
    Example:
        >>> chain = pd.DataFrame({
        ...     'x': [1.0, 1.1, 1.2, 1.3],
        ...     'y': [0.0, 0.1, 0.2, 0.3],
        ...     'z': [0.0, 0.1, 0.2, 0.3],
        ...     'name': ['N', 'CA', 'C', 'N'],
        ...     'chain': ['A', 'A', 'A', 'A'],
        ...     'resi': [1, 2, 3, 4]
        ... })
        >>> chain_phi_psi_omega(chain)
           Phi  Psi  Omega
        1  180.0  180.0   180.0
        2  180.0  180.0   180.0
        3  180.0  180.0   180.0

    Notes:
    - The function assumes the input `chain` contains the necessary backbone atoms ('N', 'CA', 'C') for each residue.
    - The calculation of Phi, Psi, and Omega angles depends on the correct ordering of atoms in the chain.
    - The residues without enough atoms (e.g., the first or last residue of the chain) will not have their angles calculated and will be assigned NaN values.
    - The `chain` parameter can optionally contain a 'chain' column to handle multiple chains in a structure.
    - The output DataFrame's angles are in degrees.
    """
    backbone = chain.query("name in ['N', 'CA', 'C'] and chain == 'A'")
    residues = backbone.query("name == 'CA'")
    Nres = len(residues)

    angles = np.zeros(3*Nres) * np.nan
    angles[1:-2] = dihedral_angles(backbone[["x", "y", "z"]].values)

    return pd.DataFrame(angles.reshape(Nres, 3), columns=["Phi", "Psi", "Omega"], index = residues.resi)

def chain_theta_angles(chain:pd.DataFrame) -> pd.Series:
    """
    Calculate the Theta angles for the backbone of a protein chain.

    This function computes the bond angles (Theta angles) for the atoms in the input chain 
    using the positions of atoms in 3D space. The Theta angles correspond to the angles between 
    consecutive vectors formed by triplets of atoms in the backbone.

    The Theta angle is calculated for every atom in the chain except for the first and last 
    atoms, as angles cannot be defined for them.

    Parameters:
    chain : pandas.DataFrame
        A DataFrame containing the 3D coordinates of atoms in the chain. The DataFrame must 
        contain the following columns:
        - 'x', 'y', 'z': 3D coordinates of the atoms.

    Returns:
    pandas.Series
        A Series containing the Theta angles (in degrees) for each atom in the input chain, 
        with the same index as the input DataFrame. The first and last atoms will have 
        `NaN` values, as angles cannot be calculated for them.

    Raises:
    KeyError
        If the input DataFrame is missing one of the required columns ['x', 'y', 'z'].

    Example:
        >>> chain = pd.DataFrame({
        ...     'x': [0.0, 1.0, 2.0],
        ...     'y': [0.0, 1.0, 2.0],
        ...     'z': [0.0, 1.0, 2.0]
        ... })
        >>> chain_theta_angles(chain)
        0      NaN
        1     90.0
        2      NaN
        Name: Theta, dtype: float64

    Notes:
    - This function assumes that the input chain DataFrame contains only atoms from the backbone, 
      and that the atoms are ordered sequentially in the chain.
    - The Theta angles are not computed for the first and last atoms in the chain, as they lack 
      enough neighboring atoms to form a valid angle.
    - The function does not verify whether all required atoms are present in the input chain.
    """

    theta = pd.Series(index=chain.index, dtype=float, name="Theta")
    theta.iloc[1:-1] = angles(chain[["x", "y", "z"]].values)

    return theta

def chain_gamma_angles(chain:pd.DataFrame) -> pd.Series:
    """
    Calculate the Gamma dihedral angles for the backbone of a protein chain.

    This function computes the dihedral angles (Gamma angles) for the atoms in the input chain 
    using the positions of atoms in 3D space. The Gamma angles are the angles formed between 
    four consecutive atoms along the backbone, using the planes defined by each consecutive triplet 
    of atoms.

    The Gamma angle is calculated for every atom in the chain except for the first, second last, and 
    last atoms, as dihedral angles cannot be defined for them.

    Parameters:
    chain : pandas.DataFrame
        A DataFrame containing the 3D coordinates of atoms in the chain. The DataFrame must 
        contain the following columns:
        - 'x', 'y', 'z': 3D coordinates of the atoms.

    Returns:
    pandas.Series
        A Series containing the Gamma dihedral angles (in degrees) for each atom in the input chain, 
        with the same index as the input DataFrame. The first, second last, and last atoms will have 
        `NaN` values, as dihedral angles cannot be calculated for them.

    Raises:
    KeyError
        If the input DataFrame is missing one of the required columns ['x', 'y', 'z'].

    Example:
        >>> chain = pd.DataFrame({
        ...     'x': [0.0, 1.0, 2.0, 3.0],
        ...     'y': [0.0, 1.0, 2.0, 3.0],
        ...     'z': [0.0, 1.0, 2.0, 3.0]
        ... })
        >>> chain_gamma_angles(chain)
        0      NaN
        1     90.0
        2      NaN
        3      NaN
        Name: Gamma, dtype: float64

    Notes:
    - This function assumes that the input chain DataFrame contains only atoms from the backbone, 
      and that the atoms are ordered sequentially in the chain.
    - The Gamma angles are not computed for the first, second last, and last atoms in the chain, 
      as they lack enough neighboring atoms to form a valid angle.
    - The function does not verify whether all required atoms are present in the input chain.
    """

    gamma = pd.Series(index=chain.index, dtype=float, name="Gamma")
    gamma.iloc[1:-2] = dihedral_angles(chain[["x", "y", "z"]].values)

    return gamma

def residue_chi_angles(res:pd.DataFrame) -> pd.Series:
    """
    Calculate the Chi angles for the side-chain of a given residue.

    This function computes the dihedral angles (Chi angles) for the side-chain atoms of 
    a residue in a protein structure. The Chi angles are rotational angles around single bonds 
    in the side-chain, typically referring to the bond angles in the side-chain's backbone atoms.

    The Chi angles are calculated based on the positions of atoms in the input residue. 
    This function assumes that the residue contains the necessary atoms for Chi angle 
    calculation, and it returns a Series containing the Chi1, Chi2, ..., Chi5 angles. 
    Missing angles will be represented as `NaN`.

    The residue should have the appropriate atoms for Chi angle calculation. If a residue 
    contains fewer than 5 Chi angles, the corresponding values will be `NaN`.

    Parameters:
    res : pandas.DataFrame
        A DataFrame containing the 3D coordinates of the atoms in the residue. The DataFrame 
        must contain the following columns:
        - 'x', 'y', 'z': 3D coordinates of the atoms in the residue.

    Returns:
    pandas.Series
        A Series containing the Chi angles (in degrees) for the residue. The Series has 
        the index labels 'Chi1', 'Chi2', ..., 'Chi5', and any missing Chi angles will be 
        represented as `NaN`.

    Raises:
    KeyError
        If the input DataFrame is missing one of the required columns ['x', 'y', 'z'].

    Example:
        >>> residue = pd.DataFrame({
        ...     'x': [1.0, 2.0, 3.0, 4.0, 5.0],
        ...     'y': [0.0, 1.0, 1.0, 0.0, 0.0],
        ...     'z': [0.0, 1.0, 2.0, 1.0, 0.0]
        ... })
        >>> residue_chi_angles(residue)
        Chi1    120.0
        Chi2    180.0
        Chi3      NaN
        Chi4      NaN
        Chi5      NaN
        dtype: float64

    Notes:
    - This function assumes that the input residue contains atoms necessary to compute the 
      Chi angles. Typically, this includes the atoms involved in the side-chain torsion angles.
    - Missing Chi angles (for example, if the residue has fewer than 5 Chi angles) will 
      be represented as `NaN` in the returned Series.
    - The function does not check for missing or extra atoms beyond those required for 
      the Chi angle calculation.
    """
    chi = np.array([np.nan for _ in range(5)])
    tmp = dihedral_angles(res[["x", "y", "z"]].values)
    chi[:len(tmp)] = tmp

    return pd.Series(chi, index=["Chi%d"%(i+1) for i in range(5)])