from ._utils import _build_df_from_atom_list, _scan_pdb_line, atom_position, _write_pdb_atom_lines

import pandas as pd
import numpy as np

import sys
from urllib.request import urlopen


### Load from Protein Data Bank ###

def fetch_protein_data_bank(pdb_code:str)->pd.DataFrame:
    """
    Returns a pandas.DataFrame associated to the structure loaded from the 'rcsb.org' website.
    """
    url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
    response = urlopen(url)

    txt = response.read()
    lines = (txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.decode("ascii")).splitlines()

    atoms = []
    for line in lines :
        if line[:6] in {"ATOM  ", "HETATM"}:
            atoms.append(_scan_pdb_line(line))

    return _build_df_from_atom_list(atoms)


### Conformational angles computations ###

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

def theta_angles(df:pd.DataFrame):
    """
    Returns a pandas.Series that contains theta angles associated to the chains in inputed df.
    """
    theta = []
    for _, chain in df.query("name == ' CA '").groupby(["chain"]):
        theta.append(np.nan)
        xyz = chain[["x", "y", "z"]].to_numpy(dtype=float)
        theta += list(angles(xyz))
        theta.append(np.nan)
    return pd.Series(theta, index=df.index, name="theta")

def gamma_angles(df:pd.DataFrame):
    """
    Returns a pandas.Series that contains the gamma angles associated to the chains in inputed df.
    """
    gamma = []
    for _, chain in df.query("name == ' CA '").groupby(["chain"]):
        gamma.append(np.nan)
        xyz = atom_position(chain)
        gamma += list(dihedral_angles(xyz))
        gamma.extend([np.nan, np.nan])
    return pd.Series(gamma, index=df.index, name="gamma")

def chi_angles(df:pd.DataFrame):
    """
    Returns a dictionary that contains the chi angles associated to the residues in inputed df.
    """
    Chis = {"%d"%i : [] for i in range(1, 6)}
    for _, residue in df.groupby(["chain", "resi"]):
        # select only atoms for the computation of chi angles
        resn = residue["resn"].unique()[0]
        residue = residue.query("name in @ATOM_NAME_SELECTION[@resn]") 

        # chi angle computation :
        xyz = atom_position(residue)
        chis = dihedral_angles(xyz)
        for i, chi in enumerate(chis):
            i += 1
            Chis["%d"%i].append(chi)
        for i in range(len(chis) + 1, 6):
            Chis["%d"%i].append(np.nan)

    return Chis

def df2pdb(df:pd.DataFrame, filename: str, title = "", remark = "") -> None:
    """
    Generates a pdb file from an inputed DataFrame.
    """
    with open(filename, "w") as file :
        if title != "":
            file.write("TITLE     " + title + "\n")
        if remark != "":
            for line in remark.splitlines():
                file.write("REMARK    " + line + "\n")
            
        file.write("MODEL        1\n")
        _write_pdb_atom_lines(file=file, df=df)
        file.write("ENDMDL\n")

def dfs2traj(dfs:list[pd.DataFrame], filename:str, title = "", remark = ""):
    """
    Generates a pdb file containing all the inputed DataFrames in the form of a trajectory file.
    """
    with open(filename, "w") as file :
        if title != "":
            file.write("TITLE     " + title + "\n")
        if remark != "":
            for line in remark.splitlines():
                file.write("REMARK    " + line + "\n")
            
        for model, df in enumerate(dfs) :
            file.write(f"MODEL    {model+1:5d}\n")
            _write_pdb_atom_lines(file=file, df=df)
            file.write("ENDMDL\n")