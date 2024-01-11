import pandas as pd
import sys
from urllib.request import urlopen
from ._PDBfile_utils import _build_df_from_atom_list, _scan_pdb_line
from .._PDBfile import PDBfile


### IMPORTED IN md_manager ###
__all__ = [
    "fetch_protein_data_bank",
    "pdb2df",
    "df2pdb",
    "df2com"
]

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

def pdb2df(filename:str) -> pd.DataFrame:
    """
    Read all the atoms in an input PDB file [filename] and returns the associated DataFrame.
    """
    with PDBfile(filename) as pdb:
        df = pdb.read2df()
    return df

def df2pdb(filename:str, df:pd.DataFrame=None, dfs:list[pd.DataFrame] = None, title:str= None):
    """
    Read all the atoms in an input DataFrame [df / dfs] and generates the associated PDB file.
    """
    NoneType = type(None)
    if type(df) == NoneType and type(dfs) == NoneType:
        raise ValueError("Please enter an input df/dfs")
    
    with PDBfile(filename, write=True) as pdb:
        if type(dfs) != NoneType:
            pdb.write_pdb(dfs=dfs)

        else :
            pdb.write_pdb(df=df)

def df2com(df:pd.DataFrame, name = " XX ", element = " X", charge = "  "):
    """
    Read all the amino-acids of an input DataFrame and generates an other DataFrame that contains their center of mass.

    This function is using the `df.groupby(groups)` method. 
    """
    # grouping by amino-acids :
    groups = ["record_name", "chain", "resi", "resn", "alt", "segi"]
    amino_acid_iterator = df.groupby(groups)

    # initialize df
    COM = []

    for (record_name, chain, resi, resn, alt, segi), aa in amino_acid_iterator:
        M = aa.m.sum()
        X = (aa.x * aa.m).sum() / M
        Y = (aa.y * aa.m).sum() / M
        Z = (aa.z * aa.m).sum() / M    

        B = aa.b.mean()

        COM.append(
            [record_name, name, alt, resn, chain, resi, " ", X, Y, Z, 1.0, B, segi, element, charge, M]
        )

    return _build_df_from_atom_list(atoms=COM)
