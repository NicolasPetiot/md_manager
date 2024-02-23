from ._PDBfile_utils import _build_df_from_atom_list, _scan_pdb_line
from .._PDBfile import PDBfile

import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial import distance_matrix

import sys
from urllib.request import urlopen



### IMPORTED IN md_manager ###
__all__ = [
    "fetch_protein_data_bank",
    "pdb2df",
    "df2pdb",
    "df2com",
    "pdb2graph",
    "df2graph"
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

def pdb2df(filename:str, atom_only = False) -> pd.DataFrame:
    """
    Read all the atoms in an input PDB file [filename] and returns the associated DataFrame.
    """
    # reading :
    pdb = PDBfile(filename)
    pdb.open()
    df = pdb.read2df()
    pdb.close()

    # query :
    if atom_only:
        df = df.query("record_name == 'ATOM'")

    return df

def df2pdb(filename:str, df:pd.DataFrame=None, dfs:list[pd.DataFrame] = None, title:str= None):
    """
    Read all the atoms in an input DataFrame [df / dfs] and generates the associated PDB file.
    """
    NoneType = type(None)
    if type(df) == NoneType and type(dfs) == NoneType:
        raise ValueError("Please enter an input df/dfs")
    
    
    pdb = PDBfile(filename, write=True)
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

def pdb2graph(filename:str, contact_threshold:float = 4.0) -> nx.Graph:
    """
    
    """
    df = pdb2df(filename)
    return df2graph(df, contact_threshold)

def df2graph(df:pd.DataFrame, contact_threshold = 4.0) -> nx.Graph:
    """
    
    """
    xyz = ["x", "y", "z"]
    distance_inter_nodes = distance_matrix(df[xyz], df[xyz])
    I, J = np.where(distance_inter_nodes < contact_threshold) 
    I = df.iloc[I].index # atom index
    J = df.iloc[J].index # atom index
    edges = [(i, j) for i, j in zip(I, J) if i < j]

    G = nx.Graph()
    G.add_edges_from(edges)
    return G
