from .._PDB import PDB
from ._params import DF_COLUMNS

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
            atoms.append(PDB.scan_pdb_line(line))

    return PDB.build_df_from_atoms(atoms)

def pdb2df(filename:str, atom_only = False) -> pd.DataFrame:
    """
    Read all the atoms in the first frame of an input PDB file [filename] and returns the associated DataFrame.
    """
    # reading :
    pdb   = PDB(filename)
    frame = iter(pdb)
    df    = next(frame)
    pdb.close()

    # query :
    if atom_only:
        query_string = f"{DF_COLUMNS[0]} == 'ATOM'"
        df = df.query(query_string)

    return df

def df2pdb(filename:str, df:pd.DataFrame=None, dfs:list[pd.DataFrame] = None):
    """
    Read all the atoms in an input DataFrame [df / dfs] and generates the associated PDB file.
    """
    new = PDB(filename, write=True)
    new.write(model=df, model_list=dfs)
   

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

    return PDB.build_df_from_atom_list(atoms=COM)

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
