
--boundary_.oOo._mrbPAYpfKmvEQ0ARb1qIFz8JXaYHJpLN
Content-Length: 3216
Content-Type: application/octet-stream
X-File-MD5: 39cd41fcd9990c83aac2d11b2b83858e
X-File-Mtime: 1709735636
X-File-Path: /md_manager/md_manager/utils/_MD_utils.py

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

--boundary_.oOo._mrbPAYpfKmvEQ0ARb1qIFz8JXaYHJpLN
Content-Length: 7645
Content-Type: application/octet-stream
X-File-MD5: 253b7ff3a49022116a2e35817eaba5ac
X-File-Mtime: 1709735636
X-File-Path: /md_manager/README.md

# MD-manager
md_manager is a set of python functions and methods that provide easy acess and manipulation of data for molecular dynamics simulations. It is developped in the research group of Physics Applied to Proteins at _université de bougrogne_. 

## Installation:

As MD-manager is still in early stage of development, it is not yet available on the [pypi](https://pypi.org) code base. To download it use 
```shell
git clone https://github.com/NicolasPetiot/md_manager
cd md_manager
pip install .
```

## AtomFrame:
The central piece of this project is the manipulation of atoms selection in the form of DataFrames provided by the pandas library. Currently, md_manager can handle pdb file format.

### DataFrame structure:
Unless youre code specifically adds, remove or rename columns in the DataFrames, they contains the following columns :
* `record_name` : Can be in the ATOM/HETATM class.
* `name` : Atom name in PDB file.
* `alt` : Alternate location indicator.
* `resn` : Residue name (coded with 3 letters).
* `chain` : Chain indicator.
* `resi` : Residue index indicator.
* `insertion` : Insertion code.
* `x` : x coordinate in Angström
* `y` : y coordinate in Angström
* `z` : z coordinate in Angström
* `occupancy` : occupancy usually measured in XRD. By default is set to `1.0`.
* `b` : Thermal B-factors in squared Angström. By default is set to `1.0`.
* `segi` : Segment identifier.
* `e` : Element symbol. Used to determine the atomic mass.
* `q` : Atom's charge.
* `m` : Atom's mass in gram per mole.

To generate such DataFrame, you can use a pdb file on youre local machine using 
```python
df = md.pdb2df(filename = "file.pdb")
```
or load structure from the protein data bank base
```python
df = md.fetch_protein_data_bank(pdb_code = "8q89")
```

### DataFrame manipulation:
Using DataFrames allow easier manipulation based on the [query](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html) and [groupby](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.groupby.html) methods.

## PDB object:

MD-manager contains a `PDB` python class that is in charge of the interaction between DataFrames and text files. 

### Reading files:

The `pdb2df` method is meant to return the first model of a pdb file. For MD trajectories, the `PDB` class allow the following syntax

```python
import md_manager as md
traj = md.PDB(filename="md_traj.pdb")
for model in traj:
    # TODO: add data processing here
```

or 

```python
import md_manager as md
traj = md.PDB(filename="md_traj.pdb")
traj = iter(traj)

# initialization :
frame = next(traj)
# TODO: add initialization code here

for frame in traj:
    # TODO: add data processing here
```

**N.B.**: The `PDB.__next__()` method uses the "ENDMDL" record name to return the atoms in the DataFrame. If none of the lines starts with "ENDMDL", the DataFrame will contains all atoms in the input file.

### Writing file:

MD-manager allows to write pdb files as follows:

```python
import md_manager as md

df = md.fetch_protein_data_bank("5f0g")

# remove hetero atoms:
df = df.query("record_name == 'ATOM'")

# save to pdb format:
new = md.PDB(filename="file.pdb", write=True)
new.write(df)
```

Trajectories can be generated as well:

```python
import md_manager as md

traj = md.PDB("full_traj.pdb")
traj = iter(traj)

# save 10 firsts models to pdb format:
new = md.PDB(filename="sub_traj.pdb", write=True)
new.write(dfs=[next(traj) for _ in range(10)]) # dfs is a list of DataFrames
```


## Normal Mode Analysis:
Normal Modes are extracted from the diagonalization of mass-weighted Hessian matricies. 
Choosing a NMA model is equivalent to chose a way to build Hessians.

### Anisotropic Network Model:
In the ANM, Hessians are computed based force constants that are identical for all connected nodes in the equivalent graph representation.
The construction of such matrix is achieved by the function `md.NMA.ANM_hessian`. 

```python
from md_manager import fetch_protein_data_bank
from md_manager.NMA import *

df = fetch_protein_data_bank("8q89")
xyz = ["x", "y", "z"] # for column selection
H = ANM_hessian(df[xyz], df.m, cutoff = 5.0, spring_constant = 1.0)
```
### Parameter-free Anisotropic Network Model :
The ANM Hessian needs a cutoff radius that ensure that only 6 eigenvalues are null from the 3 global translations and 3 global rotations. For large number of nod