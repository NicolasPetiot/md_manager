# MD-manager
md_manager is a set of python functions and methods that provide easy acess and manipulation of data for molecular dynamics simulations. It is developped in the research group of Physics Applied to Proteins at _université de bougrogne_. 

## AtomFrame
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
The ANM Hessian needs a cutoff radius that ensure that only 6 eigenvalues are null from the 3 global translations and 3 global rotations. For large number of nodes this can be problematic as it needs to perform a several diagonalization to tune the value of $r_c$. Introducing the pfANM Hessian allowing to consider interactions between all the nodes but scaling the force constant to the inversed square of the distance. The construction of such matrix is achieved by the function `md.NMA.pfANM_hessian`.

```python
from md_manager import fetch_protein_data_bank
from md_manager.NMA import *

df = fetch_protein_data_bank("8q89")
xyz = ["x", "y", "z"] # for column selection
H = pfANM_hessian(df[xyz], df.m, spring_constant = 1.0)
```
### Normal Modes & Thermal B-factors:
The normal modes are defined by the eigenvectors of the mass-weighted Hessian $\hat{H}\vec{e}_k = \tilde{\omega}^2_k\vec{e}_k$. The 6 firsts modes should not be taken into account as related to global transitions / rotations of the structures. The extraction of the modes is acieved by the function `md.NMA.collective_modes`. 

One of the useful metric that can be extracted from those modes is the Thermal B-factors of the nodes. The extraction of thermal B-factors is achieved with `md.NMA.local_MSF`. Finally the code bellow runs a Bfactor prediction and save it in a pdb file to allow easy vizualization.

```python
import md_manager as md
from md_manager.NMA import *

# Load structure
df = md.fetch_protein_data_bank("8q89")
xyz = ["x", "y", "z"] # for column selection

# Collective modes :
H = pfANM_hessian(df[xyz], df.m, spring_constant = 1.0)
eigenvals, eigenvecs = collective_modes(H)

# B-factors:
bfactors = predicted_Bfactors(eigenvals, eigenvecs, df.m)

# save prediction in pdb file:
df.b = pd.Series(b, index=df.index)
md.df2pdb("8q89.pdb", df=df)
```

Alternatively, it is possible to use the following code:
```python
import md_manager as md
from md_manager.NMA import *

# Load structure
df = md.fetch_protein_data_bank("8q89")
xyz = ["x", "y", "z"] # for column selection

# B-factors:
bfactors = predicted_Bfactors(df=df)

# save prediction in pdb file:
df.b = pd.Series(b, index=df.index)
md.df2pdb("8q89.pdb", df=df)
```

### 1D Graph representation:
Collective modes computed above are extracted from tree dimentional representation of structures but studying topological properties can be achieved with much simpler 1D graphs. (see [Einstein Model of a Graph to Characterize Protein Folded/Unfolded States](https://www.mdpi.com/1420-3049/28/18/6659)). To build such graph, one needs just the position of the set of node to consider and a contact threshold to determine wether or not two nodes are connected by an edge. Such construction is achieved using the `NMA.df2graph` function.

From such Graph, interactions are computed by a Laplacian extracted by `ANM.graph_laplacian` and collective modes are then computed from `ANM.graph_collective_modes`.

```python
import md_manager as md
import md_manager.ANM as anm

df = md.PDBfile("file.pdb").read2df()
graph  = anm.df2graph(df)
laplacian = anm.graph_laplacian(graph)
eigenvals, eigenvecs = anm.graph_collective_modes(laplacian)
```

---
Contact : Nicolas.Petiot01@u-bourgogne.fr