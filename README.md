# MD-manager
md_manager is a set of python functions and methods that provide easy acess and manipulation of data for molecular dynamics simulations. It is developped in the research group of Physics Applied to Proteins in the université de bougrogne. 

## AtomFrame
The central piece of this project is the manipulation of atoms selection in the form of DataFrames provided by the pandas library. Currently, md_manager can handle pdb file format.

### PDBfile methods:
```python
from md_manager import PDBfile
```

* `read2df` : Returns a pandas.DataFrame that contains all the atoms found in the PDB file.
```python
pdb = PDBfile(filename = "file.pdb")
df  = pdb.read2df()
```

* `df2pdb` : Generates a PDB file from an inputed DataFrame.
```python
PDBfile.df2pdb(filename = "filename.pdb", df=df)
```

* `open` & `close` : Open/Close the pdb.file wrapper.
```python
pdb.open()
for line in pdb.file:
    #TODO
pdb.close()
```

* `__len__`: Returns the number of frames in trajectory.
```python
model_number = len(pdb)
``` 

* Iterate over models using `__iter__` & `__next__`: Returns a pandas.DataFrame that corresponds to the next MODEL in PDB file.
```python
pdb.open()
for id, df in enumerate(pdb):
    # 'df' contains the atoms in MODEL n°id
    #TODO
pdb.close()
```

* With statments using `__enter__` & `__exit__` :
```python
with PDBfile("file.pdb") as pdb :
    #pdb.file is open
    #TODO
#pdb.file is close
```

### Other functions related to DataFrames:
```python
import md_manager as md
```
* `atom_position` : Returns a numpy.ndarray that contains the [x, y, z] coordinates of the atoms in df.
* `atom_mass` : Returns a numpy.ndarray that contains the masses of the atoms in df.
* `atom_bfactors` : Returns a numpy.ndarray that contains the thermal B-factors of the atoms in df.
```python
position_array = md.atom_position(df)
mass_array = md.atom_mass(df)
bfact_array = md.atom_bfactors(df)
```

### DataFrame structure:
Unless youre code specifically adds or remove columns in the DataFrames, the contains the following columns :
* `record_name` : Is the atom in the ATOM/HETATM class.
* `name` : Atom name in PDB file (caution: it's length is 4).
* `alt` : Alternate location indicator.
* `resn` : Residue name (coded with 3 letters).
* `chain` : Chain indicator.
* `resi` : Residue index indicator.
* `insertion` : Insertion code.
* `x` : x coordinate in Angström
* `y` : y coordinate in Angström
* `z` : z coordinate in Angström
* `occupancy` : occupancy usually measured in XRD. By default is set to `1.0`.
* `b` : Thermal B-factors in squared Angström.
* `segi` : Segment identifier.
* `e` : Element symbol. Used to determine the atomic mass.
* `q` : Atom's charge.
* `m` : Atom's mass in gram per mole.

## Normal Mode Analysis:
Normal Modes are extracted from the diagonalization of mass-weighted Hessian matricies. 
Choosing a NMA model is equivalent to chose a way to build Hessians.

### Anisotropic Network Model:
In the ANM, Hessians are computed based force constants that are identical for all connected nodes in the equivalent graph representation.
The construction of such matrix is achieved by the function `md.ANM_hessian`. For code efficiancy, this version of `md_manager` needs the `scipy.spatial.distance_matrix` function to compute the distance inter-nodes.

```python
import md_manager as md
import md_manager.NMA as nma
from scipy.spatial import distance_matrix

df = md.PDBfile("file.pdb").read2df()
nodes_position = df[["x", "y", "z"]].to_numpy(dtype = float)
nodes_mass = df.m.to_numpy(dtype = float)
distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

H = nma.ANM_hessian(nodes_position=nodes_position, nodes_mass=nodes_mass, distance_inter_nodes=distance_inter_nodes, cutoff_radius=4.0, spring_constant=1.0)
```
### Parameter-free Anisotropic Network Model :
The ANM Hessian needs a cutoff radius that ensure that only 6 eigenvalues are null from the 3 global translations and 3 global rotations. For large number of nodes this can be problematic as it needs to perform a several diagonalization to tune the value of $r_c$. Introducing the pfANM Hessian allowing to consider interactions between all the nodes but scaling the force constant to the inversed square of the distance. The construction of such matrix is achieved by the function `NMA.pfANM_hessian`.

```python
import md_manager as md
import md_manager.NMA as nma
from scipy.spatial import distance_matrix

df = md.PDBfile("file.pdb").read2df()
nodes_position = df[["x", "y", "z"]].to_numpy(dtype = float)
nodes_mass = df.m.to_numpy(dtype = float)
distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

H = nma.pfANM_hessian(nodes_position=nodes_position, nodes_mass=nodes_mass, distance_inter_nodes=distance_inter_nodes, spring_constant=1.0)
```
### Normal Modes & Thermal B-factors:
The normal modes are defined by the eigenvectors of the mass-weighted Hessian $\hat{H}\vec{e}_k = \tilde{\omega}^2_k\vec{e}_k$. As mensioned earlier, the 6 firsts modes should not be taken into account. The extraction of the modes is acieved by the function `NMA.collective_modes`. 

One of the useful metric that can be extracted from those modes is the Thermal B-factors of the nodes. The extraction of thermal B-factors is achieved with `NMA.local_MSF`. Finally the code bellow runs a Bfactor prediction and save it in a pdb file to allow easy vizualization.

```python
import md_manager as md
import md_manager.NMA as nma
import pandas as pd
from scipy.spatial import distance_matrix

# load structure
df = md.PDBfile("file.pdb").read2df()
nodes_position = df[["x", "y", "z"]].to_numpy(dtype = float)
nodes_mass = df.m.to_numpy(dtype = float)
distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

# compute Hessian with appropriate model
H = nma.pfANM_hessian(nodes_position=nodes_position, nodes_mass=nodes_mass, distance_inter_nodes=distance_inter_nodes, spring_constant=1.0)

# compute normal modes & thermal factors
eigenvalues, eigenvectors = nma.collective_modes(H)
B = nma.local_MSF(eigenvals=eigenvalues, eigenvecs=eigenvectors, nodes_mass=nodes_mass, convert2bfactors=True)

# save prediction in pdb file:
df.b = pd.Series(B, index=df.index)
md.PDBfile.df2pdb(df=df, filename="out.pdb")
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