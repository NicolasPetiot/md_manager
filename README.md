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
In the ANM, Hessians are computed based on the potential function
$$
    V = \frac{\gamma}{2}\sum_{i, j}(r_{ij} - r_{ij}^0)^2\theta (r_c - r_{ij}^0)
$$
* $r_{ij}$ is the distance between nodes $i$ and $j$.
* $r_{ij}^0$ is the distance between nodes $i$ and $j$ at the equilibrium.
* $\gamma$ is the force constant defining the interactions between nodes.
* $\theta$ is the heavyside step function defining if the interactions between node $i$ and $j$ are considered or not.
* $r_c$ is a cutoff radius. If $r_{ij}^0 < r_c$ the interaction between $i$ and $j$ is considered.

From this definition the Hessian is defined as 
$$
    H_{ij} = -\frac{\gamma}{(r_{ij}^0)^2}\left[ 
        \begin{array}{ccc}
            x_{ij}^0x_{ij}^0 & x_{ij}^0y_{ij}^0 & x_{ij}^0z_{ij}^0\\
            y_{ij}^0x_{ij}^0 & y_{ij}^0y_{ij}^0 & y_{ij}^0z_{ij}^0\\
            z_{ij}^0x_{ij}^0 & z_{ij}^0y_{ij}^0 & z_{ij}^0z_{ij}^0\\
        \end{array}
     \right]\theta (r_c - r_{ij}^0)
$$

For the mass-weighted Hessian, the elements $i \ne j$ are defined as
$$
    \hat{H}_{ij} = \cfrac{H_{ij}}{\sqrt{M_iM_j}}
$$
To ensure the invariance by translation, the diagonal elements are defined as 
$$
    \hat{H}_{ii} = -\sum_{j\ne i}\hat{H}_{ij}
$$

The construction of such matrix is achieved by the function `md.ANM_hessian`. For code efficiancy, this version of `md_manager` needs the `scipy.spatial.distance_matrix` function to compute the distance inter-nodes.

```python
import md_manager as md
from scipy.spatial import distance_matrix

df = md.PDBfile("file.pdb").read2df()
nodes_position = md.atom_position(df)
nodes_mass = md.atom_mass(df)
distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

H = md.ANM_hessian(nodes_position=nodes_position, nodes_mass=nodes_mass, distance_inter_nodes=distance_inter_nodes, cutoff_radius=4.0, spring_constant=1.0)
```
### Parameter-free Anisotropic Network Model :
The ANM Hessian needs a cutoff radius that ensure that only 6 eigenvalues are null from the 3 global translations and 3 global rotations. For large number of nodes this can be problematic as it needs to perform a several diagonalization to tune the value of $r_c$. Introducing the pfANM Hessian. In this representation, all nodes are interacting with each others but the intensity of the interaction depends on their relative position.
$$
    H_{ij} = -\frac{\gamma}{(r_{ij}^0)^4}\left[ 
        \begin{array}{ccc}
            x_{ij}^0x_{ij}^0 & x_{ij}^0y_{ij}^0 & x_{ij}^0z_{ij}^0\\
            y_{ij}^0x_{ij}^0 & y_{ij}^0y_{ij}^0 & y_{ij}^0z_{ij}^0\\
            z_{ij}^0x_{ij}^0 & z_{ij}^0y_{ij}^0 & z_{ij}^0z_{ij}^0\\
        \end{array}
     \right]
$$

The construction of such matrix is achieved by the function `md.pfANM_hessian`.

```python
df = md.PDBfile("file.pdb").read2df()
nodes_position = md.atom_position(df)
nodes_mass = md.atom_mass(df)
distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

H = md.pfANM_hessian(nodes_position=nodes_position, nodes_mass=nodes_mass, distance_inter_nodes=distance_inter_nodes, spring_constant=1.0)
```
### Normal Modes & Thermal B-factors:
The normal modes are defined by the eigenvectors of the mass-weighted Hessian $\hat{H}\vec{e}_k = \tilde{\omega}^2_k\vec{e}_k$. As mensioned earlier, the 6 firsts modes should not be taken into account. The extraction of the modes is acieved by the function `md.collective_modes`. 

One of the useful metric that can be extracted from those modes is the Thermal B-factors of the nodes.

$B_i = \frac{8\pi^2}{3}k_{B}T\sum_{k} \frac{1}{\tilde{\omega}_{k}^2}\left( \frac{\vec{e}_{k, i}}{\sqrt{M_{i}}} \right)^2$



The extraction of thermal B-factors is achieved with `md.local_MSF`. Finally the code bellow runs a Bfactor prediction and save it in a pdb file to allow easy vizualization.

```python
import md_manager as md
import pandas as pd
from scipy.spatial import distance_matrix

# load structure
df = md.PDBfile("PDB/AlphaFold/GstD01_AlphaFold2.pdb").read2df()
nodes_position = md.atom_position(df)
nodes_mass = md.atom_mass(df)
distance_inter_nodes = distance_matrix(nodes_position, nodes_position)

# compute Hessian with appropriate model
H = md.pfANM_hessian(nodes_position=nodes_position, nodes_mass=nodes_mass, distance_inter_nodes=distance_inter_nodes, spring_constant=1.0)

# compute normal modes & thermal factors
eigenvalues, eigenvectors = md.collective_modes(H)
B = md.local_MSF(eigenvals=eigenvalues, eigenvecs=eigenvectors, nodes_mass=nodes_mass, convert2bfactors=True)

# save prediction in pdb file:
df.b = pd.Series(B, index=df.index)
md.PDBfile.df2pdb(df=df, filename="out.pdb")
```

---
Contact : Nicolas.Petiot01@u-bourgogne.fr