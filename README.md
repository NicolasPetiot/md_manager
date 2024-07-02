# MD-manager

MD-manager is a set of python functions and methods that provide easy acess and manipulation of data for molecular dynamics simulations. It is developped in the research group of Physics Applied to Proteins at _université de bougrogne_.

## Installation:

As MD-manager is still in early stage of development, it is not yet available on the [pypi](https://pypi.org) code base. To download it use

```shell
git clone https://github.com/NicolasPetiot/md_manager
cd md_manager
pip install .
```

## PDB class:

A python class that provide methods to read and write pdb files.

Unless your code specifically adds, remove or rename columns in the DataFrames, they contains the following columns :

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

### Methods :

```python
from md_manager import PDB
```

`__init__` : Creates an instance of the `PDB` class.

args :

- `filename:str` : File name or relative path of the pdb file to read/write.
- `mode = 'r'` : Allow the creation of a the file is file not found (for write mode).

```python
pdb = PDB("inputfile.pdb")
```

`__iter__` : Initialises the iteration over models.

```python
pdb = iter(pdb)
```

`__next__` : Iterates over the lines in the file and returns the DataFrame associated for each 'ENDMDL' in the file.

```python
df = next(pdb)

# or 
for df in pdb:
    (...)
```

`open` : Set the I/O wrapper associated to the instance of PDB in open mode.

args :

- `mode:str = "r"` I/O interation mode (read by default)

```python
pdb.open()
```

`close` : Set the I/O wrapper associated to the instance of PDB in close mode.

```python
pdb.close()
```

`write` : Read the input DataFrame(s) and creates the associated pdb file.

args :

- `model:DataFrame` : Single frame to be written in the output file.
- `model_list:list[DataFrame]` : List of frames to be written in the output file.

```python
pdb.write(df=df)
```

### Reading files:

For MD trajectories, the `PDB` class allow the following syntax

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

**N.B.**: The `PDB.__next__()` method uses the "ENDMDL" record name to return the atoms in the DataFrame. If none of the lines starts with "ENDMDL", the DataFrame will contains all atoms in the input file. For files that contains a single frame, the `pdb2df` method is a much more suiable way to read pdb files.

```python
import md_manager as md
df = md.pdb2df("file.pdb")
```

The `fetch_PDB` funcion allow to create DataFrames from the pdb code of a structure deposed in the [protein data bank](https://www.rcsb.org) database.

```python
import md_manager as md
df = md.fetch_PDB("8q89")
```

### Writing file:

MD-manager allows to write pdb files as follows:

```python
import md_manager as md

df = md.fetch_PDB("5f0g")

# remove hetero atoms:
df = df.query("record_name == 'ATOM'")

# save to pdb format:
md.df2pdb("out.pdb", df)
```

Trajectories can be generated as well using a list f DataFrames:

```python
import md_manager as md

traj = md.PDB("full_traj.pdb")
traj = iter(traj)

# save 10 firsts models to pdb format:
md.df2pdb("sub_traj.pdb", [next(traj) for _ in range(10)])
```

## Conformational Angles:

MD-manager contains methods that allow to compute conformational angles for proteins.

```python
import md_manager as md
import numpy as np

xyz = np.array([
    [0, 0, 0],
    [0, 0, 1],
    [0, 1, 1],
    [1, 1, 1]
])
thetas = md.angles(xyz)
gammas = md.dihedral_angles(xyz)
```

It can be used to compute conformational landscape as did in [this](https://www.mdpi.com/2076-3417/12/16/8196) work.

## Secondary structure prediction:

MD-manager contains an implementation of the [CUTABI](https://www.frontiersin.org/articles/10.3389/fmolb.2021.786123/full) program. Such prediction can be performed as follows:

```python
import md_manager as md
from md_manager.CUTABI import(
    predict_alpha_helix, 
    predict_beta_sheets
)

CA = md.fetch_PDB("8q89").query("name == 'CA'")

alpha = predict_alpha_helix(CA, CA_only=True)
beta = predict_beta_sheets(CA, CA_only =True)

CA["helix"] = alpha
CA["sheet"] = beta

# Save structure with CUTABI predicted helices and sheets:
md.df2pdb("8q89_CA_CUTABI_SS.pdb", CA)

# Display residues identified in an alpha helix or beta sheet
CA.query("helix or sheet")
```

## Normal Mode Analysis:

Normal Modes are extracted from the diagonalization of mass-weighted Hessian matricies.
Choosing a NMA model is equivalent to chose a way to build Hessians. The computation of the modes can be performed as follows.

```python
import md_manager as md
from md_manager.NMA import *

df = md.fetch_protein_data_bank("8q89")

H = ANM_hessian(df = df, cutoff_radius=4.0, spring_constant=1.0)
# or 
H = pfANM_hessian(df = df, spring_constant=1.0)

eigenvals, modes = collective_modes(H)
```

---

Contact : Nicolas.Petiot01@u-bourgogne.fr
