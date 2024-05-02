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
- `write:bool = False` : Allow the creation of a the file is file not found.

```python
pdb = PDB("inputfile.pdb")
```

`__iter__` : Initialises the iteration over models.

```python
pdb = iter(pdb)
```

`__next__` : Iterates over the lines in the file and returns the DataFrame associated for each 'ENDMDL' in the file.

```python
pdb = iter(pdb)
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
pdb.write(model=df)
```

`build_df_from_atoms` : Creates the DataFrame associated to input `atoms` with expected column names and types.

args :

- `atoms:list[tuple]` : List created from `scan_pdb_line` method.

```python
df = PDB.build_df_from_atoms(atoms)
```

`scan_pdb_line` : Returns a tuple that contains the (record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass) informations about an atom line in the PDB file.

By default,  the occupancy is set to 1.0.

By default, the Bfactors are set to 0.0 AA².

The mass is extracted from the `ATOMIC_MASSES` dictionary (see _params.py)

args :

- `line:str` : Line of a pdb file that starts with "ATOM"/"HETATM"

```python
atom = PDB.scan_pdb_line(line)
```

`generate_atom_line` : Generates a pdb line from an input Series that contains atom's information.

args :

- `atom:Series` : Series associated to a line of a DataFrame.
- `atom_id:int` : Atom index in the pdb file.

```python
line = PDB.generate_atom_line(atom, atom_id = 1)
```

`generate_atom_lines` : Generates a list of lines to be written in a pdb file from input DataFrame and model index.

args :

- `df:DataFrame` : DataFrame to convert in lines in pdb format.
- `model_id:int = 1` : Number of the model associated to the current frame.

```python
lines = PDB.generate_atom_lines(df)
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

The `fetch_protein_data_bank` funcion allow to create DataFrames from the pdb code of a structure deposed in the [protein data bank](https://www.rcsb.org) database.

```python
import md_manager as md
df = md.fetch_protein_data_bank("8q89")
```

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

## Scodary structure prediction:

MD-manager contains an implementation of the [CUTABI](https://www.frontiersin.org/articles/10.3389/fmolb.2021.786123/full) program. Such prediction can be performed as follows:

```python
import md_manager as md
from md_manager.CUTABI import(
    predict_alpha_helix, 
    predict_beta_sheets
)

CA = md.fetch_protein_data_bank("8q89").query("name == 'CA'")

alpha = predict_alpha_helix(CA)
beta = predict_beta_sheets(CA)

CA["alpha"] = alpha
CA["beta"] = beta

# Display residues identified in an alpha helix or beta sheet
CA.query("alpha or beta")
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
