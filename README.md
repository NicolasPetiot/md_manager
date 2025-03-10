# MD-manager: Protein Conformation Analysis Toolkit

MD-manager is a set of python functions and methods that provide easy acess and manipulation of data for molecular dynamics simulations. Python package designed for analyzing protein structure and conformational properties. This toolkit allows users to calculate various structural features such as bond angles, dihedral angles, backbone and side-chain conformations, and more. It is developped in the research group of Physics Applied to Proteins at _université de bougrogne_.

## Installation:

As MD-manager is still in early stage of development, it is not yet available on the [pypi](https://pypi.org) code base.

```shell
git clone https://github.com/NicolasPetiot/md_manager
cd md_manager
pip install .
```

## PDB class:

The PDB class in MD-manager is a specialized class for handling Protein Data Bank (PDB) files, enabling efficient manipulation and analysis of molecular structures. It provides methods for reading and writing atomic data in the PDB format, extracting key information such as atom coordinates, chain IDs, and residue details.

    Reading PDB Files: The class can parse PDB files to extract and store atomic data in a structured pandas.DataFrame, making it easy to work with molecular information programmatically.
    Writing PDB Files: It also supports exporting molecular data back to the PDB format, facilitating the generation of modified or newly simulated structures.
    Frame Management: For handling molecular dynamics trajectories, the PDB class supports processing individual frames, allowing for manipulation of specific timesteps in the simulation.

The class helps streamline the management of PDB files in bioinformatics workflows and supports custom analyses like calculating angles, distances, and conformational states.

By default, the DataFrames contains the following columns :

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

**N.B.**: The `PDB.__next__()` method uses the "ENDMDL" record name to return the atoms in the DataFrame. If none of the lines starts with "ENDMDL", the DataFrame will contains all atoms in the input file. For files that contains a single frame, the `pdb2df` function is a much more suiable way to read pdb files.

```python
import md_manager as md
df = md.load("file.pdb")
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

df = md.fetch_PDB("5f0g", atom_only = True) # remove hetero atoms

# save to pdb format:
md.save("5f0g_APO.pdb", df)
```

Trajectories can be generated as well using a list of DataFrames:

```python
import md_manager as md

traj = md.PDB("full_traj.pdb")
traj = iter(traj)

# save 10 firsts models to pdb format:
md.save("sub_traj.pdb", [next(traj) for _ in range(10)])
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
