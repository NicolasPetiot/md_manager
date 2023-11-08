# MD-manager
md_manager is a set of python functions and methods that provide easy acess and manipulation of data for molecular dynamics simulations. It is developped in the research group of Physics Applied to Proteins in the université de bougrogne. 

## AtomFrame
In the current stage of the project, md_manager can handle pdb file format. Let first see how to load a pdb structure. We decided to use Pandas DataFrames to save the atom information. It allows to use pandas syntax to perform all kind of selections.

```python
import md_manager as md
df = md.pdb2df('file.pdb')
CA = df.query("name == ' CA ' and chain == 'A'")
```

Working with md trajectories often involves looping over the models in a pdb file. It is possible to do so as follows :

```python
import md_manager as md
import pandas as pd

def do_something(df:pd.DataFrame):
    # TODO
    return None

filename = "traj.pdb"
with md.MdTraj(filename) as traj :
    for model in traj :
        result = do_something(model)
```

Atom position, mass and thermal factors are features that might be often needed for a broad range of analysis. The md_manager provides functions that allows easy acess to these metrics.

```python
import md_manager as md

# load np.ndarray(s):
df    = md.pdb2df("file.pdb")
xyz   = md.atom_position(df)
mass  = md.atom_mass(df)
bfact = md.atom_bfactors(df)
```

## Anisotropic Network Models (ANM)

## Free-energy landscape analysis

---
Contact : Nicolas.Petiot01@u-bourgogne.fr