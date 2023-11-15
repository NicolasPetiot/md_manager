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

* `open` & `close` : Onpens/Close the pdb.file wrapper.
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

## Anisotropic Network Models (ANM)

## Free-energy landscape analysis

---
Contact : Nicolas.Petiot01@u-bourgogne.fr