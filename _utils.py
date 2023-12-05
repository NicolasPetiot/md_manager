from ._parameters import DF_COLUMNS, ATOMIC_MASSES

import pandas as pd
import numpy as np

def atom_position(df:pd.DataFrame) -> np.ndarray:
    """
    Returns a numpy.ndarray that contains the [x, y, z] coordinates of the atoms in df.
    """
    return df[["x", "y", "z"]].to_numpy(dtype=float)

def atom_mass(df:pd.DataFrame) -> np.ndarray:
    """
    Returns a numpy.ndarray that contains the masses of the atoms in df.
    """
    return df["m"].to_numpy(dtype=float)

def atom_bfactors(df:pd.DataFrame) -> np.ndarray:
    """
    Returns a numpy.ndarray that contains the thermal bfactors of the atoms in df.
    """
    return df["b"].to_numpy(dtype=float)

def _build_df_from_atom_list(atoms:list):
    """
    Returns a pandas.DataFrame that contains the atoms in the input list.
    """
    columns=DF_COLUMNS
    return pd.DataFrame(atoms, columns=columns)

def _scan_pdb_line(line:str):
    """
    Returns a tuple that contains the (record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass) 
    informations about an atom line in the PDB file.
    """
    if line.startswith("ATOM"):
        record_name = "ATOM  "
    elif line.startswith("HETATM"):
        record_name = "HETATM"
    else :
        raise ValueError("Input 'line' is not associated to an atom in a pdb file.")
    
    #atom_id   = line[ 6:11]
    name      = line[12:16]
    alt       = line[16:17]
    resn      = line[17:20]
    chain     = line[21:22]
    resi      = int(line[22:26])
    insertion = line[26:27]
    x         = float(line[30:38])
    y         = float(line[38:46])
    z         = float(line[46:54])
    try :
        occupancy = float(line[54:60])
    except ValueError :
        occupancy = 1.0
    try :
        b         = float(line[60:66])
    except ValueError :
        b         = 0.0
    segi      = line[72:76]
    elem      = line[76:78]
    charge    = line[78:79]
    try :
        mass      = ATOMIC_MASSES[elem]
    except KeyError:
        raise KeyError(f"Unknown element symbol '{elem}'. Please update the ATOMIC_MASSES dictionary in '_parameters.py'.")
    return record_name, name, alt, resn, chain, resi, insertion, x, y, z, occupancy, b, segi, elem, charge, mass

def _generate_pdb_line(atom:pd.Series, id:int) -> str:
    """
    Generates a pdb line from an input series that contains atom's information.
    """
    # extract strings :
    record_name = "%-6s"%atom["record_name"]
    name        = " %-3s"%atom["name"]
    alt         = "%1s"%atom["alt"]
    resn        = "%3s"%atom["resn"]
    chain       = "%1s"%atom["chain"]
    insertion   = "%1s"%atom["insertion"]
    segi        = "%-4s"%atom["segi"]
    elem        = "%2s"%atom["e"]
    charge      = "%2s"%atom["q"]

    # extract float/int
    id    = f"{id:5d}"
    resi  = f"{atom.resi:4d}"
    x     = f"{atom.x:8.3f}"
    y     = f"{atom.y:8.3f}"
    z     = f"{atom.z:8.3f}"
    occup = f"{atom.occupancy:6.2f}"
    b     = f"{atom.b:6.2f}"

    line = "%s%s%s%s%s %s%s%s   %s%s%s%s%s      %s%s%s\n"%(
        record_name, id, name, alt, resn, chain, resi, insertion, x, y, z, occup, b, segi, elem, charge
    )
    return line

def _write_pdb_atom_lines(file, df):
    """
    Write the atom lines in file object.
    """
    id = 0
    # default values :
    if df.alt.isna().all():
        df["alt"] = ""
    if df.insertion.isna().all():
        df["insertion"] = ""
    if df.occupancy.isna().all():
        df["occupancy"] = 1.0
    if df.b.isna().all():
        df["b"] = 0.0
    if df.segi.isna().all():
        df["segi"] = ""
    if df.e.isna().all():
        df["e"] = " C"
    if df.q.isna().all():
        df["q"] = ""
    
    APO = df[df.record_name == "ATOM  "]
    for _, chain in APO.groupby(["chain"]):
        for _, atom in chain.iterrows():
            id += 1
            line = _generate_pdb_line(atom, id)
            file.write(line)
        file.write("TER\n")

    
    HET = df[df.record_name == "HETATM"]
    for _, atom in HET.iterrows():
        id += 1
        line = _generate_pdb_line(atom, id)
        file.write(line)