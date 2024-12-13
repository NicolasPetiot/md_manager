from .md_files import *
from .parameters import ATOMIC_MASSES

import pandas as pd
import numpy as np

import sys
from urllib.request import urlopen
from os import path

__all__ = ["load", "save", "pdb2df", "fetch_PDB", "atomic_masses", "com", "shift_df", "rotate_df"]

######################################################################################################################################
# Load structure

def load(filename:str) -> pd.DataFrame:
    """
    Load the first frame of a molecular dynamics trajectory file into a DataFrame.

    The function determines the appropriate trajectory class (e.g., `XYZ`, `PDB`, or `GRO`)
    based on the file extension and reads the first frame of atomic coordinates.

    Args:
        filename (str): Path to the trajectory file. Must have one of the following extensions:
                        `.xyz`, `.pdb`, `.gro`.

    Returns:
        pd.DataFrame: A DataFrame containing atomic coordinates and other relevant atom-level
                      information from the first frame of the trajectory.

                      The columns and format depend on the file type:
                      - `XYZ`: ["name", "x", "y", "z"]
                      - `PDB`: ["record_name", "name", "alt", "resn", "chain", "resi",
                                "insertion", "x", "y", "z", "occupancy", "b", "segi", "e", "q"]
                      - `GRO`: ["resi", "resn", "name", "x", "y", "z", "vx", "vy", "vz"]

    Raises:
        ValueError: If the file extension is not one of `.xyz`, `.pdb`, or `.gro`.

    Example:
        >>> df = load("sample.xyz")
        >>> print(df.head())
           atom_id name      x      y      z
        1        1   H  0.000  0.000  0.000
        2        2   O  1.000  0.000  0.000
    """
    file_format = {
        ".xyz" : XYZ,
        ".pdb" : PDB,
        ".gro" : GRO
    }
    _, ext = path.splitext(filename)
    if not ext in file_format:
        raise ValueError(f"Unknown extention '{ext}'. File must be in format {list(file_format.keys())}")

    # No matter the file format, getting the next df must always be as follows:
    traj = file_format[ext](filename) # Creates traj
    traj = iter(traj)                 # Open file
    return next(traj)                 # Read coordinates

def save(filename:str, data:pd.DataFrame|list[pd.DataFrame]):
    """
    Save atomic data into a trajectory file in a specified format.

    This function writes atomic coordinates and related data into a trajectory file 
    (e.g., `.xyz`, `.pdb`, `.gro`) in the appropriate format. If a list of DataFrames is 
    provided, each DataFrame corresponds to a frame in the trajectory.

    Args:
        filename (str): The path to the output trajectory file. The file extension must be 
                        one of `.xyz`, `.pdb`, or `.gro` to determine the appropriate format.
        data (pd.DataFrame | list[pd.DataFrame]): The atomic data to be saved. This can either 
                                                  be a single DataFrame (representing one frame) 
                                                  or a list of DataFrames (representing multiple 
                                                  frames). Each DataFrame must conform to the 
                                                  column structure of the specified format:
                                                  - `XYZ`: ["name", "x", "y", "z"]
                                                  - `PDB`: ["record_name", "name", "alt", "resn", 
                                                            "chain", "resi", "insertion", "x", "y", 
                                                            "z", "occupancy", "b", "segi", "e", "q"]
                                                  - `GRO`: ["resi", "resn", "name", "x", "y", "z", 
                                                            "vx", "vy", "vz"]

    Raises:
        ValueError: If the file extension is not one of `.xyz`, `.pdb`, or `.gro`.

    Example:
        >>> df1 = pd.DataFrame({
        ...     "name": ["H", "O"],
        ...     "x": [0.0, 1.0],
        ...     "y": [0.0, 0.0],
        ...     "z": [0.0, 0.0]
        ... })
        >>> df2 = pd.DataFrame({
        ...     "name": ["H", "O"],
        ...     "x": [1.0, 2.0],
        ...     "y": [0.0, 0.0],
        ...     "z": [0.0, 0.0]
        ... })
        >>> save("output.xyz", [df1, df2])
    """
    file_format = {
        ".xyz" : XYZ,
        ".pdb" : PDB,
        ".gro" : GRO
    }
    _, ext = path.splitext(filename)
    if not ext in file_format:
        raise ValueError(f"Unknown extention '{ext}'. File must be in format {list(file_format.keys())}")
    
    new = file_format[ext](filename, "w")
    if type(data) == pd.DataFrame:
        data = [data]

    for i, df in enumerate(data):
        new.write_frame(df, model_id=i+1)


def pdb2df(filename, atom_only = False) -> pd.DataFrame:
    """
    Retruns a DataFrame associated to the atoms found in the first model of the input PDB file.
    """
    pdb = PDB(filename)
    pdb = iter(pdb)
    df = next(pdb)

    if atom_only:
        df = df.query("record_name == 'ATOM'")

    return df

def fetch_PDB(pdb_code:str, atom_only = False) -> pd.DataFrame:
    """
    Fetches a PDB structure from the RCSB PDB website and returns it as a pandas DataFrame.

    This function downloads a PDB file from the RCSB website based on the provided PDB code, 
    parses the structure, and returns it as a DataFrame. The DataFrame includes atomic data 
    from the PDB file. Optionally, only ATOM records can be returned.

    Args:
        pdb_code (str): The PDB code of the structure to fetch. This is a four-character identifier
                         assigned to each PDB entry.
        atom_only (bool, optional): If `True`, the returned DataFrame will only include 'ATOM' records
                                    (ignoring any HETATM or other record types). Defaults to `False`.

    Returns:
        pd.DataFrame: A DataFrame containing the atomic information from the fetched PDB file. 
                      The columns correspond to the fields defined by the `PDB` class.

    Raises:
        HTTPError: If the PDB code is invalid or the RCSB PDB website is unreachable.
        ValueError: If the PDB file cannot be parsed properly.

    Example:
        >>> df = fetch_PDB("1A2B", atom_only=True)
        >>> print(df.head())
           record_name name alt resn chain  resi insertion      x      y      z occupancy    b   segi e  q
        1          ATOM  CA   ALA     A    23          A  1.234  2.345  3.456     1.00  20.5  A  C  0
        2          ATOM  C    ALA     A    23          A  2.234  3.345  4.456     1.00  18.5  A  C  0
    """
    url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
    response = urlopen(url)

    txt = response.read()
    lines = (txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.decode("ascii")).splitlines()

    df = PDB.next_frame(lines)

    if atom_only:
        df = df.query("record_name == 'ATOM'")

    return df

def df2pdb(filename:str, data:pd.DataFrame|list[pd.DataFrame]):
    """Generates a pdb structure/trajctory according to the type of the input data"""
    new = PDB(filename, "w")

    if type(data) == pd.DataFrame:
        data = [data]

    for i, df in enumerate(data):
        new.write_frame(df, model_id=i+1)

######################################################################################################################################
# DataFrame manipulation

def atomic_masses(df:pd.DataFrame) -> pd.Series:
    """
    Returns a pandas Series containing atomic masses in Dalton (Da) based on the 'e' (element) column 
    of the input DataFrame.

    This function extracts the atomic element symbols from the 'e' column of the input DataFrame 
    and returns the corresponding atomic mass in Dalton (Da) for each atom. If an element is not 
    found in the `ATOMIC_MASSES` dictionary, it will be assigned a mass of 1.0 Da by default.

    Args:
        df (pd.DataFrame): A DataFrame containing at least one column named 'e' which holds the element 
                            symbols (e.g., 'H', 'C', 'O') for each atom.

    Returns:
        pd.Series: A pandas Series containing the atomic masses for each element found in the 'e' column. 
                   The index of the Series corresponds to the index of the input DataFrame.

    Warnings:
        If any unknown elements (not present in the `ATOMIC_MASSES` dictionary) are found in the 'e' column, 
        a warning will be printed indicating the number of unknown elements and the list of known elements.

    Example:
        >>> df = pd.DataFrame({'e': ['H', 'C', 'O', 'X']})
        >>> masses = atomic_masses(df)
        >>> print(masses)
        0     1.008
        1    12.011
        2    15.999
        3     1.000
        Name: e, dtype: float64
    """
    def get_atom_mass(elem:str) -> float:
        try:
            m = ATOMIC_MASSES[elem.upper()]

        except KeyError:
            m = 1.0
        return m
    
    if not "e" in df.columns:
        df["e"] = ""

    elem = df["e"]
    mass = elem.apply(get_atom_mass)

    num_unknown_elem = (mass == 1.0).sum()
    if num_unknown_elem > 0:
        print(f"Warning: Found {num_unknown_elem}/{len(mass)} unknown elements. The input DataFrame should contains 'e' (element) column. The known elements are {list(ATOMIC_MASSES.keys())}")

    return mass

def com(df:pd.DataFrame) -> pd.DataFrame:
    """
    Returns the center of mass (COM) of the atoms in the input DataFrame.

    This function calculates the center of mass of the atoms described in the input DataFrame `df`. 
    The center of mass is calculated using the atomic positions (x, y, z) and masses. If the 
    DataFrame contains a column 'm' for atomic masses, it will be used; otherwise, the function 
    will attempt to use the atomic masses based on the element symbols in the 'e' column of the 
    DataFrame.

    Args:
        df (pd.DataFrame): A DataFrame containing atomic data. The DataFrame should include 
                            at least 'x', 'y', 'z' columns for atomic positions, and optionally 
                            a 'm' column for atomic masses.

    Returns:
        pd.DataFrame: A pandas DataFrame with a single row representing the center of mass, 
                       with columns 'x', 'y', and 'z'.

    Example:
        >>> df = pd.DataFrame({
                'x': [0.0, 1.0, 2.0],
                'y': [0.0, 1.0, 2.0],
                'z': [0.0, 1.0, 2.0],
                'm': [1.0, 12.0, 16.0]  # masses in Da
            })
        >>> com(df)
        x    1.166667
        y    1.166667
        z    1.166667
        dtype: float64
    """
    xyz = ["x", "y", "z"]
    pos = df[xyz]

    if "m" in df.columns:
        mass = df["m"]
    else:
        mass = atomic_masses(df)

    return pos.apply(lambda x, m: x*m, m=mass).sum() / mass.sum()  

def radius_of_gyration(df: pd.DataFrame) -> float:
    """
    Calculate the radius of gyration (Rg) of the molecule represented by the given DataFrame.

    The radius of gyration is a measure of the compactness of a molecule, defined as the 
    square root of the average squared distance of each atom from the center of mass (COM).

    This function assumes the molecule's atomic coordinates are provided in the 'x', 'y', 
    and 'z' columns of the DataFrame. If the 'm' column (atomic mass) is not available, 
    the atomic masses will be assumed to be 1.0 Dalton by default.

    Parameters:
    df : pandas.DataFrame
        A DataFrame containing atomic coordinates of the molecule, with columns 'x', 'y', 
        'z', and optionally 'm' for atomic masses.

    Returns:
    float
        The radius of gyration (Rg) of the molecule in the same units as the atomic 
        coordinates (typically Angstroms).
    """
    com_pos = com(df)  # Calculate the center of mass
    pos = df[["x", "y", "z"]].values - com_pos.values  # Shift coordinates to COM
    rg = np.sqrt(np.sum(np.sum(pos**2, axis=1)) / len(df))
    return rg  

def shift_df(df:pd.DataFrame, translation_vector:np.ndarray|pd.Series) -> pd.DataFrame:
    """
    Shifts the coordinates of atoms in the DataFrame by the given translation vector.

    This function adds the provided `translation_vector` to the 'x', 'y', and 'z' columns in the input 
    DataFrame `df`, effectively translating the atomic positions by the specified vector.

    Args:
        df (pd.DataFrame): A DataFrame containing atomic data with 'x', 'y', and 'z' columns representing 
                            the atomic positions.
        translation_vector (np.ndarray | pd.Series): A 1D array-like object representing the translation 
                                                     vector to be added to the 'x', 'y', and 'z' coordinates.
                                                     The vector should have three components (corresponding 
                                                     to x, y, and z translations).

    Returns:
        pd.DataFrame: A new DataFrame with the 'x', 'y', and 'z' coordinates shifted by the translation vector.
        
    Example:
        >>> df = pd.DataFrame({
                'x': [0.0, 1.0, 2.0],
                'y': [0.0, 1.0, 2.0],
                'z': [0.0, 1.0, 2.0]
            })
        >>> translation_vector = np.array([1.0, -1.0, 0.5])
        >>> shifted_df = shift_df(df, translation_vector)
        >>> shifted_df
           x    y    z
        0  1.0 -1.0  0.5
        1  2.0  0.0  2.5
        2  3.0  1.0  2.5
    """
    # Get position and apply shifting using vectorized methods
    xyz = ["x", "y", "z"]
    pos = df[xyz]
    df[xyz] = pos.apply(lambda s, v: s+v, v=translation_vector, axis = 1)


def rotate_df(df: pd.DataFrame, axis: np.ndarray | pd.Series, angle: float) -> pd.DataFrame:
    """
    Rotates the coordinates of atoms in the DataFrame around a given axis by the specified angle.

    This function rotates the atomic positions (in 'x', 'y', 'z') in the input DataFrame `df` around the
    specified axis by the given angle (in radians). The rotation is performed using the 3D rotation matrix.

    Args:
        df (pd.DataFrame): A DataFrame containing atomic data with 'x', 'y', and 'z' columns representing 
                            the atomic positions.
        axis (np.ndarray | pd.Series): A 1D array-like object representing the axis of rotation.
                                       The axis should be a vector in 3D space (x, y, z).
        angle (float): The angle (in radians) by which to rotate the atoms around the axis.

    Returns:
        pd.DataFrame: A new DataFrame with the 'x', 'y', and 'z' coordinates rotated around the specified axis.

    Example:
        >>> df = pd.DataFrame({
                'x': [1.0, 0.0, 0.0],
                'y': [0.0, 1.0, 0.0],
                'z': [0.0, 0.0, 1.0]
            })
        >>> axis = np.array([0.0, 0.0, 1.0])  # Rotate around Z-axis
        >>> angle = np.pi / 2  # 90 degrees in radians
        >>> rotated_df = rotate_df(df, axis, angle)
        >>> rotated_df
           x    y    z
        0  0.0  1.0  0.0
        1 -1.0  0.0  0.0
        2  0.0 -1.0  0.0
    """
    # Normalize the axis vector
    axis = np.asarray(axis)
    axis = axis / np.linalg.norm(axis)

    # Create the rotation matrix using Rodrigues' rotation formula
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    ux, uy, uz = axis

    # Rotation matrix
    rotation_matrix = np.array([
        [cos_angle + ux**2 * (1 - cos_angle), ux * uy * (1 - cos_angle) - uz * sin_angle, ux * uz * (1 - cos_angle) + uy * sin_angle],
        [uy * ux * (1 - cos_angle) + uz * sin_angle, cos_angle + uy**2 * (1 - cos_angle), uy * uz * (1 - cos_angle) - ux * sin_angle],
        [uz * ux * (1 - cos_angle) - uy * sin_angle, uz * uy * (1 - cos_angle) + ux * sin_angle, cos_angle + uz**2 * (1 - cos_angle)]
    ])

    # Apply the rotation matrix to each atom in the DataFrame
    xyz = ["x", "y", "z"]
    pos = df[xyz].values

    # Perform the matrix multiplication to rotate the coordinates
    rotated_pos = np.dot(pos, rotation_matrix.T)

    # Update the DataFrame with the rotated positions
    df[xyz] = rotated_pos

    return df

######################################################################################################################################
# Conformational angles

class InvalidChainException(Exception):
    """Raised when a chain is missing / contains too many atoms"""
    pass

def check_chain_validity(chain:pd.DataFrame, maxwarn = 5) -> None:
    """"""
    will_raise_exception = False
    if not "resi" in chain:
        print("Warning: Input DataFrame does not contains 'resi' column. Will assume that the associated chain is not missing residue.")

    else:
        Nwarn = 0
        for win in chain.resi.rolling(2):
            if len(win) > 1:
                resi1 = win.iloc[1]
                resi0 = win.iloc[0]
                diff  = resi1 - resi0
                
                if diff == 0:
                    print(f"Warning: Found several atoms that bellongs to the same residue at index {resi0}.")
                    will_raise_exception = True
                    Nwarn += 1
                
                if diff < 0:
                    print("Warning: Found several chains in the input DataFrame.")
                    will_raise_exception = True
                    Nwarn += 1
                
                if diff > 1:
                    print(f"Warning: Missing residue(s) between resi {resi0:4d} and {resi1:4d}")
                    will_raise_exception = True
                    Nwarn += 1

                if Nwarn == maxwarn:
                    break
                    print("Maximum warning reached...")

    if will_raise_exception:
        raise InvalidChainException("Found warning(s) that does not allows conformational angles calculation.")