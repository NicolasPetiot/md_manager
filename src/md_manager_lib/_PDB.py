from ._Traj import Traj, TrajParams

import pandas as pd

import sys
from urllib.request import urlopen
from dataclasses import asdict

__all__ = [
    "PDB",
    "fetch_PDB",

]

class PDB(Traj):
    default_params = asdict(TrajParams(return_q = False, return_type = False, return_segi = False, return_icode = False))

def lines2df(lines:list[str], atom_only = False) -> pd.DataFrame:
    def read_format(line:str) -> tuple[str, str, str, str, int, str, float, float, float, float, float, float, str, str, str]:
        """
        Parses a line from a PDB file to extract atomic data.

        Args:
            line (str): A single line from the PDB file, formatted according 
                to the PDB specification.

        Returns:
            tuple[str, str, str, str, int, str, float, float, float, float, 
            float, float, str, str, str]: A tuple containing atom information, 
            including:
                - Atom name
                - Alternate location indicator
                - Residue name
                - Chain identifier
                - Residue sequence number
                - Insertion code
                - X, Y, and Z coordinates
                - Occupancy and temperature factor
                - Segment identifier, element symbol, and charge
        """
        return (
            line[:6].strip(),
            #int(line[6:11].strip()),
            line[12:16].strip(),        # Atom name
            line[16:17].strip(),        # Alternate location indicator
            line[17:20].strip() ,       # Residue name
            line[21:22].strip(),        # Chain identifier
            int(line[22:26].strip()),   # Residue sequence number
            line[26:27].strip(),        # Insertion code
            float(line[30:38].strip()), # X coordinate
            float(line[38:46].strip()), # Y coordinate
            float(line[46:54].strip()), # Z coordinate
            float(line[54:60].strip()) if line[54:60].strip() != "" else 1.0, # Occupancy
            float(line[60:66].strip()) if line[60:66].strip() != "" else 0.0, # Temperature factor
            line[72:76].strip(),        # Segment identifier
            line[76:78].strip(),        # Element symbol
            line[78:80].strip()         # Charge on the atom
        )
    
    columns = ["record_name", "name", "alt", "resn", "chain", "resi", "insertion", "x", "y", "z", "occupancy", "b", "segi", "e", "q"]
    data = [read_format(line) for line in lines if line[:6].strip() == "ATOM" or not atom_only and line[:6].strip() == "HETATM"]
    df = pd.DataFrame(data)
    df.columns = columns
    df.index.name = "atom_id"
    df.index += 1
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

    df = lines2df(lines)

    if atom_only:
        df = df.query("record_name == 'ATOM'")

    return df