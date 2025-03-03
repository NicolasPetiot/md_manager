from .Traj_ import Traj
from MDAnalysis.coordinates.PDB import PDBWriter

import pandas as pd
import sys
from urllib.request import urlopen

class PDB(Traj):
    default_params = dict(
        # Read/Write mode:
        mode = "r",

        # Topology records:
        return_topology = True,     
        return_alt = True,
        return_chain = True,
        return_occupancy = True,
        return_b = False,
        return_segi = True,
        return_e = False,
        return_q = False,
        return_m = False,

        # Data records
        return_velocity = False,
        #return_forces = False # TODO

        # Secondary structures:
        return_alpha = False,
        return_beta = False,
        ss_method = "CUTABI"
    )

    def save(self, filename:str, multiframe=False):
        with PDBWriter(filename, multiframe=multiframe) as pdb:
            if multiframe:
                self = iter(self)
                for u in self.frame:
                    pdb.write(u)

            else:
                pdb.write(self.atoms)

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
            #line[:6].strip(),
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

def fetch_PDB(pdb_code:str, atom_only = False) -> pd.DataFrame:
    url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
    response = urlopen(url)

    txt = response.read()
    lines = (txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.decode("ascii")).splitlines()

    data = [
         read_format(line) for line in lines if line.startswith("ATOM") or (not atom_only and line.startswith("HETATM"))
    ]

    df = pd.DataFrame(data, columns = ["record_name", "name", "alt", "resn", "chain", "resi", "insertion", "x", "y", "z", "occupancy", "b", "segi", "e", "q"])
    df.index.name = "atom_id"
    df.index += 1

    return df