from MDAnalysis import Universe
from urllib.request import urlopen

from pathlib import Path
import sys
from tempfile import NamedTemporaryFile

__all__ = ["fetch_PDB"]

def fetch_PDB(pdb:str) -> Universe:
    url = f"https://files.rcsb.org/download/{pdb.lower()}.pdb"
    response = urlopen(url)

    prefix = Path(__file__).parent

    with NamedTemporaryFile("w", prefix=str(prefix), suffix=".pdb") as tmp:
        txt = response.read()
        lines = (txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.decode("ascii"))

        tmp.write(lines)
        traj = Universe(prefix / tmp.name)

    return traj
