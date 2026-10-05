import sys
import warnings
from pathlib import Path
from urllib.request import urlopen
from uuid import uuid4

from .core import Traj

__all__ = ["fetch_PDB"]

def fetch_PDB(pdb:str) -> Traj:
    url = f"https://files.rcsb.org/download/{pdb.lower()}.pdb"
    response = urlopen(url)
    txt = response.read()
    lines = (txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.decode("ascii"))

    tempfile = Path(__file__).parent / f"{uuid4()}.pdb"
    tempfile.write_text(lines)

    # Note for future me:
    # in_memory=True is required to be able to delete the file after reading
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        traj = Traj(tempfile, in_memory=True)
    tempfile.unlink() # Otherwize, this raises a PermissionError

    return traj

if __name__ == "__main__":
    traj = fetch_PDB("3ein")
