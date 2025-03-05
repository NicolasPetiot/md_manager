from .Traj_ import Traj
from .PDB_ import PDB
from .GRO_ import GRO
from .XYZ_ import XYZ
from .XTC_ import XTC

import pandas as pd
from os.path import splitext

known_format = {
    ".pdb": PDB,
    ".gro": GRO,
    ".xyz": XYZ,
    ".xtc": XTC
}

def get_traj(*args, **kwargs) -> Traj:    
    filenames = [arg for arg in args if isinstance(arg, str)]# [topology, coordinate] or [topology]
    if len(filenames) == 1:
        ext = splitext(filenames[0])
    
    else:
        ext = splitext(filenames[1])

    ext = ext[-1]
    if ext in known_format:
        traj = known_format[ext](*args, **kwargs)

    else:
        traj = Traj(*args, **kwargs)

    return traj

def load(*args, **kwargs) -> pd.DataFrame:
    traj = get_traj(*args, **kwargs)
    return traj.load()

def save(*args, **kwargs) -> None:
    kwargs["mode"] = "w"
    traj = get_traj(*args, **kwargs)

    filenames = [arg for arg in args if isinstance(arg, str)]
    traj.save(*filenames)
