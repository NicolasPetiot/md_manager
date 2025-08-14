from ._Traj import Traj

import pandas as pd

__all__ = ["load", "save"]

def load(filename:str, *args, **kwargs) -> pd.DataFrame:    
    return Traj(filename, *args, **kwargs).load()

def save(df:pd.DataFrame, *filenames, frames = "all", **kwargs):
    Traj.from_df(df).save(*filenames, frames=frames, **kwargs)