import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from MDAnalysis.tests.datafiles import PDB, TPR, XTC

import md_manager as md


@pytest.mark.parametrize("top, trj", [
    (PDB, None),
    (TPR, XTC)
])
def test_can_create_Traj(top:str, trj:str|None):
    traj = md.Traj(TPR, XTC)
    assert type(traj) == md.Traj

@pytest.mark.parametrize("top, trj", [
        (PDB, None),
        (TPR, XTC)
    ])
def test_can_use_load(top:str, trj:str|None):
    path = Path(__file__).parent / "testfile.pdb"
    traj = md.load(path)
    assert isinstance(traj, pd.DataFrame)

def test_len():
    traj = md.Traj(TPR, XTC)
    assert len(traj) == 10

def test_can_get_item():
    traj = md.Traj(TPR, XTC)
    df = traj[0]
    assert isinstance(df, pd.DataFrame)

def test_item_has_position():
    traj = md.Traj(TPR, XTC)
    df = traj[0]

    pos = df[["x", "y", "z"]].values
    assert isinstance(pos, np.ndarray)
