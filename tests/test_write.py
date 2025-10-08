import pytest

import md_manager as md
import pandas as pd
import numpy as np

formats = ["pdb", "xyz"]  # + ["gro"] -> no writer...

@pytest.mark.parametrize("format", formats)
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_save(format):
    assert can_save(format)

def save(format="pdb"):
    xyz = ["x", "y", "z"]
    pos = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
    ])

    Natom, _ = pos.shape
    Nframe = 3

    df = pd.DataFrame(pos, columns=xyz)
    df["name"] = "CA"
    df["resi"] = range(1, Natom+1)
    df["chain"] = "A"
    df["resn"] = "GLY"
    df["alt"] = ""
    df["occupancy"] = 1.0
    df["b"] = 0.0
    df["e"] = "C"
    df["record_name"] = "ATOM"


    new = md.Traj.from_df(df, Nframe=Nframe)
    vec = np.array([0.5, 0.5, 0.5])
    for i in range(1, Nframe):
        pos += vec
        df[xyz] = pos
        new[i] = df

    new.save(f"test_save.{format}")

def can_save(format:str) -> bool:
    print(f"Will save in {format} format")
    try :
        save(format)
        return True

    except Exception as e:
        print(f"Got error {e}")
        return False


