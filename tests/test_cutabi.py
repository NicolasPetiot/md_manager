import numpy as np
import pandas as pd

import md_manager as md
from md_manager.cutabi import predict_secondary_structure


def test_helix_detection():
    trj = _gen_idealistic_structures("helix")
    cutabi = predict_secondary_structure(trj)

    assert cutabi.helix.iloc[1:-1].all()
    assert not cutabi.sheet.any()


def test_sheet_detection():
    trj = _gen_idealistic_structures("sheet")
    cutabi = predict_secondary_structure(trj)

    assert cutabi.groupby("chain").sheet.apply(lambda s: s.iloc[2:-1].all()).all()
    assert not cutabi.helix.any()


def test_null_detection():
    trj = _gen_idealistic_structures("none")
    cutabi = predict_secondary_structure(trj)

    assert not cutabi.helix.any()
    assert not cutabi.sheet.any()


HELIX = np.array(
    [
        [2.30, 0.00, 0.0],
        [-0.40, 2.27, 1.5],
        [-2.16, -0.79, 3.0],
        [1.15, -1.99, 4.5],
        [1.76, 1.48, 6.0],
        [-1.76, 1.48, 7.5],
        [-1.15, -1.99, 9.0],
        [2.16, -0.79, 10.5],
        [0.40, 2.27, 12.0],
        [-2.30, 0.00, 13.5],
        [0.40, -2.27, 15.0],
        [2.16, 0.79, 16.5],
    ]
)

FIRST_STAND = np.array(
    [
        [0.0, 0.0, 0.94],
        [3.3, 0.0, -0.94],
        [6.6, 0.0, 0.94],
        [9.9, 0.0, -0.94],
        [13.2, 0.0, 0.94],
        [16.5, 0.0, -0.94],
    ]
)

SECOND_STAND = np.array(
    [
        [16.5, 4.8, -0.94],
        [13.2, 4.8, 0.94],
        [9.9, 4.8, -0.94],
        [6.6, 4.8, 0.94],
        [3.3, 4.8, -0.94],
        [0.0, 4.8, 0.94],
    ]
)


def _gen_idealistic_structures(which: str) -> pd.DataFrame:
    match which.lower():
        case "helix":
            df = pd.DataFrame(HELIX, columns=["x", "y", "z"])  # pyright: ignore[reportArgumentType]

        case "sheet":
            df = pd.concat(
                [
                    pd.DataFrame(FIRST_STAND, columns=["x", "y", "z"]),  # pyright: ignore[reportArgumentType]
                    pd.DataFrame(SECOND_STAND, columns=["x", "y", "z"]),  # pyright: ignore[reportArgumentType]
                ],
                ignore_index=True,
            )
            idx = df.index
            df["resi"] = range(1, len(df) + 1)
            df.loc[idx[:6], "chain"] = "A"
            df.loc[idx[-6:], "chain"] = "B"

        case "none":
            df = pd.DataFrame(np.random.randn(6, 3), columns=["x", "y", "z"])  # pyright: ignore[reportArgumentType]

        case _:
            raise ValueError(f"Unknown structure name '{which}'")

    df["name"] = "CA"
    return md.Traj.from_df(df)  # pyright: ignore[reportReturnType]
