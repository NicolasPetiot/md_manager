import pandas as pd
import numpy as np

__all__ = ["center_of_mass", "shift", "rotate"]

def center_of_mass(df:pd.DataFrame) -> pd.Series:
    com = df.m @ df[["x", "y", "z"]] / df.m.sum()
    com.name = "com"
    return com

def shift(df:pd.DataFrame, vec:pd.Series) -> pd.DataFrame:
    xyz = ["x", "y", "z"]
    df[xyz] = df[xyz].apply(lambda s: s+vec, axis = 1)
    return df

def rotate(df:pd.DataFrame, origin: np.ndarray, axis: np.ndarray, angle: float) -> pd.DataFrame:
    xyz = ["x", "y", "z"]
    rot = rotation_matrix_from_axis_angle(axis, angle)
    def rotate_row(row):
        v = row[xyz].values - origin
        rotated = rot @ v + origin
        return pd.Series(rotated, index=xyz)
    
    df[xyz] = df.apply(rotate_row, axis=1)
    return df
    
def rotation_matrix_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """
    Create a 3D rotation matrix from an axis and angle using Rodrigues' formula.
    """
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1 - c

    R = np.array([
        [c + x*x*C,     x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s,   c + y*y*C,   y*z*C - x*s],
        [z*x*C - y*s,   z*y*C + x*s, c + z*z*C]
    ])
    return R