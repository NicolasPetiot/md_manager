
import numpy as np
import pandas as pd
import MDAnalysis as mda
from scipy.spatial import distance_matrix

from md_manager.core import universe_to_top
from .angles import angles, dihedral_angles

def predict_secondary_structure(u:mda.Universe, helix:bool = True, sheet:bool = True) -> pd.DataFrame:
    if not helix and not sheet:
        raise ValueError("At least one of helix or sheet must be True")

    df = __get_CUTABI_compatible_df(u)
    if helix:
        df["helix"] = CUTABI_helix_criterion(df)
    if sheet:
        df["sheet"] = CUTABI_sheet_criterion(df)

    expected_cols = {"chain", "resi", "helix", "sheet"}
    cols = [col for col in df.columns if col in expected_cols]
    return df[cols]



def __get_CUTABI_compatible_df(u:mda.Universe):
    ca = u.select_atoms("name CA")
    top = universe_to_top(ca)
    xyz = ["x", "y", "z"]
    top[xyz] = ca.positions

    if chain not in top.columns:
        top["chain"] = "A"

    for _, chain in top.groupby("chain"):
        idx = chain.index[1:-1]
        top.loc[idx, "theta"] = angles(chain[xyz].values)
        idx = chain.index[1:-2]
        top.loc[idx, "gamma"] = dihedral_angles(chain[xyz].values)
    return top

def CUTABI_helix_criterion(df:pd.DataFrame) -> pd.Series:
    helix = pd.Series(False, index = df.index)
    theta = df.theta.values
    gamma = df.gamma.values

    # CUTABI parameters
    theta_min, theta_max = (80.0, 105.0) # threshold for theta values
    gamma_min, gamma_max = (30.0,  80.0) # threshold for gamma values

    theta_is_valid = (theta > theta_min) & (theta < theta_max)
    gamma_is_valid = (gamma > gamma_min) & (gamma < gamma_max)

    theta_gamma_validity = pd.DataFrame({"Theta" : theta_is_valid, "Gamma": gamma_is_valid})
    for win in theta_gamma_validity.rolling(4):
        if win.Theta.all() & win.Gamma[1:-1].all():
            helix[win.index] = True

    return helix

def CUTABI_sheet_criterion(df:pd.DataFrame) -> pd.Series:
    sheet = pd.Series(False, index=df.index)
    theta = df.theta.values
    gamma = df.gamma.values
    pos   = df[["x", "y", "z"]].values

    # CUTABI parameters :
    theta_min, theta_max = (100.0, 155.0) # threshold for theta values
    gamma_lim = 80.0                      # threshold for abs(gamma) values
    contact_threshold  = 5.5              # threshold for K;I & K+1;I+-1 distances
    contact_threshold2 = 6.8              # threshold for K+1;I+-2 distances

    angle_is_valid = pd.Series(False, index = sheet.index)
    theta_is_valid = (theta > theta_min) & (theta < theta_max)
    gamma_is_valid = np.abs(gamma) > gamma_lim

    theta_gamma_validity = pd.DataFrame({"Theta" : theta_is_valid, "Gamma": gamma_is_valid})
    for win in theta_gamma_validity.rolling(2):
        if win.Theta.all() & win.Gamma[0:1].all():
            angle_is_valid[win.index] = True
    inter_atom_distance = distance_matrix(pos, pos)

    # Parallel sheet detection :
    test1 = inter_atom_distance[:-1, :-2] < contact_threshold  # K  -I   criterion
    test2 = inter_atom_distance[1:, 1:-1] < contact_threshold  # K+1-I+1 criterion
    test3 = inter_atom_distance[1:,2:]    < contact_threshold2 # K+1-I+2 criterion
    distance_is_valid = test1 & test2 & test3
    I, K = np.where(distance_is_valid)
    for i, k in zip(I, K):
        if k > i+2:
            idx = [k, k+1, i, i+1]
            if angle_is_valid.iloc[idx].all():
                sheet.iloc[idx] = True

    # Anti-parallel sheet detection :
    test1 = inter_atom_distance[:-1, 2:] < contact_threshold # K - I criterion
    #test2 = inter_atom_distance[1:, 1:-1] < contact_threshold  # K+1-I+1 criterion
    test3 = inter_atom_distance[1:,:-2] < contact_threshold2 # K+1-I-2 criterion
    distance_is_valid = test1 & test2 & test3
    I, K = np.where(distance_is_valid)
    I += 2 # because test1[0, 0] -> k = 0, i = 2
    for i, k in zip(I, K):
        if k > i+2:
            idx = [k, k+1, i, i+1]
            if angle_is_valid.iloc[idx].all():
                sheet.iloc[idx] = True

    return sheet
