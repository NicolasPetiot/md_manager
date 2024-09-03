from ..df_operations import chain_theta_angles, chain_gamma_angles

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

__all__ = ["predict_alpha_helix", "predict_beta_sheets"]

def predict_alpha_helix(df:pd.DataFrame, CA_only = False):
    """
    Returns a Series indicating if the atom bellongs to an alpha helix.
    """
    if CA_only:
        CA = df
    else :
        CA = df.query(f"name == 'CA'")

    helix = pd.Series(False, index = CA.index)

    # CUTABI parameters
    theta_min, theta_max = (80.0, 105.0) # threshold for theta values
    gamma_min, gamma_max = (30.0,  80.0) # threshold for gamma values

    for _, chain in CA.groupby("chain"):
        theta = chain_theta_angles(chain)
        gamma = chain_gamma_angles(chain)

        theta_criterion = (theta > theta_min) & (theta < theta_max)
        gamma_criterion = (gamma > gamma_min) & (gamma < gamma_max)
        tmp = pd.DataFrame([theta_criterion, gamma_criterion]).T

        for win in tmp.rolling(4):
            if win.Theta.all() & win.Gamma[1:-1].all():
                helix[win.index] = True

    return helix

def predict_beta_sheets(df:pd.DataFrame, CA_only = False):
    """
    Returns a Series indicating if the atom bellongs to an beta sheet.

    The input CA must only contains ' CA ' atoms.
    """
    if CA_only:
        CA = df
    else :
        CA = df.query(f"name == 'CA'")

    sheet = pd.Series(False, index = CA.index)

    # CUTABI parameters :
    theta_min, theta_max = (100.0, 155.0) # threshold for theta values
    gamma_lim = 80.0                      # threshold for abs(gamma) values

    contact_threshold  = 5.5              # threshold for K;I & K+1;I+-1 distances
    contact_threshold2 = 6.8              # threshold for K+1;I+-2 distances

    angle_criterion = pd.Series(False, index = CA.index)
    for _, chain in CA.groupby("chain"):
        theta = chain_theta_angles(chain)
        gamma = chain_gamma_angles(chain)

        theta_criterion = (theta > theta_min) & (theta < theta_max)
        gamma_criterion = gamma.abs() > gamma_lim
        tmp = pd.DataFrame([theta_criterion, gamma_criterion]).T

        for win in tmp.rolling(2):
            if win.Theta.all() & win.Gamma[0:1].all():
                angle_criterion[win.index] = True

    xyz = ["x", "y", "z"]
    inter_atom_distance = distance_matrix(CA[xyz], CA[xyz])

    # Parallel sheet detection :
    test1 = inter_atom_distance[:-1, :-2] < contact_threshold  # K  -I   criterion 
    test2 = inter_atom_distance[1:, 1:-1] < contact_threshold  # K+1-I+1 criterion
    test3 = inter_atom_distance[1:,2:]    < contact_threshold2 # K+1-I+2 criterion 

    distance_criterion = test1 & test2 & test3
    I, K = np.where(distance_criterion)
    for i, k in zip(I, K):
        if k > i+2:
            idx = [k, k+1, i, i+1]
            if angle_criterion.iloc[idx].all():
                sheet.iloc[idx] = True

    # Anti-parallel sheet detection :
    test1 = inter_atom_distance[:-1, 2:] < contact_threshold # K - I criterion
    test3 = inter_atom_distance[1:,:-2] < contact_threshold2 # K+1-I-2 criterion 

    distance_criterion = test1 & test2 & test3
    I, K = np.where(distance_criterion)
    I += 2 # because test1[0, 0] -> k = 0, i = 2
    for i, k in zip(I, K):
        if k > i+2:
            idx = [k, k+1, i, i+1]
            if angle_criterion.iloc[idx].all():
                sheet.iloc[idx] = True

    return sheet