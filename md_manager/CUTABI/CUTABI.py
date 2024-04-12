from .._CGA import chain_theta_angles, chain_gamma_angles

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

__all__ = ["predict_alpha_helix", "predict_beta_sheets"]

def predict_alpha_helix(CA:pd.DataFrame):
    """
    Returns a Series indicating if the atom bellongs to an alpha helix.

    The input CA must only contains ' CA ' atoms. 
    """
    # CUTABI parameters
    theta_min, theta_max = (80.0, 105.0) # threshold for theta values
    gamma_min, gamma_max = (30.0,  80.0) # threshold for gamma values

    alpha = pd.Series(False, index = CA.index, name = "alpha")
    for _, chain in CA.groupby("chain"):
        theta = chain_theta_angles(chain)
        gamma = chain_gamma_angles(chain)

        theta.index = chain.index
        gamma.index = chain.index

        # control not consecutive residues :
        # TODO 

        theta_criterion = (theta > theta_min) & (theta < theta_max)
        gamma_criterion = (gamma > gamma_min) & (gamma < gamma_max)
        df = pd.DataFrame([theta_criterion, gamma_criterion]).T

        for win in df.rolling(4):
            if win.Theta.all() & win.Gamma[1:-1].all():
                alpha[win.index] = True

    return alpha

def predict_beta_sheets(CA:pd.DataFrame):
    """
    Returns a Series indicating if the atom bellongs to an beta sheet.

    The input CA must only contains ' CA ' atoms.
    """
    # CUTABI parameters :
    theta_min, theta_max = (100.0, 155.0) # threshold for theta values
    gamma_lim = 80.0                      # threshold for abs(gamma) values

    contact_threshold  = 5.5              # threshold for K;I & K+1;I+-1 distances
    contact_threshold2 = 6.8              # threshold for K+1;I+-2 distances

    tmp = pd.Series(False, index = CA.index)
    for _, chain in CA.groupby("chain"):
        theta = chain_theta_angles(chain)
        gamma = chain_gamma_angles(chain)

        theta.index = chain.index
        gamma.index = chain.index

        # control not consecutive residues :
        # TODO 

        theta_criterion = (theta > theta_min) & (theta < theta_max)
        gamma_criterion = gamma.abs() > gamma_lim
        df = pd.DataFrame([theta_criterion, gamma_criterion]).T

        for win in df.rolling(2):
            if win.Theta.all() & win.Gamma[0:1].all():
                tmp[win.index] = True

    xyz = ["x", "y", "z"]
    inter_atom_distance = distance_matrix(CA[xyz], CA[xyz])
    beta = pd.Series(False, index = CA.index)

    # Parallel sheet detection :
    test1 = inter_atom_distance[:-1, :-2] < contact_threshold  # K  -I   criterion 
    test2 = inter_atom_distance[1:, 1:-1] < contact_threshold  # K+1-I+1 criterion
    test3 = inter_atom_distance[1:,2:]    < contact_threshold2 # K+1-I+2 criterion 

    distance_criterion = test1 & test2 & test3
    I, K = np.where(distance_criterion)
    for i, k in zip(I, K):
        if k > i+2:
            idx = [k, k+1, i, i+1]
            if tmp.iloc[idx].all():
                beta.iloc[idx] = True

    # Anti-parallel sheet detection :
    test1 = inter_atom_distance[:-1, 2:] < contact_threshold # K - I criterion
    test3 = inter_atom_distance[1:,:-2] < contact_threshold2 # K+1-I-2 criterion 

    distance_criterion = test1 & test2 & test3
    I, K = np.where(distance_criterion)
    I += 2 # because test1[0, 0] -> k = 0, i = 2
    for i, k in zip(I, K):
        if k > i+2:
            idx = [k, k+1, i, i+1]
            if tmp.iloc[idx].all():
                beta.iloc[idx] = True

    return beta

