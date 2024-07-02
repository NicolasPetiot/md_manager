from .._CGA import chain_theta_angles, chain_gamma_angles
from .._params import DF_COLUMNS

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

__all__ = ["predict_alpha_helix", "predict_beta_sheets"]

def predict_alpha_helix(df:pd.DataFrame, CA_only = False):
    """
    Returns a Series indicating if the atom bellongs to an alpha helix.

    The input CA must only contains ' CA ' atoms. 
    """
    if CA_only:
        CA = df
    else :
        CA = df.query(f"{DF_COLUMNS[1]} == 'CA'")

    helix = pd.Series(False, index = CA.index)

    # CUTABI parameters
    theta_min, theta_max = (80.0, 105.0) # threshold for theta values
    gamma_min, gamma_max = (30.0,  80.0) # threshold for gamma values

    for _, chain in CA.groupby(DF_COLUMNS[4]):
        theta = chain_theta_angles(chain)
        gamma = chain_gamma_angles(chain)

        theta.index = chain.index
        gamma.index = chain.index

        # control not consecutive residues :
        # TODO 

        theta_criterion = (theta > theta_min) & (theta < theta_max)
        gamma_criterion = (gamma > gamma_min) & (gamma < gamma_max)
        tmp = pd.DataFrame([theta_criterion, gamma_criterion]).T

        for win in tmp.rolling(4):
            if win.Theta.all() & win.Gamma[1:-1].all():
                helix[win.index] = True

    # extend CA to all atoms:
    if not CA_only:
        # get indexes and re-build sheet Series with full atom description
        idx, = np.where(helix)
        helix = pd.Series(False, index=df.index)

        # get chain and resi of the sheets:
        chains = CA.iloc[idx].chain
        res_idx = CA.iloc[idx].resi

        # dictionary of atom indices 
        groups  = df.groupby([DF_COLUMNS[4], DF_COLUMNS[5]]).indices
        for chain, resi in zip(chains, res_idx):
            atom_ids = groups[(chain, resi)]
            helix.iloc[atom_ids] = True

    return helix

def predict_beta_sheets(df:pd.DataFrame, CA_only = False):
    """
    Returns a Series indicating if the atom bellongs to an beta sheet.

    The input CA must only contains ' CA ' atoms.
    """
    if CA_only:
        CA = df
    else :
        CA = df.query(f"{DF_COLUMNS[1]} == 'CA'")

    sheet = pd.Series(False, index = CA.index)

    # CUTABI parameters :
    theta_min, theta_max = (100.0, 155.0) # threshold for theta values
    gamma_lim = 80.0                      # threshold for abs(gamma) values

    contact_threshold  = 5.5              # threshold for K;I & K+1;I+-1 distances
    contact_threshold2 = 6.8              # threshold for K+1;I+-2 distances

    angle_criterion = pd.Series(False, index = CA.index)
    for _, chain in CA.groupby(DF_COLUMNS[4]):
        theta = chain_theta_angles(chain)
        gamma = chain_gamma_angles(chain)

        theta.index = chain.index
        gamma.index = chain.index

        # control not consecutive residues :
        # TODO 

        theta_criterion = (theta > theta_min) & (theta < theta_max)
        gamma_criterion = gamma.abs() > gamma_lim
        tmp = pd.DataFrame([theta_criterion, gamma_criterion]).T

        for win in tmp.rolling(2):
            if win.Theta.all() & win.Gamma[0:1].all():
                angle_criterion[win.index] = True

    xyz = ["x", "y", "z"]
    inter_atom_distance = distance_matrix(CA[xyz], CA[xyz])
    sheet = pd.Series(False, index = CA.index)

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

    # extend CA to all atoms:
    if not CA_only:
        # get indexes and re-build sheet Series with full atom description
        idx, = np.where(sheet)
        sheet = pd.Series(False, index=df.index)

        # get chain and resi of the sheets:
        chains = CA.iloc[idx].chain
        res_idx = CA.iloc[idx].resi

        # dictionary of atom indices 
        groups  = df.groupby([DF_COLUMNS[4], DF_COLUMNS[5]]).indices
        for chain, resi in zip(chains, res_idx):
            atom_ids = groups[(chain, resi)]
            sheet.iloc[atom_ids] = True

    return sheet