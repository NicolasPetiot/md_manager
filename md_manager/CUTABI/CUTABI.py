from ..df_operations import  check_chain_validity, InvalidChainException
from ..conformation import chain_theta_angles, chain_gamma_angles

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

__all__ = ["helix_criterion", "sheet_criterion", "predict_alpha_helix", "predict_beta_sheets"]

def helix_criterion(theta:pd.Series, gamma:pd.Series) -> pd.Series:
    """
    Takes the conformational angles of a structure and returns a Series of boolean that indicates the position of alpha-helices.

    Before calling this method, please make sure that the indexes of both Series are identical.
    """
    helix = pd.Series(False, index = theta.index)

    # CUTABI parameters
    theta_min, theta_max = (80.0, 105.0) # threshold for theta values
    gamma_min, gamma_max = (30.0,  80.0) # threshold for gamma values

    theta_criterion = (theta > theta_min) & (theta < theta_max)
    gamma_criterion = (gamma > gamma_min) & (gamma < gamma_max)
    tmp = pd.DataFrame({"Theta" : theta_criterion, "Gamma": gamma_criterion})

    for win in tmp.rolling(4):
        if win.Theta.all() & win.Gamma[1:-1].all():
            helix[win.index] = True

    return helix

def sheet_criterion(theta:pd.Series, gamma:pd.Series, xyz:pd.DataFrame) -> pd.Series:
    """
    Takes the conformational angles of a structure as well as the coordinates of the CA atoms and returns a Series of boolean that indicates the position of beta-sheets.

    Before calling this method, please make sure that the indexes of both Series and DataFrame are identical.
    """
    sheet = pd.Series(False, index = theta.index)

    # CUTABI parameters :
    theta_min, theta_max = (100.0, 155.0) # threshold for theta values
    gamma_lim = 80.0                      # threshold for abs(gamma) values

    contact_threshold  = 5.5              # threshold for K;I & K+1;I+-1 distances
    contact_threshold2 = 6.8              # threshold for K+1;I+-2 distances

    angle_criterion = pd.Series(False, index = sheet.index)

    theta_criterion = (theta > theta_min) & (theta < theta_max)
    gamma_criterion = gamma.abs() > gamma_lim
    tmp = pd.DataFrame({"Theta" : theta_criterion, "Gamma": gamma_criterion})

    for win in tmp.rolling(2):
        if win.Theta.all() & win.Gamma[0:1].all():
            angle_criterion[win.index] = True

    inter_atom_distance = distance_matrix(xyz, xyz)

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

def predict_alpha_helix(df:pd.DataFrame) -> pd.Series:
    """
    Uses CUTABI criterion to predict the position of alpha-helices.
    """
    alpha = pd.Series(False, index=df.index)

    if not "name" in df:
        CA = df.copy()
    else:
        names = set(df.name)
        names.remove("CA")
        if len(names) > 0:
            CA = df.query("name == 'CA'")
        else:
            CA = df.copy()
            
    
    #check_df_infos(CA)
    if not "chain" in df:
        CA["chain"] = "A"

    for _, chain in CA.groupby("chain"):
        theta = chain_theta_angles(chain)
        gamma = chain_gamma_angles(chain)
        alpha_chain = helix_criterion(theta, gamma)

        alpha[alpha_chain.index] = alpha_chain.values

    return alpha

def predict_beta_sheets(df:pd.DataFrame) -> pd.Series:
    """
    Uses CUTABI criterion to predict the position of alpha-helices.
    """
    beta = pd.Series(False, index=df.index)

    if not "name" in df:
        CA = df.copy()
    else:
        names = set(df.name)
        names.remove("CA")
        if len(names) > 0:
            CA = df.query("name == 'CA'")
        else:
            CA = df.copy()
            
    
    #check_df_infos(CA)
    if not "chain" in df:
        CA["chain"] = "A"

    for _, chain in CA.groupby("chain"):
        theta = chain_theta_angles(chain)
        gamma = chain_gamma_angles(chain)
        xyz = chain[["x", "y", "z"]]
        beta_chain = sheet_criterion(theta, gamma, xyz)
        
        beta[beta_chain.index] = beta_chain.values

    return beta

def check_df_infos(df:pd.DataFrame) -> None:
    if not "chain" in df:
        df["chain"] = "A"

    for _, tmp in df.groupby("chain"):
        check_chain_validity(tmp)
    