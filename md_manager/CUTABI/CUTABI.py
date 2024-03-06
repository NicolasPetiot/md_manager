from .._CGA import theta_angles, gamma_angles

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

__all__ = ["predict_alpha_helix", "predict_beta_sheets"]

def predict_alpha_helix(CA:pd.DataFrame):
    """
    Returns a pandas.Series indicating the position of the predicted alpha helix.

    The input CA must be a pandas.DataFrame build he same way that the `md.PDBfile().read2df()` method and must only contains ' CA ' atoms. 
    If CA does not contains a theta column, the conformational theta angle is computed using `md.theta_angles(CA)` method. 
    If CA does not contains a gamma column, the conformational gamma angle is computed using `md.gamma_angles(CA)` method. 

    The returned pandas.Series has the same index as the input CA and contains integers 1 for CA associated to alpha helix.
    """
    # CUTABI parameters :
    theta_min, theta_max = (80.0, 105.0) # threshold for theta values
    gamma_min, gamma_max = (30.0,  80.0) # threshold for gamma values

    # conformational angles :
    columns = CA.columns
    if not "theta" in columns:
        CA["theta"] = theta_angles(df=CA)
    if not "gamma" in columns:
        CA["gamma"] = gamma_angles(df=CA)

    # CUTABI predictions alpha helix
    alpha = pd.Series(index=CA.index, dtype=int)
    for win in CA.rolling(4):
        sample_t = win.theta
        sample_g = win.gamma[:-1]

        # all the angles in the selected CA must fit the conditions on theta & gamma
        theta_test = ((sample_t > theta_min) & (sample_t < theta_max)).all()
        gamma_test = ((sample_g > gamma_min) & (sample_g < gamma_max)).all()

        if theta_test and gamma_test:
            alpha[win.index] = 1

    return alpha

def predict_beta_sheets(CA:pd.DataFrame):
    """
    Returns a tuple of pandas.Series indicating the position of the predicted beta sheets in the parallel and anti-parallel configuration.

    The input CA must be a pandas.DataFrame build he same way that the `md.PDBfile().read2df()` method and must only contains ' CA ' atoms. 
    If CA does not contains a theta column, the conformational theta angle is computed using `md.theta_angles(CA)` method. 
    If CA does not contains a gamma column, the conformational gamma angle is computed using `md.gamma_angles(CA)` method. 

    The returned pandas.Series has the same index as the input CA and contains integers 1 for CA associated to parallel beta sheets and integer -1 for anti-parallel ones.
    """
    # CUTABI parameters :
    theta_min, theta_max = (100.0, 155.0) # threshold for theta values
    gamma_lim = 80.0                      # threshold for abs(gamma) values

    contact_threshold  = 5.5              # threshold for K;I & K+1;I+-1 distances
    contact_threshold2 = 6.8              # threshold for K+1;I+-2 distances

    # conformational angles :
    columns = CA.columns
    if not "theta" in columns:
        CA["theta"] = theta_angles(df=CA)
    if not "gamma" in columns:
        CA["gamma"] = gamma_angles(df=CA)

    # CUTABI criterion on conformational angles :
    angle_criterion = pd.Series(False, index=CA.index, dtype=bool)
    for win in CA.rolling(2):
        sample_t = win.theta
        sample_g = win.gamma[:-1]

        # all the angles in the selected CA must fit the conditions on theta & gamma
        theta_test = ((sample_t > theta_min) & (sample_t < theta_max)).all()
        gamma_test = (sample_g.abs()> gamma_lim).all()

        if theta_test and gamma_test:
            angle_criterion[win.index] = True

    Nca = len(CA)
    xyz = CA[["x", "y", "z"]]
    distance_inter_CA = distance_matrix(xyz, xyz)

    # CUTABI predictions : parallel beta sheets
    beta_para = pd.Series(index = CA.index, dtype=int)

    # K-I criterion 
    test1 = distance_inter_CA[:-1, :-2] < contact_threshold
    #K+1-I+1 criterion:
    test2 = distance_inter_CA[1:, 1:-1] < contact_threshold
    # K+1-I+2 criterion :
    test3 = distance_inter_CA[1:,2:] < contact_threshold2

    distance_criterion = test1 & test2 & test3
    K, I = np.where(distance_criterion)
    for k, i in zip(K, I):
        if i > k+2:
            # test angle criterion :
            idx = [i, i+1, k, k+1]
            if angle_criterion.iloc[idx].all():
                beta_para.iloc[idx] = 1

    # CUTABI predictions : anti-parallel beta sheets :
    beta_anti = pd.Series(index = CA.index, dtype=int)

    # K-I criterion 
    test1 = distance_inter_CA[:-1, 2:] < contact_threshold
    #K+1-I+1 criterion:
    test2 = distance_inter_CA[1:, 1:-1] < contact_threshold
    # K+1-I+2 criterion :
    test3 = distance_inter_CA[1:,:-2] < contact_threshold2

    distance_criterion = test1 & test2 & test3
    K, I = np.where(distance_criterion)
    I += 2 # because test1[0, 0] -> k = 0, i = 2
    for k, i in zip(K, I):
        if i > k+2 and i+2 < Nca:
            # test angle criterion :
            idx = [i, i+1, k, k+1]
            if angle_criterion.iloc[idx].all():
                beta_anti.iloc[idx] = -1

    return beta_para, beta_anti

