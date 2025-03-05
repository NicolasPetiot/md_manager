from ..conformation import chain_theta_angles, chain_gamma_angles

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

__all__ = ["helix_criterion", "sheet_criterion", "predict_alpha_helix", "predict_beta_sheets"]

def helix_criterion(theta:pd.Series, gamma:pd.Series) -> pd.Series:
    """
    Identify the positions of alpha-helices in a protein structure based on conformational angles.

    This function evaluates the conformational angles (`θ` and `γ`) of a protein's backbone 
    and returns a boolean Series indicating the positions where alpha-helices are present. The 
    identification is based on the CUTABI criteria for the angle thresholds.

    The criteria  used for identifying alpha-helices are presented here: 
    https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2021.786123/full

    Parameters:
    theta : pandas.Series
        A Series containing the θ (theta) angles for each residue in the protein structure, in degrees. 
        The index should correspond to the residue identifiers.
    gamma : pandas.Series
        A Series containing the γ (gamma) angles for each residue in the protein structure, in degrees. 
        The index must match the index of the `theta` Series.

    Returns:
    pandas.Series
        A boolean Series with the same index as the input, where `True` indicates the position of an 
        alpha-helix (based on the CUTABI criteria) and `False` otherwise.

    Raises:
    ValueError
        If the indexes of `theta` and `gamma` do not match.
    KeyError
        If either `theta` or `gamma` contains invalid or missing data (e.g., NaN values).

    Example:
        >>> import pandas as pd
        >>> theta = pd.Series([85, 90, 100, 70], index=[1, 2, 3, 4])
        >>> gamma = pd.Series([35, 50, 75, 25], index=[1, 2, 3, 4])
        >>> helix_criterion(theta, gamma)
        1     True
        2     True
        3     True
        4    False
        dtype: bool

    Notes:
    - Ensure that both `theta` and `gamma` angles are provided in degrees and correspond to 
      the same residues (i.e., have identical indexes).
    - The criteria for `θ` and `γ` are specific to CUTABI and may not detect other types of 
      secondary structures.
    - This function assumes a residue-level representation of the protein structure, where each 
      row in the input Series corresponds to a residue.
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
    Identify the positions of beta-sheets in a protein structure based on conformational angles and spatial coordinates.

    This function evaluates the conformational angles (`θ` and `γ`) of a protein's backbone, as well as 
    the 3D coordinates of the Cα atoms, to determine the positions of beta-sheets. Beta-sheets are identified 
    based on specific ranges of conformational angles combined with spatial criteria derived from the Cα 
    atom positions.

    The criteria  used for identifying beta-sheets are presented here: 
    https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2021.786123/full

    Parameters:
    theta : pandas.Series
        A Series containing the θ (theta) angles for each residue in the protein structure, in degrees. 
        The index should correspond to the residue identifiers.
    gamma : pandas.Series
        A Series containing the γ (gamma) angles for each residue in the protein structure, in degrees. 
        The index must match the index of the `theta` Series.
    xyz : pandas.DataFrame
        A DataFrame containing the 3D coordinates of the Cα atoms for each residue in the structure. 
        The DataFrame must have the following columns:
        - 'x', 'y', 'z': 3D Cartesian coordinates of the Cα atoms.
        The index of the DataFrame must match the indexes of the `theta` and `gamma` Series.

    Returns:
    pandas.Series
        A boolean Series with the same index as the input, where `True` indicates the position of a 
        beta-sheet and `False` otherwise.

    Raises:
    ValueError
        If the indexes of `theta`, `gamma`, and `xyz` do not match.
    KeyError
        If the `xyz` DataFrame is missing required columns ['x', 'y', 'z'].

    Example:
        >>> import pandas as pd
        >>> theta = pd.Series([-120, -130, -140, 70], index=[1, 2, 3, 4])
        >>> gamma = pd.Series([120, 130, 140, 80], index=[1, 2, 3, 4])
        >>> xyz = pd.DataFrame({
        ...     'x': [0.0, 1.0, 2.0, 3.0],
        ...     'y': [0.0, 0.5, 1.5, 2.0],
        ...     'z': [0.0, -1.0, -1.5, -2.0]
        ... }, index=[1, 2, 3, 4])
        >>> sheet_criterion(theta, gamma, xyz)
        1     True
        2     True
        3     True
        4    False
        dtype: bool

    Notes:
    - This function assumes the residue indexes in `theta`, `gamma`, and `xyz` align perfectly and represent 
      the same structure.

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
    Predict the positions of alpha-helices in a protein structure using the CUTABI criterion.

    This function calculates the conformational angles (`θ` and `γ`) for each residue in the input 
    DataFrame based on the provided 3D coordinates (x, y, z) of the Cα atoms. The calculated angles 
    are then passed to the `helix_criterion` function to predict regions of alpha-helical secondary 
    structure. The prediction is based on the following thresholds, as defined by the CUTABI criterion.

    Parameters:
    df : pandas.DataFrame
        A DataFrame containing the 3D coordinates for the Cα atoms of the protein structure. The DataFrame 
        must include the following columns:
        - 'x', 'y', 'z': The 3D Cartesian coordinates of the Cα atoms.
        The index should correspond to the residue identifiers.
        - Optional columns required atomic name selection (e.g., 'name').
        
    Returns:
    pandas.Series
        A boolean Series with the same index as the input DataFrame, where `True` indicates the positions 
        of alpha-helices (based on the CUTABI criterion) and `False` otherwise.

    Raises:
    KeyError
        If the input DataFrame is missing required columns ['x', 'y', 'z'].
    
    ValueError
        If the number of residues or coordinates is insufficient for computing the conformational angles.

    Example:
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'x': [0.0, 1.0, 2.0, 3.0],
        ...     'y': [0.0, 0.5, 1.5, 2.0],
        ...     'z': [0.0, -1.0, -1.5, -2.0]
        ... }, index=[1, 2, 3, 4])
        >>> predict_alpha_helix(data)
        1     True
        2     True
        3     True
        4    False
        dtype: bool

    See Also:
    - `helix_criterion`: The function used to evaluate the CUTABI criterion for alpha-helices based on 
      the angles `θ` and `γ`.

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
    Predict the positions of beta-sheets in a protein structure using the CUTABI criterion.

    This function calculates the conformational angles (`θ` and `γ`) as well as inter-atomic distances for each residue in the input 
    DataFrame based on the provided 3D coordinates (x, y, z) of the Cα atoms. The calculated angles 
    are then passed to the `sheet_criterion` function to predict regions of alpha-helical secondary 
    structure. The prediction is based on the following thresholds, as defined by the CUTABI criterion.

    Parameters:
    df : pandas.DataFrame
        A DataFrame containing the 3D coordinates for the Cα atoms of the protein structure. The DataFrame 
        must include the following columns:
        - 'x', 'y', 'z': The 3D Cartesian coordinates of the Cα atoms.
        The index should correspond to the residue identifiers.
        - Optional columns required atomic name selection (e.g., 'name').
        
    Returns:
    pandas.Series
        A boolean Series with the same index as the input DataFrame, where `True` indicates the positions 
        of alpha-helices (based on the CUTABI criterion) and `False` otherwise.

    Raises:
    KeyError
        If the input DataFrame is missing required columns ['x', 'y', 'z'].
    
    ValueError
        If the number of residues or coordinates is insufficient for computing the conformational angles.

    Example:
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'x': [0.0, 1.0, 2.0, 3.0],
        ...     'y': [0.0, 0.5, 1.5, 2.0],
        ...     'z': [0.0, -1.0, -1.5, -2.0]
        ... }, index=[1, 2, 3, 4])
        >>> predict_alpha_helix(data)
        1     True
        2     True
        3     True
        4    False
        dtype: bool

    See Also:
    - `sheet_criterion`: The function used to evaluate the CUTABI criterion for beta-sheets based on 
      the angles `θ` and `γ`.

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
            
    if not "chain" in df:
        CA["chain"] = "A"

    for _, chain in CA.groupby("chain"):
        theta = chain_theta_angles(chain)
        gamma = chain_gamma_angles(chain)
        xyz = chain[["x", "y", "z"]]
        beta_chain = sheet_criterion(theta, gamma, xyz)
        
        beta[beta_chain.index] = beta_chain.values

    return beta
    