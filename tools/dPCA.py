import md_manager as md
import pandas as pd
import numpy as np
from scipy.linalg import eigh, ishermitian

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--pkl", help="Theta Gamma angles in PKL format")
    args = parser.parse_args()

    df = pd.read_pickle(args.pkl)
    cov, res_idx = get_covariance_matrix(df.theta, df.gamma)

    if not ishermitian(cov):
        raise ValueError("Non-symmetrical covariance matrix")    
    
    eigenvals, eigenvecs = eigh(cov)

    save_eigenvalues(eigenvals=eigenvals)
    save_eigenvectors(eigenvecs, res_idx)

def save_eigenvectors(eigenvecs:np.ndarray, res_idx:pd.Index, csv_out="eigenvecs.csv"):
    """
    
    """
    Nres = len(res_idx)
    _, Nmode = eigenvecs.shape

    eigenvecs = eigenvecs.T.reshape(Nmode * Nres, 3)
    idx = [(mode, chain, resi) for mode in reversed(range(Nmode)) for chain, resi in res_idx.values]
    idx = pd.MultiIndex.from_tuples(idx, names=["mode", "chain", "resi"])

    df = pd.DataFrame(eigenvecs, index=idx, columns=["x", "y", "z"])
    df = df.sort_index()
    df.to_csv(csv_out)

def save_eigenvalues(eigenvals, csv_out = "eigenvals.csv"):
    """
    
    """
    Nmode = len(eigenvals)
    weights = eigenvals / eigenvals.sum()

    df = pd.DataFrame(dict(
        eigenvals = eigenvals,
        weights = weights
    ))

    df.index = pd.Index(reversed(range(Nmode)), name="mode")
    df = df.sort_index()
    df.to_csv(csv_out)


def get_covariance_matrix(theta:pd.Series, gamma:pd.Series, deg = True) -> tuple[np.ndarray, pd.Index]:
    """
    
    """
    if deg:
        theta = np.deg2rad(theta)
        gamma = np.deg2rad(gamma)

    # Compute conformation vectors associated to each CA
    x = np.sin(theta) * np.cos(gamma)
    y = np.sin(theta) * np.sin(gamma)
    z = np.cos(theta)
    df = pd.DataFrame(dict(x=x, y=y, z=z)).dropna()
    res_index = df.loc[0].index

    Nframe = len(df.groupby("frame"))
    Nres = len(res_index)

    # Compute covariance matrix:
    data = df.values.reshape(3*Nres, Nframe)
    cov = np.cov(data)

    return cov, res_index

if __name__ == "__main__":
    main()