import md_manager as md
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Input trajectory file in PDB format")
    parser.add_argument("--theta-gamma", help = "Output conformation file in FEATHER format", required=False, default=None)
    parser.add_argument("--chi"        , help = "Output conformation file in FEATHER format", required=False, default=None)
    args = parser.parse_args()

    if args.theta_gamma is not None:
        traj = md.PDB(args.input)
        confs = traj_theta_gamma(traj)
        confs.to_feather(args.theta_gamma)

    if args.chi is not None:
        traj = md.PDB(args.input)
        confs = traj_theta_gamma(traj)
        confs.to_feather(args.chi)

def traj_theta_gamma(traj:md.Traj) -> pd.DataFrame:
    confs = []
    for df in traj:
        conf = md.backbone_conformation(df)
        confs.append(conf)

    Nframe = len(confs)
    # Generate MultiIndexed DataFrame:
    chains = df.chain.unique()
    prod = [
        [i for i in range(Nframe)],
        chains,
        [i+1 for i in range(len(conf)//(len(chains)))]
    ]
    idx = pd.MultiIndex.from_product(prod, names = ["frame", "chain", "resi"])

    confs = pd.concat(confs, ignore_index=True)
    confs.index = idx

    return confs

if __name__ == "__main__":
    main()