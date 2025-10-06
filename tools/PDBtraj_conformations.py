import md_manager as md
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help = "Input trajectory file in PDB format")
    parser.add_argument("--theta-gamma", help = "Output conformation file in PKL format", required=False, default=None)
    parser.add_argument("--chi"        , help = "Output conformation file in PKL format", required=False, default=None)
    args = parser.parse_args()

    if args.theta_gamma is not None:
        traj = md.PDB(args.input)
        confs = traj_confs(traj, conf_function=md.backbone_conformation)
        confs.to_pickle(args.theta_gamma)

    if args.chi is not None:
        traj = md.PDB(args.input)
        confs = traj_confs(traj, conf_function=md.side_chain_conformation)
        confs.to_pickle(args.chi)

def traj_confs(traj:md.Traj, conf_function) -> pd.DataFrame:
    confs = []
    for i, df in enumerate(traj):
        conf = conf_function(df)
        conf["frame"] = i
        conf = conf.set_index("frame", append=True).reorder_levels([2, 0, 1])
        confs.append(conf)
    
    return pd.concat(confs)

if __name__ == "__main__":
    main()