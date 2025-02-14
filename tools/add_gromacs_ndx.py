import md_manager as md
import numpy as np

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--gro")
    parser.add_argument("-n", "--ndx")
    parser.add_argument("-o", "--out")
    args = parser.parse_args()

    make_ndx(args.gro, args.ndx, args.out)

def make_ndx(gro:str, ndx:str, out:str):
    """
    Creates a new NDX file containing all chain and residue groups.
    """
    # Load GRO:
    box = md.load(gro)
    protein_resn = list(md.ONE_LETTER_CODE.keys())
    Protein = box.query("resn in @protein_resn").copy()

    # Add chains:
    chain = chain_iter()
    Protein["chain"] = next(chain)
    for id in Protein[Protein.resi.diff() < 0].index:
        Protein.loc[id:, "chain"] = next(chain)
    chains = Protein.chain.unique()

    # Load NDX
    lines = open(ndx, "r").readlines()
    groups = [line.split()[1] for line in lines if line.startswith("[")]

    # Add chains:
    for chain in chains:
        if not f"Chain{chain}" in groups:
            idx = Protein.query(f"chain == '{chain}'").index.values
            lines += gen_ndx_lines(name = f"Chain{chain}", idx=idx)

    # One group per residue:
    for (chain, resi), group in Protein.groupby(["chain", "resi"]):
        name = f"r_{resi}_{chain}"
        if not name in groups:
            idx = group.index.values
            lines += gen_ndx_lines(name, idx)

    with open(out, "w") as file:
        file.writelines(lines)

def gen_ndx_lines(name:str, idx:np.ndarray) -> list[str]:
    """
    Generates la list of lines in a Gromacs NDX compatible format based on a group name and an array of index.
    """
    ### To be added in a future md.gromacs_tools module ? ###
    
    print("Generating", name, "indexes")
    lines = ["[ %s ]\n"%name]
    Nlines_full = len(idx) // 15
    for line_idx in idx[:15*Nlines_full].reshape(Nlines_full, 15):
        line_idx = [f"{id:4d}" for id in line_idx]
        lines.append(" ".join(line_idx) + "\n")
    
    # Possible non-full line
    line_idx = [f"{id:4d}" for id in idx[15*Nlines_full:]]
    if len(line_idx) > 0:
        lines.append(" ".join(line_idx) + "\n")

    return lines

def chain_iter():
    for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        yield letter

if __name__ == "__main__":
    main()