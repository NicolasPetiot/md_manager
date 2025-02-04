# Tools:

The scripts presented here are meant to display simple examples of the way I am using MD-manager.

### PDBtraj_conformations:

```shell
python PDBtraj_conformations.py --input traj.pdb --theta-gamma Theta_Gamma.pkl -chi Chis.pkl
```

It is used to iterate over the frames of a PDB trajectory and compute the $\theta$; $\gamma$ and $\chi$ angles.
