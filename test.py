from MDAnalysisTests.datafiles import PDB

import md_manager as md
from md_manager.angles import backbone_theta_gamma

traj = md.Traj(PDB)

print(len(traj))
df = backbone_theta_gamma(traj)
print(df)
