import md_manager as md
from md_manager.angles import backbone_theta_gamma

from MDAnalysisTests.datafiles import PDB

traj = md.Traj(PDB)

print(len(traj))
df = backbone_theta_gamma(traj)
print(df)
