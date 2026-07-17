# %%
from md_manager import NMA

import MDAnalysis as mda
from MDAnalysisTests.datafiles import TPR, XTC

u = mda.Universe(TPR, XTC)
ag = u.select_atoms("name CA")

debey = NMA.predict_thermal_factors(ag)
