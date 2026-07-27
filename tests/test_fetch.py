from MDAnalysis import Universe

import md_manager as md


def test_fetch_PDB(pdb_code = "3ein"):
    u = md.fetch_PDB(pdb_code)
    assert isinstance(u, Universe)
