from pathlib import Path

import md_manager as md


def test_backbone_theta_gamma():
    path = Path(__file__).parent / "testfile.pdb"
    traj = md.Traj(path)

    df = md.backbone_theta_gamma(traj)
    assert df is not None

def test_side_chain_chis():
    path = Path(__file__).parent / "testfile.pdb"
    traj = md.Traj(path)

    df = md.side_chain_chis(traj)
    assert df is not None
