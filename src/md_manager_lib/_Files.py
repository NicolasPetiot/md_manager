from ._Traj import Traj, TrajParams

from dataclasses import asdict

__all__ = [
     "GRO",
     "XYZ",
     "XTC",
]
class GRO(Traj):
    default_params = asdict(TrajParams(
        return_record_name = False,
        return_alt = False,
        return_chain = False,
        return_icode = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False,
        return_q = False,
        return_type = False
    ))


class XYZ(Traj):
        default_params = asdict(TrajParams(
        # Topology Records:
        return_record_name = False,
        return_alt = False,
        return_resn = False,
        return_chain = False,
        return_resi = False,
        return_icode = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False, 
        return_e = False,
        return_q = False,
        return_m = False,
        return_type = False,

        # Data Record:
        return_v = False,
    ))
        
class XTC(Traj):
    default_params = asdict(TrajParams(
        # Topology Records:
        return_record_name = False,
        return_alt = False,
        return_chain = False,
        return_icode = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False, 
        return_e = False,
        return_m = False,
        return_type = False,

        # Data Record:
        return_v = False
    ))
