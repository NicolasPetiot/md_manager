from .Traj_ import Traj
from MDAnalysis.coordinates.XTC import XTCWriter

class XTC(Traj):    
    default_params = dict(
        # Read/Write mode:
        mode = "r",

        # Topology records:
        return_topology = True,     
        return_alt = False,
        return_chain = False,
        return_occupancy = False,
        return_b = False,
        return_segi = False,
        return_e = False,
        return_q = True,
        return_m = False,

        # Data records
        return_velocity = False,
        #return_forces = False # TODO

        # Secondary structures:
        return_alpha = False,
        return_beta = False,
        ss_method = "CUTABI"
    )

    def save(self, filename:str, multiframe=False):
        with XTCWriter(filename, multiframe=multiframe) as xyz:
            if multiframe:
                self = iter(self)
                for u in self.frame:
                    xyz.write(u)

            else:
                xyz.write(self.atoms)