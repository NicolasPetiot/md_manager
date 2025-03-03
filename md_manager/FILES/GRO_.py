from .Traj_ import Traj
from MDAnalysis.coordinates.GRO import GROWriter

class GRO(Traj):
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
        return_q = False,
        return_m = False,

        # Data records
        return_velocity = False,
        #return_forces = False # TODO

        # Secondary structures:
        return_alpha = False,
        return_beta = False,
        ss_method = "CUTABI"
    )

    def save(self, filename:str):
        with GROWriter(filename) as gro:
            gro.write(self.atoms)