"""Loading-driven model with film coupling and a Danckwerts inlet."""
from importlib import import_module

import numpy as np

PoreModel = import_module("Column_Simulation.src.kinetics.3_PoreDiffusion_ceq_Framework").PoreModel


class Model(PoreModel):
    profile_names = ("salt",)

    def inlet(self, adjacent, feed):
        # v*(cm_in - feed) = Dl*(adjacent - cm_in)/dz
        return self.space.boundary_value(adjacent, feed, self.p.v, self.Dl, self.p.dz)

    def adsorption(self, t, c, q):
        physical_c = np.maximum(c, 0)
        if hasattr(self.iso, "pore_q_eq"):
            equilibrium = self.iso.pore_q_eq(physical_c, **self.environment(t))
        elif self.p.IsothermModel ==1:
            equilibrium = np.copysign(self.iso.q_eq(np.abs(c)), c)
        else:
            equilibrium = self.iso.q_eq(physical_c)
        return self.kb * (equilibrium - q)

    def coupling(self, c, cm):
        return 3 * self.kf / self.radius * (cm - c[:, 0])

