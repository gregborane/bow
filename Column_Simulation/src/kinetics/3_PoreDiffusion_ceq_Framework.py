"""Concentration-driven pore kinetics and shared pore state/time integration.

Python uses species-last arrays and C-order packing, not MATLAB linear
indices. Independent states are c[:,1:], q (all radii), and cm[1:-1].
The surface c and column endpoints are algebraic and reconstructed for every
RHS evaluation and output. Surface q is a diagnostic at a zero-storage node.
"""
from importlib import import_module
from types import SimpleNamespace

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse import lil_matrix

from Column_Simulation.src.misc.isotherm import _species_parameter


class PoreModel:
    """Common balances. An injected backend implements the spatial operators
    documented in solvers/1_finite_difference_method.py, including
    particle_interior(c, De, inward_flux) for algebraic surface elimination.
    """
    def __init__(self, inputs, isotherm, space=None):
        self.p, self.iso = inputs, isotherm
        self.nz, self.nr, self.nc = inputs.Nz, inputs.Nr, len(inputs.species)
        if min(self.nz, self.nr) < 3 or self.nc != isotherm.nc:
            raise ValueError("Need Nz and Nr >= 3 and matching species counts")
        # Optional measured/simulated gradient tables, time by axial position.
        # Never infer a spatial gradient from just the feed endpoints.
        self.positions = np.linspace(0, inputs.L, self.nz)
        self.profiles = {}
        if inputs.ElutionType != 1:
            for name in self.profile_names:
                table = getattr(inputs, name + "_profile", None)
                if table is None:
                    raise ValueError(f"Gradient elution requires {name}_profile=(times, z, values)")
                times, z, values = (np.asarray(a, dtype=float) for a in table)
                if (not all(np.isfinite(a).all() for a in (times, z, values))
                        or times.ndim != 1 or z.ndim != 1 or values.shape != (times.size, z.size)
                        or times.size < 2 or z.size < 2 or times.min() > 0
                        or times.max() < inputs.tfinal or z.min() > 0 or z.max() < inputs.L):
                    raise ValueError(f"{name}_profile must span the experiment with finite data")
                self.profiles[name] = RegularGridInterpolator((times, z), values)

        missing = [k for k in ('De', 'kb', 'Dl', 'Cfeed') if getattr(inputs, k, None) is None]
        if missing:
            raise ValueError("Missing pore experiment parameters: " + ", ".join(missing))
        for name in ('De', 'kb', 'Dl', 'Cfeed'):
            setattr(self, name, _species_parameter(getattr(inputs, name), self.nc, name,
                                                    allow_zero=name != 'De'))
        self.kf = _species_parameter([inputs.GetSpeciesResults(s).kf for s in inputs.species],
                                     self.nc, "kf", allow_zero=True)
        self.radius = inputs.Rp  # StoreInput.Rp in µm, diffusivities in cm²/s
        if not (0 < inputs.ep_p <= 1 and inputs.v > 0 and self.radius > 0):
            raise ValueError("Expected positive flow, radius and particle porosity <= 1")
        if space is None:
            backend = import_module("solvers.1_finite_difference_method")
            space = backend.SpatialOperators(np.linspace(0, inputs.L, self.nz),
                                             np.linspace(self.radius, 0, self.nr))
        self.space = space
        self.c_shape = (self.nz, self.nr - 1, self.nc)
        self.q_shape = (self.nz, self.nr, self.nc)
        self.cm_shape = (self.nz - 2, self.nc)

    @staticmethod
    def pack(c, q, cm):
        return np.concatenate((c.ravel(), q.ravel(), cm.ravel()))

    def unpack(self, state):
        split_c = np.prod(self.c_shape)
        split_q = split_c + np.prod(self.q_shape)
        return (state[:split_c].reshape(self.c_shape),
                state[split_c:split_q].reshape(self.q_shape),
                state[split_q:].reshape(self.cm_shape))

    def fields(self, state, feed):
        interior, q, bulk = self.unpack(state)
        cm = np.concatenate((self.inlet(bulk[0], feed)[None], bulk, bulk[-1:]))
        surface = self.space.boundary_value(interior[:, 0], cm, self.kf, self.De, self.p.dr)
        c = np.concatenate((surface[:, None], interior), axis=1)
        return c, q, cm

    def environment(self, t):
        """Return species-broadcastable (Nz, 1, 1) mobile-phase fields."""
        points = np.column_stack((np.full(self.nz, t), self.positions))
        initial = {"salt": self.p.Cso, "pH": self.p.pHo}
        return {name: (self.profiles[name](points).reshape(self.nz, 1, 1)
                       if name in self.profiles else initial[name])
                for name in self.profile_names}

    def rhs(self, t, state, feed):
        c, q, cm = self.fields(state, feed)
        dq = self.adsorption(t, c, q)
        flux = self.kf * (cm - c[:, 0])
        diffusion = self.space.particle_interior(c, self.De, flux)
        dc = (
            diffusion
            - (1.0 - self.p.ep_p) * dq[:, 1:]
            ) / self.p.ep_p        
        dcm = (-self.p.v * self.space.first(cm, self.p.v) + self.Dl * self.space.second(cm)
               - self.p.theta * self.coupling(c, cm))
        return self.pack(dc, dq, dcm[1:-1])

    def sparsity(self):
        """BDF coloring for the local stencil; include competitive species
        and algebraic endpoint dependencies at both ends of the column.
        """
        n = sum(np.prod(s) for s in (self.c_shape, self.q_shape, self.cm_shape))
        ci, qi, mi = self.unpack(np.arange(n))
        pattern = lil_matrix((n, n), dtype=bool)
        for z in range(self.nz):
            bulk = mi[np.clip(z - 1, 0, self.nz - 3)]
            for r in range(self.nr):
                j = max(0, r - 1)
                rows = qi[z, r] if r == 0 else np.r_[qi[z, r], ci[z, j]]
                cols = np.r_[qi[z, r], ci[z, max(0, j-1):j+2].ravel(), bulk]
                pattern[np.ix_(rows, cols)] = True
            if 0 < z < self.nz - 1:
                cols = np.r_[mi[max(0, z-2):z+1].ravel(), ci[z, 0]]
                pattern[np.ix_(bulk, cols)] = True
        return pattern.tocsc()

    def simulate(self):
        """Empty column at t=0. Split BDF integration at the pulse edge."""
        p = self.p
        times = np.asarray(p.tspan)
        initial = self.pack(*(np.zeros(s) for s in (self.c_shape, self.q_shape, self.cm_shape)))
        pattern = self.sparsity() if getattr(self.space, "local_three_point", False) else None
        states = np.empty((len(times), initial.size))
        edges = sorted({0., float(np.clip(p.tinj, 0, p.tfinal)), float(p.tfinal)})
        for start, end in zip(edges, edges[1:]):
            feed = self.Cfeed if start < p.tinj else np.zeros(self.nc)
            mask = (times >= start) & (times <= end)
            evaluation = np.union1d(times[mask], [end])
            solution = solve_ivp(self.rhs, (start, end), initial, method="BDF",
                                 t_eval=evaluation, args=(feed,), jac_sparsity=pattern,
                                 rtol=1e-6, atol=1e-10)
            if not solution.success:
                raise RuntimeError(solution.message)
            states[mask] = solution.y[:, np.searchsorted(evaluation, times[mask])].T
            initial = solution.y[:, -1]
        fields = [self.fields(y, self.Cfeed if (t < p.tinj or p.InputType == 1)
                              else np.zeros(self.nc)) for t, y in zip(times, states)]
        c, q, cm = (np.stack(a) for a in zip(*fields))
        return SimpleNamespace(t=times, c=c, q=q, cm=cm, breakthrough=cm[:, -1])


class Model(PoreModel):
    profile_names = ("salt", "pH")

    def inlet(self, adjacent, feed):
        return feed

    def adsorption(self, t, c, q):
        # Evaluate equilibrium on the physical domain during BDF trial steps.
        physical_q = np.maximum(q, 0)
        if hasattr(self.iso, "pore_c_eq"):
            equilibrium = self.iso.pore_c_eq(c, physical_q, **self.environment(t))
        elif self.p.IsothermModel == 2:
            equilibrium = self.iso.c_eq(physical_q, c_context=np.maximum(c, 0))
        else:  # Linear equilibrium is independent of salt and pH.
            equilibrium = self.iso.c_eq(physical_q)
        return self.kb * (c - equilibrium)

    def coupling(self, c, cm):
        return 3 * self.De / self.radius * self.space.outward_gradient(c)

