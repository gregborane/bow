"""
Boundary concentrations are supplied by the kinetic model. Reuse BOW's
central axial and conservative spherical diffusion operators (SI units).
"""
from importlib import import_module

import numpy as np

FiniteDifferences = import_module(".1_finite_difference_method", __package__).SpatialOperators


class SpatialOperators(FiniteDifferences):
    local_three_point = False
    axial_radius = 3

    def __init__(self, z, r, epsilon=1e-12):
        super().__init__(z, r)
        if not np.isscalar(epsilon) or not np.isfinite(epsilon) or epsilon <= 0:
            raise ValueError("WENO epsilon must be positive and finite")
        self.epsilon = float(epsilon)

    def _positive_faces(self, u):
        """Faces between nodes, ordered upstream to downstream.

        Reduce WENO5 to WENO3, then to upwind, as boundary stencils shrink.
        Each interior control volume shares both face values with its neighbors.
        """
        face = u[:-1].copy()
        for i in (1, len(u) - 2):
            left = u[i] - u[i-1]
            right = u[i+1] - u[i]
            wl = (1/3) / (self.epsilon + left**2)**2
            wr = (2/3) / (self.epsilon + right**2)**2
            face[i] = u[i] + 0.5 * (wl*left + wr*right) / (wl + wr)
        a, b, c, d, e = u[:-4], u[1:-3], u[2:-2], u[3:-1], u[4:]
        p0 = (2*a - 7*b + 11*c) / 6
        p1 = (-b + 5*c + 2*d) / 6
        p2 = (2*c + 5*d - e) / 6
        beta0 = (13/12)*(a - 2*b + c)**2 + 0.25*(a - 4*b + 3*c)**2
        beta1 = (13/12)*(b - 2*c + d)**2 + 0.25*(b - d)**2
        beta2 = (13/12)*(c - 2*d + e)**2 + 0.25*(3*c - 4*d + e)**2
        w0 = 0.1 / (self.epsilon + beta0)**2
        w1 = 0.6 / (self.epsilon + beta1)**2
        w2 = 0.3 / (self.epsilon + beta2)**2
        face[2:-1] = (w0*p0 + w1*p1 + w2*p2) / (w0 + w1 + w2)
        return face

    def first(self, u, velocity):
        """Return dc/dz for a finite scalar velocity, along axis 0 only."""
        u = np.asarray(u, dtype=float)
        if u.ndim == 0 or u.shape[0] != self.z.size:
            raise ValueError("Expected one concentration per axial node on axis 0")
        if np.ndim(velocity) or not np.isfinite(velocity):
            raise ValueError("Expected a spatially constant finite scalar velocity")
        faces = self._positive_faces(u) if velocity >= 0 else self._positive_faces(u[::-1])[::-1]
        interior = np.diff(faces, axis=0) / self.dz
        return np.concatenate(((u[1:2]-u[:1])/self.dz, interior,
                               (u[-1:]-u[-2:-1])/self.dz))

    def second(self, u):
        """Second-order central interiors and one-sided four-node endpoints."""
        u = np.asarray(u, dtype=float)
        interior = np.diff(u, n=2, axis=0)
        left = 2*u[:1] - 5*u[1:2] + 4*u[2:3] - u[3:4]
        right = 2*u[-1:] - 5*u[-2:-1] + 4*u[-3:-2] - u[-4:-3]
        return np.concatenate((left, interior, right)) / self.dz**2

    def flux_derivative(self, u, velocity):
        """Return d(velocity*c)/dz for spatially constant velocity."""
        return velocity * self.first(u, velocity)
