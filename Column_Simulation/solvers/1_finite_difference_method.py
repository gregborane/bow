"""Spatial operators on uniform axial and descending spherical grids.

Backend contract: first/second differentiate along axis 0, radial_diffusion
acts on axis 1 of (Nz, Nr, Nc). boundary_value eliminates a linear flux
boundary. outward_gradient differentiates at the outer radial endpoint.
All returned derivatives use physical coordinates (cm).
"""
import numpy as np


class SpatialOperators:
    axial_radius: int = 1
    local_three_point = True  # permits the local BDF Jacobian sparsity optimization

    def __init__(self, z, r):
        self.z, self.r = (np.asarray(grid, dtype=float) for grid in (z, r))
        for grid, direction in ((self.z, 1), (self.r, -1)):
            if grid.ndim != 1 or grid.size < 3 or not np.isfinite(grid).all():
                raise ValueError("Each spatial grid needs at least three finite nodes")
            steps = direction * np.diff(grid)
            if np.any(steps <= 0) or not np.allclose(steps, steps[0], atol=0):
                raise ValueError("Expected uniform z increasing and r decreasing")
        if self.r[-1] != 0:
            raise ValueError("Last radial node must be the particle centre")
        self.dz = self.z[1] - self.z[0]
        self.dr = self.r[0] - self.r[1]
        self.faces = (self.r[:-1] + self.r[1:]) / 2
        bounds = np.r_[self.r[0], self.faces, 0.]
        self.volumes = -np.diff(bounds**3) / 3

    def _positive_faces(self, u):
        return u[:-1]

    def column_transport(self, bulk, velocity, dispersion, inlet_flux):
        """Positive-flow convection and centered dispersion on shared faces."""
        interior = (velocity * self._positive_faces(bulk)
                    - dispersion * np.diff(bulk, axis=0) / self.dz)
        faces = np.concatenate((inlet_flux[None], interior, velocity * bulk[-1:]))
        return -np.diff(faces, axis=0) / self.dz

    def first(self, u, velocity):
        """Upwind derivative, with a one-sided closure at the inflow end."""
        slopes = np.diff(u, axis=0) / self.dz
        if velocity >= 0:
            return np.concatenate((slopes[:1], slopes))
        return np.concatenate((slopes, slopes[-1:]))

    def second(self, u):
        """Central interior stencil and three-point one-sided end stencils."""
        interior = np.diff(u, n=2, axis=0) / self.dz**2
        return np.concatenate((interior[:1], interior, interior[-1:]))

    @staticmethod
    def boundary_value(adjacent, external, transfer, diffusivity, spacing):
        """Solve D*(boundary-adjacent)/h = transfer*(external-boundary)."""
        conductance = diffusivity / spacing
        return adjacent + transfer / (transfer + conductance) * (external - adjacent)

    def outward_gradient(self, u):
        return (u[:, 0] - u[:, 1]) / self.dr

    def radial_diffusion(self, u, diffusivity, surface_flux):
        """Shell balance of r² D dc/dr. surface_flux is D dc/dr at Rp.

        The centre shell has zero inner face area and positive volume,
        yielding 6 D (c_neighbor - c_center) / dr² at r=0.
        """
        flux = np.empty((u.shape[0], self.r.size + 1, u.shape[-1]))
        flux[:, 0] = self.r[0]**2 * surface_flux
        flux[:, 1:-1] = (-np.diff(u, axis=1) / self.dr * diffusivity
                         * self.faces[None, :, None]**2)
        flux[:, -1] = 0
        return -np.diff(flux, axis=1) / self.volumes[None, :, None]

    def particle_interior(self, u, diffusivity, surface_flux):
        """Algebraic surface has no independent storage. Extend the first
        dynamic shell to Rp so that the retained shell balances conserve flux.
        """
        rates = self.radial_diffusion(u, diffusivity, surface_flux)
        rates[:, 1] = ((rates[:, 0] * self.volumes[0] + rates[:, 1] * self.volumes[1])
                       / (self.volumes[0] + self.volumes[1]))
        return rates[:, 1:]

