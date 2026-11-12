from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bow.Col_Simulation.misc.isotherm import _affinities, _species_parameter, _state

FloatArray = NDArray[np.float64]

# %% 
@dataclass(frozen=True, eq=False)
class Isotherm:
    """Competitive multicomponent Langmuir isotherm.

    q_i = qm_i*K_i*c_i / (1 + sum_j K_j*c_j)
          + (kinf_i/theta)*c_i  [optional LDF-only term]

    For pore-diffusion Langmuir, omit `kinf` and `theta`.
    Inverse methods are available for the pure competitive Langmuir case
    (kinf absent or zero). Deliberately rejecting nonzero `kinf`.
    """

    K: ArrayLike
    qm: ArrayLike
    kinf: ArrayLike | None = None
    theta: float | None = None

    def __post_init__(self) -> None:
        k = _affinities(self.K, allow_zero=True)
        object.__setattr__(self, "K", k)
        object.__setattr__(self, "qm", _species_parameter(self.qm, k.size, "qm"))

    @property
    def nc(self) -> int:
        return self.K.size

    def q_eq(self, c: ArrayLike) -> FloatArray:
        """Forward equilibrium q(c), competitive along the last axis."""
        c = _state(c, self.nc, "c")
        kc = self.K * c
        denom = 1.0 + np.sum(kc, axis=-1, keepdims=True)
        q = self.qm * kc / denom
        if self.kinf is not None and np.any(self.kinf > 0):
            q = q + (self.kinf / self.theta) * c
        return q

    def loading_jacobian(self, c):
        """dq_i/dc_j for pure competitive pore Langmuir, components last."""
        active = c >= 0
        c = np.maximum(c, 0)
        denominator = 1 + np.sum(self.K * c, axis=-1, keepdims=True)
        slope = self.qm * self.K / denominator
        q = slope * c
        return (slope[..., :, None] * np.eye(self.nc)
                - q[..., :, None] * (self.K / denominator)[..., None, :]) * active[..., None, :]

    def c_eq(
        self, q: ArrayLike, c_context: ArrayLike | None = None
    ) -> FloatArray:
        """Invert the *pure* competitive Langmuir model.

        When c_context is None, invert all components **simultaneously**:
            S = sum_i (q_i / qm_i)
            c_i = q_i / [qm_i*K_i*(1 - S)]
        This requires S < 1 at every grid location.

        When c_context is supplied, 
        inverse, which treats competing c_j as given (j != i):
            c_i = q_i*(1 + sum_{j!=i} K_j*c_j)/(K_i*(qm_i - q_i))

        """
        q = _state(q, self.nc, "q")

        if c_context is None:
            occupancy = np.sum(q / self.qm, axis=-1, keepdims=True)
            if np.any(occupancy >= 1.0):
                raise ValueError("Langmuir inverse requires sum(q_i/qm_i) < 1")
            return q / (self.qm * self.K * (1.0 - occupancy))

        c = _state(c_context, self.nc, "c_context")
        try:
            q, c = np.broadcast_arrays(q, c)
        except ValueError as exc:
            raise ValueError("q and c_context have incompatible grid shapes") from exc
        if np.any(q >= self.qm):
            raise ValueError("Requires q_i < qm_i")
        other = np.sum(self.K * c, axis=-1, keepdims=True) - self.K * c
        return q * (1.0 + other) / (self.K * (self.qm - q))
