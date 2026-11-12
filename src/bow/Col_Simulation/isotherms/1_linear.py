"""Linear adsorption isotherms for BOW.
"""
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bow.Col_Simulation.misc.isotherm import _affinities, _species_parameter, _state

FloatArray = NDArray[np.float64]

# %%

@dataclass(frozen=True, eq=False)
class Isotherm:
    """LDF linear isotherm: q_i = qm_i * K_i * c_i.

    `K` is one affinity per component; `qm` may be a scalar common to all
    components or one capacity per component.
    """

    K: ArrayLike
    qm: ArrayLike

    def __post_init__(self) -> None:
        k = _affinities(self.K, allow_zero=True)
        object.__setattr__(self, "K", k)
        object.__setattr__(self, "qm", _species_parameter(self.qm, k.size, "qm"))

    @property
    def nc(self) -> int:
        return self.K.size

    def q_eq(self, c: ArrayLike) -> FloatArray:
        """Forward equilibrium q(c), preserving the input grid dimensions."""
        c = _state(c, self.nc, "c")
        return self.qm * self.K * c

    def loading_jacobian(self, c):
        return np.diag(self.qm * self.K) * (c >= 0)[..., None, :]

    def c_eq(self, q: ArrayLike) -> FloatArray:
        """Inverse equilibrium c(q)."""
        q = _state(q, self.nc, "q")
        if np.any(self.K == 0):
            raise ValueError("Zero affinity has no unique inverse concentration")
        return q / (self.qm * self.K)

# %%

