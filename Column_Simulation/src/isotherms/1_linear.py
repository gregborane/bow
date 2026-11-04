"""Linear adsorption isotherms for BOW.
"""
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from misc.isotherm import _affinities, _species_parameter, _state
from numpy.typing import ArrayLike, NDArray

SRC_DIR = Path.cwd().resolve().parents[1]
sys.path.insert(0, str(SRC_DIR))


FloatArray = NDArray[np.float64]

# %%

@dataclass(frozen=True, eq=False)
class LinearIsotherm:
    """LDF linear isotherm: q_i = qm_i * K_i * c_i.

    `K` is one affinity per component; `qm` may be a scalar common to all
    components or one capacity per component.
    """

    K: ArrayLike
    qm: ArrayLike

    def __post_init__(self) -> None:
        k = _affinities(self.K)
        object.__setattr__(self, "K", k)
        object.__setattr__(self, "qm", _species_parameter(self.qm, k.size, "qm"))

    @property
    def nc(self) -> int:
        return self.K.size

    def q_eq(self, c: ArrayLike) -> FloatArray:
        """Forward equilibrium q(c), preserving the input grid dimensions."""
        c = _state(c, self.nc, "c")
        return self.qm * self.K * c

    def c_eq(self, q: ArrayLike) -> FloatArray:
        """Inverse equilibrium c(q)."""
        q = _state(q, self.nc, "q")
        return q / (self.qm * self.K)

# %%

