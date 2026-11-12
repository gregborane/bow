from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]

def _affinities(values: ArrayLike, *, allow_zero=False) -> FloatArray:
    """Normalise K to a nonempty, finite (Nc,) array on the requested domain."""
    array = np.array(values, dtype=np.float64, copy=True)
    if array.ndim == 0:
        array = array.reshape(1)
    if (array.ndim != 1 or array.size == 0 or
            not np.all(np.isfinite(array)) or np.any(array < 0 if allow_zero else array <= 0)):
        raise ValueError("Invalid affinity for this isotherm")
    array.flags.writeable = False
    return array


def _species_parameter(
    values: ArrayLike, nc: int, name: str, *, allow_zero: bool = False
) -> FloatArray:
    """Normalise a scalar or Nc-vector parameter to shape (Nc,)."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim == 0:
        arr = np.full(nc, float(arr))
    elif arr.ndim == 1 and arr.size == nc:
        arr = np.array(arr, dtype=np.float64, copy=True)
    else:
        raise ValueError(f"{name} must be a scalar or a vector of length {nc}")

    if not np.all(np.isfinite(arr)) or np.any(arr < 0 if allow_zero else arr <= 0):
        condition = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{name} must contain finite, {condition} values")
    arr.flags.writeable = False
    return arr


def _state(values: ArrayLike, nc: int, name: str) -> FloatArray:
    """Validate a nonnegative concentration/loading array, with components last."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim == 0 and nc == 1:
        arr = arr.reshape(1)
    if arr.ndim == 0 or arr.shape[-1] != nc:
        raise ValueError(f"{name} must have final dimension Nc={nc}; got {arr.shape}")
    if not np.all(np.isfinite(arr)) or np.any(arr < 0):
        raise ValueError(f"{name} must contain finite, nonnegative values")
    return arr
