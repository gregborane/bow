import sys
from dataclasses import dataclass
from pathlib import Path

ROOT_DIR=Path.cwd().resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from isotherms import *
from kinetics import *
from solvers import *


@dataclass
class StoreParameters:

    K  : float | None = None
    qm : float | None = None
    kinf : float | None = None
    theta : float | None = None

# %%

