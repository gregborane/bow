"""Set up, validate, and inspect one chromatography simulation with N species.

All editable inputs live here. Column, solution, operating, and numerical
settings are defined once. Add any number of species to the same StoreInput.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as scp

from bow.Col_Simulation.inputs.CalculateAllVariables import (
    ColumnParameters,
    SolutionParameters,
    SpeciesParameters,
    StoreInput,
)
from bow.Col_Simulation.inputs.VarError import ErrorDetection
from bow.Col_Simulation.main import DEFAULT_EXPERIMENT_ID, run


def build_inputs() -> StoreInput:
    """Build one shared simulation containing all configured species."""
    inputs = StoreInput()

    # ==========================================================
    # 1. COLUMN PARAMETERS: shared geometry and stationary phase
    # ==========================================================

    inputs.SetColumn(ColumnParameters(
        dp=9.0e-5,
        dc=0.0115,
        L=0.6,
        e=0.37,
        ep_p=0.33,
        ep_s=0.7,
    ))

    # =========================================================
    # 2. SOLUTION PARAMETERS: shared mobile phase
    # =========================================================
    
    inputs.SetSolution(SolutionParameters(
        Cso=0,              
        Csf=0,              
        pHo=7,              
        pHf=7,              
        nu=0.001,            
    ))

    # ============================================================
    # 3. SPECIES: add N proteins/solutes to this ONE simulation
    #    Each species has its own transport/adsorption parameters.
    # ============================================================
    
    inputs.AddSpecies("Protein_A", SpeciesParameters(
        Do = 1e-7,
        kc = 0,
        K = 2,
        qm = 1.0,
        kb = 0.33,
        Dl = 2e-7,
        De = 1e-7,
        Cfeed = 1.0,
        kf = 1e-4,
    ))

    # Add more species using another inputs.AddSpecies("name", SpeciesParameters(...)).

    # =============================================================
    # 4. SIMULATION: shared operating conditions and numerical grid
    # =============================================================

    inputs.Q = 1.0e-7
    column_volume = np.pi * (inputs.dc / 2)**2 * inputs.L
    inputs.CVfinal = 1200 * inputs.Q / column_volume
    inputs.Nt = 601
    inputs.Nz = 40
    inputs.Nr = 5

    # Pulse mode only (InputType = 2)
    inputs.Vinj = 3.0e-6
    inputs.CVwash = inputs.CVfinal - inputs.Vinj / column_volume

    # Gradient mode only (ElutionType != 1); unused in this isocratic reference
    inputs.CVG = 5

    # Model selections
    inputs.InputType = 2
    inputs.ElutionType = 1
    inputs.KineticModel = 4
    inputs.binding_equilibrium = True
    inputs.SolverModel = 5
    inputs.IsothermModel = 2

    # =========================================================
    # 5. VALIDATION AND CALCULATION
    # =========================================================

    ErrorDetection().ValidateAll(inputs)
    inputs.CalculateAllVariables()
    return inputs


# Run configuration
inputs = build_inputs()
EXPERIMENT_ID = DEFAULT_EXPERIMENT_ID
OUTPUT_DIR = None
OPEN_EXTERNAL_DISPLAY = True
EXPORT_CONFIG = True

results, names = run(inputs=inputs,
                     experiment_id=EXPERIMENT_ID,
                     output_dir=None,
                     config_file=None,
                     )

time = results[:, 0] 
c = results[:, 1]
cs = scp.interpolate.CubicSpline(time, c)
intergration = scp.integrate.quad(cs, 0, time[-1])[0]
max_index = np.where(c == np.max(c))[0][0]
maxy = c[max_index]
maxx = time[max_index]

plt.figure()
plt.plot(time, c)
plt.plot([maxx, maxx], [0, maxy], label=f'retention time = {maxx.round(1)} s')
plt.plot([0, maxx], [maxy, maxy], label=f'concentration = {maxy.round(2)} mol/m³')
plt.fill(time, c, hatch='/', label=f"integration value = {np.float64(intergration).round(2)}")
plt.title('Column (Outlet)')
plt.xlabel('$time~/~s$')
plt.ylabel('$concentration~/~mM$')
plt.legend()
plt.savefig('./run_1/run.png')
plt.show() 


