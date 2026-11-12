"""Set up, validate, and inspect one chromatography simulation with N species.

All editable inputs live here. Column, solution, operating, and numerical
settings are defined once. Add any number of species to the same StoreInput.
"""

# This configuration is numerically aligned to the uploaded CADET langmuir.py
# reference. See inline comments for parameters that CADET does not expose in
# the same form (notably film diffusion, solution viscosity, and ep_s).
# CADET reference values: film_diffusion=1e-4 m/s, is_kinetic=0,
# mclldf_keq=0, mclldf_kkin=0.33 1/s, mclldf_qmax=1.
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy as scp

from Column_Simulation import main
from Column_Simulation.src.inputs.CalculateAllVariables import (
    ColumnParameters,
    SolutionParameters,
    SpeciesParameters,
    StoreInput,
)
from Column_Simulation.src.inputs.VarError import ErrorDetection
from Column_Simulation.src.misc.results_output import prepare_table


def build_inputs() -> StoreInput:
    """Build one shared simulation containing all configured species."""
    inputs = StoreInput()

    # ==========================================================
    # 1. COLUMN PARAMETERS: shared geometry and stationary phase
    # ==========================================================

    inputs.SetColumn(ColumnParameters(
        dp=9.0e-5,          # Particle DIAMETER (m); CADET par_radius = 4.5e-5 m
        dc=0.0115,          # Column diameter (m), from A = pi*dc^2/4 = 1.0386890711e-4 m^2
        L=0.6,              # Column length (m)
        e=0.37,             # Column porosity (-); CADET col_porosity
        ep_p=0.33,          # Particle porosity (-); CADET par_porosity
        ep_s=0.7,           # No direct CADET counterpart; retained as a VarInputs-only parameter
    ))

    # =========================================================
    # 2. SOLUTION PARAMETERS: shared mobile phase
    # =========================================================
    
    inputs.SetSolution(SolutionParameters(
        Cso=0,              # No salt modifier in the CADET Langmuir reference
        Csf=0,              # Keep salt concentration constant at zero (isocratic case)
        pHo=7,              # Fixed pH placeholder; CADET reference has no pH dependence
        pHf=7,              # Fixed pH placeholder; no pH gradient
        nu=0.001,            # Not specified by CADET; retained pending film-transfer calibration
    ))

    # ============================================================
    # 3. SPECIES: add N proteins/solutes to this ONE simulation
    #    Each species has its own transport/adsorption parameters.
    # ============================================================
    
    inputs.AddSpecies("Protein_A", SpeciesParameters(
        Do = 1e-7,           # VarInputs molecular diffusivity; CADET has no separate Do field
        kc = 0.333,          # Same numeric value as CADET mclldf_kkin (inactive: is_kinetic=0)
        K = 2,        # CADET mclldf_keq = 0; preserves the reference's zero-affinity case
        qm = 1.0,            # CADET mclldf_qmax
        kb = 10000,           # No separate counterpart in CADET's quasi-stationary LDF setup
        Dl = 2e-7,           # CADET col_dispersion (m^2/s)
        De = 1e-7,           # CADET par_diffusion (m^2/s)
        Cfeed = 1.0,          # CADET inlet concentration during the pulse
        kf = 1e-4
    ))

    # Add more species using another inputs.AddSpecies("name", SpeciesParameters(...)).

    # =============================================================
    # 4. SIMULATION: shared operating conditions and numerical grid
    # =============================================================

    # CADET reference: Q = 1e-7 m^3/s, t_inj = 30 s, t_final = 1200 s.
    # Vc = A*L = 6.2321344266e-5 m^3, so CV = Q*t/Vc.
    inputs.Q = 1.0e-7                  # CADET connection flow rate (m^3/s)
    inputs.CVfinal = 1.9255040374      # Q*1200/Vc: total duration using geometric bed volume
    # If this project defines 1 CV as the interstitial volume e*Vc instead,
    # use CVfinal=5.204064966 and CVwash=5.073963342 below.
    inputs.Nt = 601                    # CADET user_solution_times has 601 points
    inputs.Nz = 40                     # CADET set_discretization(..., n_col=40)
    inputs.Nr = 5                      # Assumes the common set_discretization helper default n_par=5

    # Pulse mode only (InputType = 2)
    inputs.Vinj = 3.0e-6               # Q*30 s: gives the same 30 s inlet pulse
    inputs.CVwash = 1.8773664365       # Q*(1200-30)/Vc: 1170 s wash after injection

    # Gradient mode only (ElutionType != 1); unused in this isocratic reference
    inputs.CVG = 5                     # Retained positive because the field is unused in isocratic mode

    # Model selections
    inputs.InputType = 2               # Pulse: inlet concentration 1 until 30 s, then 0
    inputs.ElutionType = 1             # Isocratic: no salt or pH gradient
    inputs.KineticModel = 4            # Pore-diffusion model, matching CADET GENERAL_RATE_MODEL
    inputs.SolverModel = 1             # WENO, matching CADET's WENO reconstruction
    inputs.IsothermModel = 2           # Langmuir

    # =========================================================
    # 5. VALIDATION AND CALCULATION
    # =========================================================

    ErrorDetection().ValidateAll(inputs)
    inputs.CalculateAllVariables()
    return inputs


# Run configuration
inputs = build_inputs()
EXPERIMENT_ID = main.DEFAULT_EXPERIMENT_ID  # Change to a name or number for this run.
OUTPUT_DIR = None  # None selects ./run_<experiment_id>/ in the working directory.
OPEN_EXTERNAL_DISPLAY = True
EXPORT_CONFIG = True


def display_results(inputs: StoreInput):
    """Display shared results once, followed by each species' results."""
    print("\n========== SHARED SIMULATION VARIABLES ==========\n")
    print(f"Particle radius Rp       = {inputs.Rp} m")
    print(f"Column area Ac           = {inputs.Ac} m²")
    print(f"Column volume Vc         = {inputs.Vc} m³")
    print(f"Superficial velocity u   = {inputs.u} m/s")
    print(f"Interstitial velocity v  = {inputs.v} m/s")
    print(f"Final simulation time    = {inputs.tfinal} s")
    print(f"Time points tspan        = {inputs.tspan}")
    print(f"Phase ratio theta        = {inputs.theta}")
    print(f"Salt retention km        = {inputs.km}")
    print(f"Non-binding kinf         = {inputs.kinf}")
    print(f"Axial step dz            = {inputs.dz} m")
    print(f"Radial step dr           = {inputs.dr} m")
    print(f"Injection time tinj      = {inputs.tinj} s")

    if inputs.InputType == 2:
        print(f"Wash end time twash      = {inputs.twash} s")

    if inputs.ElutionType != 1:
        print(f"Gradient time tG         = {inputs.tG} s")
        print(f"Salt gradient Bs         = {inputs.Bs} M/s")
        print(f"pH gradient BpH          = {inputs.BpH} pH units/s")
        print(f"Normalized slope gamma   = {inputs.gamma}")

    if inputs.SolverModel in (3, 4, 5):
        print(f"Reynolds number Re       = {inputs.Re}")

    print("\n========== SPECIES VARIABLES ==========\n")
    for name, parameters in inputs.species.items():
        results = inputs.GetSpeciesResults(name)
        print(f"{name}:")
        print(f"  Molecular diffusivity Do = {parameters.Do}")
        print(f"  LDF coefficient kc       = {parameters.kc}")
        print(f"  Henry constant           = {parameters.K}")
        print(f"  Qm                       = {parameters.qm}")

        if inputs.SolverModel == 2:
            print(f"  LDF parameter n          = {results.n}")
        elif inputs.SolverModel in (3, 4, 5):
            print(f"  Schmidt number Sc        = {results.Sc}")
            print(f"  Sherwood number Sh       = {results.Sh}")
            print(f"  Film transfer coeff kf   = {results.kf} m/s")
            print("  CADET reference film kf  = 0.0001 m/s (explicit CADET input)")
        print()

    print("==========================================\n")

def _config_source():
    # Editors can execute the buffer in a kernel namespace without __file__.
    if "__file__" in globals():
        return Path(__file__).resolve()
    from Column_Simulation.experiments import VarInputs
    return Path(VarInputs.__file__).resolve()

def run_simulation(custom_inputs=None, *, experiment_id=None):
    """Launch the configured experiment or use a caller-supplied input object."""
    selected_inputs = inputs if custom_inputs is None else custom_inputs
    display_results(selected_inputs)
    return main.run(
        selected_inputs,
        experiment_id=EXPERIMENT_ID if experiment_id is None else experiment_id,
        output_dir=OUTPUT_DIR,
        config_file=_config_source() if EXPORT_CONFIG else None,
    )

results_raw = run_simulation()
results, names = prepare_table(results_raw, inputs.species)

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
plt.savefig('./cadet.png')
plt.show() 

# %%

