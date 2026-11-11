"""Set up, validate, and inspect one chromatography simulation with N species.

All editable inputs live here. Column, solution, operating, and numerical
settings are defined once. Add any number of species to the same StoreInput.
"""
from pathlib import Path

from Column_Simulation import main
from Column_Simulation.src.inputs.CalculateAllVariables import (
    ColumnParameters,
    SolutionParameters,
    SpeciesParameters,
    StoreInput,
)
from Column_Simulation.src.inputs.VarError import ErrorDetection


def build_inputs() -> StoreInput:
    """Build one shared simulation containing all configured species."""
    inputs = StoreInput()

    # ==========================================================
    # 1. COLUMN PARAMETERS: shared geometry and stationary phase
    # ==========================================================

    inputs.SetColumn(ColumnParameters(
        dp=4.5e-5,          # Particle diameter (m)
        dc=0.011,           # Column diameter (m)
        L=0.6,              # Column length (m)
        e=0.4,              # Column porosity (-)
        ep_p=0.4,           # Particle phase parameter (-)
        ep_s=0.7,           # Salt phase parameter (-)
    ))

    # =========================================================
    # 2. SOLUTION PARAMETERS: shared mobile phase
    # =========================================================
    
    inputs.SetSolution(SolutionParameters(
        Cso=0,              # Initial salt concentration (M)
        Csf=100,            # Final salt concentration (M)
        pHo=None,           # Initial pH
        pHf=7,              # Final pH
        nu=0.01,            # Kinematic viscosity (m²/s)
    ))

    # ============================================================
    # 3. SPECIES: add N proteins/solutes to this ONE simulation
    #    Each species has its own transport/adsorption parameters.
    # ============================================================
    
    inputs.AddSpecies("Protein_A", SpeciesParameters(
        Do = 1e-7,            # Molecular diffusivity (m²/s)
        kc = 0.33,            # LDF mass-transfer coefficient (1/s)
        K = 0.0001,         # Henry Constant
        qm = 1,             # Maximum quantity adsorbed 
        kb = 0.33,          # Ads rate constant
        Dl = 2e-7,          # Axial dispersion coefficient
        De = 1e-7,          # Effective Pore diffusivity
        Cfeed = 1           # Molecular concentration
    ))

    # Add more species using another inputs.AddSpecies("name", SpeciesParameters(...)).

    # =============================================================
    # 4. SIMULATION: shared operating conditions and numerical grid
    # =============================================================

    inputs.Q = 1            # Volumetric flow rate (m³/s)
    inputs.CVfinal = 10     # Total simulation duration (column volumes)
    inputs.Nt = 100         # Number of time points
    inputs.Nz = 50          # Axial grid points
    inputs.Nr = 20          # Radial particle grid points

    # Pulse mode only (InputType = 2)
    inputs.Vinj = 1         # Injection volume (m³)
    inputs.CVwash = 1       # Wash duration (column volumes)

    # Gradient mode only (ElutionType != 1)
    inputs.CVG = 5          # Gradient duration (column volumes)

    # Model selections (original numeric codes retained)
    inputs.InputType = 2        # 1 step | 2 pulse
    inputs.ElutionType = 1      # 1 isocratic | 2-6 gradient types
    inputs.KineticModel = 3     # 2 LDF c | 3,4,5 pore diffusion
    inputs.SolverModel = 1      # 1 FDM | 5 WENO 
    inputs.IsothermModel = 2    # 2 Langmuir

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
        open_external_display=OPEN_EXTERNAL_DISPLAY,
        config_file=_config_source() if EXPORT_CONFIG else None,
    )

results = run_simulation()
