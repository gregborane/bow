"""Set up, validate, and inspect one chromatography simulation with N species.

All editable inputs live here. Column, solution, operating, and numerical
settings are defined once. Add any number of species to the same StoreInput.
"""
from argparse import ArgumentParser
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

    # =========================================================
    # 1. COLUMN PARAMETERS: shared geometry and stationary phase
    # =========================================================
    inputs.SetColumn(ColumnParameters(
        dp=1,               # Particle diameter (µm)
        dc=1,               # Column diameter (cm)
        L=1,                # Column length (cm)
        e=0.4,              # Column porosity (-)
        ep_p=0.7,           # Particle phase parameter (-)
        ep_s=0.7,           # Salt phase parameter (-)
    ))

    # =========================================================
    # 2. SOLUTION PARAMETERS: shared mobile phase
    # =========================================================
    inputs.SetSolution(SolutionParameters(
        Cso=0,              # Initial salt concentration (mM)
        Csf=100,            # Final salt concentration (mM)
        pHo=7,              # Initial pH
        pHf=7,              # Final pH
        nu=0.01,            # Kinematic viscosity (cm²/s)
    ))

    # =========================================================
    # 3. SPECIES: add N proteins/solutes to this ONE simulation
    #    Each species has its own transport/adsorption parameters.
    # =========================================================
    inputs.AddSpecies("Protein_A", SpeciesParameters(
        Do=1e-6,            # Molecular diffusivity (cm²/s)
        kc=0.01,            # LDF mass-transfer coefficient (1/s)
        K = 1,              # Henry Constant
        qm = 100,          # Maximum quantity adsorbed 
        kb = 0.01,          # Ads rate constant
        Dl = 1e-4,          # Axial dispersion coefficient
        De = 1e-7,          # Effective Pore diffusivity
        Cfeed = 1           # Molecular concentration
    ))

    inputs.AddSpecies("Protein_B", SpeciesParameters(
        Do=4e-6,            # Molecular diffusivity (cm²/s)
        kc=0.02,            # LDF mass-transfer coefficient (1/s)
        K = 20,              # Henry Constant
        qm = 10,            # Maximum quantity adsorbed 
        kb = 1,            # Ads rate constant
        Dl = 1e-3,          # Axial dispersion coefficient
        De = 1e-6,          # Effective Pore diffusivity
        Cfeed = 1           # Molecular concentration
    ))
    # Add more species using another inputs.AddSpecies("name", SpeciesParameters(...)).

    # =========================================================
    # 4. SIMULATION: shared operating conditions and numerical grid
    # =========================================================
    inputs.Q = 1            # Volumetric flow rate (cm³/min)
    inputs.CVfinal = 10     # Total simulation duration (column volumes)
    inputs.Nt = 100         # Number of time points
    inputs.Nz = 50          # Axial grid points
    inputs.Nr = 20          # Radial particle grid points

    # Pulse mode only (InputType = 2)
    inputs.Vinj = 1         # Injection volume (cm³)
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
OPEN_EXTERNAL_DISPLAY = False
EXPORT_CONFIG = True


def display_results(inputs: StoreInput):
    """Display shared results once, followed by each species' results."""
    print("\n========== SHARED SIMULATION VARIABLES ==========\n")
    print(f"Particle radius Rp       = {inputs.Rp} µm")
    print(f"Column area Ac           = {inputs.Ac} cm²")
    print(f"Column volume Vc         = {inputs.Vc} cm³")
    print(f"Superficial velocity u   = {inputs.u} cm/s")
    print(f"Interstitial velocity v  = {inputs.v} cm/s")
    print(f"Final simulation time    = {inputs.tfinal} s")
    print(f"Time points tspan        = {inputs.tspan}")
    print(f"Phase ratio theta        = {inputs.theta}")
    print(f"Salt retention km        = {inputs.km}")
    print(f"Non-binding kinf         = {inputs.kinf}")
    print(f"Axial step dz            = {inputs.dz} cm")
    print(f"Radial step dr           = {inputs.dr} cm")
    print(f"Injection time tinj      = {inputs.tinj} s")

    if inputs.InputType == 2:
        print(f"Wash end time twash      = {inputs.twash} s")

    if inputs.ElutionType != 1:
        print(f"Gradient time tG         = {inputs.tG} s")
        print(f"Salt gradient Bs         = {inputs.Bs} mM/s")
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
            print(f"  Film transfer coeff kf   = {results.kf} cm/s")
        print()

    print("==========================================\n")


def run_simulation(custom_inputs=None, *, experiment_id=None):
    """Launch the configured experiment or use a caller-supplied input object."""
    selected_inputs = inputs if custom_inputs is None else custom_inputs
    display_results(selected_inputs)
    return main.run(
        selected_inputs,
        experiment_id=EXPERIMENT_ID if experiment_id is None else experiment_id,
        output_dir=OUTPUT_DIR,
        open_external_display=OPEN_EXTERNAL_DISPLAY,
        config_file=Path(__file__) if EXPORT_CONFIG else None,
    )


def console(argv=None):
    parser = ArgumentParser(description="Run the BOW column simulation")
    parser.add_argument("--experiment-id", help="Override the configured experiment ID")
    parser.add_argument("--output", help="Export all result fields to an additional NPZ file")
    args = parser.parse_args(argv)
    results = run_simulation(experiment_id=args.experiment_id)
    if args.output:
        import numpy as np

        from Column_Simulation.paths import output_path

        if results is None:
            parser.error("The selected kinetic model produced no results")
        target = output_path(args.output)
        with target.open("wb") as stream:
            np.savez_compressed(stream, **vars(results))
        print(f"Saved {target}")


if __name__ == "__main__":
    console()
