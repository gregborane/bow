"""Set up, validate, and inspect one chromatography simulation with N species.

All editable inputs live here. Column, solution, operating, and numerical
settings are defined once. Add any number of species to the same StoreInput.
"""

from CalculateAllVariables import (
    ColumnParameters,
    SolutionParameters,
    SpeciesParameters,
    StoreInput,
)
from VarError import ErrorDetection


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
    ))
    inputs.AddSpecies("Protein_B", SpeciesParameters(
        Do=5e-7,            # Replace example value with measured data
        kc=0.02,            # Replace example value with measured data
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
    inputs.InputType = 1        # 1 step | 2 pulse
    inputs.ElutionType = 1      # 1 isocratic | 2-6 gradient types
    inputs.SolverModel = 3      # 2 LDF c | 3,4,5 pore diffusion
    inputs.IsothermModel = 2    # 2 Langmuir

    # =========================================================
    # 5. VALIDATION AND CALCULATION
    # =========================================================
    ErrorDetection().ValidateAll(inputs)
    inputs.CalculateAllVariables()
    return inputs


# One StoreInput instance represents the entire physical simulation.
inputs = build_inputs()


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

        if inputs.SolverModel == 2:
            print(f"  LDF parameter n          = {results.n}")
        elif inputs.SolverModel in (3, 4, 5):
            print(f"  Schmidt number Sc        = {results.Sc}")
            print(f"  Sherwood number Sh       = {results.Sh}")
            print(f"  Film transfer coeff kf   = {results.kf} cm/s")
        print()

    print("==========================================\n")


if __name__ == "__main__":
    display_results(inputs)

# %%

