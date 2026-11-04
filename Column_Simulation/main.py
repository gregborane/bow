import sys
from importlib import import_module
from pathlib import Path

ROOT_DIR = Path.cwd().resolve()
sys.path.insert(0, str(ROOT_DIR))

SRC_DIR = ROOT_DIR / "src"
sys.path.insert(0, str(SRC_DIR))

INPUTS_DIR = ROOT_DIR / "experiments"
sys.path.insert(0, str(SRC_DIR))

# %%
# Import the input file
from experiments.VarInputs import display_results, inputs

display_results(inputs)

# %%

def display_model(isotherm, kinetic, solver):
    print('\n========= SIMULATION FUNCTION AND ALGORITHM SELECTED ==========\n')
    print(f'Isotherm Model selected is : {isotherm.stem}')
    print(f'Kinetic Model selected is : {kinetic.stem}') 
    print(f'Solver selected is : {solver.stem}')
    print('\n===============================================================\n')


# Import the kinetics
kinetics_dir = SRC_DIR / "kinetics"
kinetic_id = inputs.KineticModel

for kinetic in kinetics_dir.glob("*.py"):
    if kinetic.stem.partition("_")[0] == str(kinetic_id):
        kinetics_model = import_module(f"kinetics.{kinetic.stem}")
        break

# Importing the isotherm
isotherms_dir = SRC_DIR / "isotherms"
model_id = inputs.IsothermModel
check_m = None

for isotherm in isotherms_dir.glob("*.py"):
    if isotherm.stem.partition("_")[0] == str(model_id):
        model = import_module(f"isotherms.{isotherm.stem}")
        break

# Import the solver 
solver_dir = SRC_DIR / "solvers"
solver_id = inputs.SolverModel
check_s = None

for solver in solver_dir.glob("*.py"):
    if solver.stem.partition("_")[0] == str(model_id):
        model = import_module(f"solvers.{solver.stem}")
        break

display_model(isotherm, kinetic, solver) 

# %%

