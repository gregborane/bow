import sys
from importlib import import_module
from pathlib import Path

ROOT_DIR = Path.cwd().resolve().parents[0]
sys.path.insert(0, str(ROOT_DIR))

print(ROOT_DIR)

from experiments.VarInputs import display_results, inputs

display_results(inputs)

# %%
# Import the kinetics

kinetics_dir = Path("kinetics")
kinetic_id = inputs.KineticModel

for file in kinetics_dir.glob("*.py"):
    if file.stem.split("_")[0] == str(kinetic_id):
        model = import_module(f"{kinetics_dir}.{file.stem}")
        break

# %%
# Importing the isotherm
isotherms_dir = Path("isotherms")
model_id = inputs.IsothermModel

for file in isotherms_dir.glob("*.py"):
    if file.stem.split("_")[0] == str(model_id):
        model = import_module(f"{isotherms_dir}.{file.stem}")
        break

# %%
# Import the solver 
solver_dir = Path("solver")
solver_id = inputs.SolverModel
for file in solver_dir.glob("*.py"):
    if file.stem.split("_")[0] == str(model_id):
        model = import_module(f"{solver_dir}.{file.stem}")
        break


# %%

