from importlib import import_module
from pathlib import Path
from pkgutil import iter_modules

from Column_Simulation.src.misc.results_output import (
    DEFAULT_EXPERIMENT_ID,
    save_simulation_outputs,
)


def display_model(isotherm, kinetic, solver):
    print('\n========= SIMULATION FUNCTION AND ALGORITHM SELECTED ==========\n')
    print(f'Isotherm Model selected is : {isotherm.stem}')
    print(f'Kinetic Model selected is : {kinetic.stem}') 
    print(f'Solver selected is : {solver.stem}')
    print('\n===============================================================\n')


def selected(group, identifier):
    package = import_module(f"Column_Simulation.src.{group}")
    prefix = f"{identifier}_"
    names = [item.name for item in iter_modules(package.__path__) if item.name.startswith(prefix)]
    if len(names) != 1:
        raise ValueError(f"Expected one {group} module for ID {identifier}")
    module = import_module(f".{names[0]}", package.__name__)
    return Path(module.__file__), module


def run(inputs, experiment_id=DEFAULT_EXPERIMENT_ID, output_dir=None,
        config_file=None,):
    """Run the selected pore experiment and return time-indexed fields."""

    kinetic, kinetics_model = selected("kinetics", inputs.KineticModel)
    isotherm, equilibrium = selected("isotherms", inputs.IsothermModel)
    solver, backend = selected("solvers", inputs.SolverModel)

    display_model(isotherm, kinetic, solver)
    
    if inputs.KineticModel not in (3, 4):
        raise NotImplementedError(f"Selectec inputs.KineticModel in (3,4) \n now is {inputs.KineticModel}") 
    if inputs.IsothermModel in (1,2):
        isotherm_class = equilibrium.Isotherm
    else:
        raise NotImplementedError(f"No pore model adapter for isotherm {inputs.IsothermModel}")
    
    if not hasattr(backend, "SpatialOperators"):
        raise NotImplementedError(f"Solver {inputs.SolverModel} lacks spatial operators")
    
    species = list(inputs.species.values())
    iso = isotherm_class(
        [s.K for s in species], [s.qm for s in species])
    space = backend.SpatialOperators(inputs.z, inputs.r)
    
    results = kinetics_model.Model(inputs, iso, space).simulate()
    
    csv_path = save_simulation_outputs(
        results, 
        experiment_id, 
        output_dir, 
        component_names=inputs.species,
        config_file=config_file, )
    return results
