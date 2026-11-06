import sys
from argparse import ArgumentParser
from importlib import import_module
from pathlib import Path
from pkgutil import iter_modules

from Column_Simulation.experiments import VarInputs
from Column_Simulation.experiments.VarInputs import display_results, inputs
from Column_Simulation.paths import output_path
from Column_Simulation.src.misc.results_output import save_simulation_outputs


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


def run(inputs, *, experiment_id=1, output_dir=None):
    """Run the selected pore experiment and return time-indexed fields."""
    import numpy as np

    kinetic, kinetics_model = selected("kinetics", inputs.KineticModel)
    isotherm, equilibrium = selected("isotherms", inputs.IsothermModel)
    solver, backend = selected("solvers", inputs.SolverModel)
    display_model(isotherm, kinetic, solver)
    if inputs.KineticModel not in (3, 4):
        return  # Other kinetic modules are placeholders in this repository.
    if inputs.IsothermModel == 1:
        isotherm_class = equilibrium.LinearIsotherm
    elif inputs.IsothermModel == 2:
        isotherm_class = equilibrium.LangmuirIsotherm
    else:
        raise NotImplementedError(f"No pore model adapter for isotherm {inputs.IsothermModel}")
    if not hasattr(backend, "SpatialOperators"):
        raise NotImplementedError(f"Solver {inputs.SolverModel} lacks spatial operators")
    species = list(inputs.species.values())
    iso = isotherm_class(
        [s.K for s in species], [s.qm for s in species])
    space = backend.SpatialOperators(np.linspace(0, inputs.L, inputs.Nz),
                                     np.linspace(inputs.Rp * 1e-4, 0, inputs.Nr))
    results = kinetics_model.Model(inputs, iso, space).simulate()
    csv_path, png_path = save_simulation_outputs(
        results, experiment_id, output_dir, component_names=inputs.species,
        config_file=VarInputs.__file__)
    print(f"Saved {csv_path}\nSaved {png_path}")
    return results


def cli(argv=()):
    """Run the configured experiment and optionally export all result fields."""
    import numpy as np

    parser = ArgumentParser(description="Run the BOW column simulation")
    parser.add_argument("--output", help="Save results as NPZ under the output directory")
    parser.add_argument("--experiment-id", default="1", help="Identifier for the run directory")
    args = parser.parse_args(argv)
    display_results(inputs)
    results = run(inputs, experiment_id=args.experiment_id)
    if args.output:
        if results is None:
            parser.error("The selected kinetic model produced no results")
        target = output_path(args.output)
        with target.open("wb") as stream:
            np.savez_compressed(stream, **vars(results))
        print(f"Saved {target}")
    return results


def console():
    cli(sys.argv[1:])


if __name__ == "__main__":
    # An editor can execute this file in the kernel's __main__ namespace.
    # Only a launch of this script owns the process command-line arguments.
    script = globals().get("__file__")
    launched_as_script = bool(script and sys.argv and Path(sys.argv[0]).resolve() == Path(script).resolve()
                              and Path(script).name == "main.py")
    results = cli(sys.argv[1:] if launched_as_script else ())

# %%

