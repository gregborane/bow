from importlib import import_module
from pathlib import Path
from pkgutil import iter_modules

from bow.Col_Simulation.inputs.CalculateAllVariables import StoreInput
from bow.Col_Simulation.misc.results_output import (
    DEFAULT_EXPERIMENT_ID,
    prepare_table,
    save_simulation_outputs,
)


def display_model(isotherm, kinetic, solver):
    print('\n========= SIMULATION FUNCTION AND ALGORITHM SELECTED ==========\n')
    print(f'Isotherm Model selected is : {isotherm.stem}')
    print(f'Kinetic Model selected is : {kinetic.stem}') 
    print(f'Solver selected is : {solver.stem}')
    print('\n===============================================================\n')

def display_results(inputs: StoreInput):
    """Display shared results once, followed by each species' results."""
    print("\n========== SHARED SIMULATION VARIABLES ==========\n")
    print(f"Particle radius Rp       = {inputs.Rp} m")
    print(f"Column area Ac           = {inputs.Ac} m²")
    print(f"Column volume Vc         = {inputs.Vc} m³")
    print(f"Superficial velocity u   = {inputs.u} m/s")
    print(f"Interstitial velocity v  = {inputs.v} m/s")
    print(f"Final simulation time    = {inputs.tfinal} s")
    print(f"Phase ratio theta        = {inputs.theta}")
    print(f"Salt retention km        = {inputs.km}")
    print(f"Non-binding kinf         = {inputs.kinf}")
    print(f"Axial step dz            = {inputs.dz} m")
    print(f"Radial step dr           = {inputs.dr} m")
    print(f"Time step                = {inputs.tfinal / inputs.Nt} s")
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

    print("\n==========================================\n")

def selected(group, identifier):
    package = import_module(f"bow.Col_Simulation.{group}")
    prefix = f"{identifier}_"
    names = [item.name for item in iter_modules(package.__path__) if item.name.startswith(prefix)]
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

    display_results(inputs)
    results_clean, names = prepare_table(results, inputs.species)
    return results_clean, names

