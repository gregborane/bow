# Path and import audit

The entry point imports the experiment, which builds and validates the input store.
The runner imports the numbered kinetic, isotherm and solver modules from the installed
package. The q-equilibrium model imports the shared pore model. Both isotherms use
the shared equilibrium helpers. The path module handles input and output filenames
separately from Python package discovery. The existing checkout contained no external
data reader, result writer, automated tests, packaging or IDE configuration. This
refactor adds an optional NPZ export and path resolution for external input files.

The table covers the path refactor. The pore model, finite-difference solver and species
inputs already included the numerical changes from the preceding implementation.

| Repository path | Finding | Action | Status |
|---|---|---|---|
| `Column_Simulation/__init__.py` | Missing package boundary | Mark the directory as an importable package | Added |
| `Column_Simulation/experiments/VarInputs.py` | Experiment used working directory to find inputs | Import parameter builders and validator by package name | Modified |
| `Column_Simulation/experiments/__init__.py` | Missing package boundary | Mark the directory as an importable package | Added |
| `Column_Simulation/main.py` | Path injection and directory glob for numbered imports | Import experiments and discover models by package; expose CLI and explicit NPZ export | Modified |
| `Column_Simulation/paths.py` | No resource location policy | Provide validated roots and file resolution helpers | Added |
| `Column_Simulation/src/__init__.py` | Missing package boundary | Mark the directory as an importable package | Added |
| `Column_Simulation/src/inputs/CalculateAllVariables.py` | No path dependency | Keep existing implementation | No path edit |
| `Column_Simulation/src/inputs/VarError.py` | Ambiguous top-level inputs import | Import StoreInput and SpeciesParameters from the installed package | Modified |
| `Column_Simulation/src/inputs/VarInputsExemple.py` | Sibling-only imports, missing commas and missing kinetic selection blocked importing the example | Qualify imports, fix syntax and separate FDM from kinetic selection | Modified |
| `Column_Simulation/src/inputs/__init__.py` | Missing package boundary | Mark the directory as an importable package | Added |
| `Column_Simulation/src/isotherms/10_modified_logarithmic_model.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/11_steric_mass_action_linear_pH.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/12_hydrophobic_interaction_chromatography.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/13_hic_polynomial.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/1_linear.py` | Path injection and ambiguous helper import | Remove injection and qualify the isotherm helper import | Modified |
| `Column_Simulation/src/isotherms/2_langmuir.py` | Path injection and ambiguous helper import | Remove injection and qualify the isotherm helper import | Modified |
| `Column_Simulation/src/isotherms/3_ion_exchange_mass_action.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/4_steric_action.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/5_steric_mass_action_linear_equation.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/6_interpolated_isotherm.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/7_modified_bet.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/8_bimodal_isotherm.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/9_modified_langmuir.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/isotherms/__init__.py` | Missing package boundary | Mark the directory as an importable package | Added |
| `Column_Simulation/src/kinetics/1_LDF_qeq_Framework.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/kinetics/2_LDF_ceq_Framework.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/kinetics/3_PoreDiffusion_ceq_Framework.py` | Unqualified model, solver or helper imports | Use the installed package namespace | Modified |
| `Column_Simulation/src/kinetics/4_PoreDiffusion_qeq_Framework.py` | Unqualified model, solver or helper imports | Use the installed package namespace | Modified |
| `Column_Simulation/src/kinetics/5_PoreDiffusion_ceq_Framework_twosite.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/kinetics/6_LDF_ceq_Framework_reaction.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/kinetics/__init__.py` | Missing package boundary | Mark the directory as an importable package | Added |
| `Column_Simulation/src/misc/__init__.py` | Missing package boundary | Mark the directory as an importable package | Added |
| `Column_Simulation/src/misc/isotherm.py` | No path dependency | Keep existing implementation | No path edit |
| `Column_Simulation/src/solvers/1_finite_difference_method.py` | No path dependency | Keep existing implementation | No path edit |
| `Column_Simulation/src/solvers/2_finite_volume_method.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/solvers/3_finite_element_method.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/solvers/4_orthogonal_collocation_on_finite_elements.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/solvers/5_weighted_essentially_non_oscillatory.py` | Unused working-directory mutation in a placeholder | Remove the bootstrap and retain the module | Modified |
| `Column_Simulation/src/solvers/__init__.py` | Missing package boundary | Mark the directory as an importable package | Added |
| `pyproject.toml` | Missing package metadata and entry point | Declare dependencies, packages, marker and console command | Added |
| `tests/test_paths.py` | Missing path-independence tests | Exercise subprocess launches, imports, simulation and I/O | Added |
| `readme.md` | Missing install, run and path policy | Document package and filesystem interfaces | Modified |
| `.gitignore` | Local installation and test runs produce generated files | Exclude generated files | Added |
| `docs/PathAudit.md` | Missing path audit | Record dependencies and file changes | Added |
