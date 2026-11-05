# Choosing a BOW column experiment

Edit `Column_Simulation/experiments/VarInputs.py` to define the column, solution, species and operating conditions. The supplied experiment injects a step of Protein_A under isocratic conditions. It uses concentration-driven pore kinetics, a Langmuir isotherm and finite differences. These are five independent choices.

| Setting | Default | What it selects |
|---|---:|---|
| `InputType` | 1 | Protein feed schedule. 1 is a step and 2 is a pulse. |
| `ElutionType` | 1 | Salt and pH program. 1 holds the initial conditions constant. |
| `KineticModel` | 3 | Transport and adsorption balances. Models 3 and 4 run in Python. |
| `IsothermModel` | 2 | Binding equilibrium. The runner integrates 1 (linear) and 2 (competitive Langmuir). |
| `SolverModel` | 1 | Spatial discretization. The implemented backend is finite differences. Both kinetic models use BDF for time integration. |

## Select an injection

| InputType | Feed supplied to the inlet condition | Timing |
|---|---|---|
| 1 · Step | Each species receives its `Cfeed` throughout the run. | `tinj = tfinal` |
| 2 · Pulse | Each species receives its `Cfeed` before `tinj`, then zero feed. | `tinj = 60 * Vinj / Q` |

`Q` is in cm³/min, `Vinj` and column volume `Vc` are in cm³, and simulation time is in seconds. The run ends at `tfinal = 60 * CVfinal * Vc / Q`. For a pulse the builder also calculates `twash = tinj + 60 * CVwash * Vc / Q`. The switch to zero protein feed occurs at `tinj`. Model 3 sets the inlet concentration directly. In model 4 the specified feed enters a flux boundary, so the computed concentration at the inlet node can differ from the feed.

## Select elution conditions

| ElutionType | Program identified by the input code |
|---:|---|
| 1 | Isocratic |
| 2 | Salt gradient |
| 3 | pH gradient |
| 4 | Linear salt and pH gradients |
| 5 | Salt gradient with a nonlinear pH relationship |
| 6 | Salt gradient with an experimental pH profile |

For codes 2–6 the builder computes `tG = 60 * CVG * Vc / Q`, `Bs = (Csf - Cso) / tG` and `BpH = (pHf - pHo) / tG`. The Python pore models require spatial profiles supplied by the caller for **every** non-isocratic code. Model 3 requires both `inputs.salt_profile` and `inputs.pH_profile`. Model 4 requires `inputs.salt_profile`. To hold a required field fixed, supply a constant table.

Each profile is a tuple `(times_s, positions_cm, values)`. Use strictly monotonic time and position axes with at least two points each. The finite table must have shape `(len(times_s), len(positions_cm))` and span the full run and column. The pore framework supports interpolation for environment-aware equilibrium functions. The two currently wired isotherms use fixed parameters, so changing the profile alone does not change their binding equilibrium. The MATLAB program names in the table do not imply that Python generates these profiles from their endpoints.

## Select kinetics and equilibrium

| KineticModel | 3 · Concentration-driven | 4 · Loading-driven |
|---|---|---|
| Equilibrium query | Liquid concentration consistent with the current loading, `ceq` | Loading consistent with the current pore concentration, `qeq` |
| Adsorption rate | `dq/dt = kb * (c - ceq)` | `dq/dt = kb * (qeq - q)` |
| Particle uptake in column balance | `3 * De / Rp * (c_surface - c_adjacent) / dr` | `3 * kf / Rp * (cm - c_surface)` |
| Inlet | Prescribed feed | Dispersion-aware (Danckwerts) inlet |

The linear isotherm (ID 1) uses $q_{eq,i}=q_{m,i}K_i c_i$. Competitive Langmuir (ID 2) uses $q_{eq,i}=q_{m,i}K_i c_i/(1+\sum_j K_j c_j)$. Model 3 calls the inverse isotherm. With Langmuir it passes the current pore concentrations as context for the componentwise inverse.

### Species inputs

Define these on each `SpeciesParameters` passed to `inputs.AddSpecies`. Registration order controls the species axis of all results, legend entries and CSV columns.

| Field | Meaning and units |
|---|---|
| `Cfeed` | Concentration of this species in the feed, in the same units as liquid concentration `c`. |
| `K`, `qm` | Isotherm affinity and capacity or loading scale. Together they determine equilibrium. |
| `kb` | Adsorption rate coefficient. Model 4 uses s⁻¹. Model 3 uses loading per liquid concentration per second. A fitted value depends on the driving force. |
| `Do` | Molecular diffusivity in free solution (cm²/s), used to calculate film transfer. |
| `De` | Effective diffusivity within the particle (cm²/s). Must be positive. |
| `Dl` | Axial dispersion coefficient in the column (cm²/s). |
| `kc` | LDF rate parameter. The pore adsorption balances use `kb`. |

`Cfeed`, `Dl` and `kb` may be zero. `kf` (cm/s) is calculated for each species by the input builder. The values in the supplied experiment are example inputs.

### Other numbered modules

The repository retains the following IDs from the larger model catalogue. Selecting a numbered file alone does not supply a working simulation. For kinetic IDs outside 3 and 4, the current runner returns without a result. For unsupported pore isotherm or spatial adapters, validation or dispatch raises an error.

| Family | IDs and reserved implementations |
|---|---|
| Kinetics | 1 LDF loading-driven, 2 LDF concentration-driven, 5 two-site pore diffusion, 6 LDF with reaction |
| Isotherms | 3 ion-exchange mass action, 4 steric mass action, 5 linear steric mass action, 6 interpolated, 7 modified BET, 8 bimodal, 9 modified Langmuir |
| Isotherms | 10 modified logarithmic, 11 linear steric mass action with pH, 12 hydrophobic interaction chromatography, 13 HIC polynomial |
| Spatial methods | 2 finite volume, 3 finite element, 4 orthogonal collocation on finite elements, 5 WENO |

## Run a linear pulse from the REPL

After installing the checkout with `python -m pip install -e .`, run:

```python
from Column_Simulation.experiments.VarInputs import build_inputs
from Column_Simulation.src.inputs.VarError import ErrorDetection
from Column_Simulation.main import run

p = build_inputs()
p.InputType = 2
p.KineticModel = 4
p.IsothermModel = 1
ErrorDetection().ValidateAll(p)
p.CalculateAllVariables()
result = run(p, experiment_id="linear_pulse")
```

Recalculate after changing flow, grid, species or timing because `build_inputs()` has already derived the original experiment's values. This example writes to `experiments/run_linear_pulse/` by default. For a reproducible experiment source, make the edits in `VarInputs.py` before calling its builder. The exported source copy records that file, while changes made in a live REPL affect only the run.

See [the simulation algorithm](AlgorithmicLoop.md) for state dimensions, boundary conditions and output fields.
