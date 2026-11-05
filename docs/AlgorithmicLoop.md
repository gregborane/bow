# How BOW runs a column simulation

`Column_Simulation.main.run(inputs)` takes a prepared experiment, integrates one column trajectory, exports its breakthrough curve, and returns the numerical fields. The repeated numerical work happens inside the stiff ODE integrator. The `Nt` requested output times specify when to report the solution; BDF chooses its own internal steps.

## From experiment to model

| Stage | Owner | Result |
|---|---|---|
| 1. Define and validate | `experiments/VarInputs.py`, `src/inputs/VarError.py` | One `StoreInput` with shared column, solution and operating conditions plus a named entry per species. |
| 2. Derive parameters | `src/inputs/CalculateAllVariables.py` | Column volume, velocity, phase ratio, time grid, injection schedule, axial and radial spacing, and species film coefficients. |
| 3. Select modules | `main.selected()`, `main.run()` | Resolve the numeric prefix within each installed package. Construct the isotherm, spatial operators and selected kinetic model. |
| 4. Integrate | `PoreModel.simulate()` and the selected model's rate and boundary methods | Start from an empty column, integrate each feed interval and reconstruct the fields at the requested times. |
| 5. Export | `src/misc/results_output.py` | Use the returned time and outlet arrays for both CSV and PNG. With terminal `--output`, the CLI also writes the full fields to NPZ. |

The axial grid goes from the inlet at 0 to the outlet at `L`. The radial grid goes from the particle surface at `Rp` to its centre at 0. Both grids include their endpoints and require at least three nodes. The input store expresses `Rp` in µm. The spatial operators use cm for both axes, so the runner multiplies `Rp` by $10^{-4}$.

## Inside one BDF evaluation

1. **Unpack the trial state.** Split the one-dimensional vector into independent pore-liquid concentrations, adsorbed loadings and bulk concentrations. Species occupy the final axis.
2. **Reconstruct the algebraic boundaries.** Model 3 sets the bulk inlet to the current feed. Model 4 solves $v(c_{in}-c_{feed})=D_l(c_{next}-c_{in})/\Delta z$. Both copy the adjacent bulk value to the outlet for zero axial gradient. At each axial position the surface concentration satisfies $D_e(c_s-c_a)/\Delta r=k_f(c_m-c_s)$, where $c_a$ is the first interior particle concentration.
3. **Evaluate adsorption.** Model 3 computes the concentration driving force from inverse equilibrium. Model 4 computes the loading driving force from forward equilibrium. Each returns $\dot q$ for all radial nodes.
4. **Evaluate transport.** The FDM backend uses upwind axial convection, centred interior axial dispersion, and a spherical shell flux balance. At the centre the inner face has zero area and the limiting diffusion term is $6D_e(c_{neighbor}-c_{center})/\Delta r^2$. The first dynamic shell includes the surface shell's volume.
5. **Assemble the derivatives.** For independent particle liquid nodes, $\dot c=(\mathcal D_r(c)-\dot q)/\epsilon_p$. For independent bulk nodes, $\dot c_m=-v\,\partial_z c_m+D_l\,\partial_{zz}c_m-\theta S$, with $\theta=(1-e)/e$. Model 3 obtains $S$ from the surface gradient; model 4 uses film flux. Pack the derivatives in the same order as the state and return them to BDF.

### Packed state versus reported fields

| Block | Stored for integration | Returned at each output time |
|---|---|---|
| Pore liquid `c` | `(Nz, Nr-1, Nc)` | `(Nz, Nr, Nc)`, including the reconstructed surface at index 0 |
| Adsorbed loading `q` | `(Nz, Nr, Nc)` | `(Nz, Nr, Nc)`; the surface loading is tracked at a node with no independent liquid storage |
| Bulk liquid `cm` | `(Nz-2, Nc)` | `(Nz, Nc)`, including both reconstructed column endpoints |

The blocks are flattened in NumPy C order and concatenated as `c, q, cm`. Here `Nc = len(inputs.species)`. These shapes matter when adding a model because integrating an algebraic boundary as an independent unknown changes the system.

## Injection switching and time integration

The independent initial state is zero. A step has one feed interval. A pulse ending strictly inside the run has two intervals separated at `tinj`. The first solve advances exactly to the switch, even when it falls between requested output times. Its final state initializes the zero-feed interval. When output is requested at the pulse edge, BOW reconstructs its boundary using the new feed. For model 4, zero feed can still give a nonzero inlet-node concentration through the dispersion term.

Each interval uses SciPy `solve_ivp(method="BDF", rtol=1e-6, atol=1e-10)`. The finite-difference backend supplies local stencil information, which the model uses to construct a sparse Jacobian pattern. On integration failure the model raises an error before export.

## Read the result

| Attribute | Shape | Interpretation |
|---|---|---|
| `result.t` | `(Nt,)` | Requested time samples, seconds |
| `result.c`, `result.q` | `(Nt, Nz, Nr, Nc)` | Particle liquid concentration and adsorbed loading |
| `result.cm` | `(Nt, Nz, Nc)` | Bulk liquid concentration along the column |
| `result.breakthrough` | `(Nt, Nc)` | Outlet concentration, exactly `result.cm[:, -1, :]` |

The exporter validates time, species names and finite signal values, then uses one table for both files. The CSV starts with `time` and has one `outlet_concentration_<species>` column per protein in registration order. The PNG shows a line per protein against time in seconds. BOW reports concentration directly; detector response would need a separate calibration.

The default destination is `Column_Simulation/experiments/run_1/`. An experiment ID chooses `run_<id>/`. An explicit `output_dir` chooses the exact directory, while `BOW_OUTPUT_DIR` overrides the parent of the run folders. Repeating an ID replaces the CSV, PNG and experiment source copy. For optional NPZ export, a relative `--output` filename follows the separate full-field output policy in [the README](../readme.md).

## Reproduce the workflow

From the checkout root, after editable installation:

```sh
python -m Column_Simulation.main --experiment-id baseline
python -m unittest discover -s tests -v
```

For an interactive run, start with [the model selection guide](InputsMethod.md). It includes a complete pulse example and explains which numbered models currently run.
