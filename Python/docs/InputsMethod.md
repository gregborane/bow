# Model selection for simulations

## Injection method

To select an injection method InputType variable must be defined. 

Allowed values are : 

InputType == 
        
        1 : Step elution
        2 : Pulse Elution    

2 are implemented : 

### Step is defined as is: 

InputType == 1 Constant injection for a certain amount of time

```
y (concentration)
 
|
|
|
|
|            -----------------------   
|            |
|            |
------------------------------------ x (time)
```

For a step input, the feed concentration is maintained throughout the chosen simulation horizon. The current python code implements this by setting `tinj = tfinal` when `InputType == 1`. At the inlet, the protein concentration is set to `Cfeed` until `tinj`.

**Relevant parameters:** `Cfeed` (feed protein concentration), `tinj` (injection end time) and `tfinal` (simulation end time).

### Pulse is defined as is : 

InputType == 2 Constant injection for a certain amount of time

```
y (concentration)
 
|
|
|
|
|            ----
|            |  |
|            |  |
------------------------------------ x (time)
```

A pulse is a finite injection: protein enters at the feed concentration for a limited time, then the inlet protein concentration returns to zero. The pulse duration is computed from the injection volume and flow rate:

$$
t_{inj}=60\,\frac{V_{inj}}{Q}
$$

The factor `60` converts minutes to seconds when `Q` is expressed in volume/minute. The code also calculates `twash`, the nominal end of the wash step. The injection profile itself is applied inside the transport framework as `cm(1,:) = Cfeed` before `tinj` and `cm(1,:) = 0` afterward.

**Relevant parameters:** `Vinj` (injected volume), `Q` (flow rate), `CVwash` (wash in column volumes), `Vc` (column volume), `tinj`, and `twash`.

**Note:** `InputType` describes how protein enters the column, not how it is subsequently eluted. Injection and elution are independent settings.

## Elution method

To select an injection method ElutionType variable must be defined. 

Allowed values are :

ElutionType == 

        1 : Isocratic
        2 : Salt gradient only
        3 : pH gradient only
        4 : Linear salt and pH gradient
        5 : Salt gradient with a non linear pH relationship
        6 : Linear salt gradient with an experimental pH profile

### Isocratic : 

ElutionType == 1 : Every parameters are constant

In an isocratic run, the elution conditions are held fixed rather than applying a programmed salt or pH gradient. The main controller assigns a constant salt value `Cs = Cso`. This does **not** imply that protein concentration stays constant: the injected protein can still move, adsorb, and elute over time.

### Gradient : 

ElutionType == 2,3,4 : pH, Salt concentration or both vary linearly

ElutionType == 5,6 : pH, Salt concentration, one vary linearly and the other non linearly

In `CalculateAllVariables.py`, the programmed gradients use the duration `tG`, the initial/final salt concentrations `Cso`, `Csf`, and initial/final pH values `pHo`, `pHf`:

$$
t_G=60\,\frac{CV_G\,V_c}{Q},\qquad
B_s=\frac{C_{sf}-C_{so}}{t_G},\qquad
B_{pH}=\frac{pH_f-pH_o}{t_G}
$$

The controller then computes salt and/or pH profiles along the column using `LDF_NonBinding_Framework` (or `LDF_NonBinding_Framework_ExppH` for the experimental-pH case). These profiles are used by adsorption models that depend on salt or pH.

- **Type 2:** salt gradient; no varying pH profile is passed by the controller.
- **Type 3:** pH gradient; salt remains at `Cso` in the controller.
- **Type 4:** both gradients are computed from the linear-gradient setup.
- **Type 5:** the salt profile is computed, while the pH state is converted using the supplied piecewise-polynomial relationship `pH_calc`.
- **Type 6:** the salt profile is computed and pH uses an experimental-profile framework.

## Kinetic model for transport inside the chromatograpghic system

To select a transport method SolverType variable must be defined. 

Allowed values are : 

SolverType ==

        1 : LDF q driving force
        2 : LDF c driving force
        3 : Pore diffusion c equilibrium
        4 : Pore diffusion q equilibrium
        5 : Two site pore diffusion
        6 : LDF with kinetics

### Linear Driving Force (LDF)

SolverType == 1,2,6 : 

The Linear Driving Force model approximates mass transfer using a difference between the current state and its equilibrium state. It does not explicitly discretize the interior of each particle in the radial direction.

Two common formulations explain the names used by the code:

- **q-driving force (solver 1):** the rate is expressed using the difference between equilibrium and actual adsorbed loading, conceptually `dq/dt = k(qeq - q)`.
- **c-driving force (solver 2):** the rate is expressed using a liquid-concentration driving force, conceptually `dq/dt = kb(c - ceq)`.
- **LDF with kinetics (solver 6):** the original controller calls `LDF_ceq_Framework_reaction`, suggesting an additional reaction term or reaction-related kinetic behavior. The exact equations cannot be established without that missing file.

Here, `c` is protein concentration in the liquid phase, `q` is adsorbed loading, `ceq` is the equilibrium liquid concentration for the current loading, and `qeq` is the equilibrium adsorbed loading for the current liquid concentration.

The uploaded main controller sends solvers `1`, `2`, and `6` to different functions: `LDF_qeq_Framework`, `LDF_ceq_Framework`, and `LDF_ceq_Framework_reaction`. These three framework files were not provided, so the equations above are explanatory formulations rather than verified transcriptions of those implementations.

### Pore diffusion model

SolverType == 3,4 : 

The pore diffusion model explicitly represents concentration gradients inside a porous particle. The column is discretized along its length (`Nz` grid points), and each particle is discretized radially (`Nr` grid points). If `Nc` protein components are modeled, the intraparticle arrays have shape `Nz × Nr × Nc`.

- **Solver 3 (`PoreDiffusion_ceq_Framework`):** uses a liquid-phase equilibrium concentration `ceq` obtained from `IsothermLibraryPore_c`.
- **Solver 4 (`PoreDiffusion_qeq_Framework`):** selects a q-equilibrium formulation, but its framework file was not provided, so its exact dynamic equations remain unverified.

For solver 3, the uploaded framework calculates the adsorption rate as:

$$
\frac{dq}{dt}=k_b\,(c-c_{eq})
$$

It also computes intraparticle transport, transfer across the external film, convection along the column, and axial dispersion. The external liquid concentration is `cm`; `c` and `q` refer to particle states. MATLAB's `ode15s` integrates the coupled differential equations after all state arrays are packed into a one-dimensional vector.

**Principal parameters:** `Rp`, `dr`, `dz`, `De`, `kf`, `kb`, `Dl`, `theta`, `Nz`, `Nr`, and `Nc`.

### Two site pore diffusion

SolverType == 5 :

This option selects the specialized `PoreDiffusion_ceq_Framework_twosite` solver. Unlike the ordinary pore formulation, its output-handling code stores two distinct intraparticle adsorbed-state arrays, `q1_intra` and `q2_intra`, in addition to `c_intra`.

A two-site equilibrium model often describes two classes of adsorption sites with different capacities or affinities. For a simple independent-site Langmuir example:

$$
q_{eq}(c)=\frac{q_{m1}K_1c}{1+K_1c}+\frac{q_{m2}K_2c}{1+K_2c}
$$

**Important:** that expression illustrates the scientific idea, but it is not verified as the actual two-site solver equation. The function `PoreDiffusion_ceq_Framework_twosite.m` was not included. In this codebase, solver `5` is therefore best understood as a **specialized solver selection**, not as an isotherm ID.

## Isotherm for binding interaction

To select an isotherm IsothermModel variable must be defined.

Allowed values are : 

IsothermModel ==  
       
        1 : Linear
        2 : Langmuir
        3 : Ion Exchange mass action
        4 : Steric mass action
        5 : Steric mass action linear equation
        6 : Interpolated isotherm
        7 : Modified BET
        8 : Bimodal isotherm
        9 : Modified Langmuir
        10 : Modified logarithmic model
        11 : Steric mass action linear pH
        12 : hydrophobic interaction chromatography
        13 : HIC polynomial

### 1 : Linear


The linear isotherm assumes that the amount bound to the stationary phase is proportional to protein concentration in the liquid phase:

$$q_{eq}=Hc$$

`H` is a partition/affinity coefficient; it does not represent saturation at high protein concentrations. In the uploaded LDF libraries, the linear relationship uses `qm` and `K` (`ceq = q/(qm*K)` and the corresponding q-formulation). **Caution:** the pore libraries use different special-case definitions for `Linear`, so the actual coefficient is not consistent across all frameworks.

### 2 : Langmuir


The Langmuir isotherm models saturable adsorption onto a finite population of binding sites. The conventional one-component form is:

$$q_{eq}=\frac{q_m Kc}{1+Kc}$$

`qm` is the maximum adsorbed capacity and `K` is an affinity parameter. The uploaded MATLAB code also contains multi-component variants that account for competing proteins through a denominator containing the sum of their `K*c` contributions. Some LDF expressions additionally include a nonbinding contribution involving `kinf`.

### 3 : Ion Exchange mass action


The mass-action model describes adsorption coupled to displacement of counterions in ion-exchange chromatography. Adsorption depends on the salt concentration, the charge-related parameter `z`, the binding constant `K`, and the available ligand/counterion capacity `qo`.

In `IsothermLibraryLDF_c.m`, the model computes equilibrium liquid concentration from loading using a term of the form:

$$c_{eq}=\frac{q}{K\left(\frac{q_o-\sum_i z_i q_i}{C_{salt}}\right)^z+k_{inf}/\theta}$$

The q-formulation uses a numerical nonlinear solve. Only the LDF c/q libraries have an explicit `MassAction` switch branch in the provided files.

### 4 : Steric mass action


Steric Mass Action (SMA) extends the ion-exchange mass-action idea by including steric exclusion: an adsorbed protein can prevent neighboring ligand sites from being used by other proteins.

In addition to `K`, `z`, `qo`, and `Csalt`, it uses the steric shielding parameter `sig`. The free-capacity term depends on sums involving `(z + sig) * q` across protein components. This introduces competition for available adsorption sites. The pore q library explicitly rejects this isotherm and asks the user to use the c-equilibrium formulation instead.

### 5 : Steric mass action linear equation


This option uses a simplified, salt-dependent relationship that is **linear in protein concentration** while still changing with salt concentration. A representative expression from `IsothermLibraryLDF_q.m` is:

$$q_{eq}=\frac{c}{\theta}\left(A C_{salt}^{-z}+k_{inf}\right)$$

`A` is a fitted coefficient, `z` controls salt sensitivity, `theta` is the phase ratio, and `kinf` is the nonbinding retention contribution. The pore formulation implements a closely related relationship without the extra `kinf` term in the q branch. This is not identical to the full competitive SMA model.

### 6 : Interpolated isotherm


Instead of enforcing a fixed analytical isotherm formula, the model evaluates a supplied interpolation function built from experimental or fitted adsorption data.

The implementation uses function handles such as `F_c(Csalt, q)`, `F_q(Csalt, c)`, and, for pH-dependent cases, `F_c(Csalt, q, pH)`. In some multi-component cases it uses `F_c1` and `F_c2` to represent component coupling. The quality and valid range of the model depend on the interpolation data and the interpolation method defined elsewhere.

### 7 : Modified BET


This is a modified Brunauer–Emmett–Teller-type adsorption relationship. Unlike simple Langmuir adsorption, it can capture more complex concentration-dependent behavior and uses a nonlinear rational expression.

The implementation computes salt-dependent coefficients `qm`, `K`, and `K2` from parameters `A1`, `B1`, `C1`, `A2`, `B2`, `C2`, `A3`, `B3`, and `C3`. The expression also includes a pore-liquid/nonbinding contribution proportional to `ep_p * c`.

In the LDF c library, the equilibrium liquid concentration is found with `lsqnonlin`, rather than by a direct closed-form inversion. The model is implemented in LDF c and Pore q, but not in every other library.

### 8 : Bimodal isotherm


The bimodal isotherm is intended to describe adsorption behavior involving more than one concentration-dependent mechanism. In this code it combines an ion-exchange-like term `alpha` with a salt-dependent term `Beta`.

The LDF c implementation defines `alpha` using `K`, `qo`, `z`, `sig`, `Csalt`, and `ep_p`. It defines `Beta` using `A`, `B`, `C`, and `Csalt`, then solves a quadratic expression for equilibrium liquid concentration. A branch is present for this model in `IsothermLibraryLDF_c.m`; the other provided isotherm switches do not implement it.

### 9 : Modified Langmuir


This variant makes the effective Langmuir capacity and affinity depend on salt concentration. In the LDF c implementation:

$$q_m(C_{salt})=\frac{A_1}{1+\exp[B_1(C_{salt}-C_1)]}$$

$$K(C_{salt})=A_2\exp(-B_2 C_{salt})$$

It also defines a salt-dependent term `A` and solves a quadratic for `ceq`.


### 10 : Modified logarithmic model


The implementation calls this model `MLog`. It uses salt-dependent parameters `beta` and `qm` with a nonlinear, fractional-power adsorption expression:

$$q = q_m\frac{c^{\beta}}{1+c^{\beta}}+\varepsilon_p c$$

The model finds the equilibrium concentration `ceq` using `lsqnonlin`, because the inverse relationship is solved numerically. The fitted parameters include `A1`, `B1`, `A2`, `B2`, and `C2`, and the nonbinding contribution uses `ep_p`.


### 11 : Steric mass action linear pH


This option is intended to make a simplified SMA-like isotherm sensitive to pH as well as salt concentration. pH can affect protein net charge, so a charge-related parameter such as `z` may vary with pH.

The available code contains a commented candidate relation for `z2D` as a function of `pH_i` and an expression using `K`, `qo`, `z2D`, and `cs_i`. **However, the assignments to `K`, `qo`, and `z2D` are commented out.** As supplied, this local function does not define those required variables and should be considered incomplete until corrected and tested.

### 12 : hydrophobic interaction chromatography


The HIC model represents adsorption associated with hydrophobic interactions, with salt concentration affecting the effective retention or partition coefficient. It uses an empirical salt-dependent expression involving parameters `ko`, `A`, `B`, and `C`.

In the LDF c library, the code calculates:

$$k=k_m+\exp\left[\ln(k_o)-\frac{B\sqrt{C_{salt}/1000}}{1+C\sqrt{C_{salt}/1000}}+A C_{salt}/1000\right]$$

It then uses `m = k/theta` and `ceq = q/m`. In the pore c library the expression omits the additive `km`, so those two implementations are not exactly identical.

### 13 : HIC polynomial


This extends the HIC-related equilibrium calculation by introducing an additional polynomial/rational relationship. The implementation first computes a salt-dependent coefficient, then constructs coefficients from parameters `pa`, `pb`, `pc`, and `pd`.

The equilibrium concentration is recovered from a quadratic formula in `IsothermLibraryLDF_c.m`. Its intended empirical operating range and the physical interpretation of each polynomial parameter are not documented in the supplied files; those should be established from the original fitting workflow before using it for extrapolation.

---

## Implementation reference and important compatibility notes

| Draft ID | Documented model | Current Python `IsothermModel` ID |
|---:|----|
| 1 | Linear |  
| 2 | Langmuir | 
| 3 | Ion Exchange mass action |
| 4 | Steric mass action |
| 5 | Steric mass action linear equation |
| 6 | Interpolated isotherm |
| 7 | Modified BET |
| 8 | Bimodal isotherm | 
| 9 | Modified Langmuir (`MLang`) |
| 10 | Modified logarithmic model |
| 11 | Steric mass action linear pH |
| 12 | hydrophobic interaction chromatography |
| 13 | HIC polynomial |


### Actual solver IDs in the uploaded controller

`MainColumnSimulation_Input.m` uses `s.RateModel`, not `s.SolverType`, and has a **seventh** branch:

| Runtime ID (`RateModel`) | Called function |
|---:|---|
| 1 | `LDF_qeq_Framework` |
| 2 | `LDF_ceq_Framework` |
| 3 | `PoreDiffusion_ceq_Framework` |
| 4 | `PoreDiffusion_qeq_Framework` |
| 5 | `PoreDiffusion_ceq_Framework_twosite` |
| 6 | `LDF_ceq_Framework_reaction` |
| 7 | `PoreDiffusion_ceq_DualFlow_Framework` |


### Availability of isotherms across libraries

The following table records `switch` branches, not full verification of scientific or numerical correctness.

| MATLAB isotherm ID | LDF c | LDF q | Pore c | Pore q |
|---:|:---:|:---:|:---:|:---:|
| 1 | Yes | Yes | Yes | Yes |
| 2 | Yes | Yes | Yes | Yes |
| 3 | Yes | Yes | No | No |
| 4 | Yes | Yes | Yes | Explicit error |
| 5 | Yes | Yes | Yes | Yes |
| 6 | Yes | Yes | Yes | Yes |
| 7 | Yes | No | No | Yes |
| 8 | Yes | No | No | No |
| 9 | Yes (`MLang`) | Yes (`BiLangmuir`) | Yes (`BiLangmuir`) | Yes (`BiLangmuir`) |
| 10 | Yes | No | No | No |
| 11 | Yes (incomplete function) | No | No | No |
| 12 | Yes | No | Yes | No |
| 13 | Yes | No | No | No |

### Numerical and physical parameter glossary

| Variable | Meaning |
|---|---|
| `s` / `si` | Configuration structure; the name changes between functions |
| `Cfeed` | Protein concentration at the column inlet during injection |
| `Cso`, `Csf` | Initial and final salt concentrations |
| `pHo`, `pHf` | Initial and final pH values |
| `Vinj`, `Q` | Protein injection volume and volumetric flow rate |
| `L`, `dc`, `dp`, `Rp` | Column length, column diameter, particle diameter, particle radius |
| `Ac`, `Vc` | Column cross-sectional area and geometric volume |
| `e`, `theta` | Bed void fraction and solid-to-interstitial phase ratio |
| `u`, `v` | Superficial and interstitial velocities |
| `tinj`, `tfinal`, `tspan` | Injection end, simulation end, requested output times |
| `Nz`, `Nr`, `Nc`, `Nt` | Axial grid points, radial grid points, proteins, requested times |
| `dz`, `dr` | Axial and radial grid spacing |
| `c`, `q`, `cm` | Pore-liquid concentration, adsorbed loading, bulk-mobile concentration |
| `ceq`, `qeq` | Equilibrium pore-liquid or adsorbed concentration |
| `K`, `qm` | Affinity parameter and saturation capacity in relevant isotherms |
| `De`, `Do`, `kf`, `kb`, `Dl` | Effective pore diffusion, molecular diffusion, film transfer, adsorption kinetic coefficient, axial dispersion |
| `ep_p`, `kinf`, `km` | Porosity/nonbinding and retention-related parameters |
| `t`, `Yp`, `Cp` | Solver times, packed state trajectory, protein outlet history |

### Recommended simulation sequence

1. Define the numeric model options and the required physical parameters.
2. Calculate derived variables using `CalculateAllVariables(si)`.
3. Compute salt and/or pH profiles as required by `ElutionType`.
4. Select exactly one transport solver with `RateModel` (or a renamed `SolverType`).
5. Inside that solver, evaluate the chosen isotherm in the needed direction (`ceq` or `qeq`).
6. Integrate the coupled equations and extract outlet and intraparticle results.
