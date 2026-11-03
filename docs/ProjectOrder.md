Structure of the projects

```

bow --
     | - results
     |     | - folder per project
     |                | - inputs.py
     |                | - inputs.json
     |                | - ouputs.csv
     |                | - chromatograms.png
     |
     | - src
          | - inputs : 
          |     | -  CalculateAllVariables.py
          |     | -  VarError.py 
          |     | -  VarInputs.py
          |
          | - isotherms 
          |      | - hic_polynomial.py
          |      | - hydrophobic_interaction_chromatography.py
          |      | - interpolated_isotherm.py
          |      | - ion_exchange_mass_action.py
          |      | - langmuir.py
          |      | - linear.py
          |      | - modified_bet.py
          |      | - modified_langmuir.py
          |      | - modified_logarithmic_model.py
          |      | - steric_action.py
          |      | - steric_mass_action_linear_equation.py
          |      | - steric_mass_action_linear_pH.py
          |
          | - kinetics 
          |       | - LDF_ceq_Framework.py
          |       | - LDF_ceq_Framework_reaction.py
          |       | - LDF_qeq_Framework.py
          |       | - PoreDiffusion_ceq_Framework.py
          |       | - PoreDiffusion_ceq_Framework_twosite.py
          |       | - PoreDiffusion_qeq_Framework.py
          |
          | - misc : additionnal functions (mainly data science stuffs) to treat, visualize computations results
          | - solvers
                | - FineElement.py 
                | - Weno.py

```
