# BOW software 

A simulation chromatographic software. 

Very similar to [CADET](https://github.com/cadet/) in most aspects.
But way more inferior in terms of speed, optimisation and modularity.

The main objective is to have a quick and easily modifiable project. 

Developed in python to quickly implement isotherms and kinetics models.

Eventual port to golang will be envisaged for faster computation while keeping
easiness of implementation.

## Installation and use

Install from the checkout root with `python -m pip install -e .`. Edit
`Column_Simulation/experiments/VarInputs.py` to select the inputs and set
`EXPERIMENT_ID`, `OUTPUT_DIR`, `EXPORT_CONFIG`, and `OPEN_EXTERNAL_DISPLAY`.
From the checkout root, run:

```sh
python Column_Simulation/experiments/VarInputs.py
```

The default run writes `run.csv`, `run.png`, and a copy of `VarInputs.py` to
`./run_1/`. Relative output directories start at the working directory.
The installed `bow-sim` command also accepts `--experiment-id` and `--output`.
For the optional NPZ export, relative `--output` filenames start in
`BOW_OUTPUT_DIR` (default: the checkout's `results/` directory).

## Goals : 

Next in line is incorporation of : 

- Linear & Langmuir Isotherm
- Pore diffusion c\q equilibrium
- Finite difference method

### 1. Inputs / Ouputs Builders

#### 1. Inputs

| Inputs Builders               | Implemented | GoLang |
|:------------------------------|:-----------:|:------:|
| Experimental parameters setup |      Y      |   N    |
| Dependent parameters          |      Y      |   N    |
| Error checks                  |      Y      |   N    |
| Template builders             |      Y      |   N    |

#### 2. Outputs

| Outputs Builders              | Implemented | GoLang |
|:------------------------------|:-----------:|:------:|
| Save InputFile                |      Y      |   N    |
| Save chromatogram             |      Y      |   N    |
| Save results table            |      Y      |   N    |

### 2. Isotherms

| Isotherms                                   | Implemented | GoLang |
|:--------------------------------------------|:-----------:|:------:|
| 1\. Linear                                  |      Y      |   N    |
| 2\. Langmuir                                |      Y      |   N    |
| 3\. Ion Exchange mass action                |      N      |   N    |
| 4\. Steric mass action                      |      N      |   N    |
| 5\. Steric mass action linear equation      |      N      |   N    |
| 6\. Interpolated isotherm                   |      N      |   N    |
| 7\. Modified BET                            |      N      |   N    |
| 8\. Bimodal isotherm                        |      N      |   N    |
| 9\. Modified Langmuir                       |      N      |   N    |
| 10\. Modified logarithmic model             |      N      |   N    |
| 11\. Steric mass action linear pH           |      N      |   N    |
| 12\. Hydrophobic interaction chromatography |      N      |   N    |
| 13\. HIC polynomial                         |      N      |   N    |

### 3. Kinetics

| Kinetics                         | Implemented | GoLang |
|:---------------------------------|:-----------:|:------:|
| 1\. LDF q driving force          |      N      |   N    |
| 2\. LDF c driving force          |      N      |   N    |
| 3\. Pore diffusion c equilibrium |      Y      |   N    |
| 4\. Pore diffusion q equilibrium |      Y      |   N    |
| 5\. Two site pore diffusion      |      N      |   N    |
| 6\. LDF with kinetics            |      N      |   N    |

### 4. Solvers

| Solvers                                   | Implemented | GoLang |
|:------------------------------------------|:-----------:|:------:|
| Finite Difference Method                  |      Y      |   N    |
| Finite Volume Method                      |      N      |   N    |
| Finite Element Method                     |      N      |   N    |
| Orthogonal Collocation on Finite Elements |      N      |   N    |
| Weighted Essentially Non-Oscillatory      |      N      |   N    |

### 5. Software wise

| User Experience                           | Implemented | GoLang |
|:------------------------------------------|:-----------:|:------:|
| GUI                                       |      N      |   N    |
| Use a chosen file for inputs              |      N      |   N    |
| Use fail safe for inputs                  |      Y      |   N    |
| Use .json file for inputs                 |      N      |   N    |
| Custom save location for ouputs           |      Y      |   N    |
| Use fail safe location for ouputs         |      Y      |   N    |
| Custom format for output image            |      N      |   N    |
