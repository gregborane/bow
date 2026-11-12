"""Store shared chromatography inputs and calculate variables for N species.

There is one StoreInput per physical simulation. Column, solution, operating,
and numerical settings are shared. Species-specific inputs and derived
mass-transfer variables are stored in dictionaries keyed by species name.
"""

from collections.abc import Mapping
from dataclasses import dataclass, fields
from numbers import Real

import numpy as np


@dataclass(frozen=True)
class ColumnParameters:
    """Column geometry, particle geometry, and stationary-phase parameters."""

    dp:  float | None = None              # Particle diameter (m)
    dc:  float | None = None              # Column diameter (m)
    L:  float | None = None               # Column length (m)
    e:  float | None = None               # Interstitial/column porosity (-)
    ep_s:  float | None = None            # Salt-related phase parameter (-)
    ep_p:  float | None = None            # Particle-related phase parameter (-)


@dataclass(frozen=True)
class SolutionParameters:
    """Mobile-phase composition and physical properties."""

    Cso:  float | None = None    # Initial salt concentration (mol/m³)
    Csf:  float | None = None    # Final salt concentration (mol/m³)
    pHo:  float | None = None    # Initial pH
    pHf:  float | None = None    # Final pH
    nu:  float | None = None     # Kinematic viscosity (m²/s)


@dataclass(frozen=True)
class SpeciesParameters:
    """Input properties belonging to one protein or other solute."""

    Do:  float | None = None     # Molecular diffusivity (m²/s)
    kc:  float | None = None     # LDF mass-transfer rate constant (1/s)
    K: float | None = None       # Thermodynamic equilibrium constant ka / kb
    qm : float | None = None     # Maximum quantity adsorbed on materials

    Cfeed: float | None = None   # Feed concentration, same units as mobile-phase c
    kb: float | None = None      # Adsorption rate coefficient (driving-force dependent)
    Dl: float | None = None      # Axial dispersion coefficient (m²/s)
    De: float | None = None      # Effective pore diffusivity (m²/s)
    kf: float | None = None      # film diffusion (m/s)


@dataclass
class SpeciesResults:
    """Derived mass-transfer variables for one species."""

    n:  float | None = None     # Dimensionless LDF parameter
    Sc:  float | None = None    # Schmidt number
    Sh:  float | None = None    # Sherwood number
    kf:  float | None = None    # External film mass-transfer coefficient (m/s)


class StoreInput:
    """One simulation input store containing shared settings and N species.

    Shared values retain the original flat interface, for example
    ``inputs.dp``, ``inputs.Q``, and ``inputs.tfinal``.

    Species-specific values are accessed by name:

        inputs.species["Protein_A"].Do
        inputs.species_results["Protein_A"].kf

    Use ``AddSpecies`` to add a species or ``SetSpecies`` to replace the whole
    species collection. Calculations of shared variables happen once, while
    species-dependent transport variables are calculated for every species.
    """

    # Shared input groups. Species inputs are deliberately NOT part of INPUT_NAMES.
    COLUMN_INPUTS = ("dp", "dc", "L", "e", "ep_s", "ep_p")
    SOLUTION_INPUTS = ("Cso", "Csf", "pHo", "pHf", "nu")
    SIMULATION_INPUTS = (
        "Q", "CVfinal", "Nt", "Nz", "Nr", "InputType", "ElutionType",
        "SolverModel", "KineticModel", "IsothermModel", "Vinj", "CVwash", "CVG",
    )
    INPUT_NAMES = COLUMN_INPUTS + SOLUTION_INPUTS + SIMULATION_INPUTS

    # Variables calculated once per simulation (not once per species).
    DERIVED_NAMES = (
        "Rp", "Ac", "Vc", "u", "v", "tfinal", "tspan", "tinj", "twash",
        "theta", "km", "kinf", "dz", "dr", "tG", "Bs", "BpH", "gamma",
        "Re", "z", "r"
    )

    species: dict[str, SpeciesParameters]
    species_results: dict[str, SpeciesResults]

    def __init__(self):
        # Shared inputs and derived values start unset. The lists are also used
        # by validation, so each shared attribute has one central declaration.
        # Class-level annotations make shared attributes visible to static analysis.
        # Column and particle parameters
        self.dp: float       # m
        self.dc: float       # m
        self.L: float        # m
        self.e: float        # dimensionless
        self.ep_s: float     # dimensionless
        self.ep_p: float     # dimensionless

        # Solution parameters
        self.Cso: float      # mol/m³
        self.Csf: float      # mol/m³
        self.pHo: float      # dimensionless (pH)
        self.pHf: float      # dimensionless (pH)
        self.nu: float       # m²/s

        # Operating parameters
        self.Q: float           # m³/s
        self.CVfinal: float     # dimensionless
        self.Nt: int            # number of time points
        self.Nz: int            # number of axial points
        self.Nr: int            # number of radial points
        self.InputType: int     # model selection code
        self.ElutionType: int   # model selection code
        self.KineticModel: int  # model selection code
        self.SolverModel: int   # model selection code
        self.IsothermModel: int # model selection code
        self.Vinj: float        # m³
        self.CVwash: float      # dimensionless
        self.CVG: float         # dimensionless

        # Geometry and calculated variables
        self.Rp: float          # m
        self.Ac: float          # m²
        self.Vc: float          # m³
        self.u: float           # m/s
        self.v: float           # m/s
        self.tfinal: float      # s
        self.tspan: np.ndarray  # s
        self.tinj: float        # s
        self.twash: float       # s

        # Phase parameters
        self.theta: float       # dimensionless
        self.km: float          # dimensionless
        self.kinf: float        # dimensionless

        # Numerical grid
        self.dz: float          # m
        self.dr: float          # m
        self.z: float           # m
        self.r: float           # m

        # Elution parameters
        self.tG: float          # s
        self.Bs: float          # mol/(m³.s)
        self.BpH: float         # pH/s (formally s⁻¹)
        self.gamma: float       # mol/m³

        # Dimensionless transport parameter
        self.Re: float          # dimensionless

        for name in self.INPUT_NAMES + self.DERIVED_NAMES:
            setattr(self, name, None)

        # Species are stored once in this simulation, not as duplicated models.
        self.species = {}
        self.species_results = {}

    def __getattr__(self, name):
        """Expose pore inputs as (Nc,) vectors in species registration order.

        Called only when an ordinary attribute is absent, so older callers
        that explicitly assigned a scalar or vector on the store still work.
        """
        if name not in ("Cfeed", "kb", "Dl", "De"):
            raise AttributeError(name)
        values = []
        for species_name, parameters in self.species.items():
            value = getattr(parameters, name)
            if (isinstance(value, bool) or not isinstance(value, Real)
                    or not np.isfinite(value) or value < 0 or (name == "De" and value == 0)):
                requirement = "positive" if name == "De" else "nonnegative"
                raise ValueError(f"{species_name}.{name} must be a finite {requirement} number")
            values.append(value)
        return np.asarray(values, dtype=float)

    def _set_group(self, parameters):
        """Copy fields from a column or solution parameter dataclass."""
        for field in fields(parameters):
            setattr(self, field.name, getattr(parameters, field.name))
        return self

    def SetColumn(self, parameters: ColumnParameters):
        """Set the shared column and stationary-phase inputs."""
        if not isinstance(parameters, ColumnParameters):
            raise TypeError("parameters must be a ColumnParameters instance")
        return self._set_group(parameters)

    def SetSolution(self, parameters: SolutionParameters):
        """Set the shared mobile-phase inputs."""
        if not isinstance(parameters, SolutionParameters):
            raise TypeError("parameters must be a SolutionParameters instance")
        return self._set_group(parameters)

    def AddSpecies(self, name: str, parameters: SpeciesParameters):
        """Add one named species to this simulation.

        Raises an error for duplicate names to prevent accidental replacement.
        Use ``SetSpecies`` when intentionally replacing the full collection.
        """
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Species name must be a non-empty string")
        if name in self.species:
            raise ValueError(f"Species {name!r} already exists")
        if not isinstance(parameters, SpeciesParameters):
            raise TypeError("parameters must be a SpeciesParameters instance")

        self.species[name] = parameters
        self.species_results.pop(name, None)
        return self

    def SetSpecies(self, species: Mapping[str, SpeciesParameters]):
        """Replace the complete species collection with a name-to-input mapping."""
        if not isinstance(species, Mapping):
            raise TypeError(
                "SetSpecies expects a mapping such as "
                "{'Protein_A': SpeciesParameters(...)}; use AddSpecies for one species"
            )

        # Build a temporary mapping first so a bad entry cannot partially replace it.
        new_species: dict[str, SpeciesParameters] = {}
        for name, parameters in species.items():
            if not isinstance(name, str) or not name.strip():
                raise ValueError("Every species name must be a non-empty string")
            if not isinstance(parameters, SpeciesParameters):
                raise TypeError(
                    f"Parameters for species {name!r} must be SpeciesParameters"
                )
            new_species[name] = parameters

        self.species = new_species
        self.species_results = {}
        return self

    def RemoveSpecies(self, name: str):
        """Remove a species and its calculated variables."""
        if name not in self.species:
            raise KeyError(f"Unknown species {name!r}")
        del self.species[name]
        self.species_results.pop(name, None)
        return self

    def CalculateAllVariables(self):
        """Calculate shared variables once and transport variables per species."""
        if not self.species:
            raise ValueError("Add at least one species before calculating variables")
 
        # Reset previous results so values from an earlier mode cannot survive
        # Fail safe before recalling any new inputs.
        self._reset_derived_variables()
        
        # Create system
        self._calculate_geometry()
        self._calculate_flow_and_time()
        self._calculate_phase_parameters()
        self._calculate_spatial_grid()
        self._calculate_injection()
        self._calculate_elution()
        self._calculate_mass_transfer()
        return self

    def _reset_derived_variables(self):
        """Reset shared and species-specific results before recalculation."""
        
        for name in self.DERIVED_NAMES:
            setattr(self, name, None)
        
        self.species_results.clear()

    def _calculate_geometry(self):
        self.Rp = self.dp / 2                        # m
        self.Ac = np.pi * (self.dc / 2) ** 2         # m²
        self.Vc = self.L * self.Ac                   # m³

    def _calculate_flow_and_time(self):
        self.u = self.Q / self.Ac             # m/s
        self.tfinal = self.CVfinal * self.Vc  / self.Q
        self.tspan = np.linspace(0, self.tfinal, self.Nt)
        self.v = self.u / self.e                     # m/s

    def _calculate_phase_parameters(self):
        self.theta = (1 - self.e) / self.e
        self.km = (1 - self.e) * self.ep_s / self.e
        self.kinf = self.theta * self.ep_p

    def _calculate_spatial_grid(self):
        self.z = np.linspace(0, self.L, self.Nz)
        self.r = np.linspace(self.Rp, 0, self.Nr)
        self.dz = self.L / (self.Nz - 1)             # m
        # Rp is stored in m; dr is required in m.
        self.dr = (self.Rp) / (self.Nr - 1)          # m

    def _calculate_injection(self):
        if self.InputType == 1:                      # Step input
            self.tinj = self.tfinal
        else:                                        # Pulse input
            self.tinj = self.Vinj / self.Q
            self.twash = self.tinj + self.CVwash * self.Vc  / self.Q

    def _calculate_elution(self):
        if self.ElutionType != 1:                    # Gradient elution
            self.tG = self.CVG * self.Vc / self.Q
            self.Bs = (self.Csf - self.Cso) / self.tG     # mol/(m³.s)
            self.BpH = (self.pHf - self.pHo) / self.tG    # pH units/s
            self.gamma = self.Bs * self.L / self.v

    def _calculate_mass_transfer(self):
        """Calculate transport results once for each registered species."""
        dp = self.dp  # Particle diameter converted from m to m.

        # Reynolds number depends on shared column, flow, and solution data.
        if self.KineticModel in (3, 4, 5):
            self.Re = self.u * dp / self.nu

        for name, parameters in self.species.items():
            # Create the result object here, where its values are calculated.
            results = SpeciesResults()

            if self.KineticModel == 2:                # Linear Driving Force
                results.n = self.theta * parameters.kc * self.L / self.v

            elif self.KineticModel in (3, 4, 5):      # Pore diffusion
                results.Sc = self.nu / parameters.Do
                results.Sh = (
                    1.09 / self.e
                    * self.Re ** 0.33
                    * results.Sc ** 0.33
                )
                results.kf = (results.Sh * parameters.Do / dp if parameters.kf is None
                              else parameters.kf)

            self.species_results[name] = results

    def GetSpeciesResults(self, name: str) -> SpeciesResults:
        """Return calculated transport variables for a named species."""
        if name not in self.species_results:
            raise KeyError(
                f"No calculated results for species {name!r}; "
                "call CalculateAllVariables() first"
            )
        return self.species_results[name]


