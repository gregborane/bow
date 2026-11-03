"""Validation of shared chromatography inputs and N species parameters."""

from numbers import Real

from CalculateAllVariables import SpeciesParameters, StoreInput


class ErrorDetection:
    """Validate numeric inputs, selected models, and every registered species."""

    def CheckIsNumber(self, input: StoreInput):
        """Check all supplied shared and species-specific input values."""
        invalid = {}

        for name in input.INPUT_NAMES:
            value = getattr(input, name)
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, Real)
            ):
                invalid[name] = value

        for species_name, parameters in input.species.items():
            if not isinstance(parameters, SpeciesParameters):
                invalid[f"species[{species_name!r}]"] = parameters
                continue
            for parameter_name in ("Do", "kc"):
                value = getattr(parameters, parameter_name)
                if value is not None and (
                    isinstance(value, bool) or not isinstance(value, Real)
                ):
                    invalid[f"species[{species_name!r}].{parameter_name}"] = value

        if invalid:
            raise ValueError(f"Input parameters must be numeric: {invalid}")

    def ValidateInput(self, input: int):
        """1 = step input; 2 = pulse input."""
        if input not in (1, 2):
            raise ValueError(f"InputType must be 1 or 2; got {input}")

    def ValidateElution(self, elution: int):
        """1 = isocratic; 2-6 = gradient variants."""
        if elution not in range(1, 7):
            raise ValueError(f"ElutionType must be between 1 and 6; got {elution}")

    def ValidateRate(self, solver: int):
        """1-6 = documented rate-model selections."""
        if solver not in range(1, 7):
            raise ValueError(f"SolverModel must be between 1 and 6; got {solver}")

    def ValidateIsotherm(self, isotherm: int, solver: int):
        """Keep the original allowed isotherms for each solver."""
        allowed_by_solver = {
            1: set(range(1, 15)),
            2: {1, 2, 3, 4, 5, 6, 9},
            3: {1, 2, 4, 5, 6, 9, 12},
            4: {1, 2, 5, 6, 7, 9},
        }
        allowed = allowed_by_solver.get(solver)
        if allowed is not None and isotherm not in allowed:
            raise ValueError(
                f"IsothermModel {isotherm} is not supported by SolverModel "
                f"{solver}; allowed values: {sorted(allowed)}"
            )
        # As in the original project, no combinations are specified for solvers 5 or 6.

    def CheckRequiredInputs(self, input: StoreInput):
        """Require shared and species-specific inputs for the selected models."""
        if not input.species:
            raise ValueError("At least one species must be added to the simulation")

        required = [
            *input.COLUMN_INPUTS,
            "Q", "CVfinal", "Nt", "Nz", "Nr", "InputType",
            "ElutionType", "SolverModel", "IsothermModel",
        ]
        if input.InputType == 2:
            required.extend(("Vinj", "CVwash"))
        if input.ElutionType != 1:
            required.extend(("CVG", "Cso", "Csf", "pHo", "pHf"))
        if input.SolverModel in (3, 4, 5):
            required.append("nu")

        missing = [name for name in required if getattr(input, name) is None]
        if missing:
            raise ValueError(f"Missing required shared inputs: {', '.join(missing)}")

        for species_name, parameters in input.species.items():
            if not isinstance(parameters, SpeciesParameters):
                raise TypeError(
                    f"Species {species_name!r} must use SpeciesParameters"
                )
            species_required = []
            if input.SolverModel == 2:
                species_required.append("kc")
            elif input.SolverModel in (3, 4, 5):
                species_required.append("Do")

            species_missing = [
                name for name in species_required
                if getattr(parameters, name) is None
            ]
            if species_missing:
                raise ValueError(
                    f"Species {species_name!r} is missing required inputs: "
                    f"{', '.join(species_missing)}"
                )

    def CheckRanges(self, input: StoreInput):
        """Check physical and numerical ranges for shared and species inputs."""
        for name in ("dp", "dc", "L", "Q"):
            if getattr(input, name) <= 0:
                raise ValueError(f"{name} must be positive")
        if not 0 < input.e < 1:
            raise ValueError("e must be between 0 and 1 (exclusive)")
        for name in ("ep_s", "ep_p", "CVfinal"):
            if getattr(input, name) < 0:
                raise ValueError(f"{name} cannot be negative")
        for name in ("Nt", "Nz", "Nr"):
            value = getattr(input, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 2:
                raise ValueError(f"{name} must be an integer >= 2")

        if input.InputType == 2 and (input.Vinj < 0 or input.CVwash < 0):
            raise ValueError("Vinj and CVwash cannot be negative")
        if input.ElutionType != 1 and input.CVG <= 0:
            raise ValueError("CVG must be positive for gradient elution")
        if input.SolverModel in (3, 4, 5) and input.nu <= 0:
            raise ValueError("nu must be positive")

        for species_name, parameters in input.species.items():
            if input.SolverModel == 2 and parameters.kc < 0:
                raise ValueError(f"kc for species {species_name!r} cannot be negative")
            if input.SolverModel in (3, 4, 5) and parameters.Do <= 0:
                raise ValueError(f"Do for species {species_name!r} must be positive")

    def ValidateAll(self, input: StoreInput):
        """Run all validation checks and return the validated StoreInput."""
        self.CheckIsNumber(input)
        self.CheckRequiredInputs(input)
        self.ValidateInput(input.InputType)
        self.ValidateElution(input.ElutionType)
        self.ValidateRate(input.SolverModel)
        self.ValidateIsotherm(input.IsothermModel, input.SolverModel)
        self.CheckRanges(input)
        return input

