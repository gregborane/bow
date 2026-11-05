"""Filesystem policy for BOW. Module imports always use the installed package."""
import os
import tomllib
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
SOURCE_ROOT = PACKAGE_ROOT / "src"
EXPERIMENTS_DIR = PACKAGE_ROOT / "experiments"
ISOTHERMS_DIR = SOURCE_ROOT / "isotherms"
KINETICS_DIR = SOURCE_ROOT / "kinetics"
SOLVERS_DIR = SOURCE_ROOT / "solvers"


def _project_root():
    override = os.environ.get("BOW_PROJECT_ROOT")
    root = Path(override).expanduser() if override is not None else PACKAGE_ROOT.parent
    if not root.is_absolute():
        raise ValueError("BOW_PROJECT_ROOT must be an absolute directory")
    root = root.resolve()
    marker = root / "pyproject.toml"
    if marker.is_file():
        with marker.open("rb") as stream:
            if tomllib.load(stream).get("project", {}).get("name") == "bow-chromatography":
                return root
    raise ValueError(f"Expected a BOW pyproject.toml at {marker}. Set BOW_PROJECT_ROOT to a BOW checkout.")


PROJECT_ROOT = _project_root()


def project_path(path):
    """Absolute paths pass through; relative paths start at the project root."""
    return (PROJECT_ROOT / Path(path).expanduser()).resolve()


def _directory(variable, default):
    value = os.environ.get(variable, default)
    if not value.strip():
        raise ValueError(f"{variable} must be a directory path")
    directory = project_path(value)
    if directory.exists() and not directory.is_dir():
        raise NotADirectoryError(directory)
    if variable == "BOW_DATA_DIR" and variable in os.environ and not directory.is_dir():
        raise FileNotFoundError(directory)
    return directory


DATA_DIR = _directory("BOW_DATA_DIR", "data")
OUTPUT_DIR = _directory("BOW_OUTPUT_DIR", "results")


def config_path(path, config_file):
    """Relative config references start beside the config file."""
    return (project_path(config_file).parent / Path(path).expanduser()).resolve()


def input_path(path, *, config_file=None):
    """Require an existing input file, relative to DATA_DIR or an explicit config."""
    target = (config_path(path, config_file) if config_file is not None
              else (DATA_DIR / Path(path).expanduser()).resolve())
    if not target.is_file():
        raise FileNotFoundError(f"Expected input file at {target}")
    return target


def output_path(path):
    """Relative exports start in OUTPUT_DIR. Create parent directories on demand."""
    target = (OUTPUT_DIR / Path(path).expanduser()).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    return target
