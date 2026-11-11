"""Export BOW's time-by-species outlet concentrations without resampling."""
import csv
import os
import re
from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory

import numpy as np

DEFAULT_EXPERIMENT_ID = 1

def prepare_table(results, component_names=None):

    """Return one table with time first, plus species labels in simulation order."""
    
    for field in ("t", "breakthrough"):
        if getattr(results, field, None) is None:
            raise ValueError(f"Simulation results require {field!r}")
    
    t = np.asarray(results.t)
    outlet = np.asarray(results.breakthrough)
    
    if t.ndim != 1 or not t.size:
        raise ValueError("Expected a nonempty one-dimensional time vector")
    
    if outlet.ndim == 1:
        outlet = outlet[:, None]
    
    if outlet.ndim != 2 or outlet.shape[0] != t.size or outlet.shape[1] == 0:
        raise ValueError("Expected breakthrough shape (number of times, number of species)")
    
    for array in (t, outlet):
        if array.dtype.kind not in "fiu" or not np.isfinite(array).all():
            raise ValueError("Time and outlet concentrations must be finite real numbers")
    
    if np.any(t[1:] <= t[:-1]):
        raise ValueError("Time must be strictly increasing")
    
    names = (tuple(f"component_{i+1}" for i in range(outlet.shape[1]))
             if component_names is None else tuple(component_names))
    
    if (len(names) != outlet.shape[1] or any(not isinstance(n, str) or not n.strip() for n in names)
            or len(set(names)) != len(names)):
        raise ValueError("Provide one distinct, nonempty name for each outlet component")
    
    return np.column_stack((t, outlet)), names


def write_csv(path, table, names):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time", *(f"outlet_concentration_{name}" for name in names)])
        writer.writerows(table.tolist())

def _publish(stage, directory, filenames):
    """Serialize replacements and restore the previous set on a write failure."""
    lock = directory / ".bow-export.lock"
    try:
        lock.mkdir()
    except FileExistsError as exc:
        raise FileExistsError(f"Export already locked at {lock}") from exc
    replaced = []
    try:
        for name in filenames:
            target = directory / name
            if target.exists():
                copyfile(target, stage / (name + ".previous"))
        try:
            for name in filenames:
                (stage / name).replace(directory / name)
                replaced.append(name)
        except OSError:
            for name in reversed(replaced):
                backup = stage / (name + ".previous")
                if backup.exists():
                    backup.replace(directory / name)
                else:
                    (directory / name).unlink()
            raise
    finally:
        lock.rmdir()


def save_simulation_outputs(
    results,
    experiment_id=DEFAULT_EXPERIMENT_ID,
    output_dir=None,
    component_names=None,
    config_file=None,
):
    """Publish the breakthrough CSV andoptional source copy. Return (csv_path)."""

    if (
        isinstance(experiment_id, bool)
        or not isinstance(experiment_id, (int, str))
        or not re.fullmatch(
            r"[A-Za-z0-9][A-Za-z0-9_-]*",
            str(experiment_id)
        )
    ):
        raise ValueError(
            "experiment_id must be a nonempty alphanumeric ID "
            "(internal '-' and '_' allowed)"
        )

    table, names = prepare_table(results, component_names)
    if output_dir is None:
        output_dir = f"run_{experiment_id}"

    try:
        if not os.fspath(output_dir).strip():
            raise ValueError("empty directory path")
    
        directory = Path(output_dir).resolve()
        directory.mkdir(parents=True, exist_ok=True)
    
    except (OSError, ValueError, TypeError, RuntimeError) as exc:
        raise ValueError(f"Cannot use output_dir {output_dir!r}: {exc}") from exc
    
    csv_path = directory / "run.csv"
    
    with TemporaryDirectory(dir=directory) as staging:
        stage = Path(staging)
        write_csv(stage / "run.csv", table, names)
        files = ["run.csv"]
        
        if config_file is not None:
            copyfile(
                    Path(config_file),
                    stage / "VarInputs.py"
                    )
            files.append("VarInputs.py")
        _publish(stage, directory, files)
    return csv_path
