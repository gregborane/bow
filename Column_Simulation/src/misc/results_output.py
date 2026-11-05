"""Export BOW's time-by-species outlet concentrations without resampling."""
import csv
import os
import re
from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory

import numpy as np

from Column_Simulation import paths


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


def breakthrough_figure(table, names, signal_name=None, time_unit="s", signal_unit=None):
    from matplotlib.figure import Figure
    from matplotlib import colormaps

    fig = Figure(figsize=(7.2, 4.5), layout="constrained")
    ax = fig.subplots()
    colors = colormaps["tab10"].colors if len(names) <= 10 else colormaps["turbo"](np.linspace(.05, .95, len(names)))
    for index, name in enumerate(names):
        ax.plot(table[:, 0], table[:, index+1], label=name, color=colors[index],
                linewidth=1.9, marker="o" if len(table) == 1 else None)
    ax.set_xlabel(f"Time ({time_unit or 's'})", fontsize=12)
    label = signal_name or "Outlet concentration"
    ax.set_ylabel(f"{label} ({signal_unit})" if signal_unit else label, fontsize=12)
    ax.set_title("Column breakthrough", loc="left", fontsize=15, pad=12)
    ax.grid(color="#dce2e6", linewidth=.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    if len(names) > 1:
        ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    return fig


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


def save_simulation_outputs(results, experiment_id=1, output_dir=None, signal_name=None,
                            time_unit="s", signal_unit=None, *, component_names=None,
                            config_file=None):
    """Return (csv_path, png_path). output_dir overrides the exact run directory.

    Relative overrides resolve against the project root. Otherwise BOW_OUTPUT_DIR,
    when configured, replaces experiments/ as the parent of run_<id>. config_file
    optionally archives the experiment source; it is not a serialization of inputs.
    Axis unit labels describe the supplied values; no conversion is performed.
    """
    if (isinstance(experiment_id, bool) or not isinstance(experiment_id, (int, str))
            or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", str(experiment_id))):
        raise ValueError("experiment_id must be a nonempty alphanumeric ID (internal '-' and '_' allowed)")
    table, names = prepare_table(results, component_names)
    base = paths.OUTPUT_DIR if "BOW_OUTPUT_DIR" in os.environ else paths.EXPERIMENTS_DIR
    directory = (paths.project_path(output_dir) if output_dir is not None
                 else base / f"run_{experiment_id}")
    directory.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=directory) as staging:
        stage = Path(staging)
        write_csv(stage / "run.csv", table, names)
        fig = breakthrough_figure(table, names, signal_name, time_unit, signal_unit)
        fig.savefig(stage / "run.png", dpi=300, format="png")
        files = ["run.csv", "run.png"]
        if config_file is not None:
            copyfile(paths.input_path(config_file), stage / "VarInputs.py")
            files.append("VarInputs.py")
        _publish(stage, directory, files)
    return directory / "run.csv", directory / "run.png"
