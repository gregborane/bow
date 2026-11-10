"""Integration tests. Run after editable installation with unittest discovery."""
import os
import pkgutil
import subprocess
import sys
import sysconfig
import tempfile
import unittest
from importlib import import_module
from pathlib import Path

import numpy as np
import Column_Simulation
from Column_Simulation import paths


class PathsTest(unittest.TestCase):
    def setUp(self):
        base = Path(tempfile.gettempdir()).resolve()
        if base.is_relative_to(paths.PROJECT_ROOT):
            base = paths.PROJECT_ROOT.parent
        temp = tempfile.TemporaryDirectory(prefix="BOW path é ", dir=base)
        self.addCleanup(temp.cleanup)
        self.temp = Path(temp.name)
        self.env = os.environ.copy()
        for key in ("PYTHONPATH", "BOW_PROJECT_ROOT", "BOW_DATA_DIR", "BOW_OUTPUT_DIR"):
            self.env.pop(key, None)
        self.env.update(BOW_OUTPUT_DIR=str(self.temp / "résults dir"), OPENBLAS_NUM_THREADS="1")

    def launch(self, args, cwd, *, success=True, **env):
        proc = subprocess.run(args, cwd=cwd, env=self.env | env, capture_output=True, text=True)
        if success:
            self.assertEqual(proc.returncode, 0, proc.stderr)
        else:
            self.assertNotEqual(proc.returncode, 0)
        return proc

    def test_run_from_four_directories(self):
        expected = None
        script = paths.EXPERIMENTS_DIR / "VarInputs.py"
        for i, cwd in enumerate((paths.PROJECT_ROOT, paths.SOURCE_ROOT, paths.KINETICS_DIR, self.temp)):
            with self.subTest(cwd=cwd):
                # Keep exports isolated even when the launch directory is inside the checkout.
                run_dir = self.temp / f"run_{i}"
                script_copy = self.temp / "VarInputs.py"
                script_copy.write_text(script.read_text(encoding="utf-8").replace(
                    "OUTPUT_DIR = None", f"OUTPUT_DIR = {str(run_dir)!r}"), encoding="utf-8")
                self.launch([sys.executable, str(script_copy), "--output", f"{i}.npz"], cwd)
                self.assertEqual((run_dir / "VarInputs.py").read_bytes(), script_copy.read_bytes())
                with np.load(Path(self.env["BOW_OUTPUT_DIR"]) / f"{i}.npz") as result:
                    current = result["breakthrough"]
                    self.assertTrue(np.isfinite(current).all())
                    if expected is not None:
                        np.testing.assert_array_equal(current, expected)
                    expected = current

    def test_package_and_console(self):
        console = Path(sysconfig.get_path("scripts")) / ("bow-sim.exe" if os.name == "nt" else "bow-sim")
        self.launch([sys.executable, str(paths.EXPERIMENTS_DIR / "VarInputs.py")], self.temp)
        self.assertEqual((self.temp / "run_1/VarInputs.py").read_bytes(),
                         (paths.EXPERIMENTS_DIR / "VarInputs.py").read_bytes())
        self.assertTrue((self.temp / "run_1/run.png").is_file())
        self.assertTrue((self.temp / "run_1/run.csv").is_file())
        self.launch([str(console)], self.temp)

    def test_all_imports(self):
        before = sys.path.copy()
        for module in pkgutil.walk_packages(Column_Simulation.__path__, "Column_Simulation."):
            import_module(module.name)
        self.assertEqual(sys.path, before)

    def test_linear_experiment(self):
        from Column_Simulation.experiments.VarInputs import build_inputs
        from Column_Simulation.main import run, selected
        inputs = build_inputs()
        inputs.IsothermModel = 1
        self.assertEqual(selected("isotherms", 1)[1].__name__, "Column_Simulation.src.isotherms.1_linear")
        res = run(inputs, output_dir=self.temp / "linear")
        self.assertEqual(res.cm.shape, (inputs.Nt, inputs.Nz, len(inputs.species)))
        self.assertTrue(np.isfinite(res.cm).all())
        self.assertAlmostEqual(float(res.breakthrough[-1, 0]), 0.0012081088866, places=6)

    def test_file_resolution(self):
        sample = self.temp / "data" / "échantillon 1.csv"
        sample.parent.mkdir()
        sample.write_text("1,2\n", encoding="utf-8")
        config = self.temp / "experiment.py"
        code = """
from Column_Simulation.paths import input_path, output_path
from pathlib import Path
import sys
sample, config, expected = map(Path, sys.argv[1:])
assert input_path(sample) == sample
assert input_path(sample.name) == sample
assert input_path('data/' + sample.name, config_file=config) == sample
if sys.platform == 'win32':
    assert input_path('data\\\\' + sample.name, config_file=config) == sample
assert output_path('deep/curve.npz') == expected
assert expected.parent.is_dir()
try:
    input_path(sample.name, config_file=config)
except FileNotFoundError:
    pass
else:
    raise AssertionError('unexpected fallback to data directory')
"""
        self.launch([sys.executable, "-c", code, str(sample), str(config),
                     str(Path(self.env["BOW_OUTPUT_DIR"]) / "deep/curve.npz")], self.temp,
                    BOW_DATA_DIR=str(sample.parent))

    def test_overrides(self):
        code = "from Column_Simulation.paths import PROJECT_ROOT, OUTPUT_DIR; print(OUTPUT_DIR)"
        (self.temp / "pyproject.toml").write_text('[project]\nname="bow-chromatography"\n')
        proc = self.launch([sys.executable, "-c", code], self.temp,
                           BOW_PROJECT_ROOT=str(self.temp), BOW_OUTPUT_DIR="nested/out")
        self.assertEqual(Path(proc.stdout.strip()), self.temp / "nested/out")
        for key, value in (("BOW_PROJECT_ROOT", "relative"), ("BOW_PROJECT_ROOT", str(self.temp / "missing")),
                           ("BOW_DATA_DIR", str(self.temp / "missing")),
                           ("BOW_OUTPUT_DIR", str(self.temp / "pyproject.toml"))):
            self.launch([sys.executable, "-c", code], self.temp, success=False, **{key: value})
