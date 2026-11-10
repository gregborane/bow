"""Regression tests for configured launches and transactional run outputs."""
from contextlib import chdir, redirect_stdout
from io import StringIO
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from Column_Simulation.experiments import VarInputs
from Column_Simulation.src.misc import results_output as output


class SimulationEntryTest(unittest.TestCase):
    def setUp(self):
        temp = TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.directory = Path(temp.name).resolve()
        self.enterContext(chdir(self.directory))
        self.results = SimpleNamespace(t=[0, 1], breakthrough=[0, 0.125])

    def test_imports_do_not_run_or_load_configuration_into_engine(self):
        code = """
import sys
from Column_Simulation import main
name = 'Column_Simulation.experiments.VarInputs'
assert name not in sys.modules
def unexpected(*args, **kwargs):
    raise AssertionError('simulation started during import')
main.run = unexpected
from Column_Simulation.experiments import VarInputs
assert VarInputs.main is main
"""
        proc = subprocess.run([sys.executable, "-c", code], check=True,
                              capture_output=True, text=True)
        self.assertEqual(proc.stdout, "")
        self.assertEqual(list(self.directory.iterdir()), [])

    def test_custom_input_and_relative_output(self):
        custom = VarInputs.build_inputs()
        custom.Nt = 7
        custom.CalculateAllVariables()
        with patch.object(VarInputs, "OUTPUT_DIR", "nested/run"), redirect_stdout(StringIO()):
            result = VarInputs.run_simulation(custom)
        self.assertEqual(result.breakthrough.shape, (7, len(custom.species)))
        self.assertEqual(Path("nested/run/VarInputs.py").read_bytes(),
                         Path(VarInputs.__file__).read_bytes())

    def test_launcher_reads_current_configuration(self):
        configured = VarInputs.build_inputs()
        with (patch.multiple(VarInputs, inputs=configured, EXPERIMENT_ID="trial",
                             OUTPUT_DIR=self.directory, EXPORT_CONFIG=False,
                             OPEN_EXTERNAL_DISPLAY=True),
              patch.object(VarInputs.main, "run", return_value=self.results) as run,
              patch.object(VarInputs, "display_results") as display):
            self.assertIs(VarInputs.run_simulation(), self.results)
            display.assert_called_once_with(configured)
            run.assert_called_once_with(configured, experiment_id="trial",
                                        output_dir=self.directory, config_file=None,
                                        open_external_display=True)

    def test_repeat_export_and_rollback(self):
        source = Path("config.py")
        source.write_text("# first run\n")
        target = self.directory / "new" / "experiment"

        def save():
            return output.save_simulation_outputs(self.results, output_dir=target,
                                                  config_file=source)

        self.assertEqual(save(), (target / "run.csv", target / "run.png"))
        source.write_text("# second run\n")
        self.results.breakthrough = [0, 0.25]
        save()
        self.assertIn("0.25", (target / "run.csv").read_text())
        before = {p.name: p.read_bytes() for p in target.iterdir()}
        self.assertEqual(before["VarInputs.py"], source.read_bytes())
        source.write_text("# third run\n")
        self.results.breakthrough = [0, 0.5]
        replace = Path.replace

        def fail_last(path, destination):
            if path.name == "VarInputs.py":
                raise OSError("interrupted export")
            return replace(path, destination)

        with patch.object(Path, "replace", fail_last), self.assertRaisesRegex(OSError, "interrupted"):
            save()
        self.assertEqual({p.name: p.read_bytes() for p in target.iterdir()}, before)
        lock = target / ".bow-export.lock"
        lock.mkdir()
        with self.assertRaisesRegex(FileExistsError, "locked"):
            save()
        lock.rmdir()
        self.assertEqual({p.name: p.read_bytes() for p in target.iterdir()}, before)

    def test_invalid_destination_or_id(self):
        Path("file").touch()
        for directory in ("", "  ", "file", "file/nested", 123):
            with self.subTest(directory=directory), self.assertRaisesRegex(ValueError, "output_dir"):
                output.save_simulation_outputs(self.results, output_dir=directory)
        for identifier in ("", "..", "a/b", True):
            with self.subTest(identifier=identifier), self.assertRaisesRegex(ValueError, "experiment_id"):
                output.save_simulation_outputs(self.results, experiment_id=identifier)
        self.assertEqual(list(Path.cwd().iterdir()), [self.directory / "file"])

    def test_default_directory_and_external_display_option(self):
        with patch.object(output, "open_image") as viewer:
            for enabled in (False, True):
                csv, png = output.save_simulation_outputs(
                    self.results, "named-run", open_external_display=enabled)
                self.assertEqual(csv, self.directory / "run_named-run/run.csv")
                self.assertEqual(png, self.directory / "run_named-run/run.png")
                self.assertEqual(viewer.call_count, int(enabled))
            viewer.assert_called_once_with(png)
