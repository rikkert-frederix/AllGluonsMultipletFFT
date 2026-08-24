from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from su3wigner import read_table


PROJECT_DIRECTORY = Path(__file__).resolve().parents[1]
GENERATOR = PROJECT_DIRECTORY / "generate_table.py"


class GenerateTableTests(unittest.TestCase):
    def test_cutoff_three_cross_checks_all_backends(self) -> None:
        with tempfile.TemporaryDirectory() as directory_text:
            output = Path(directory_text) / "swap.tbl"
            result = subprocess.run(
                (
                    sys.executable,
                    str(GENERATOR),
                    "--method",
                    "recursive-reduction",
                    "--cross-check",
                    "--max-prefix-gluons",
                    "3",
                    "--braid-check-depth",
                    "2",
                    "--output",
                    str(output),
                ),
                cwd=PROJECT_DIRECTORY,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("cross-check[explicit-index]", result.stdout)
            self.assertIn("cross-check[fundamental-split]", result.stdout)
            self.assertIn("cross-check[dimension-only]", result.stdout)
            self.assertIn("cross-check[direct-specht]", result.stdout)
            self.assertNotIn("direct-specht coverage=", result.stdout)
            self.assertEqual(len(read_table(output)), 89)

    def test_dimension_only_primary_cross_check(self) -> None:
        with tempfile.TemporaryDirectory() as directory_text:
            output = Path(directory_text) / "swap.tbl"
            result = subprocess.run(
                (
                    sys.executable,
                    str(GENERATOR),
                    "--method",
                    "dimension-only",
                    "--cross-check",
                    "--max-prefix-gluons",
                    "1",
                    "--braid-check-depth",
                    "-1",
                    "--output",
                    str(output),
                ),
                cwd=PROJECT_DIRECTORY,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("method=dimension-only", result.stdout)
            self.assertIn("cross-check[explicit-index]", result.stdout)
            self.assertIn("cross-check[fundamental-split]", result.stdout)
            self.assertIn("cross-check[recursive-reduction]", result.stdout)
            self.assertIn("cross-check[direct-specht]", result.stdout)
            self.assertNotIn("direct-specht coverage=", result.stdout)
            self.assertEqual(len(read_table(output)), 13)

    def test_direct_specht_primary(self) -> None:
        with tempfile.TemporaryDirectory() as directory_text:
            output = Path(directory_text) / "swap.tbl"
            result = subprocess.run(
                (sys.executable, str(GENERATOR), "--method", "direct-specht",
                 "--max-prefix-gluons", "1", "--braid-check-depth", "1",
                 "--output", str(output)),
                cwd=PROJECT_DIRECTORY, capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("method=direct-specht", result.stdout)
            self.assertEqual(len(read_table(output)), 13)

    def test_failed_validation_preserves_existing_output(self) -> None:
        sentinel = "existing table must survive\n"
        with tempfile.TemporaryDirectory() as directory_text:
            directory = Path(directory_text)
            output = directory / "swap.tbl"
            output.write_text(sentinel, encoding="ascii")

            result = subprocess.run(
                (
                    sys.executable,
                    str(GENERATOR),
                    "--max-prefix-gluons",
                    "0",
                    "--output",
                    str(output),
                    "--zero-tolerance",
                    "2",
                    "--braid-check-depth",
                    "-1",
                ),
                cwd=PROJECT_DIRECTORY,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("invalid swap block", result.stderr)
            self.assertEqual(output.read_text(encoding="ascii"), sentinel)
            self.assertEqual(set(directory.iterdir()), {output})

    def test_nonfinite_zero_tolerance_is_rejected_before_generation(self) -> None:
        with tempfile.TemporaryDirectory() as directory_text:
            directory = Path(directory_text)
            for value in ("nan", "inf", "-inf"):
                with self.subTest(value=value):
                    output = directory / f"{value}.tbl"
                    result = subprocess.run(
                        (
                            sys.executable,
                            str(GENERATOR),
                            "--max-prefix-gluons",
                            "0",
                            f"--zero-tolerance={value}",
                            "--output",
                            str(output),
                        ),
                        cwd=PROJECT_DIRECTORY,
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("finite and non-negative", result.stderr)
                    self.assertFalse(output.exists())

    def test_braid_depth_below_minus_one_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory_text:
            output = Path(directory_text) / "swap.tbl"
            result = subprocess.run(
                (
                    sys.executable,
                    str(GENERATOR),
                    "--max-prefix-gluons",
                    "0",
                    "--braid-check-depth",
                    "-2",
                    "--output",
                    str(output),
                ),
                cwd=PROJECT_DIRECTORY,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("must be -1 or non-negative", result.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
