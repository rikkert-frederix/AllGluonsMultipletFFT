import argparse
from itertools import product
import subprocess
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from run_benchmark import (
    ADJOINT_BACKENDS,
    BACKENDS,
    DEFAULT_BG_BACKENDS,
    DIRECT_COLOUR_BACKENDS,
    MULTIPLET_BACKENDS,
    TRACE_BACKENDS,
    DriverRun,
    estimate_adjoint_footprint_gib,
    format_memory,
    calculate_weighted_timings,
    helicity_path,
    is_analytic_zero_helicity,
    parse_arguments,
    parse_driver_output,
    parse_proxy_output,
    read_event_helicities,
    run_driver,
    select_helicity_samples,
    validate_arguments,
    validate_table_multiplicity,
    preflight_backend,
    write_single_helicity_event,
    BenchmarkError,
)


def _parse_args(*extra_args: str) -> argparse.Namespace:
    previous = sys.argv
    sys.argv = ["benchmark.py", *extra_args]
    try:
        return parse_arguments()
    finally:
        sys.argv = previous


def _validated_args(*extra_args: str) -> argparse.Namespace:
    arguments = _parse_args(*extra_args)
    validate_arguments(arguments)
    return arguments


def _write_event(path: Path, rows: list[tuple[int, ...]], marker: str = "") -> None:
    total_gluons = len(rows[0])
    contents = [
        "AMPLIGLUON_EVENT_V1",
        f"FINAL_GLUONS {total_gluons - 2}",
        "STRONG_COUPLING 1.0",
        "BEGIN_MOMENTA",
        *(f"{leg}.0 0.0 0.0 0.0" for leg in range(1, total_gluons + 1)),
        "END_MOMENTA",
        f"NHELICITIES {len(rows)}",
        "BEGIN_HELICITIES",
        *(" ".join(f"{value:+d}" for value in row) for row in rows),
        "END_HELICITIES",
    ]
    if marker:
        contents.append(marker)
    path.write_text("\n".join(contents) + "\n", encoding="utf-8")


class TableDepthValidationTest(unittest.TestCase):
    def test_prefix_depth_accepts_both_parities_at_the_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            table = Path(temporary_directory) / "table.tbl"
            for prefix_depth in (2, 3, 6):
                with self.subTest(prefix_depth=prefix_depth):
                    table.write_text(
                        f"MAX_PREFIX_GLUONS {prefix_depth}\n", encoding="utf-8"
                    )
                    self.assertEqual(
                        validate_table_multiplicity(table, 2 * prefix_depth),
                        prefix_depth,
                    )
                    self.assertEqual(
                        validate_table_multiplicity(table, 2 * prefix_depth + 1),
                        prefix_depth,
                    )
                    with self.assertRaisesRegex(
                        BenchmarkError,
                        rf"supports at most {2 * prefix_depth + 1} total gluons",
                    ):
                        validate_table_multiplicity(table, 2 * prefix_depth + 2)


class BackendSelectionTest(unittest.TestCase):
    def test_defaults_select_both_mhv_and_colour_modes_and_ten_points(self) -> None:
        arguments = _validated_args(
            "--min-gluons",
            "4",
            "--max-gluons",
            "4",
        )

        self.assertEqual(arguments.points, 10)
        self.assertEqual(arguments.mhv_samples, 3)
        self.assertEqual(arguments.non_mhv_samples, 1)
        self.assertEqual(arguments.backends, BACKENDS)
        expected_modes = {
            "AmpliGluonMultipletOptimizedMHV": "optimized-mhv",
            "AmpliGluonMultipletDefaultBG": "default-bg",
            "AmpliGluonTraceOptimizedMHV": "optimized-mhv",
            "AmpliGluonTraceDefaultBG": "default-bg",
            "AmpliGluonAdjointOptimizedMHV": "optimized-mhv",
            "AmpliGluonAdjointDefaultBG": "default-bg",
            "AmpliGluonTraceOptimizedMHVDirectColour": "optimized-mhv",
            "AmpliGluonTraceDefaultBGDirectColour": "default-bg",
            "AmpliGluonAdjointOptimizedMHVDirectColour": "optimized-mhv",
            "AmpliGluonAdjointDefaultBGDirectColour": "default-bg",
        }
        actual_modes = {
            backend: modes[backend]
            for modes in (MULTIPLET_BACKENDS, TRACE_BACKENDS, ADJOINT_BACKENDS)
            for backend in arguments.backends
            if backend in modes
        }
        self.assertEqual(actual_modes, expected_modes)
        self.assertEqual(
            DIRECT_COLOUR_BACKENDS,
            {
                "AmpliGluonTraceOptimizedMHVDirectColour",
                "AmpliGluonTraceDefaultBGDirectColour",
                "AmpliGluonAdjointOptimizedMHVDirectColour",
                "AmpliGluonAdjointDefaultBGDirectColour",
            },
        )
        self.assertEqual(
            DEFAULT_BG_BACKENDS,
            {
                "AmpliGluonMultipletDefaultBG",
                "AmpliGluonTraceDefaultBG",
                "AmpliGluonAdjointDefaultBG",
                "AmpliGluonTraceDefaultBGDirectColour",
                "AmpliGluonAdjointDefaultBGDirectColour",
            },
        )

    def test_driver_routes_fft_and_direct_colour_modes(self) -> None:
        cases = (
            ("AmpliGluonTraceOptimizedMHV", "optimized-mhv", "fft"),
            (
                "AmpliGluonTraceDefaultBGDirectColour",
                "default-bg",
                "direct",
            ),
            ("AmpliGluonAdjointDefaultBG", "default-bg", "fft"),
            (
                "AmpliGluonAdjointOptimizedMHVDirectColour",
                "optimized-mhv",
                "direct",
            ),
        )
        for backend, mhv_mode, colour_mode in cases:
            with self.subTest(backend=backend):
                arguments = _validated_args(
                    "--backend",
                    backend,
                    "--points",
                    "1",
                    "--min-gluons",
                    "4",
                    "--max-gluons",
                    "4",
                )
                process = subprocess.CompletedProcess(
                    [],
                    0,
                    stdout=(
                        f"BACKEND {backend}\n"
                        "DIMENSION 1\n"
                        "INITIALIZATION_SECONDS 0.0\n"
                    ),
                    stderr="BENCHMARK_MAX_RSS_KIB 1024\n",
                )
                executable = Path("/tmp/benchmark-driver")
                event = Path("gg_to_2g_000001.event")
                with patch("run_benchmark.run_command", return_value=process) as run:
                    run_driver(backend, executable, [event], arguments, True)

                measured_command = run.call_args.args[0]
                separator = measured_command.index("--")
                backend_command = measured_command[separator + 1 :]
                self.assertEqual(
                    backend_command[:3],
                    [str(executable), mhv_mode, colour_mode],
                )
                self.assertEqual(
                    backend_command[-2:], [str(event), "--initialization-only"]
                )

    def test_repetition_quantum_is_opt_in_and_trace_only(self) -> None:
        arguments = _validated_args(
            "--backend",
            "AmpliGluonTraceDefaultBG",
            "--min-gluons",
            "4",
            "--max-gluons",
            "4",
            "--repetition-quantum",
            "128",
            "--max-memory-gib",
            "20",
            "--batches",
            "2",
        )
        self.assertEqual(arguments.repetition_quantum, 128)
        process = subprocess.CompletedProcess(
            [],
            0,
            stdout=(
                "BACKEND AmpliGluonTraceDefaultBG\n"
                "DIMENSION 1\n"
                "INITIALIZATION_SECONDS 0.01\n"
                "FIRST_SAMPLE_PASS_SECONDS 0.02\n"
                "MATRIX_ELEMENT 1 1 3.0\n"
                "CALIBRATION_CELL_TOTAL_SECONDS 1 1 0.25\n"
                "EVALUATION_CELL_SECONDS 1 1 1 0.001\n"
                "EVALUATION_CELL_REPETITIONS 1 1 1 256\n"
                "EVALUATION_CELL_SECONDS 2 1 1 0.001\n"
                "EVALUATION_CELL_REPETITIONS 2 1 1 256\n"
            ),
            stderr="BENCHMARK_MAX_RSS_KIB 1024\n",
        )
        with patch("run_benchmark.run_command", return_value=process) as run:
            measured = run_driver(
                "AmpliGluonTraceDefaultBG",
                Path("/tmp/benchmark-driver"),
                [Path("gg_to_2g_000001.event")],
                arguments,
                False,
            )
        self.assertEqual(set(measured.cell_repetitions.values()), {256})
        self.assertIn("--repetition-quantum=128", run.call_args.args[0])
        undercalibrated = subprocess.CompletedProcess(
            [],
            0,
            stdout=process.stdout.replace(
                "CALIBRATION_CELL_TOTAL_SECONDS 1 1 0.25",
                "CALIBRATION_CELL_TOTAL_SECONDS 1 1 0.24",
            ),
            stderr=process.stderr,
        )
        with patch("run_benchmark.run_command", return_value=undercalibrated):
            with self.assertRaisesRegex(BenchmarkError, "repetition quantum"):
                run_driver(
                    "AmpliGluonTraceDefaultBG",
                    Path("/tmp/benchmark-driver"),
                    [Path("gg_to_2g_000001.event")],
                    arguments,
                    False,
                )

        invalid = _parse_args(
            "--backend",
            "AmpliGluonAdjointDefaultBG",
            "--min-gluons",
            "4",
            "--max-gluons",
            "4",
            "--repetition-quantum",
            "128",
        )
        with self.assertRaisesRegex(BenchmarkError, "trace backends"):
            validate_arguments(invalid)

    def test_trace_default_bg_can_request_one_complete_helicity_sum(self) -> None:
        arguments = _validated_args(
            "--backend",
            "AmpliGluonTraceDefaultBG",
            "--min-gluons",
            "4",
            "--max-gluons",
            "4",
            "--repetition-quantum",
            "1",
            "--max-memory-gib",
            "8",
            "--batches",
            "2",
        )
        process = subprocess.CompletedProcess(
            [],
            0,
            stdout=(
                "BACKEND AmpliGluonTraceDefaultBG\n"
                "DIMENSION 6\n"
                "INITIALIZATION_SECONDS 0.01\n"
                "FIRST_SAMPLE_PASS_SECONDS 0.02\n"
                "MATRIX_ELEMENT 1 1 3.0\n"
                "MATRIX_ELEMENT 1 2 0.0\n"
                "CALIBRATION_CELL_TOTAL_SECONDS 1 1 0.25\n"
                "EVALUATION_CELL_SECONDS 1 1 1 0.001\n"
                "EVALUATION_CELL_REPETITIONS 1 1 1 4\n"
                "EVALUATION_CELL_SECONDS 2 1 1 0.0011\n"
                "EVALUATION_CELL_REPETITIONS 2 1 1 4\n"
            ),
            stderr="BENCHMARK_MAX_RSS_KIB 1024\n",
        )
        with patch("run_benchmark.run_command", return_value=process) as run:
            measured = run_driver(
                "AmpliGluonTraceDefaultBG",
                Path("/tmp/benchmark-driver"),
                [Path("exhaustive.event")],
                arguments,
                False,
                sum_helicities=True,
            )

        measured_command = run.call_args.args[0]
        separator = measured_command.index("--")
        backend_command = measured_command[separator + 1 :]
        self.assertIn("--sum-helicities", backend_command)
        self.assertLess(
            backend_command.index("--sum-helicities"),
            backend_command.index("--repetition-quantum=1"),
        )
        self.assertEqual(set(measured.cell_timings), {(1, 1, 1), (2, 1, 1)})

    def test_helicity_sample_counts_must_be_positive(self) -> None:
        for option in ("--mhv-samples", "--non-mhv-samples"):
            with self.subTest(option=option):
                arguments = _parse_args(
                    "--backend",
                    "AmpliGluonTraceOptimizedMHV",
                    "--min-gluons",
                    "4",
                    "--max-gluons",
                    "4",
                    option,
                    "0",
                )
                with self.assertRaisesRegex(
                    BenchmarkError, "helicity sample counts must be positive"
                ):
                    validate_arguments(arguments)

    def test_non_multiplet_selection_ignores_missing_table(self) -> None:
        selected = (*TRACE_BACKENDS, *ADJOINT_BACKENDS)
        backend_arguments = [item for name in selected for item in ("--backend", name)]
        arguments = _validated_args(
            *backend_arguments,
            "--table",
            "/path/does/not/exist.tbl",
            "--points",
            "1",
            "--min-gluons",
            "4",
            "--max-gluons",
            "4",
        )
        self.assertEqual(arguments.backends, selected)
        self.assertFalse(any(name in MULTIPLET_BACKENDS for name in arguments.backends))
        self.assertIsNone(arguments.table)
        self.assertIsNone(arguments.table_depth)

    def test_adjoint_memory_preflight_rejects_small_cap(self) -> None:
        for backend in ADJOINT_BACKENDS:
            with self.subTest(backend=backend):
                arguments = _validated_args(
                    "--backend",
                    backend,
                    "--max-memory-gib",
                    "0.001",
                    "--points",
                    "1",
                    "--min-gluons",
                    "4",
                    "--max-gluons",
                    "4",
                )
                event = Path("gg_to_7g_000001.event")
                feasible, reason = preflight_backend(
                    backend,
                    Path("/tmp/does-not-matter"),
                    [event],
                    arguments,
                )
                self.assertFalse(feasible)
                self.assertIn("memory estimate", reason)

    def test_adjoint_memory_estimate_no_longer_assumes_packed_gram(self) -> None:
        event = Path("gg_to_8g_000001.event")
        estimate = estimate_adjoint_footprint_gib(event)
        self.assertGreater(estimate, 0.1)
        self.assertLess(estimate, 0.5)

    def test_single_shot_skips_runtime_preflight(self) -> None:
        for backend in (*TRACE_BACKENDS, *ADJOINT_BACKENDS):
            with self.subTest(backend=backend):
                arguments = _validated_args(
                    "--backend",
                    backend,
                    "--skip-initialization-preflight",
                    "--points",
                    "1",
                    "--min-gluons",
                    "4",
                    "--max-gluons",
                    "4",
                )
                feasible, reason = preflight_backend(
                    backend,
                    Path("/tmp/does-not-exist"),
                    [Path("gg_to_2g_000001.event")],
                    arguments,
                )
                self.assertTrue(feasible)
                self.assertIsNone(reason)


class FormattingTest(unittest.TestCase):
    def test_peak_memory_uses_binary_units(self) -> None:
        self.assertEqual(format_memory(None), "N/A")
        self.assertEqual(format_memory(1536), "1.5 MiB")
        self.assertEqual(format_memory(2 * 1024 * 1024), "2 GiB")


class ExhaustiveHelicityTest(unittest.TestCase):
    def test_reader_accepts_all_canonical_four_gluon_configurations(self) -> None:
        rows = list(product((-1, 1), repeat=4))
        contents = [
            "AMPLIGLUON_EVENT_V1",
            "FINAL_GLUONS 2",
            "NHELICITIES 16",
            "BEGIN_HELICITIES",
            *(" ".join(str(value) for value in row) for row in rows),
            "END_HELICITIES",
        ]
        with tempfile.TemporaryDirectory() as temporary_directory:
            event = Path(temporary_directory) / "gg_to_2g_000001.event"
            event.write_text("\n".join(contents) + "\n", encoding="utf-8")

            parsed = read_event_helicities(event)

        self.assertEqual(parsed, rows)
        self.assertEqual(len(set(parsed)), 1 << 4)
        self.assertEqual(parsed[1], (-1, -1, -1, 1))
        self.assertEqual(parsed[8], (1, -1, -1, -1))

    def test_reader_rejects_a_duplicated_helicity_row(self) -> None:
        rows = list(product((-1, 1), repeat=4))
        rows[7] = rows[6]
        contents = [
            "AMPLIGLUON_EVENT_V1",
            "FINAL_GLUONS 2",
            "NHELICITIES 16",
            "BEGIN_HELICITIES",
            *(" ".join(str(value) for value in row) for row in rows),
            "END_HELICITIES",
        ]
        with tempfile.TemporaryDirectory() as temporary_directory:
            event = Path(temporary_directory) / "gg_to_2g_000001.event"
            event.write_text("\n".join(contents) + "\n", encoding="utf-8")

            with self.assertRaisesRegex(BenchmarkError, "missing or duplicated"):
                read_event_helicities(event)

    def test_analytic_zero_classifier_accounts_for_incoming_crossing(self) -> None:
        nonzero = (
            (-1, -1, -1, -1),
            (1, 1, 1, 1),
            (-1, 1, 1, -1),
        )
        zero = (
            (1, 1, -1, -1),
            (-1, -1, 1, 1),
            (-1, 1, -1, -1),
        )

        for helicities in nonzero:
            with self.subTest(helicities=helicities):
                self.assertFalse(is_analytic_zero_helicity(helicities))
        for helicities in zero:
            with self.subTest(helicities=helicities):
                self.assertTrue(is_analytic_zero_helicity(helicities))

    def test_path_classifier_distinguishes_zero_mhv_and_general_bg(self) -> None:
        cases = {
            (1, 1, -1, -1, -1, -1): "zero",
            (-1, -1, -1, -1, -1, -1): "mhv",
            (1, 1, 1, 1, 1, 1): "mhv",
            (-1, -1, 1, -1, -1, -1): "bg",
            (1, 1, 1, 1, 1, -1): "bg",
        }

        for helicities, expected in cases.items():
            with self.subTest(helicities=helicities):
                self.assertEqual(helicity_path(helicities), expected)


class ProxyOutputParsingTest(unittest.TestCase):
    helicities = {
        (1, 1): (1, 1, -1, -1, -1, -1),
        (1, 2): (-1, -1, -1, -1, -1, -1),
        (1, 3): (-1, -1, 1, -1, -1, -1),
        (2, 1): (-1, -1, 1, 1, 1, 1),
        (2, 2): (1, 1, 1, 1, 1, 1),
        (2, 3): (1, 1, 1, 1, 1, -1),
    }
    output = """\
LC_PROXY CanonicalTraceOrderBG
LC_WEIGHT 1 1 0.0
LC_WEIGHT 1 2 2.0
LC_WEIGHT 1 3 6.0
LC_HELICITY_SUM 1 8.0
LC_MHV_FRACTION 1 0.25
LC_WEIGHT 2 1 1.0e-30
LC_WEIGHT 2 2 3.0
LC_WEIGHT 2 3 1.0
LC_HELICITY_SUM 2 4.0
LC_MHV_FRACTION 2 0.75
"""

    def test_parser_reconstructs_point_and_global_proxy_mixtures(self) -> None:
        weights, totals, fractions, global_fraction = parse_proxy_output(
            self.output, self.helicities, 1.0e-24
        )

        self.assertEqual(weights[(1, 3)], 6.0)
        self.assertEqual(weights[(2, 1)], 1.0e-30)
        self.assertEqual(totals, {1: 8.0, 2: 4.0})
        self.assertEqual(fractions, {1: 0.25, 2: 0.75})
        self.assertAlmostEqual(global_fraction, 5.0 / 12.0)

    def test_parser_rejects_proxy_weight_in_an_analytic_zero_sector(self) -> None:
        output = self.output.replace("LC_WEIGHT 1 1 0.0", "LC_WEIGHT 1 1 1.0")

        with self.assertRaisesRegex(
            BenchmarkError, "nonzero leading-colour proxy in analytic-zero sector"
        ):
            parse_proxy_output(output, self.helicities, 1.0e-24)


class HelicitySamplingTest(unittest.TestCase):
    zero = (1, 1, -1, -1, -1, -1)
    mhv = (-1, -1, -1, -1, -1, -1)
    bg = (-1, -1, 1, -1, -1, -1)
    anti_zero = (-1, -1, 1, 1, 1, 1)
    anti_mhv = (1, 1, 1, 1, 1, 1)
    other_bg = (1, 1, 1, 1, 1, -1)

    def test_single_helicity_event_replaces_only_the_helicity_block(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            source = directory / "source.event"
            destination = directory / "sampled" / "sample.event"
            _write_event(source, [self.zero, self.mhv, self.bg], "SOURCE_MARKER")

            write_single_helicity_event(source, destination, self.other_bg)

            lines = destination.read_text(encoding="utf-8").splitlines()
            begin = lines.index("BEGIN_HELICITIES")
            self.assertIn("NHELICITIES 1", lines)
            self.assertNotIn("NHELICITIES 3", lines)
            self.assertEqual(
                lines[begin + 1], " ".join(f"{value:+d}" for value in self.other_bg)
            )
            self.assertEqual(lines[begin + 2], "END_HELICITIES")
            self.assertEqual(lines[-1], "SOURCE_MARKER")

    def test_selector_is_deterministic_bounded_and_writes_one_row_per_sample(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            events = [directory / "point1.event", directory / "point2.event"]
            _write_event(
                events[0], [self.zero, self.mhv, self.bg], "SOURCE_EVENT_1"
            )
            _write_event(
                events[1],
                [self.anti_zero, self.anti_mhv, self.other_bg],
                "SOURCE_EVENT_2",
            )
            helicities = {
                (1, 1): self.zero,
                (1, 2): self.mhv,
                (1, 3): self.bg,
                (2, 1): self.anti_zero,
                (2, 2): self.anti_mhv,
                (2, 3): self.other_bg,
            }
            weights = {
                (1, 1): 0.0,
                (1, 2): 1.0,
                (1, 3): 1.0,
                (2, 1): 0.0,
                (2, 2): 9.0,
                (2, 3): 9.0,
            }
            arguments = _validated_args(
                "--backend",
                "AmpliGluonTraceOptimizedMHV",
                "--min-gluons",
                "6",
                "--max-gluons",
                "6",
                "--points",
                "2",
                "--seed",
                "5",
                "--build-dir",
                str(directory / "build"),
                "--mhv-samples",
                "2",
                "--non-mhv-samples",
                "2",
            )

            first = select_helicity_samples(
                events, helicities, weights, 6, arguments
            )
            second = select_helicity_samples(
                events, helicities, weights, 6, arguments
            )

            sampled_events, sampled_helicities, samples = first
            self.assertEqual(first, second)
            self.assertEqual(len(sampled_events), 5)
            sample_paths = [sample.path for sample in samples]
            # The first row is a deterministic draw from the complete nonzero
            # proxy distribution for DefaultBG; the second is the zero sentinel.
            self.assertEqual(
                (samples[0].source_event, samples[0].source_configuration),
                (2, 2),
            )
            self.assertEqual(sample_paths[:2], ["mhv", "zero"])
            self.assertEqual(sample_paths.count("bg"), 2)
            self.assertEqual(sample_paths.count("mhv"), 2)
            self.assertEqual(sample_paths.count("zero"), 1)
            self.assertEqual(
                set(sampled_helicities), {(index, 1) for index in range(1, 6)}
            )
            for index, sample in enumerate(samples, start=1):
                self.assertEqual(sampled_helicities[(index, 1)], sample.helicities)
                self.assertEqual(
                    sample.proxy_weight,
                    weights[(sample.source_event, sample.source_configuration)],
                )
                lines = sample.event_file.read_text(encoding="utf-8").splitlines()
                begin = lines.index("BEGIN_HELICITIES")
                self.assertIn("NHELICITIES 1", lines)
                self.assertEqual(
                    lines[begin + 1],
                    " ".join(f"{value:+d}" for value in sample.helicities),
                )
                self.assertEqual(lines[begin + 2], "END_HELICITIES")
                self.assertIn(f"SOURCE_EVENT_{sample.source_event}", lines)


class SampledTimingMixtureTest(unittest.TestCase):
    @staticmethod
    def _run(
        matrix_elements: dict[tuple[int, int], float],
        cell_timings: dict[tuple[int, int, int], float],
    ) -> DriverRun:
        return DriverRun(
            initialization=0.0,
            first_helicity_sweep=0.0,
            cell_timings=cell_timings,
            matrix_elements=matrix_elements,
            dimension=1,
        )

    def test_optimized_mixture_and_default_bg_representative(self) -> None:
        helicities = {
            (1, 1): (-1, -1, 1, -1, -1, -1),
            (2, 1): (-1, -1, -1, -1, -1, -1),
            (3, 1): (1, 1, 1, 1, 1, 1),
            (4, 1): (1, 1, -1, -1, -1, -1),
        }
        matrix_elements = {
            (1, 1): 10.0,
            (2, 1): 20.0,
            (3, 1): 30.0,
            (4, 1): 0.0,
        }
        optimized = self._run(
            matrix_elements,
            {
                (1, 1, 1): 10.0,
                (1, 2, 1): 2.0,
                (1, 3, 1): 4.0,
                (2, 1, 1): 20.0,
                (2, 2, 1): 6.0,
                (2, 3, 1): 8.0,
            },
        )
        default_bg = self._run(
            matrix_elements,
            {
                (1, 1, 1): 12.0,
                (2, 1, 1): 14.0,
            },
        )
        runs = {
            "AmpliGluonTraceOptimizedMHV": optimized,
            "AmpliGluonTraceDefaultBG": default_bg,
        }

        calculate_weighted_timings(runs, helicities, tuple(runs), 0.25)

        self.assertEqual(optimized.mhv_batches, [3.0, 7.0])
        self.assertEqual(optimized.bg_batches, [10.0, 20.0])
        self.assertEqual(optimized.weighted_batches, [8.25, 16.75])
        self.assertIsNone(default_bg.mhv_batches)
        self.assertEqual(default_bg.bg_batches, [12.0, 14.0])
        self.assertEqual(default_bg.weighted_batches, [12.0, 14.0])


class DriverOutputParsingTest(unittest.TestCase):
    output = """\
BACKEND TestBackend
DIMENSION 17
INITIALIZATION_SECONDS 1.25
FIRST_HELICITY_SWEEP_SECONDS 2.5
MATRIX_ELEMENT 1 1 3.0
MATRIX_ELEMENT 1 2 4.0
EVALUATION_CELL_SECONDS 1 1 1 0.1
EVALUATION_CELL_SECONDS 1 1 2 0.2
EVALUATION_CELL_SECONDS 2 1 1 0.3
EVALUATION_CELL_SECONDS 2 1 2 0.4
"""
    sparse_output = """\
BACKEND TestBackend
DIMENSION 17
INITIALIZATION_SECONDS 1.25
FIRST_SAMPLE_PASS_SECONDS 2.5
MATRIX_ELEMENT 1 1 3.0
MATRIX_ELEMENT 2 1 4.0
MATRIX_ELEMENT 3 1 0.0
CALIBRATION_CELL_TOTAL_SECONDS 1 1 0.25
EVALUATION_CELL_SECONDS 1 1 1 0.1
EVALUATION_CELL_REPETITIONS 1 1 1 128
EVALUATION_CELL_SECONDS 2 1 1 0.2
EVALUATION_CELL_REPETITIONS 2 1 1 128
"""

    def test_parser_covers_every_batch_event_helicity_cell(self) -> None:
        run = parse_driver_output(self.output, "TestBackend", False)

        self.assertEqual(run.dimension, 17)
        self.assertEqual(run.initialization, 1.25)
        self.assertEqual(run.first_helicity_sweep, 2.5)
        self.assertEqual(run.matrix_elements, {(1, 1): 3.0, (1, 2): 4.0})
        self.assertEqual(
            run.cell_timings,
            {
                (1, 1, 1): 0.1,
                (1, 1, 2): 0.2,
                (2, 1, 1): 0.3,
                (2, 1, 2): 0.4,
            },
        )

    def test_parser_accepts_sparse_default_bg_timing_cells(self) -> None:
        run = parse_driver_output(self.sparse_output, "TestBackend", False)

        self.assertEqual(run.first_helicity_sweep, 2.5)
        self.assertEqual(
            run.matrix_elements, {(1, 1): 3.0, (2, 1): 4.0, (3, 1): 0.0}
        )
        self.assertEqual(
            run.cell_timings, {(1, 1, 1): 0.1, (2, 1, 1): 0.2}
        )
        self.assertEqual(
            run.cell_repetitions, {(1, 1, 1): 128, (2, 1, 1): 128}
        )
        self.assertEqual(run.cell_calibration_seconds, {(1, 1): 0.25})

    def test_parser_rejects_duplicate_matrix_elements_and_timings(self) -> None:
        cases = {
            "matrix element": self.output + "MATRIX_ELEMENT 1 1 3.0\n",
            "cell timing": self.output + "EVALUATION_CELL_SECONDS 1 1 1 0.1\n",
        }
        for label, output in cases.items():
            with self.subTest(label=label):
                with self.assertRaisesRegex(BenchmarkError, f"duplicate {label}"):
                    parse_driver_output(output, "TestBackend", False)


if __name__ == "__main__":
    unittest.main()
