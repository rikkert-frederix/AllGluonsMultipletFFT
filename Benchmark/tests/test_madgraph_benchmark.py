from itertools import product
import hashlib
import json
from pathlib import Path
import re
import shlex
import sys
import tempfile
from types import SimpleNamespace
import unittest


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
if str(BENCHMARK_DIR) not in sys.path:
    sys.path.insert(0, str(BENCHMARK_DIR))

from run_madgraph_benchmark import (  # noqa: E402
    CompilerIdentity,
    MadGraphBenchmarkError,
    MadGraphResult,
    build_signature,
    bundled_reference,
    canonical_helicity_index,
    direct_matrix_body,
    event_path,
    parse_driver_output,
    read_single_helicity,
    render_report,
)
from vendor_madgraph_processes import (  # noqa: E402
    process_command,
    sanitize_process_card,
)


VENDORED_ROOT = BENCHMARK_DIR / "MadGraph5"


def _write_helicity_event(path: Path, rows: list[tuple[int, ...]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            (
                "AMPLIGLUON_EVENT_V1",
                f"FINAL_GLUONS {len(rows[0]) - 2}",
                f"NHELICITIES {len(rows)}",
                "BEGIN_HELICITIES",
                *(" ".join(str(value) for value in row) for row in rows),
                "END_HELICITIES",
                "",
            )
        ),
        encoding="utf-8",
    )


def _valid_driver_output() -> str:
    return """\
BACKEND MadGraph5_aMCatNLOFixedHelicity
HELICITY_EVALUATOR MATRIX_DIRECT_VECTOR
TOTAL_GLUONS 6
INITIALIZATION_SECONDS 1.0e-6
FIRST_SAMPLE_PASS_SECONDS 2.0e-6
MATRIX_ELEMENT 1 1 12.5
EVALUATIONS_PER_SWEEP 1
EVALUATION_CELL_SECONDS 2 1 1 4.0e-6
EVALUATION_CELL_SECONDS 1 1 1 3.0e-6
CHECKSUM 25.0
"""


class HelicitySelectionTest(unittest.TestCase):
    def test_canonical_index_matches_generated_all_gluon_order(self) -> None:
        for expected, helicities in enumerate(product((-1, 1), repeat=6), start=1):
            with self.subTest(expected=expected):
                self.assertEqual(canonical_helicity_index(helicities), expected)

    def test_canonical_index_rejects_nonphysical_helicity(self) -> None:
        with self.assertRaisesRegex(
            MadGraphBenchmarkError, r"helicities must be -1 or \+1"
        ):
            canonical_helicity_index((-1, 0, 1, -1))

    def test_single_row_reader_accepts_exactly_one_physical_row(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            event = Path(temporary_directory) / "N4.event"
            row = (1, -1, -1, 1)
            _write_helicity_event(event, [row])

            self.assertEqual(read_single_helicity(event, 4), row)

    def test_single_row_reader_rejects_multiple_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            event = Path(temporary_directory) / "N4.event"
            _write_helicity_event(
                event,
                [(-1, -1, -1, -1), (1, 1, 1, 1)],
            )

            with self.assertRaisesRegex(
                MadGraphBenchmarkError, "exactly one helicity row"
            ):
                read_single_helicity(event, 4)

    def test_single_row_reader_rejects_declared_multiplicity_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            event = Path(temporary_directory) / "N6.event"
            _write_helicity_event(event, [(-1, -1, -1, -1, -1, -1)])
            event.write_text(
                event.read_text(encoding="utf-8").replace(
                    "FINAL_GLUONS 4", "FINAL_GLUONS 3"
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(MadGraphBenchmarkError, "event multiplicity"):
                read_single_helicity(event, 6)

    def test_event_locator_prefers_vendored_layout_then_native_layout(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            vendored = root / "N6.event"
            native = root / "N6" / "sampled" / "gg_to_4g_sample_000001.event"
            native.parent.mkdir(parents=True)
            native.touch()

            self.assertEqual(event_path(root, 6), native)
            vendored.touch()
            self.assertEqual(event_path(root, 6), vendored)


class DriverOutputTest(unittest.TestCase):
    def test_parser_orders_batches_and_extracts_fixed_helicity_result(self) -> None:
        initialization, first_pass, batches, matrix_element = parse_driver_output(
            _valid_driver_output(), 2
        )

        self.assertEqual(initialization, 1.0e-6)
        self.assertEqual(first_pass, 2.0e-6)
        self.assertEqual(batches, (3.0e-6, 4.0e-6))
        self.assertEqual(matrix_element, 12.5)

    def test_parser_rejects_duplicate_timing_cell(self) -> None:
        output = _valid_driver_output() + "EVALUATION_CELL_SECONDS 1 1 1 9e-6\n"

        with self.assertRaisesRegex(
            MadGraphBenchmarkError, "unexpected MadGraph timing cell"
        ):
            parse_driver_output(output, 2)

    def test_parser_reports_missing_required_scalar_as_benchmark_error(self) -> None:
        output = _valid_driver_output().replace(
            "FIRST_SAMPLE_PASS_SECONDS 2.0e-6\n", ""
        )

        with self.assertRaisesRegex(
            MadGraphBenchmarkError, "incomplete MadGraph timing output"
        ):
            parse_driver_output(output, 2)

    def test_parser_rejects_duplicate_scalar(self) -> None:
        output = _valid_driver_output() + "INITIALIZATION_SECONDS 8.0e-6\n"

        with self.assertRaisesRegex(
            MadGraphBenchmarkError, "duplicate MadGraph scalar"
        ):
            parse_driver_output(output, 2)

    def test_parser_rejects_wrong_total_gluon_count(self) -> None:
        with self.assertRaisesRegex(MadGraphBenchmarkError, "total-gluon count"):
            parse_driver_output(_valid_driver_output(), 2, expected_total_gluons=5)

    def test_parser_rejects_more_than_one_evaluation_per_sweep(self) -> None:
        output = _valid_driver_output().replace(
            "EVALUATIONS_PER_SWEEP 1", "EVALUATIONS_PER_SWEEP 2"
        )

        with self.assertRaisesRegex(MadGraphBenchmarkError, "one event per sweep"):
            parse_driver_output(output, 2, expected_total_gluons=6)

    def test_parser_rejects_nonfinite_checksum(self) -> None:
        output = _valid_driver_output().replace("CHECKSUM 25.0", "CHECKSUM nan")

        with self.assertRaisesRegex(MadGraphBenchmarkError, "invalid MadGraph"):
            parse_driver_output(output, 2, expected_total_gluons=6)

    def test_parser_rejects_missing_checksum(self) -> None:
        output = _valid_driver_output().replace("CHECKSUM 25.0\n", "")

        with self.assertRaisesRegex(
            MadGraphBenchmarkError, "incomplete MadGraph timing output"
        ):
            parse_driver_output(output, 2, expected_total_gluons=6)

    def test_zero_results_have_zero_relative_difference(self) -> None:
        result = MadGraphResult(
            total_gluons=4,
            process="g g > g g",
            feynman_diagrams=4,
            generated_matrix_graphs=6,
            colour_flows=6,
            helicities=(1, 1, -1, -1),
            helicity_index=13,
            event_path=Path("N4.event"),
            event_hash="0" * 64,
            initialization_seconds=0.0,
            first_pass_seconds=0.0,
            batch_seconds=(1.0,),
            matrix_element=0.0,
            reference_matrix_element=0.0,
            peak_rss_kib=1,
        )

        self.assertEqual(result.relative_difference, 0.0)


class CacheSignatureTest(unittest.TestCase):
    def test_compiler_identity_participates_in_build_signature(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source = Path(temporary_directory) / "source.f90"
            source.write_text("end\n", encoding="ascii")
            first = CompilerIdentity(Path("/compiler"), "compiler 1", 10, 20)
            second = CompilerIdentity(Path("/compiler"), "compiler 2", 10, 20)

            self.assertNotEqual(
                build_signature([source], ["-O3"], first),
                build_signature([source], ["-O3"], second),
            )


class BundledReferenceTest(unittest.TestCase):
    def test_reference_is_keyed_by_total_gluons_and_exact_event_hash(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            processes = Path(temporary_directory)
            event = processes / "custom.event"
            event.write_bytes(b"event contents\n")
            event_hash = hashlib.sha256(event.read_bytes()).hexdigest()
            reference = processes / "events" / "reference.json"
            reference.parent.mkdir()
            reference.write_text(
                json.dumps({"6": {event_hash: 12.5}}), encoding="utf-8"
            )

            self.assertEqual(bundled_reference(event, 6, processes), 12.5)
            self.assertIsNone(bundled_reference(event, 5, processes))
            event.write_bytes(b"different event\n")
            self.assertIsNone(bundled_reference(event, 6, processes))

    def test_malformed_matching_reference_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            processes = Path(temporary_directory)
            event = processes / "custom.event"
            event.write_bytes(b"event contents\n")
            event_hash = hashlib.sha256(event.read_bytes()).hexdigest()
            reference = processes / "events" / "reference.json"
            reference.parent.mkdir()
            reference.write_text(
                json.dumps({"6": {event_hash: "not a number"}}), encoding="utf-8"
            )

            with self.assertRaisesRegex(
                MadGraphBenchmarkError, "invalid bundled reference value"
            ):
                bundled_reference(event, 6, processes)


class ReportProvenanceTest(unittest.TestCase):
    def test_report_names_custom_event_source_and_exact_hash(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            processes = root / "processes"
            processes.mkdir()
            (processes / "VERSION").write_text("version = 3.6.0\n", encoding="ascii")
            custom_events = root / "my-events"
            event = custom_events / "N4.event"
            event.parent.mkdir()
            event.write_bytes(b"custom event\n")
            event_hash = hashlib.sha256(event.read_bytes()).hexdigest()
            result = MadGraphResult(
                total_gluons=4,
                process="g g > g g",
                feynman_diagrams=4,
                generated_matrix_graphs=6,
                colour_flows=6,
                helicities=(1, 1, -1, -1),
                helicity_index=13,
                event_path=event.resolve(),
                event_hash=event_hash,
                initialization_seconds=0.0,
                first_pass_seconds=0.0,
                batch_seconds=(1.0e-6,),
                matrix_element=1.0,
                reference_matrix_element=1.0,
                peak_rss_kib=1024,
            )
            arguments = SimpleNamespace(
                processes=processes,
                events_dir=custom_events.resolve(),
                batches=1,
                fflags="-O3",
            )

            report = render_report([result], arguments, "compiler 1")

            self.assertIn(str(custom_events.resolve()), report)
            self.assertIn(str(event.resolve()), report)
            self.assertIn(event_hash, report)


class GeneratedKernelTest(unittest.TestCase):
    def test_direct_matrix_body_stops_at_function_end(self) -> None:
        source = """\
      REAL*8 FUNCTION MATRIX(P,NHEL,IC)
      INTEGER NHEL(4), IC(4)
      CALL VXXXXX(P,NHEL(1),IC(1))
      MATRIX = 1D0
      END
      SUBROUTINE SMATRIX(P,ANS)
      INTEGER IHEL, NCOMB
      END
"""
        with tempfile.TemporaryDirectory() as temporary_directory:
            matrix = Path(temporary_directory) / "matrix.f"
            matrix.write_text(source, encoding="utf-8")

            body = direct_matrix_body(matrix)

        self.assertIn("CALL VXXXXX", body)
        self.assertNotIn("SMATRIX", body)
        self.assertNotRegex(body, r"\b(?:IHEL|NCOMB)\b")

    def test_direct_matrix_body_rejects_missing_function(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            matrix = Path(temporary_directory) / "matrix.f"
            matrix.write_text("      SUBROUTINE SMATRIX(P,ANS)\n      END\n")

            with self.assertRaisesRegex(MadGraphBenchmarkError, "MATRIX function"):
                direct_matrix_body(matrix)

    def test_coupling_adapter_initializes_every_direct_kernel_coupling(self) -> None:
        required_couplings: set[str] = set()
        for total_gluons in range(4, 7):
            matrix = VENDORED_ROOT / "Processes" / f"N{total_gluons}" / "matrix.f"
            required_couplings.update(
                re.findall(r"\bGC_\d+\b", direct_matrix_body(matrix))
            )

        adapter = (BENCHMARK_DIR / "src" / "benchmark_madgraph_coupling.f").read_text(
            encoding="utf-8"
        )
        initialized_couplings = set(re.findall(r"(?im)^\s*(GC_\d+)\s*=", adapter))
        self.assertEqual(required_couplings, {"GC_10", "GC_12"})
        self.assertEqual(required_couplings - initialized_couplings, set())


class VendoringTest(unittest.TestCase):
    def test_generation_card_requests_one_pure_gluon_standalone(self) -> None:
        output = Path("standalone_N6")
        command = process_command(6, output)

        self.assertIn("generate g g > g g g g QED=0", command)
        self.assertIn("output standalone standalone_N6 -f", command)
        self.assertTrue(command.endswith("quit\n"))

    def test_process_card_sanitizer_removes_transient_output_and_ansi(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            card = Path(temporary_directory) / "proc_card_mg5.dat"
            card.write_text(
                "\x1b[31mdevelopment build\x1b[0m\n"
                "generate g g > g g QED=0\n"
                "output standalone /tmp/transient/standalone_N4 -f\n",
                encoding="utf-8",
            )

            sanitize_process_card(card, "gg_to_2g")
            sanitized = card.read_text(encoding="utf-8")

        self.assertNotIn("\x1b", sanitized)
        self.assertNotIn("/tmp/transient", sanitized)
        self.assertIn("output standalone standalone_gg_to_2g -f", sanitized)

    def test_checked_in_metadata_and_matrix_hashes_agree(self) -> None:
        version = (
            re.search(
                r"version\s*=\s*([^\n]+)",
                (VENDORED_ROOT / "VERSION").read_text(encoding="utf-8"),
            )
            .group(1)
            .strip()
        )
        for total_gluons in range(4, 7):
            with self.subTest(total_gluons=total_gluons):
                process_dir = VENDORED_ROOT / "Processes" / f"N{total_gluons}"
                matrix = process_dir / "matrix.f"
                metadata = json.loads(
                    (process_dir / "process.json").read_text(encoding="utf-8")
                )
                digest = hashlib.sha256(matrix.read_bytes()).hexdigest()

                self.assertEqual(metadata["total_gluons"], total_gluons)
                self.assertEqual(metadata["madgraph_version"], version)
                self.assertEqual(metadata["matrix_sha256"], digest)
                self.assertGreater(metadata["feynman_diagrams"], 0)
                self.assertGreater(metadata["generated_matrix_graphs"], 0)
                self.assertGreater(metadata["colour_flows"], 0)

    def test_generated_dimensions_and_helicity_tables_match_canonical_order(
        self,
    ) -> None:
        row_pattern = re.compile(
            r"(?im)^\s*DATA \(NHEL\(I,\s*(\d+)\),I=1,(\d+)\)\s*" r"/([^/]+)/"
        )
        for total_gluons in range(4, 7):
            with self.subTest(total_gluons=total_gluons):
                process_dir = VENDORED_ROOT / "Processes" / f"N{total_gluons}"
                nexternal = process_dir.joinpath("nexternal.inc").read_text(
                    encoding="utf-8"
                )
                self.assertRegex(
                    nexternal,
                    rf"(?i)PARAMETER\s*\(\s*NEXTERNAL\s*=\s*{total_gluons}\s*\)",
                )

                matrix = process_dir.joinpath("matrix.f").read_text(encoding="utf-8")
                rows = row_pattern.findall(matrix)
                expected = list(product((-1, 1), repeat=total_gluons))
                self.assertEqual(len(rows), len(expected))
                for expected_index, ((index, width, values), helicities) in enumerate(
                    zip(rows, expected, strict=True), start=1
                ):
                    self.assertEqual(int(index), expected_index)
                    self.assertEqual(int(width), total_gluons)
                    parsed = tuple(int(value) for value in values.split(","))
                    self.assertEqual(parsed, helicities)
                    self.assertEqual(canonical_helicity_index(parsed), expected_index)

    def test_direct_kernels_resolve_all_generated_helas_calls(self) -> None:
        definitions: set[str] = set()
        definition_pattern = re.compile(
            r"(?im)^\s*(?:SUBROUTINE|(?:[A-Z0-9*]+\s+)?FUNCTION)\s+([A-Z0-9_]+)"
        )
        for source in (VENDORED_ROOT / "Source" / "DHELAS").glob("*.f"):
            definitions.update(
                name.upper() for name in definition_pattern.findall(source.read_text())
            )

        for total_gluons in range(4, 7):
            with self.subTest(total_gluons=total_gluons):
                matrix = VENDORED_ROOT / "Processes" / f"N{total_gluons}" / "matrix.f"
                body = direct_matrix_body(matrix)
                calls = set(re.findall(r"(?im)^\s*CALL\s+([A-Z0-9_]+)", body))

                self.assertTrue(calls)
                self.assertEqual(calls - definitions, set())
                self.assertNotRegex(body, r"(?i)\b(?:IHEL|NCOMB)\b")

    def test_vendored_events_are_the_reported_sample_one_rows(self) -> None:
        expected = {
            4: (
                16,
                "c43bc4fd882517442025aaef2432618744cfd0cf005813596a8aa3f13f21c41d",
            ),
            5: (
                32,
                "45891c0b556ef87ad995bb12e5bf7adb91e2290f9fbc297dda99f0bc2530b94a",
            ),
            6: (
                1,
                "411e7a6e37e4db895338834b748c98cfb9f53905d586a282db7b74d7a4cbb024",
            ),
        }
        for total_gluons, (expected_index, expected_hash) in expected.items():
            with self.subTest(total_gluons=total_gluons):
                event = VENDORED_ROOT / "events" / f"N{total_gluons}.event"
                helicities = read_single_helicity(event, total_gluons)

                self.assertEqual(canonical_helicity_index(helicities), expected_index)
                self.assertEqual(
                    hashlib.sha256(event.read_bytes()).hexdigest(), expected_hash
                )

    def test_process_cards_do_not_retain_transient_paths_or_terminal_escapes(
        self,
    ) -> None:
        for total_gluons in range(4, 7):
            with self.subTest(total_gluons=total_gluons):
                card = (
                    VENDORED_ROOT
                    / "Processes"
                    / f"N{total_gluons}"
                    / "proc_card_mg5.dat"
                )
                text = card.read_text(encoding="utf-8")
                output_lines = [
                    line for line in text.splitlines() if line.startswith("output ")
                ]

                self.assertNotIn("\x1b", text)
                self.assertEqual(len(output_lines), 1)
                fields = shlex.split(output_lines[0])
                self.assertEqual(fields[:2], ["output", "standalone"])
                self.assertFalse(Path(fields[2]).is_absolute())


if __name__ == "__main__":
    unittest.main()
