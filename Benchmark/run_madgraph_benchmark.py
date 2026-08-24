#!/usr/bin/env python3
"""Time vendored MadGraph all-gluon kernels at one assigned helicity."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import statistics
import subprocess
import sys
from typing import Sequence


BENCHMARK_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = BENCHMARK_DIR.parent
SOURCE_DIR = BENCHMARK_DIR / "src"
DEFAULT_PROCESSES = BENCHMARK_DIR / "MadGraph5"


class MadGraphBenchmarkError(RuntimeError):
    """A vendored-process build, validation, or timing step failed."""


@dataclass(frozen=True)
class CompilerIdentity:
    """Compiler details that must participate in an object-cache key."""

    executable: Path
    version: str
    binary_size: int
    binary_mtime_ns: int

    @property
    def cache_key(self) -> str:
        return "\n".join(
            (
                str(self.executable),
                str(self.binary_size),
                str(self.binary_mtime_ns),
                self.version,
            )
        )

    @property
    def display_name(self) -> str:
        return self.version.splitlines()[0]


@dataclass(frozen=True)
class MadGraphResult:
    total_gluons: int
    process: str
    feynman_diagrams: int
    generated_matrix_graphs: int
    colour_flows: int
    helicities: tuple[int, ...]
    helicity_index: int
    event_path: Path
    event_hash: str
    initialization_seconds: float
    first_pass_seconds: float
    batch_seconds: tuple[float, ...]
    matrix_element: float
    reference_matrix_element: float | None
    peak_rss_kib: int

    @property
    def warm_seconds(self) -> float:
        return statistics.median(self.batch_seconds)

    @property
    def relative_difference(self) -> float | None:
        if self.reference_matrix_element is None:
            return None
        denominator = max(abs(self.matrix_element), abs(self.reference_matrix_element))
        if denominator == 0.0:
            return 0.0
        return abs(self.matrix_element - self.reference_matrix_element) / denominator


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-gluons", type=int, default=4)
    parser.add_argument("--max-gluons", type=int, default=6)
    parser.add_argument("--processes", type=Path, default=DEFAULT_PROCESSES)
    parser.add_argument("--events-dir", type=Path, default=DEFAULT_PROCESSES / "events")
    parser.add_argument(
        "--build-dir", type=Path, default=BENCHMARK_DIR / "build" / "madgraph"
    )
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    parser.add_argument(
        "--fflags",
        default="-O3",
        help="optimization flags shared by the generated kernel and timing driver",
    )
    parser.add_argument("--target-seconds", type=float, default=0.25)
    parser.add_argument("--batches", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument(
        "--reference-executable",
        type=Path,
        default=(
            BENCHMARK_DIR
            / "build"
            / "benchmark_ampligluon_adjoint"
            / "benchmark_ampligluon_adjoint"
        ),
        help="optional native DDM driver used only to validate normalization",
    )
    parser.add_argument("--skip-reference-check", action="store_true")
    parser.add_argument("--relative-tolerance", type=float, default=5.0e-10)
    parser.add_argument(
        "--output",
        type=Path,
        default=BENCHMARK_DIR / "build" / "results-madgraph-fixed.md",
    )
    return parser.parse_args()


def validate_arguments(arguments: argparse.Namespace) -> None:
    if arguments.min_gluons < 4 or arguments.max_gluons < arguments.min_gluons:
        raise MadGraphBenchmarkError("invalid total-gluon range")
    if arguments.batches < 1:
        raise MadGraphBenchmarkError("--batches must be positive")
    for name in ("target_seconds", "timeout", "relative_tolerance"):
        value = getattr(arguments, name)
        if not math.isfinite(value) or value <= 0.0:
            raise MadGraphBenchmarkError(f"--{name.replace('_', '-')} must be positive")
    for name in ("processes", "events_dir", "build_dir", "output"):
        setattr(arguments, name, getattr(arguments, name).expanduser().resolve())
    arguments.reference_executable = (
        arguments.reference_executable.expanduser().resolve()
    )
    if not arguments.processes.is_dir():
        raise MadGraphBenchmarkError(
            f"vendored MadGraph process directory not found: {arguments.processes}"
        )
    if (
        not arguments.skip_reference_check
        and not arguments.reference_executable.is_file()
    ):
        raise MadGraphBenchmarkError(
            "native reference executable is missing; run the main benchmark or "
            "pass --skip-reference-check"
        )


def run_command(
    command: Sequence[str], description: str, timeout: float | None = None
) -> subprocess.CompletedProcess[str]:
    try:
        process = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=timeout,
            env={**os.environ, "LC_ALL": "C"},
        )
    except subprocess.TimeoutExpired as error:
        raise MadGraphBenchmarkError(
            f"{description} timed out after {error.timeout:.1f}s"
        ) from error
    except OSError as error:
        raise MadGraphBenchmarkError(
            f"{description} could not be started: {error}"
        ) from error
    if process.returncode != 0:
        details = "\n".join(
            value for value in (process.stdout.strip(), process.stderr.strip()) if value
        )
        raise MadGraphBenchmarkError(
            f"{description} failed:\n{shlex.join(command)}\n{details}"
        )
    return process


def identify_compiler(compiler: str) -> CompilerIdentity:
    executable_name = shutil.which(compiler)
    if executable_name is None:
        raise MadGraphBenchmarkError(f"Fortran compiler not found: {compiler}")
    executable = Path(executable_name).resolve()
    process = run_command([str(executable), "--version"], "compiler query")
    version = "\n".join(
        value.strip() for value in (process.stdout, process.stderr) if value.strip()
    )
    if not version:
        raise MadGraphBenchmarkError(
            f"Fortran compiler reported no version: {executable}"
        )
    status = executable.stat()
    return CompilerIdentity(
        executable=executable,
        version=version,
        binary_size=status.st_size,
        binary_mtime_ns=status.st_mtime_ns,
    )


def compiler_version(compiler: str) -> str:
    """Return the compiler's display version (kept for callers of the old helper)."""

    return identify_compiler(compiler).display_name


def hash_sources(paths: Sequence[Path], command: Sequence[str]) -> str:
    digest = hashlib.sha256("\0".join(command).encode())
    for path in paths:
        digest.update(str(path.resolve()).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def build_signature(
    paths: Sequence[Path], command: Sequence[str], compiler: CompilerIdentity
) -> str:
    return hash_sources(paths, [compiler.cache_key, *command])


def build_common(
    arguments: argparse.Namespace, compiler: CompilerIdentity
) -> tuple[list[Path], list[str]]:
    output = arguments.build_dir / "common"
    output.mkdir(parents=True, exist_ok=True)
    module_dir = output / "modules"
    module_dir.mkdir(parents=True, exist_ok=True)
    common_sources = [
        SOURCE_DIR / "benchmark_events.f90",
        SOURCE_DIR / "benchmark_madgraph_coupling.f",
        SOURCE_DIR / "benchmark_madgraph_fixed_helicity.f90",
        *sorted((arguments.processes / "Source" / "DHELAS").glob("*.f")),
        arguments.processes / "Source" / "coupl.inc",
    ]
    if not all(path.is_file() for path in common_sources):
        missing = [str(path) for path in common_sources if not path.is_file()]
        raise MadGraphBenchmarkError(
            "incomplete shared MadGraph source: " + ", ".join(missing)
        )
    optimization = shlex.split(arguments.fflags)
    signature = build_signature(
        common_sources,
        [arguments.fc, *optimization, "fixed-helicity-common-v1"],
        compiler,
    )
    signature_file = output / "build.sha256"
    object_names = [
        "benchmark_events.o",
        "benchmark_madgraph_coupling.o",
        "benchmark_madgraph_fixed_helicity.o",
        *(f"dhelas_{path.stem}.o" for path in common_sources[3:-1]),
    ]
    objects = [output / name for name in object_names]
    if (
        signature_file.is_file()
        and signature_file.read_text(encoding="ascii").strip() == signature
        and all(path.is_file() for path in objects)
    ):
        return objects, optimization

    fixed_flags = [*optimization, "-std=legacy", "-ffixed-line-length-132"]
    free_flags = [*optimization, "-std=f2018", "-Wall", "-Wextra"]
    coupl_include = arguments.processes / "Source"
    run_command(
        [
            arguments.fc,
            *free_flags,
            f"-J{module_dir}",
            "-c",
            str(common_sources[0]),
            "-o",
            str(objects[0]),
        ],
        "compilation of the benchmark event reader",
    )
    run_command(
        [
            arguments.fc,
            *fixed_flags,
            f"-I{coupl_include}",
            "-c",
            str(common_sources[1]),
            "-o",
            str(objects[1]),
        ],
        "compilation of the lightweight MadGraph coupling adapter",
    )
    run_command(
        [
            arguments.fc,
            *free_flags,
            f"-I{module_dir}",
            "-c",
            str(common_sources[2]),
            "-o",
            str(objects[2]),
        ],
        "compilation of the MadGraph timing driver",
    )
    for source, object_file in zip(common_sources[3:-1], objects[3:], strict=True):
        run_command(
            [
                arguments.fc,
                *fixed_flags,
                "-c",
                str(source),
                "-o",
                str(object_file),
            ],
            f"compilation of {source.name}",
        )
    signature_file.write_text(signature + "\n", encoding="ascii")
    return objects, optimization


def direct_matrix_body(matrix_source: Path) -> str:
    text = matrix_source.read_text(encoding="utf-8")
    match = re.search(r"(?ims)^\s*REAL\*8 FUNCTION MATRIX\b(.*?)(?=^\s*END\s*$)", text)
    if match is None:
        raise MadGraphBenchmarkError(f"MATRIX function missing from {matrix_source}")
    return match.group(1)


def build_process(
    total_gluons: int,
    common_objects: Sequence[Path],
    optimization: Sequence[str],
    arguments: argparse.Namespace,
    compiler: CompilerIdentity,
) -> tuple[Path, dict[str, object]]:
    process_dir = arguments.processes / "Processes" / f"N{total_gluons}"
    matrix_source = process_dir / "matrix.f"
    metadata_path = process_dir / "process.json"
    if not matrix_source.is_file() or not metadata_path.is_file():
        raise MadGraphBenchmarkError(
            f"no vendored MadGraph process for N={total_gluons}"
        )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("total_gluons") != total_gluons:
        raise MadGraphBenchmarkError(f"invalid process metadata at N={total_gluons}")
    version_text = (arguments.processes / "VERSION").read_text(encoding="utf-8")
    version_match = re.search(r"version\s*=\s*([^\n]+)", version_text)
    required_positive_integers = (
        "feynman_diagrams",
        "generated_matrix_graphs",
        "colour_flows",
    )
    if (
        version_match is None
        or metadata.get("madgraph_version") != version_match.group(1).strip()
        or not isinstance(metadata.get("process"), str)
        or any(
            isinstance(metadata.get(name), bool)
            or not isinstance(metadata.get(name), int)
            or int(metadata[name]) < 1
            for name in required_positive_integers
        )
    ):
        raise MadGraphBenchmarkError(f"invalid process metadata at N={total_gluons}")
    if hashlib.sha256(matrix_source.read_bytes()).hexdigest() != metadata.get(
        "matrix_sha256"
    ):
        raise MadGraphBenchmarkError(f"matrix source hash mismatch at N={total_gluons}")
    body = direct_matrix_body(matrix_source)
    if re.search(r"\bIHEL\b|\bNCOMB\b", body, re.IGNORECASE):
        raise MadGraphBenchmarkError(
            f"direct MATRIX kernel unexpectedly loops over helicities at N={total_gluons}"
        )

    output = arguments.build_dir / f"N{total_gluons}"
    output.mkdir(parents=True, exist_ok=True)
    matrix_object = output / "matrix.o"
    executable = output / "benchmark_madgraph_fixed_helicity"
    coupl_include = arguments.processes / "Source"
    fixed_flags = [*optimization, "-std=legacy", "-ffixed-line-length-132"]
    signature = build_signature(
        [matrix_source, *common_objects],
        [arguments.fc, *fixed_flags, "fixed-helicity-process-v1"],
        compiler,
    )
    signature_file = output / "build.sha256"
    if not (
        executable.is_file()
        and signature_file.is_file()
        and signature_file.read_text(encoding="ascii").strip() == signature
    ):
        run_command(
            [
                arguments.fc,
                *fixed_flags,
                f"-I{coupl_include}",
                "-c",
                str(matrix_source),
                "-o",
                str(matrix_object),
            ],
            f"compilation of MadGraph N={total_gluons} matrix",
            arguments.timeout,
        )
        run_command(
            [
                arguments.fc,
                *optimization,
                str(matrix_object),
                *(str(path) for path in common_objects),
                "-o",
                str(executable),
            ],
            f"linking of MadGraph N={total_gluons} benchmark",
            arguments.timeout,
        )
        signature_file.write_text(signature + "\n", encoding="ascii")
    return executable, metadata


def read_single_helicity(event: Path, total_gluons: int) -> tuple[int, ...]:
    lines = event.read_text(encoding="utf-8").splitlines()
    final_gluon_lines = [
        line.split() for line in lines if line.startswith("FINAL_GLUONS ")
    ]
    if len(final_gluon_lines) != 1 or final_gluon_lines[0] != [
        "FINAL_GLUONS",
        str(total_gluons - 2),
    ]:
        raise MadGraphBenchmarkError(
            f"event multiplicity does not match N={total_gluons}: {event}"
        )
    try:
        number_line = next(line for line in lines if line.startswith("NHELICITIES "))
        begin = lines.index("BEGIN_HELICITIES")
        end = lines.index("END_HELICITIES")
    except (StopIteration, ValueError) as error:
        raise MadGraphBenchmarkError(f"invalid sampled event: {event}") from error
    if number_line.split() != ["NHELICITIES", "1"] or end != begin + 2:
        raise MadGraphBenchmarkError(
            f"MadGraph event must contain exactly one helicity row: {event}"
        )
    row = tuple(int(value) for value in lines[begin + 1].split())
    if len(row) != total_gluons or any(abs(value) != 1 for value in row):
        raise MadGraphBenchmarkError(f"invalid helicity row in {event}")
    return row


def event_path(events_dir: Path, total_gluons: int) -> Path:
    """Locate either a vendored row or the main harness's sampled row."""

    vendored = events_dir / f"N{total_gluons}.event"
    if vendored.is_file():
        return vendored
    return (
        events_dir
        / f"N{total_gluons}"
        / "sampled"
        / f"gg_to_{total_gluons-2}g_sample_000001.event"
    )


def canonical_helicity_index(helicities: Sequence[int]) -> int:
    index = 1
    for offset, helicity in enumerate(reversed(helicities)):
        if helicity == 1:
            index += 1 << offset
        elif helicity != -1:
            raise MadGraphBenchmarkError("helicities must be -1 or +1")
    return index


def parse_driver_output(
    output: str, batches: int, expected_total_gluons: int | None = None
) -> tuple[float, float, tuple[float, ...], float]:
    scalars: dict[str, str] = {}
    batch_values: dict[int, float] = {}
    matrix_element: float | None = None
    try:
        for line in output.splitlines():
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "MATRIX_ELEMENT" and len(fields) == 4:
                if matrix_element is not None or fields[1:3] != ["1", "1"]:
                    raise MadGraphBenchmarkError("unexpected MadGraph matrix output")
                matrix_element = float(fields[3])
            elif fields[0] == "EVALUATION_CELL_SECONDS" and len(fields) == 5:
                batch = int(fields[1])
                if batch in batch_values or fields[2:4] != ["1", "1"]:
                    raise MadGraphBenchmarkError("unexpected MadGraph timing cell")
                batch_values[batch] = float(fields[4])
            elif len(fields) == 2:
                if fields[0] in scalars:
                    raise MadGraphBenchmarkError(
                        f"duplicate MadGraph scalar: {fields[0]}"
                    )
                scalars[fields[0]] = fields[1]
    except ValueError as error:
        raise MadGraphBenchmarkError("invalid numeric MadGraph output") from error
    if scalars.get("BACKEND") != "MadGraph5_aMCatNLOFixedHelicity":
        raise MadGraphBenchmarkError("unexpected MadGraph backend identity")
    if scalars.get("HELICITY_EVALUATOR") != "MATRIX_DIRECT_VECTOR":
        raise MadGraphBenchmarkError("MadGraph did not use the direct-helicity kernel")
    if set(batch_values) != set(range(1, batches + 1)) or matrix_element is None:
        raise MadGraphBenchmarkError("incomplete MadGraph timing output")
    required_scalars = {
        "TOTAL_GLUONS",
        "INITIALIZATION_SECONDS",
        "FIRST_SAMPLE_PASS_SECONDS",
        "EVALUATIONS_PER_SWEEP",
        "CHECKSUM",
    }
    if not required_scalars.issubset(scalars):
        raise MadGraphBenchmarkError("incomplete MadGraph timing output")
    try:
        total_gluons = int(scalars["TOTAL_GLUONS"])
        evaluations_per_sweep = int(scalars["EVALUATIONS_PER_SWEEP"])
        checksum = float(scalars["CHECKSUM"])
        values = (
            float(scalars["INITIALIZATION_SECONDS"]),
            float(scalars["FIRST_SAMPLE_PASS_SECONDS"]),
            tuple(batch_values[index] for index in range(1, batches + 1)),
            matrix_element,
        )
    except ValueError as error:
        raise MadGraphBenchmarkError("invalid numeric MadGraph output") from error
    if total_gluons < 4 or (
        expected_total_gluons is not None and total_gluons != expected_total_gluons
    ):
        raise MadGraphBenchmarkError("unexpected MadGraph total-gluon count")
    if evaluations_per_sweep != 1:
        raise MadGraphBenchmarkError(
            "MadGraph fixed-helicity timing must evaluate one event per sweep"
        )
    flattened = [values[0], values[1], *values[2], values[3], checksum]
    if any(not math.isfinite(value) or value < 0.0 for value in flattened):
        raise MadGraphBenchmarkError("invalid MadGraph timing output")
    return values


def native_reference(event: Path, arguments: argparse.Namespace) -> float | None:
    process = run_command(
        [
            str(arguments.reference_executable),
            "default-bg",
            "0.001",
            "1",
            str(event),
        ],
        "native DDM normalization check",
        arguments.timeout,
    )
    values = []
    for line in process.stdout.splitlines():
        fields = line.split()
        if fields[:3] == ["MATRIX_ELEMENT", "1", "1"] and len(fields) == 4:
            values.append(float(fields[3]))
    if len(values) != 1 or not math.isfinite(values[0]):
        raise MadGraphBenchmarkError("native DDM normalization check was incomplete")
    return values[0]


def bundled_reference(event: Path, total_gluons: int, processes: Path) -> float | None:
    """Look up an optional native-DDM reference by multiplicity and event hash."""

    reference_path = processes / "events" / "reference.json"
    if not reference_path.is_file():
        return None
    try:
        manifest = json.loads(reference_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise MadGraphBenchmarkError(
            f"invalid bundled reference manifest: {reference_path}"
        ) from error
    entries = manifest.get(str(total_gluons)) if isinstance(manifest, dict) else None
    if entries is None:
        return None
    if not isinstance(entries, dict):
        raise MadGraphBenchmarkError(
            f"invalid bundled reference entries for N={total_gluons}"
        )
    event_hash = hashlib.sha256(event.read_bytes()).hexdigest()
    raw_value = entries.get(event_hash)
    if raw_value is None:
        return None
    if isinstance(raw_value, bool) or not isinstance(raw_value, (int, float)):
        raise MadGraphBenchmarkError(
            f"invalid bundled reference value for N={total_gluons}"
        )
    value = float(raw_value)
    if not math.isfinite(value) or value < 0.0:
        raise MadGraphBenchmarkError(
            f"invalid bundled reference value for N={total_gluons}"
        )
    return value


def matrix_element_reference(
    event: Path, total_gluons: int, arguments: argparse.Namespace
) -> float | None:
    if arguments.skip_reference_check:
        return bundled_reference(event, total_gluons, arguments.processes)
    return native_reference(event, arguments)


def format_duration(seconds: float) -> str:
    for scale, unit in ((1.0e-9, "ns"), (1.0e-6, "us"), (1.0e-3, "ms")):
        if seconds < 1000.0 * scale:
            return f"{seconds/scale:.3g} {unit}"
    return f"{seconds:.3g} s"


def render_report(
    results: Sequence[MadGraphResult], arguments: argparse.Namespace, compiler: str
) -> str:
    version_text = (arguments.processes / "VERSION").read_text(encoding="utf-8")
    version = re.search(r"version\s*=\s*([^\n]+)", version_text).group(1).strip()
    lines = [
        "# MadGraph5_aMC@NLO direct fixed-helicity benchmark",
        "",
        (
            "Each entry calls the generated `MATRIX(P,NHEL,IC)` kernel directly "
            "with one explicit physical-helicity vector. The helicity-summing "
            "`SMATRIX` routine and the index-scanning `SMATRIXHEL` wrapper are "
            "not called."
        ),
        "",
        (
            "Each process is timed on the one assigned-helicity event listed "
            "in the provenance table below. Timings are warm medians "
            f"of {arguments.batches} calibrated batches and exclude setup."
        ),
        "",
        "| Total gluons | Feynman diagrams | Generated amplitudes (`NGRAPHS`) | Colour structures (`NCOLOR`) | Physical helicity | Canonical index | Warm fixed-helicity time | Matrix element | Relative difference vs DDM | Peak RSS |",
        "|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|",
    ]
    for result in results:
        difference = (
            f"{result.relative_difference:.3e}"
            if result.relative_difference is not None
            else "N/A"
        )
        lines.append(
            "| "
            + " | ".join(
                (
                    str(result.total_gluons),
                    str(result.feynman_diagrams),
                    str(result.generated_matrix_graphs),
                    str(result.colour_flows),
                    " ".join(f"{value:+d}" for value in result.helicities),
                    str(result.helicity_index),
                    format_duration(result.warm_seconds),
                    f"{result.matrix_element:.12e}",
                    difference,
                    f"{result.peak_rss_kib/1024.0:.3g} MiB",
                )
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Event provenance",
            "",
            f"Event input directory: `{arguments.events_dir}`",
            "",
            "| Total gluons | Event file | SHA-256 |",
            "|---:|:---|:---|",
            *(
                f"| {result.total_gluons} | `{result.event_path}` | "
                f"`{result.event_hash}` |"
                for result in results
            ),
            "",
            f"MadGraph5_aMC@NLO version: `{version}`",
            "",
            f"Compiler: `{compiler}`",
            "",
            f"Optimization flags: `{arguments.fflags}`",
            "",
            (
                "Generated fixed-form sources additionally use "
                "`-std=legacy -ffixed-line-length-132`. The pure-gluon coupling "
                "adapter sets `g_s` from the event file; no parameter-card parsing "
                "occurs inside a timed call. `MATRIX` returns the full colour sum "
                "without initial-colour/helicity averages or a final-state symmetry "
                "factor, matching the native benchmark convention."
            ),
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    arguments = parse_arguments()
    validate_arguments(arguments)
    compiler = identify_compiler(arguments.fc)
    common_objects, optimization = build_common(arguments, compiler)
    results: list[MadGraphResult] = []
    for total_gluons in range(arguments.min_gluons, arguments.max_gluons + 1):
        print(f"Building and timing MadGraph N={total_gluons}...", file=sys.stderr)
        executable, metadata = build_process(
            total_gluons, common_objects, optimization, arguments, compiler
        )
        event = event_path(arguments.events_dir, total_gluons)
        if not event.is_file():
            raise MadGraphBenchmarkError(
                f"sampled native benchmark event not found: {event}"
            )
        helicities = read_single_helicity(event, total_gluons)
        timed_command = [
            "/usr/bin/time",
            "-f",
            "MADGRAPH_PEAK_RSS_KIB %M",
            str(executable),
            f"{arguments.target_seconds:.17g}",
            str(arguments.batches),
            str(event),
        ]
        command = [
            "/usr/bin/timeout",
            "--signal=TERM",
            "--kill-after=5s",
            f"{arguments.timeout:.17g}s",
            *timed_command,
        ]
        process = run_command(
            command,
            f"MadGraph fixed-helicity timing at N={total_gluons}",
            arguments.timeout + 15.0,
        )
        rss_matches = re.findall(r"MADGRAPH_PEAK_RSS_KIB (\d+)", process.stderr)
        if len(rss_matches) != 1:
            raise MadGraphBenchmarkError("GNU time did not report MadGraph peak RSS")
        initialization, first_pass, batches, matrix_element = parse_driver_output(
            process.stdout, arguments.batches, total_gluons
        )
        reference = matrix_element_reference(event, total_gluons, arguments)
        result = MadGraphResult(
            total_gluons=total_gluons,
            process=str(metadata["process"]),
            feynman_diagrams=int(metadata["feynman_diagrams"]),
            generated_matrix_graphs=int(metadata["generated_matrix_graphs"]),
            colour_flows=int(metadata["colour_flows"]),
            helicities=helicities,
            helicity_index=canonical_helicity_index(helicities),
            event_path=event.resolve(),
            event_hash=hashlib.sha256(event.read_bytes()).hexdigest(),
            initialization_seconds=initialization,
            first_pass_seconds=first_pass,
            batch_seconds=batches,
            matrix_element=matrix_element,
            reference_matrix_element=reference,
            peak_rss_kib=int(rss_matches[0]),
        )
        if (
            result.relative_difference is not None
            and result.relative_difference > arguments.relative_tolerance
        ):
            raise MadGraphBenchmarkError(
                f"MadGraph disagrees with native DDM at N={total_gluons}: "
                f"{result.relative_difference:.3e}"
            )
        results.append(result)
    report = render_report(results, arguments, compiler.display_name)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MadGraphBenchmarkError as error:
        print(f"MadGraph benchmark error: {error}", file=sys.stderr)
        raise SystemExit(1)
