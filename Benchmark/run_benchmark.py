#!/usr/bin/env python3
"""Estimate all-gluon single-helicity evaluation times from bounded samples."""

from __future__ import annotations

import argparse
import bisect
from dataclasses import dataclass
import hashlib
import math
import os
from pathlib import Path
import random
import shlex
import statistics
import subprocess
import sys
from typing import Iterable, Mapping, Sequence


BENCHMARK_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = BENCHMARK_DIR.parent
SOURCE_DIR = BENCHMARK_DIR / "src"
BACKENDS = (
    "AmpliGluonMultipletOptimizedMHV",
    "AmpliGluonMultipletDefaultBG",
    "AmpliGluonTraceOptimizedMHV",
    "AmpliGluonTraceDefaultBG",
    "AmpliGluonAdjointOptimizedMHV",
    "AmpliGluonAdjointDefaultBG",
    "AmpliGluonTraceOptimizedMHVDirectColour",
    "AmpliGluonTraceDefaultBGDirectColour",
    "AmpliGluonAdjointOptimizedMHVDirectColour",
    "AmpliGluonAdjointDefaultBGDirectColour",
)
MULTIPLET_BACKENDS = {
    "AmpliGluonMultipletOptimizedMHV": "optimized-mhv",
    "AmpliGluonMultipletDefaultBG": "default-bg",
}
TRACE_BACKENDS = {
    "AmpliGluonTraceOptimizedMHV": "optimized-mhv",
    "AmpliGluonTraceDefaultBG": "default-bg",
    "AmpliGluonTraceOptimizedMHVDirectColour": "optimized-mhv",
    "AmpliGluonTraceDefaultBGDirectColour": "default-bg",
}
ADJOINT_BACKENDS = {
    "AmpliGluonAdjointOptimizedMHV": "optimized-mhv",
    "AmpliGluonAdjointDefaultBG": "default-bg",
    "AmpliGluonAdjointOptimizedMHVDirectColour": "optimized-mhv",
    "AmpliGluonAdjointDefaultBGDirectColour": "default-bg",
}
DIRECT_COLOUR_BACKENDS = {
    "AmpliGluonTraceOptimizedMHVDirectColour",
    "AmpliGluonTraceDefaultBGDirectColour",
    "AmpliGluonAdjointOptimizedMHVDirectColour",
    "AmpliGluonAdjointDefaultBGDirectColour",
}
DEFAULT_BG_BACKENDS = {
    name
    for modes in (MULTIPLET_BACKENDS, TRACE_BACKENDS, ADJOINT_BACKENDS)
    for name, mode in modes.items()
    if mode == "default-bg"
}
BACKEND_LABELS = {
    "AmpliGluonMultipletOptimizedMHV": "Multiplet (optimized MHV)",
    "AmpliGluonMultipletDefaultBG": "Multiplet (default BG)",
    "AmpliGluonTraceOptimizedMHV": "Trace (optimized MHV, FFT colour)",
    "AmpliGluonTraceDefaultBG": "Trace (default BG, FFT colour)",
    "AmpliGluonAdjointOptimizedMHV": "Adjoint (optimized MHV, FFT colour)",
    "AmpliGluonAdjointDefaultBG": "Adjoint (default BG, FFT colour)",
    "AmpliGluonTraceOptimizedMHVDirectColour": (
        "Trace (optimized MHV, direct colour)"
    ),
    "AmpliGluonTraceDefaultBGDirectColour": "Trace (default BG, direct colour)",
    "AmpliGluonAdjointOptimizedMHVDirectColour": (
        "Adjoint (optimized MHV, direct colour)"
    ),
    "AmpliGluonAdjointDefaultBGDirectColour": (
        "Adjoint (default BG, direct colour)"
    ),
}


class BenchmarkError(RuntimeError):
    """A build, generation, validation, or timing step failed."""


@dataclass
class DriverRun:
    initialization: float
    first_helicity_sweep: float | None
    cell_timings: dict[tuple[int, int, int], float]
    matrix_elements: dict[tuple[int, int], float]
    dimension: int
    peak_rss_kib: int | None = None
    weighted_batches: list[float] | None = None
    mhv_batches: list[float] | None = None
    bg_batches: list[float] | None = None
    cell_repetitions: dict[tuple[int, int, int], int] | None = None
    cell_calibration_seconds: dict[tuple[int, int], float] | None = None


@dataclass(frozen=True)
class HelicitySample:
    source_event: int
    source_configuration: int
    path: str
    proxy_weight: float
    helicities: tuple[int, ...]
    event_file: Path


@dataclass
class MultiplicityResult:
    total_gluons: int
    event_hash: str
    event_files: list[Path]
    helicities: dict[tuple[int, int], tuple[int, ...]]
    samples: list[HelicitySample]
    proxy_weights: dict[tuple[int, int], float]
    proxy_point_totals: dict[int, float]
    proxy_point_mhv_fractions: dict[int, float]
    proxy_mhv_fraction: float
    initialization_samples: dict[str, list[float]]
    runs: dict[str, DriverRun]
    skipped_backends: dict[str, str]
    compared_backends: tuple[str, ...]
    maximum_relative_difference: float
    maximum_zero_value: float


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate RAMBO-on-diet events, estimate their leading-colour "
            "helicity mixture, benchmark bounded MHV and non-MHV samples, "
            "validate the all-gluon generators, and emit Markdown tables."
        )
    )
    parser.add_argument("--min-gluons", type=int, default=4)
    parser.add_argument("--max-gluons", type=int, default=8)
    parser.add_argument("--points", type=int, default=10)
    parser.add_argument("--sqrt-s", type=float, default=1000.0)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument(
        "--mhv-samples",
        type=int,
        default=3,
        help=(
            "number of leading-colour-weighted MHV/anti-MHV timing samples "
            "per multiplicity"
        ),
    )
    parser.add_argument(
        "--non-mhv-samples",
        type=int,
        default=1,
        help=(
            "number of leading-colour-weighted general-helicity timing samples "
            "per multiplicity"
        ),
    )
    parser.add_argument(
        "--table",
        type=Path,
        default=REPOSITORY_ROOT
        / "Wigner6j"
        / "data"
        / "su3_adjoint_swap_prefix_6.tbl",
        help=(
            "multiplet Wigner table; a prefix-depth p table supports total "
            "multiplicities through 2p+1"
        ),
    )
    parser.add_argument(
        "--build-dir", type=Path, default=BENCHMARK_DIR / "build"
    )
    parser.add_argument(
        "--fc", default=os.environ.get("FC", "gfortran"), help="Fortran compiler"
    )
    parser.add_argument(
        "--fflags",
        default="-O3 -std=f2018 -Wall -Wextra",
        help="identical compiler flags used for all generators and drivers",
    )
    parser.add_argument(
        "--backend",
        action="append",
        choices=BACKENDS,
        help=(
            "select one backend; repeat for multiple. "
            "default: all ten native backend variants"
        ),
    )
    parser.add_argument(
        "--backend-timeout",
        type=float,
        default=600.0,
        help=(
            "wall-time cap (seconds) per backend invocation; "
            "timed out invocations are skipped for that multiplicity"
        ),
    )
    parser.add_argument(
        "--max-memory-gib",
        type=float,
        default=None,
        help=(
            "maximum estimated runtime memory per backend; defaults to available "
            "system memory times --memory-fraction"
        ),
    )
    parser.add_argument(
        "--memory-fraction",
        type=float,
        default=0.85,
        help="default memory cap fraction of system available memory",
    )
    parser.add_argument(
        "--max-target-seconds",
        type=float,
        default=None,
        help="cap the calibration target per backend process",
    )
    parser.add_argument(
        "--target-seconds",
        type=float,
        default=0.25,
        help="minimum CPU time used to calibrate an evaluation batch",
    )
    parser.add_argument("--batches", type=int, default=3)
    parser.add_argument(
        "--repetition-quantum",
        type=int,
        default=None,
        help=(
            "calibrate trace-backend timing samples in indivisible groups of "
            "this many scalar evaluations"
        ),
    )
    parser.add_argument("--initialization-runs", type=int, default=3)
    parser.add_argument(
        "--skip-initialization-preflight",
        action="store_true",
        help=(
            "do not launch a separate initialization-only feasibility run; "
            "useful for one-shot high-multiplicity measurements"
        ),
    )
    parser.add_argument("--relative-tolerance", type=float, default=5.0e-10)
    parser.add_argument("--zero-tolerance", type=float, default=1.0e-24)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def validate_arguments(arguments: argparse.Namespace) -> None:
    if arguments.min_gluons < 4 or arguments.min_gluons > arguments.max_gluons:
        raise BenchmarkError("the total-gluon range must start at four or above")
    if arguments.points < 1 or arguments.points > 100:
        raise BenchmarkError("--points must be between 1 and 100")
    backend_names = list(BACKENDS) if arguments.backend is None else arguments.backend
    if not backend_names:
        raise BenchmarkError("select at least one backend")
    arguments.backends = tuple(dict.fromkeys(backend_names))
    if not any(name in MULTIPLET_BACKENDS for name in arguments.backends):
        arguments.table = None
        arguments.table_depth = None
    for name in ("sqrt_s", "target_seconds", "relative_tolerance", "zero_tolerance"):
        value = getattr(arguments, name)
        if not math.isfinite(value) or value <= 0.0:
            raise BenchmarkError(f"--{name.replace('_', '-')} must be positive")
    if arguments.batches < 1 or arguments.initialization_runs < 1:
        raise BenchmarkError("batch and initialization-run counts must be positive")
    if arguments.repetition_quantum is not None:
        if arguments.repetition_quantum < 1:
            raise BenchmarkError("--repetition-quantum must be positive")
        if any(name not in TRACE_BACKENDS for name in arguments.backends):
            raise BenchmarkError(
                "--repetition-quantum is currently supported only by trace backends"
            )
    if arguments.mhv_samples < 1 or arguments.non_mhv_samples < 1:
        raise BenchmarkError("helicity sample counts must be positive")
    for timeout_option in ("backend_timeout", "memory_fraction"):
        timeout_or_fraction = getattr(arguments, timeout_option)
        if timeout_or_fraction is None or not math.isfinite(timeout_or_fraction):
            raise BenchmarkError(f"--{timeout_option.replace('_', '-')} must be finite")
        if timeout_option == "backend_timeout" and timeout_or_fraction <= 0.0:
            raise BenchmarkError(f"--{timeout_option.replace('_', '-')} must be positive")
        if timeout_option == "memory_fraction" and not 0.05 <= timeout_or_fraction <= 0.99:
            raise BenchmarkError("--memory-fraction must be in [0.05, 0.99]")
    if arguments.max_target_seconds is not None:
        if not math.isfinite(arguments.max_target_seconds) or arguments.max_target_seconds <= 0.0:
            raise BenchmarkError("--max-target-seconds must be positive if provided")
        arguments.target_seconds = min(
            arguments.target_seconds, arguments.max_target_seconds
        )
    arguments.build_dir = arguments.build_dir.expanduser().resolve()
    if arguments.output is not None:
        arguments.output = arguments.output.expanduser().resolve()
    if arguments.table is not None:
        arguments.table = arguments.table.expanduser().resolve()
        if not arguments.table.is_file():
            raise BenchmarkError(f"Wigner table not found: {arguments.table}")
        arguments.table_depth = validate_table_multiplicity(
            arguments.table, arguments.max_gluons
        )
    else:
        arguments.table_depth = None

    if arguments.max_memory_gib is None:
        arguments.max_memory_gib = arguments.memory_fraction * available_memory_gib()
    if not math.isfinite(arguments.max_memory_gib) or arguments.max_memory_gib <= 0.0:
        raise BenchmarkError("--max-memory-gib must be positive")


def available_memory_gib() -> float:
    meminfo = Path("/proc/meminfo")
    if meminfo.is_file():
        for line in meminfo.read_text(encoding="utf-8").splitlines():
            if not line.startswith("MemAvailable:"):
                continue
            _, value, unit = line.split()
            if unit != "kB":
                continue
            return float(int(value)) / (1024.0 ** 2)
    return float("inf")


def read_table_depth(table: Path) -> int:
    with table.open(encoding="utf-8") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) == 2 and fields[0] == "MAX_PREFIX_GLUONS":
                return int(fields[1])
    raise BenchmarkError(f"MAX_PREFIX_GLUONS is missing from {table}")


def validate_table_multiplicity(table: Path, maximum_gluons: int) -> int:
    """Check the current multiplet-consumer cutoff and return table depth."""

    table_depth = read_table_depth(table)
    if table_depth < 0:
        raise BenchmarkError(f"invalid MAX_PREFIX_GLUONS in {table}")
    maximum_supported_gluons = 2 * table_depth + 1
    if maximum_gluons > maximum_supported_gluons:
        raise BenchmarkError(
            f"{table} supports at most {maximum_supported_gluons} total gluons "
            f"(prefix depth {table_depth})"
        )
    return table_depth


def run_command(
    command: Sequence[str],
    description: str,
    timeout_seconds: float | None = None,
) -> subprocess.CompletedProcess[str]:
    try:
        process = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            env={**os.environ, "LC_ALL": "C"},
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        raise BenchmarkError(
            f"{description} timed out after {error.timeout:.1f}s:\n"
            f"{shlex.join(command)}"
        ) from error
    if process.returncode != 0:
        details = "\n".join(
            part for part in (process.stdout.strip(), process.stderr.strip()) if part
        )
        raise BenchmarkError(
            f"{description} failed with status {process.returncode}:\n"
            f"{shlex.join(command)}\n{details}"
        )
    return process


def compiler_version(compiler: Sequence[str]) -> str:
    return run_command([*compiler, "--version"], "compiler query").stdout.splitlines()[0]


def build_signature(
    sources: Iterable[Path],
    compiler: Sequence[str],
    compiler_identity: str,
    flags: Sequence[str],
) -> str:
    digest = hashlib.sha256()
    digest.update("\0".join([*compiler, compiler_identity, *flags]).encode())
    for source in sources:
        digest.update(str(source.resolve()).encode())
        digest.update(source.read_bytes())
    return digest.hexdigest()


def build_program(
    name: str,
    sources: Sequence[Path],
    compiler: Sequence[str],
    compiler_identity: str,
    flags: Sequence[str],
    build_dir: Path,
) -> Path:
    output_dir = build_dir / name
    module_dir = output_dir / "modules"
    executable = output_dir / name
    signature_file = output_dir / "build.sha256"
    signature = build_signature(sources, compiler, compiler_identity, flags)
    if (
        executable.is_file()
        and signature_file.is_file()
        and signature_file.read_text(encoding="ascii").strip() == signature
    ):
        return executable

    print(f"Building {name}...", file=sys.stderr, flush=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    module_dir.mkdir(parents=True, exist_ok=True)
    objects: list[Path] = []
    for index, source in enumerate(sources[:-1]):
        output = output_dir / f"{index:02d}_{source.stem}.o"
        run_command(
            [
                *compiler,
                *flags,
                f"-J{module_dir}",
                f"-I{module_dir}",
                "-c",
                str(source),
                "-o",
                str(output),
            ],
            f"compilation of {source}",
        )
        objects.append(output)
    driver = sources[-1]
    run_command(
        [
            *compiler,
            *flags,
            f"-J{module_dir}",
            f"-I{module_dir}",
            str(driver),
            *(str(path) for path in objects),
            "-o",
            str(executable),
        ],
        f"linking of {name}",
    )
    signature_file.write_text(signature + "\n", encoding="ascii")
    return executable


def build_executables(
    arguments: argparse.Namespace,
    compiler: Sequence[str],
    compiler_identity: str,
    flags: Sequence[str],
) -> tuple[Path, Path, dict[str, Path]]:
    multiplet = REPOSITORY_ROOT / "AmpliGluonMultiplet" / "src"
    trace = REPOSITORY_ROOT / "AmpliGluonTrace" / "src"
    adjoint = REPOSITORY_ROOT / "AmpliGluonAdjoint" / "src"
    rambo = REPOSITORY_ROOT / "RamboOnDiet"
    common_event = SOURCE_DIR / "benchmark_events.f90"

    rambo_executable = build_program(
        "generate_ampligluon_events",
        [
            rambo / "find_zero.f90",
            rambo / "rambo_on_diet.f90",
            rambo / "generate_ampligluon_events.f90",
        ],
        compiler,
        compiler_identity,
        flags,
        arguments.build_dir,
    )
    proxy_executable = build_program(
        "benchmark_helicity_proxy",
        [
            common_event,
            trace / "ampligluon_common.f90",
            trace / "gluon_kinematics.f90",
            trace / "trace_order_recursion.f90",
            SOURCE_DIR / "benchmark_helicity_proxy.f90",
        ],
        compiler,
        compiler_identity,
        flags,
        arguments.build_dir,
    )
    executables: dict[str, Path] = {}
    selected_multiplet_backends = [
        name for name in arguments.backends if name in MULTIPLET_BACKENDS
    ]
    if selected_multiplet_backends:
        multiplet_executable = build_program(
            "benchmark_ampligluon_multiplet",
            [
                common_event,
                multiplet / "ampligluon_multiplet_kinds.f90",
                multiplet / "wigner_table.f90",
                multiplet / "multiplet_paths.f90",
                multiplet / "recoupling_plan.f90",
                multiplet / "multiplet_radiation.f90",
                multiplet / "qcd_kinematics.f90",
                multiplet / "multiplet_mhv.f90",
                multiplet / "ampligluon_multiplet.f90",
                SOURCE_DIR / "benchmark_ampligluon_multiplet.f90",
            ],
            compiler,
            compiler_identity,
            flags,
            arguments.build_dir,
        )
        for backend in selected_multiplet_backends:
            executables[backend] = multiplet_executable
    selected_trace_backends = [
        name for name in arguments.backends if name in TRACE_BACKENDS
    ]
    if selected_trace_backends:
        trace_executable = build_program(
            "benchmark_ampligluon_trace",
            [
                common_event,
                trace / "ampligluon_common.f90",
                trace / "symmetric_group_fft.f90",
                trace / "trace_colour_kernel.f90",
                trace / "gluon_kinematics.f90",
                trace / "trace_colour_matrix.f90",
                trace / "trace_current_dag.f90",
                trace / "trace_order_recursion.f90",
                trace / "trace_mhv.f90",
                trace / "ampligluon_trace.f90",
                SOURCE_DIR / "benchmark_ampligluon_trace.f90",
            ],
            compiler,
            compiler_identity,
            flags,
            arguments.build_dir,
        )
        for backend in selected_trace_backends:
            executables[backend] = trace_executable
    selected_adjoint_backends = [
        name for name in arguments.backends if name in ADJOINT_BACKENDS
    ]
    if selected_adjoint_backends:
        adjoint_executable = build_program(
            "benchmark_ampligluon_adjoint",
            [
                common_event,
                trace / "ampligluon_common.f90",
                trace / "symmetric_group_fft.f90",
                trace / "trace_colour_kernel.f90",
                trace / "gluon_kinematics.f90",
                trace / "trace_current_dag.f90",
                adjoint / "adjoint_colour_matrix.f90",
                adjoint / "adjoint_current_dag.f90",
                adjoint / "ampligluon_adjoint.f90",
                SOURCE_DIR / "benchmark_ampligluon_adjoint.f90",
            ],
            compiler,
            compiler_identity,
            flags,
            arguments.build_dir,
        )
        for backend in selected_adjoint_backends:
            executables[backend] = adjoint_executable
    return rambo_executable, proxy_executable, executables


def generate_events(
    executable: Path, total_gluons: int, arguments: argparse.Namespace
) -> tuple[list[Path], str, dict[tuple[int, int], tuple[int, ...]]]:
    event_dir = arguments.build_dir / "events" / f"N{total_gluons}"
    event_dir.mkdir(parents=True, exist_ok=True)
    prefix = event_dir / f"gg_to_{total_gluons - 2}g"
    seed = arguments.seed + total_gluons
    run_command(
        [
            str(executable),
            str(total_gluons - 2),
            str(arguments.points),
            str(prefix),
            f"{arguments.sqrt_s:.17g}",
            str(seed),
        ],
        f"RAMBO generation at {total_gluons} total gluons",
    )
    events = [
        Path(f"{prefix}_{point:06d}.event")
        for point in range(1, arguments.points + 1)
    ]
    if not all(path.is_file() for path in events):
        raise BenchmarkError("RAMBO did not create every requested event")
    digest = hashlib.sha256()
    helicities: dict[tuple[int, int], tuple[int, ...]] = {}
    for event_index, event in enumerate(events, start=1):
        digest.update(event.read_bytes())
        for configuration, row in enumerate(read_event_helicities(event), start=1):
            helicities[(event_index, configuration)] = row
    return events, digest.hexdigest(), helicities


def read_event_helicities(event: Path) -> list[tuple[int, ...]]:
    """Read and validate the exhaustive physical-helicity rows in one event."""

    lines = event.read_text(encoding="utf-8").splitlines()
    final_gluons: int | None = None
    number_of_helicities: int | None = None
    begin_index: int | None = None
    for index, line in enumerate(lines):
        fields = line.split()
        if len(fields) == 2 and fields[0] == "FINAL_GLUONS":
            final_gluons = int(fields[1])
        elif len(fields) == 2 and fields[0] == "NHELICITIES":
            number_of_helicities = int(fields[1])
        elif fields == ["BEGIN_HELICITIES"]:
            begin_index = index + 1
    if final_gluons is None or number_of_helicities is None or begin_index is None:
        raise BenchmarkError(f"incomplete helicity metadata in {event}")

    total_gluons = final_gluons + 2
    expected_configurations = 1 << total_gluons
    if number_of_helicities != expected_configurations:
        raise BenchmarkError(
            f"{event} contains {number_of_helicities} helicities; "
            f"expected all {expected_configurations} configurations"
        )
    rows: list[tuple[int, ...]] = []
    for line in lines[begin_index : begin_index + number_of_helicities]:
        row = tuple(int(value) for value in line.split())
        if len(row) != total_gluons or any(value not in (-1, 1) for value in row):
            raise BenchmarkError(f"invalid helicity row in {event}")
        rows.append(row)
    if len(rows) != number_of_helicities or len(set(rows)) != number_of_helicities:
        raise BenchmarkError(f"helicity rows in {event} are missing or duplicated")
    expected_rows = [
        tuple(
            1 if configuration & (1 << (total_gluons - leg - 1)) else -1
            for leg in range(total_gluons)
        )
        for configuration in range(number_of_helicities)
    ]
    if rows != expected_rows:
        raise BenchmarkError(f"helicity ordering in {event} is not canonical")
    end_index = begin_index + number_of_helicities
    if end_index >= len(lines) or lines[end_index].split() != ["END_HELICITIES"]:
        raise BenchmarkError(f"missing END_HELICITIES row in {event}")
    return rows


def is_analytic_zero_helicity(helicities: Sequence[int]) -> bool:
    """Return whether a physical-input helicity row vanishes after crossing."""

    outgoing_positive = outgoing_positive_helicities(helicities)
    return outgoing_positive < 2 or outgoing_positive > len(helicities) - 2


def outgoing_positive_helicities(helicities: Sequence[int]) -> int:
    """Count positive helicities after crossing both incoming gluons."""

    return sum(value < 0 for value in helicities[:2]) + sum(
        value > 0 for value in helicities[2:]
    )


def helicity_path(helicities: Sequence[int]) -> str:
    """Classify a physical helicity row as zero, MHV, or general BG."""

    outgoing_positive = outgoing_positive_helicities(helicities)
    total_gluons = len(helicities)
    if outgoing_positive < 2 or outgoing_positive > total_gluons - 2:
        return "zero"
    if outgoing_positive == 2 or outgoing_positive == total_gluons - 2:
        return "mhv"
    return "bg"


def parse_proxy_output(
    output: str,
    helicities: Mapping[tuple[int, int], Sequence[int]],
    zero_tolerance: float,
) -> tuple[
    dict[tuple[int, int], float],
    dict[int, float],
    dict[int, float],
    float,
]:
    """Parse and validate exhaustive canonical leading-colour proxy weights."""

    proxy_name: str | None = None
    weights: dict[tuple[int, int], float] = {}
    reported_totals: dict[int, float] = {}
    reported_fractions: dict[int, float] = {}
    for line in output.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "LC_PROXY" and len(fields) == 2:
            proxy_name = fields[1]
        elif fields[0] == "LC_WEIGHT" and len(fields) == 4:
            key = (int(fields[1]), int(fields[2]))
            if key in weights:
                raise BenchmarkError(f"duplicate leading-colour proxy weight: {key}")
            weights[key] = float(fields[3])
        elif fields[0] == "LC_HELICITY_SUM" and len(fields) == 3:
            event = int(fields[1])
            if event in reported_totals:
                raise BenchmarkError(f"duplicate proxy helicity sum for event {event}")
            reported_totals[event] = float(fields[2])
        elif fields[0] == "LC_MHV_FRACTION" and len(fields) == 3:
            event = int(fields[1])
            if event in reported_fractions:
                raise BenchmarkError(f"duplicate proxy MHV fraction for event {event}")
            reported_fractions[event] = float(fields[2])

    if proxy_name != "CanonicalTraceOrderBG":
        raise BenchmarkError("unexpected or missing leading-colour proxy identity")
    if set(weights) != set(helicities):
        raise BenchmarkError("leading-colour proxy does not cover every generated helicity")
    if any(not math.isfinite(value) or value < 0.0 for value in weights.values()):
        raise BenchmarkError("leading-colour proxy returned an invalid weight")

    events = sorted({event for event, _ in helicities})
    totals: dict[int, float] = {}
    fractions: dict[int, float] = {}
    global_total = 0.0
    global_mhv = 0.0
    for event in events:
        total = 0.0
        mhv_total = 0.0
        for key, row in helicities.items():
            if key[0] != event:
                continue
            value = weights[key]
            path = helicity_path(row)
            if path == "zero":
                if value > zero_tolerance:
                    raise BenchmarkError(
                        f"nonzero leading-colour proxy in analytic-zero sector at {key}"
                    )
                continue
            total += value
            if path == "mhv":
                mhv_total += value
        if not math.isfinite(total) or total <= 0.0:
            raise BenchmarkError(
                f"non-positive leading-colour proxy helicity sum at event {event}"
            )
        fraction = mhv_total / total
        totals[event] = total
        fractions[event] = fraction
        global_total += total
        global_mhv += mhv_total
        if event not in reported_totals or not math.isclose(
            reported_totals[event], total, rel_tol=5.0e-13, abs_tol=zero_tolerance
        ):
            raise BenchmarkError(f"inconsistent proxy helicity sum at event {event}")
        if event not in reported_fractions or not math.isclose(
            reported_fractions[event], fraction, rel_tol=5.0e-13, abs_tol=5.0e-15
        ):
            raise BenchmarkError(f"inconsistent proxy MHV fraction at event {event}")
    if set(reported_totals) != set(events) or set(reported_fractions) != set(events):
        raise BenchmarkError("leading-colour proxy event metadata is incomplete")
    if not math.isfinite(global_total) or global_total <= 0.0:
        raise BenchmarkError("non-positive global leading-colour proxy sum")
    return weights, totals, fractions, global_mhv / global_total


def run_helicity_proxy(
    executable: Path,
    events: Sequence[Path],
    helicities: Mapping[tuple[int, int], Sequence[int]],
    arguments: argparse.Namespace,
) -> tuple[
    dict[tuple[int, int], float],
    dict[int, float],
    dict[int, float],
    float,
]:
    process = run_command(
        [str(executable), *(str(path) for path in events)],
        "canonical leading-colour helicity proxy",
        timeout_seconds=arguments.backend_timeout,
    )
    return parse_proxy_output(process.stdout, helicities, arguments.zero_tolerance)


def _weighted_draws(
    keys: Sequence[tuple[int, int]],
    weights: Mapping[tuple[int, int], float],
    count: int,
    generator: random.Random,
) -> list[tuple[int, int]]:
    positive = [(key, weights[key]) for key in sorted(keys) if weights[key] > 0.0]
    if not positive:
        raise BenchmarkError("cannot sample a helicity path with zero proxy weight")
    cumulative: list[float] = []
    running = 0.0
    for _, value in positive:
        running += value
        cumulative.append(running)
    if not math.isfinite(running) or running <= 0.0:
        raise BenchmarkError("invalid cumulative leading-colour proxy weight")
    draws: list[tuple[int, int]] = []
    for _ in range(count):
        target = generator.random() * running
        index = min(bisect.bisect_right(cumulative, target), len(positive) - 1)
        draws.append(positive[index][0])
    return draws


def write_single_helicity_event(
    source: Path,
    destination: Path,
    helicities: Sequence[int],
) -> None:
    """Copy one generated phase-space point with exactly one helicity row."""

    lines = source.read_text(encoding="utf-8").splitlines()
    number_index: int | None = None
    begin_index: int | None = None
    end_index: int | None = None
    for index, line in enumerate(lines):
        fields = line.split()
        if fields and fields[0] == "NHELICITIES":
            number_index = index
        elif fields == ["BEGIN_HELICITIES"]:
            begin_index = index
        elif fields == ["END_HELICITIES"]:
            end_index = index
    if (
        number_index is None
        or begin_index is None
        or end_index is None
        or not number_index < begin_index < end_index
    ):
        raise BenchmarkError(f"cannot derive sampled event from {source}")
    replacement = [
        "NHELICITIES 1",
        "BEGIN_HELICITIES",
        " ".join(f"{value:+d}" for value in helicities),
        "END_HELICITIES",
    ]
    sampled_lines = [*lines[:number_index], *replacement, *lines[end_index + 1 :]]
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(sampled_lines) + "\n", encoding="utf-8")


def select_helicity_samples(
    events: Sequence[Path],
    helicities: Mapping[tuple[int, int], tuple[int, ...]],
    weights: Mapping[tuple[int, int], float],
    total_gluons: int,
    arguments: argparse.Namespace,
) -> tuple[list[Path], dict[tuple[int, int], tuple[int, ...]], list[HelicitySample]]:
    """Draw a bounded joint phase-space/helicity sample from the LC proxy."""

    by_path = {
        path: [key for key, row in helicities.items() if helicity_path(row) == path]
        for path in ("zero", "mhv", "bg")
    }
    generator = random.Random(arguments.seed + 104729 * total_gluons)
    # Draw the timing strata separately so the relatively rare MHV path cannot
    # disappear from a small sample.  Also draw the default-BG representative
    # from the complete nonzero population.  Replacing one same-stratum draw
    # with that representative preserves the requested stratum counts without
    # adding another expensive optimized-backend evaluation.  Draw it first so
    # changing either stratum sample count does not change the DefaultBG row.
    nonzero_population = [*by_path["mhv"], *by_path["bg"]]
    representative = _weighted_draws(
        nonzero_population, weights, 1, generator
    )[0]
    representative_path = helicity_path(helicities[representative])
    bg_sample_count = arguments.non_mhv_samples if by_path["bg"] else 0
    bg_draws = (
        _weighted_draws(by_path["bg"], weights, bg_sample_count, generator)
        if by_path["bg"]
        else []
    )
    mhv_draws = _weighted_draws(
        by_path["mhv"], weights, arguments.mhv_samples, generator
    )
    if representative_path == "bg":
        bg_draws[0] = representative
    elif representative_path == "mhv":
        mhv_draws[0] = representative
    else:
        raise BenchmarkError("default-BG representative is analytically zero")

    remaining_path_samples: list[tuple[tuple[int, int], str]] = []
    representative_removed = False
    for path, draws in (("bg", bg_draws), ("mhv", mhv_draws)):
        for key in draws:
            if (
                not representative_removed
                and path == representative_path
                and key == representative
            ):
                representative_removed = True
                continue
            remaining_path_samples.append((key, path))
    if not representative_removed:
        raise BenchmarkError("internal default-BG representative mismatch")

    # Keep one representative nonzero row first and the zero sentinel second.
    # Default-BG invocations receive just this prefix, while optimized
    # invocations receive the complete path-stratified sample pool.
    selected = [(representative, representative_path)]
    if by_path["zero"]:
        selected.append((generator.choice(sorted(by_path["zero"])), "zero"))
    selected.extend(remaining_path_samples)

    expected_samples = bg_sample_count + arguments.mhv_samples + int(bool(by_path["zero"]))
    if len(selected) != expected_samples:
        raise BenchmarkError("internal helicity sample-count mismatch")

    sample_dir = arguments.build_dir / "events" / f"N{total_gluons}" / "sampled"
    sampled_events: list[Path] = []
    sampled_helicities: dict[tuple[int, int], tuple[int, ...]] = {}
    samples: list[HelicitySample] = []
    for sample_index, (source_key, path) in enumerate(selected, start=1):
        source_event, source_configuration = source_key
        row = helicities[source_key]
        destination = sample_dir / (
            f"gg_to_{total_gluons - 2}g_sample_{sample_index:06d}.event"
        )
        write_single_helicity_event(events[source_event - 1], destination, row)
        sampled_events.append(destination)
        sampled_helicities[(sample_index, 1)] = row
        samples.append(
            HelicitySample(
                source_event=source_event,
                source_configuration=source_configuration,
                path=path,
                proxy_weight=weights[source_key],
                helicities=row,
                event_file=destination,
            )
        )
    return sampled_events, sampled_helicities, samples


def parse_driver_output(output: str, backend: str, initialization_only: bool) -> DriverRun:
    scalars: dict[str, str] = {}
    cell_timings: dict[tuple[int, int, int], float] = {}
    cell_repetitions: dict[tuple[int, int, int], int] = {}
    cell_calibration_seconds: dict[tuple[int, int], float] = {}
    matrix_elements: dict[tuple[int, int], float] = {}
    for line in output.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "MATRIX_ELEMENT" and len(fields) == 4:
            key = (int(fields[1]), int(fields[2]))
            if key in matrix_elements:
                raise BenchmarkError(f"duplicate matrix element from {backend}: {key}")
            matrix_elements[key] = float(fields[3])
        elif fields[0] == "EVALUATION_CELL_SECONDS" and len(fields) == 5:
            key = (int(fields[1]), int(fields[2]), int(fields[3]))
            if key in cell_timings:
                raise BenchmarkError(f"duplicate cell timing from {backend}: {key}")
            cell_timings[key] = float(fields[4])
        elif fields[0] == "EVALUATION_CELL_REPETITIONS" and len(fields) == 5:
            key = (int(fields[1]), int(fields[2]), int(fields[3]))
            if key in cell_repetitions:
                raise BenchmarkError(
                    f"duplicate cell repetition count from {backend}: {key}"
                )
            cell_repetitions[key] = int(fields[4])
        elif fields[0] == "CALIBRATION_CELL_TOTAL_SECONDS" and len(fields) == 4:
            key = (int(fields[1]), int(fields[2]))
            if key in cell_calibration_seconds:
                raise BenchmarkError(
                    f"duplicate cell calibration timing from {backend}: {key}"
                )
            cell_calibration_seconds[key] = float(fields[3])
        elif len(fields) == 2:
            scalars[fields[0]] = fields[1]
    if scalars.get("BACKEND") != backend:
        raise BenchmarkError(f"unexpected {backend} output:\n{output}")
    for key in ("DIMENSION", "INITIALIZATION_SECONDS"):
        if key not in scalars:
            raise BenchmarkError(f"missing {key} in {backend} output")
    first_sample_key = (
        "FIRST_SAMPLE_PASS_SECONDS"
        if "FIRST_SAMPLE_PASS_SECONDS" in scalars
        else "FIRST_HELICITY_SWEEP_SECONDS"
    )
    if not initialization_only and (
        first_sample_key not in scalars or not cell_timings or not matrix_elements
    ):
        raise BenchmarkError(f"incomplete full output from {backend}")
    run = DriverRun(
        initialization=float(scalars["INITIALIZATION_SECONDS"]),
        first_helicity_sweep=(
            float(scalars[first_sample_key])
            if first_sample_key in scalars
            else None
        ),
        cell_timings=cell_timings,
        matrix_elements=matrix_elements,
        dimension=int(scalars["DIMENSION"]),
        cell_repetitions=cell_repetitions,
        cell_calibration_seconds=cell_calibration_seconds,
    )
    numbers = [run.initialization, *run.cell_timings.values()]
    if run.first_helicity_sweep is not None:
        numbers.append(run.first_helicity_sweep)
    numbers.extend(run.matrix_elements.values())
    numbers.extend(cell_calibration_seconds.values())
    if any(not math.isfinite(value) for value in numbers):
        raise BenchmarkError(f"non-finite result from {backend}")
    if any(value < 0.0 for value in [run.initialization, *run.cell_timings.values()]):
        raise BenchmarkError(f"negative timing from {backend}")
    if run.first_helicity_sweep is not None and run.first_helicity_sweep < 0.0:
        raise BenchmarkError(f"negative first-sweep timing from {backend}")
    if not initialization_only:
        batch_ids = {batch for batch, _, _ in run.cell_timings}
        timing_sets = {
            batch: {
                (event, configuration)
                for timing_batch, event, configuration in run.cell_timings
                if timing_batch == batch
            }
            for batch in batch_ids
        }
        first_timing_set = next(iter(timing_sets.values()))
        if not first_timing_set or any(
            cells != first_timing_set for cells in timing_sets.values()
        ):
            raise BenchmarkError(f"inconsistent event/helicity timing cells from {backend}")
        if not first_timing_set <= set(run.matrix_elements):
            raise BenchmarkError(f"timing cell without a matrix element from {backend}")
        if cell_repetitions and (
            set(cell_repetitions) != set(cell_timings)
            or any(value < 1 for value in cell_repetitions.values())
        ):
            raise BenchmarkError(
                f"invalid cell repetition counts from {backend}"
            )
        if cell_calibration_seconds and (
            set(cell_calibration_seconds) != first_timing_set
            or any(value <= 0.0 for value in cell_calibration_seconds.values())
        ):
            raise BenchmarkError(f"invalid cell calibration timings from {backend}")
    return run


def run_driver(
    backend: str,
    executable: Path,
    events: Sequence[Path],
    arguments: argparse.Namespace,
    initialization_only: bool,
    *,
    sum_helicities: bool = False,
) -> DriverRun:
    """Run one backend executable and return its parsed metrics."""
    command = [str(executable)]
    if backend in MULTIPLET_BACKENDS:
        command.extend([str(arguments.table), MULTIPLET_BACKENDS[backend]])
    elif backend in TRACE_BACKENDS:
        command.extend(
            [
                TRACE_BACKENDS[backend],
                "direct" if backend in DIRECT_COLOUR_BACKENDS else "fft",
            ]
        )
    elif backend in ADJOINT_BACKENDS:
        command.extend(
            [
                ADJOINT_BACKENDS[backend],
                "direct" if backend in DIRECT_COLOUR_BACKENDS else "fft",
            ]
        )
    command.extend(
        [
            f"{arguments.target_seconds:.17g}",
            str(arguments.batches),
            *(str(path) for path in events),
        ]
    )
    if sum_helicities:
        if backend != "AmpliGluonTraceDefaultBG":
            raise BenchmarkError(
                "the exhaustive helicity-sum timing mode is implemented only "
                "for AmpliGluonTraceDefaultBG"
            )
        command.append("--sum-helicities")
    if arguments.repetition_quantum is not None:
        command.append(f"--repetition-quantum={arguments.repetition_quantum}")
    if initialization_only:
        command.append("--initialization-only")

    # GNU time measures the complete fresh-process run, including
    # initialization and evaluation.  The address-space limit is enforced on
    # the backend itself; GNU timeout owns the process group so a timed-out
    # backend cannot remain behind as an orphaned high-memory process.
    memory_bytes = max(1, int(arguments.max_memory_gib * 1024.0**3))
    measured_command = [
        "/usr/bin/timeout",
        "--signal=TERM",
        "--kill-after=5s",
        f"{arguments.backend_timeout:.17g}s",
        "/usr/bin/time",
        "--format=BENCHMARK_MAX_RSS_KIB %M",
        "/usr/bin/prlimit",
        f"--as={memory_bytes}",
        "--",
        *command,
    ]
    process = run_command(
        measured_command,
        f"{backend} benchmark",
        timeout_seconds=arguments.backend_timeout + 15.0,
    )
    run = parse_driver_output(process.stdout, backend, initialization_only)
    if not initialization_only:
        batch_ids = {batch for batch, _, _ in run.cell_timings}
        expected_batch_ids = set(range(1, arguments.batches + 1))
        if batch_ids != expected_batch_ids:
            raise BenchmarkError(
                f"{backend} returned timing batches {sorted(batch_ids)}; "
                f"expected {sorted(expected_batch_ids)}"
            )
        if arguments.repetition_quantum is not None and (
            run.cell_repetitions is None
            or set(run.cell_repetitions) != set(run.cell_timings)
            or any(
                count < arguments.repetition_quantum
                or count % arguments.repetition_quantum != 0
                for count in run.cell_repetitions.values()
            )
            or run.cell_calibration_seconds is None
            or set(run.cell_calibration_seconds)
            != {(event, configuration) for _, event, configuration in run.cell_timings}
            or any(
                seconds
                < arguments.target_seconds / len(run.cell_calibration_seconds)
                for seconds in run.cell_calibration_seconds.values()
            )
        ):
            raise BenchmarkError(
                f"{backend} did not use the requested repetition quantum"
            )
    for line in process.stderr.splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[0] == "BENCHMARK_MAX_RSS_KIB":
            run.peak_rss_kib = int(fields[1])
            break
    if run.peak_rss_kib is None:
        raise BenchmarkError(f"missing peak-RSS measurement from {backend}")
    return run


def parse_total_gluons_from_event(path: Path) -> int:
    """Return total gluons from event filename ``gg_to_<n>g_....event``."""

    stem = path.name
    if not stem.startswith("gg_to_"):
        raise BenchmarkError(f"cannot parse final multiplicity from event path: {path}")
    suffix = stem[len("gg_to_") :]
    final_tag = suffix.split("_", 1)[0]
    if not final_tag.endswith("g"):
        raise BenchmarkError(f"cannot parse final multiplicity from event filename: {path}")
    final_gluons = int(final_tag[:-1])
    if final_gluons < 1:
        raise BenchmarkError(f"invalid final multiplicity in event filename: {path}")
    return final_gluons + 2


def estimate_adjoint_footprint_gib(event: Path) -> float:
    """Conservative persistent and initialization memory estimate for Adjoint."""

    total_gluons = parse_total_gluons_from_event(event)
    adjoint_orders = math.factorial(total_gluons - 2)
    trace_orders = math.factorial(total_gluons - 1)
    # The colour kernel, amplitudes, and either Fourier or direct work buffers
    # use O((N-2)!) storage.  Kernel initialization additionally keeps the trace
    # kernel, two integer finite-difference vectors, and the lexicographic
    # S_(N-1) permutation table.
    # The additional 512 bytes per trace order conservatively cover the
    # selected-current DAG and its construction workspaces at supported
    # multiplicities.  This deliberately overestimates the measured N=10
    # peak while preserving the factorial preflight guard at larger N.
    bytes_estimate = (
        64 * adjoint_orders
        + (512 + 24 + 4 * (total_gluons - 1)) * trace_orders
    )
    return float(bytes_estimate) / (1024.0**3)


def preflight_backend(
    backend: str,
    executable: Path,
    events: Sequence[Path],
    arguments: argparse.Namespace,
) -> tuple[bool, str | None]:
    """Return ``(is_feasible, skip_reason)`` for a backend at this multiplicity."""

    if backend in ADJOINT_BACKENDS:
        estimated_gib = estimate_adjoint_footprint_gib(events[0])
        if estimated_gib > arguments.max_memory_gib:
            return False, (
                f"Adjoint preflight memory estimate {estimated_gib:.2f} GiB "
                f"exceeds cap {arguments.max_memory_gib:.2f} GiB"
            )
    if arguments.skip_initialization_preflight:
        return True, None
    try:
        run_driver(backend, executable, events[:1], arguments, True)
    except BenchmarkError as error:
        return False, str(error)
    return True, None


def compare_matrix_elements(
    runs: dict[str, DriverRun],
    helicities: Mapping[tuple[int, int], Sequence[int]],
    relative_tolerance: float,
    zero_tolerance: float,
    backends: Sequence[str],
) -> tuple[float, float]:
    """Validate matrix-element agreement and return (max relative spread, max zero abs)."""

    if not backends:
        return float("nan"), float("nan")
    key_sets = {backend: set(runs[backend].matrix_elements) for backend in backends}
    for backend, keys in key_sets.items():
        if not keys or not keys <= set(helicities):
            raise BenchmarkError(
                f"{backend} returned an invalid sampled event/configuration set"
            )
    keys = set().union(*key_sets.values())

    maximum_relative = 0.0
    maximum_zero_value = 0.0
    for key in keys:
        if key not in helicities:
            raise BenchmarkError(f"missing helicities for event/configuration {key}")
        available_backends = [name for name in backends if key in key_sets[name]]
        values = [runs[name].matrix_elements[key] for name in available_backends]
        if is_analytic_zero_helicity(helicities[key]):
            largest_absolute_value = max(abs(value) for value in values)
            maximum_zero_value = max(maximum_zero_value, largest_absolute_value)
            if largest_absolute_value > zero_tolerance:
                rendered = ", ".join(
                    f"{name}={runs[name].matrix_elements[key]:.16e}"
                    for name in available_backends
                )
                raise BenchmarkError(
                    f"nonzero result in an analytic-zero sector at {key}: "
                    f"{rendered}"
                )
            continue

        # Optimized/default variants of one colour basis should agree tightly.
        # Distinct bases can have larger convention-level numerical spread at
        # high multiplicity, so retain the existing all-backend comparison.
        if len(available_backends) < 2:
            continue
        reference = values[0]
        scale = max(abs(value) for value in values)
        for value in values[1:]:
            difference = abs(value - reference) / scale if scale > 0.0 else 0.0
            maximum_relative = max(maximum_relative, difference)
            if difference > relative_tolerance:
                rendered = ", ".join(
                    f"{name}={runs[name].matrix_elements[key]:.16e}"
                    for name in available_backends
                )
                raise BenchmarkError(
                    f"matrix-element mismatch at event/configuration {key}: {rendered}"
                )
    return (
        maximum_relative if len(backends) >= 2 else float("nan"),
        maximum_zero_value,
    )


def calculate_weighted_timings(
    runs: Mapping[str, DriverRun],
    helicities: Mapping[tuple[int, int], Sequence[int]],
    backends: Sequence[str],
    proxy_mhv_fraction: float,
) -> None:
    """Aggregate sampled timings into MHV, BG, and LC-mixture estimates."""

    if not 0.0 <= proxy_mhv_fraction <= 1.0:
        raise BenchmarkError("leading-colour proxy MHV fraction is outside [0, 1]")
    nonzero_keys = {
        key for key, row in helicities.items() if helicity_path(row) != "zero"
    }
    mhv_keys = {key for key in nonzero_keys if helicity_path(helicities[key]) == "mhv"}
    bg_keys = nonzero_keys - mhv_keys
    for backend in backends:
        run = runs[backend]
        for key, value in run.matrix_elements.items():
            if helicity_path(helicities[key]) != "zero" and value < 0.0:
                raise BenchmarkError(
                    f"negative sampled nonzero matrix element from {backend} "
                    f"at {key}: {value:.16e}"
                )
        batch_ids = sorted({batch for batch, _, _ in run.cell_timings})
        if not batch_ids:
            raise BenchmarkError(f"{backend} returned no timing batches")
        timed_keys = {
            (event, configuration)
            for timing_batch, event, configuration in run.cell_timings
            if timing_batch == batch_ids[0]
        }
        is_default_bg = backend in DEFAULT_BG_BACKENDS
        if is_default_bg:
            if len(timed_keys) != 1 or not timed_keys <= nonzero_keys:
                raise BenchmarkError(
                    f"{backend} must time exactly one representative nonzero cell"
                )
        elif timed_keys != nonzero_keys:
            raise BenchmarkError(
                f"{backend} did not time every selected nonzero helicity"
            )

        mhv_batches: list[float] = []
        bg_batches: list[float] = []
        combined_batches: list[float] = []
        for batch in batch_ids:
            batch_keys = {
                (event, configuration)
                for timing_batch, event, configuration in run.cell_timings
                if timing_batch == batch
            }
            if batch_keys != timed_keys:
                raise BenchmarkError(
                    f"{backend} batch {batch} has a different event/helicity set"
                )
            if is_default_bg:
                bg_value = run.cell_timings[(batch, *next(iter(timed_keys)))]
                bg_batches.append(bg_value)
                combined_batches.append(bg_value)
                continue

            mhv_values = [
                run.cell_timings[(batch, *key)] for key in timed_keys & mhv_keys
            ]
            bg_values = [
                run.cell_timings[(batch, *key)] for key in timed_keys & bg_keys
            ]
            if not mhv_values:
                raise BenchmarkError(f"{backend} has no timed MHV sample")
            mhv_value = statistics.fmean(mhv_values)
            mhv_batches.append(mhv_value)
            if bg_values:
                bg_value = statistics.fmean(bg_values)
                bg_batches.append(bg_value)
                combined_batches.append(
                    proxy_mhv_fraction * mhv_value
                    + (1.0 - proxy_mhv_fraction) * bg_value
                )
            else:
                if not math.isclose(proxy_mhv_fraction, 1.0, abs_tol=5.0e-14):
                    raise BenchmarkError(
                        f"{backend} has no BG sample but the proxy has BG weight"
                    )
                combined_batches.append(mhv_value)
        run.mhv_batches = mhv_batches or None
        run.bg_batches = bg_batches or None
        run.weighted_batches = combined_batches


def rotate_backend_order(backends: Sequence[str], total_gluons: int) -> list[str]:
    if not backends:
        return []
    offset = total_gluons % len(backends)
    return [*backends[offset:], *backends[:offset]]


def collect_results(
    arguments: argparse.Namespace,
    rambo_executable: Path,
    proxy_executable: Path,
    executables: dict[str, Path],
    compiler_name: str,
    table_hash: str,
) -> list[MultiplicityResult]:
    results: list[MultiplicityResult] = []
    selected_backends = list(arguments.backends)

    for total_gluons in range(arguments.min_gluons, arguments.max_gluons + 1):
        print(
            f"Generating and benchmarking {total_gluons} total gluons...",
            file=sys.stderr,
            flush=True,
        )
        events, event_hash, exhaustive_helicities = generate_events(
            rambo_executable, total_gluons, arguments
        )
        (
            proxy_weights,
            proxy_point_totals,
            proxy_point_mhv_fractions,
            proxy_mhv_fraction,
        ) = run_helicity_proxy(
            proxy_executable, events, exhaustive_helicities, arguments
        )
        sampled_events, helicities, samples = select_helicity_samples(
            events,
            exhaustive_helicities,
            proxy_weights,
            total_gluons,
            arguments,
        )
        order = rotate_backend_order(selected_backends, total_gluons)

        runs: dict[str, DriverRun] = {}
        initialization_samples: dict[str, list[float]] = {}
        skipped_backends: dict[str, str] = {}

        for backend in order:
            backend_events = (
                sampled_events[:2]
                if backend in DEFAULT_BG_BACKENDS
                else sampled_events
            )
            feasible, reason = preflight_backend(
                backend, executables[backend], backend_events, arguments
            )
            if not feasible:
                skipped_backends[backend] = reason
                continue
            try:
                run = run_driver(
                    backend, executables[backend], backend_events, arguments, False
                )
            except BenchmarkError as error:
                skipped_backends[backend] = str(error)
                continue
            runs[backend] = run
            initialization_samples[backend] = [run.initialization]

        for sample in range(1, arguments.initialization_runs):
            if not runs:
                break
            sample_order = order if sample % 2 == 0 else list(reversed(order))
            for backend in sample_order:
                if backend not in runs:
                    continue
                try:
                    initialization_samples[backend].append(
                        run_driver(
                            backend,
                            executables[backend],
                            sampled_events[:1],
                            arguments,
                            True,
                        ).initialization
                    )
                except BenchmarkError as error:
                    print(
                        f"{backend} initialization-sample preflight failed for {total_gluons} gluons: {error}",
                        file=sys.stderr,
                    )

        compared_backends = tuple(
            backend for backend in selected_backends if backend in runs
        )
        maximum_relative, maximum_zero_value = compare_matrix_elements(
            runs,
            helicities,
            arguments.relative_tolerance,
            arguments.zero_tolerance,
            compared_backends,
        )
        calculate_weighted_timings(
            runs, helicities, compared_backends, proxy_mhv_fraction
        )

        result = MultiplicityResult(
            total_gluons=total_gluons,
            event_hash=event_hash,
            event_files=events,
            helicities=helicities,
            samples=samples,
            proxy_weights=proxy_weights,
            proxy_point_totals=proxy_point_totals,
            proxy_point_mhv_fractions=proxy_point_mhv_fractions,
            proxy_mhv_fraction=proxy_mhv_fraction,
            initialization_samples=initialization_samples,
            runs=runs,
            skipped_backends=skipped_backends,
            compared_backends=compared_backends,
            maximum_relative_difference=maximum_relative,
            maximum_zero_value=maximum_zero_value,
        )
        results.append(result)

        if arguments.output is not None:
            arguments.output.parent.mkdir(parents=True, exist_ok=True)
            report = render_report(
                arguments=arguments,
                results=results,
                compiler_name=compiler_name,
                table_hash=table_hash,
                completed=total_gluons,
            )
            arguments.output.write_text(report, encoding="utf-8")

    return results


def format_duration(seconds: float) -> str:
    if not math.isfinite(seconds):
        return "N/A"
    if seconds >= 1.0:
        return f"{seconds:.3g} s"
    if seconds >= 1.0e-3:
        return f"{seconds * 1.0e3:.3g} ms"
    if seconds >= 1.0e-6:
        return f"{seconds * 1.0e6:.3g} µs"
    return f"{seconds * 1.0e9:.3g} ns"


def format_memory(kibibytes: int | None) -> str:
    if kibibytes is None:
        return "N/A"
    mebibytes = kibibytes / 1024.0
    if mebibytes >= 1024.0:
        return f"{mebibytes / 1024.0:.3g} GiB"
    return f"{mebibytes:.3g} MiB"


def fastest(timings: dict[str, float]) -> str:
    valid = [(name, value) for name, value in timings.items() if math.isfinite(value)]
    if len(valid) < 2:
        return "N/A"
    ordered = sorted(valid, key=lambda item: item[1])
    name, best = ordered[0]
    label = BACKEND_LABELS.get(name, name)
    next_best = ordered[1][1]
    if math.isclose(best, next_best, rel_tol=5.0e-3):
        return "tie"
    if best == 0.0:
        return f"{label} (timer resolution)"
    return f"{label} {next_best / best:.2f}×"


def _backend_header_columns(backends: Sequence[str]) -> tuple[str, str]:
    return (
        "| Total gluons | "
        + " | ".join(BACKEND_LABELS.get(name, name) for name in backends)
        + " | Fastest vs next |",
        "|---:|" + "|".join("---:" for _ in backends) + "|---:|",
    )


def _backend_timing_cells(
    result: MultiplicityResult, backends: Sequence[str], timings: Mapping[str, float]
) -> list[str]:
    cells = [str(result.total_gluons)]
    for backend in backends:
        cells.append(format_duration(timings.get(backend, float("nan"))))
    cells.append(fastest(timings))
    return cells


def _backend_labels(backends: Sequence[str]) -> str:
    return " | ".join(BACKEND_LABELS.get(name, name) for name in backends)


def _markdown_table_text(value: object) -> str:
    text = str(value)
    maximum_length = 500
    if len(text) > maximum_length:
        text = text[:maximum_length].rstrip() + "…"
    return text.replace("|", "\\|").replace("\n", "<br>")


def _median_batch_value(values: Sequence[float] | None) -> float:
    return statistics.median(values) if values else float("nan")


def render_report(
    arguments: argparse.Namespace,
    results: Sequence[MultiplicityResult],
    compiler_name: str,
    table_hash: str,
    completed: int | None = None,
) -> str:
    selected_backends = list(arguments.backends)
    if not selected_backends:
        raise BenchmarkError("no backends selected")

    lines = [
        "# Proxy-weighted sampled-helicity all-gluon matrix-element benchmark",
        "",
        (
            f"RAMBO-on-diet phase-space points at `sqrt(s) = {arguments.sqrt_s:g}`, "
            f"{arguments.points} point(s) per multiplicity. A canonical "
            "leading-colour proxy weight is produced for all `2^N` physical "
            "helicities (with analytic zeros assigned directly), while the "
            "expensive full-colour backends use a bounded proxy-weighted sample "
            "pool. Optimized backends use every selected nonzero row; default-BG "
            "backends use its representative row and the common zero sentinel."
        ),
        "",
        "Selected backends: `"
        + "`, `".join(BACKEND_LABELS.get(name, name) for name in selected_backends)
        + "`.",
    ]
    if completed is not None:
        lines.append(f"Progress: completed up to N={completed} total gluons.")

    lines.extend(
        [
            "",
            "## Problem dimensions",
            "",
            "| Total gluons | " + _backend_labels(selected_backends) + " |",
            "|---:|" + "|".join("---:" for _ in selected_backends) + "|",
        ]
    )
    for result in results:
        row = [str(result.total_gluons)]
        for backend in selected_backends:
            run = result.runs.get(backend)
            row.append(str(run.dimension) if run is not None else "N/A")
        lines.append("| " + " | ".join(row) + " |")

    lines.extend(
        [
            "",
            "## Peak resident memory",
            "",
            (
                "Maximum resident set size of the full fresh-process run, including "
                "initialization and all requested evaluations."
            ),
            "",
            "| Total gluons | " + _backend_labels(selected_backends) + " |",
            "|---:|" + "|".join("---:" for _ in selected_backends) + "|",
        ]
    )
    for result in results:
        row = [str(result.total_gluons)]
        for backend in selected_backends:
            run = result.runs.get(backend)
            row.append(format_memory(run.peak_rss_kib if run is not None else None))
        lines.append("| " + " | ".join(row) + " |")

    header, divider = _backend_header_columns(selected_backends)
    sampled_pass_header = (
        "| Total gluons | " + _backend_labels(selected_backends) + " |"
    )
    sampled_pass_divider = (
        "|---:|" + "|".join("---:" for _ in selected_backends) + "|"
    )
    lines.extend(
        [
            "",
            "## Initialization",
            "",
            f"Median of {arguments.initialization_runs} fresh-process runs.",
            "",
            header,
            divider,
        ]
    )
    for result in results:
        timings: dict[str, float] = {}
        for backend in selected_backends:
            samples = result.initialization_samples.get(backend)
            timings[backend] = statistics.median(samples) if samples else float("nan")
        lines.append("| " + " | ".join(_backend_timing_cells(result, selected_backends, timings)) + " |")

    lines.extend(
        [
            "",
            "## First sampled pass",
            "",
            (
                "CPU time for one pass over the selected validation rows. This "
                "includes any lazy MHV or BG setup triggered by the sample."
            ),
            "",
            sampled_pass_header,
            sampled_pass_divider,
        ]
    )
    for result in results:
        timings = {backend: float("nan") for backend in selected_backends}
        for backend in result.runs:
            timings[backend] = result.runs[backend].first_helicity_sweep or 0.0
        cells = [str(result.total_gluons)]
        cells.extend(format_duration(timings[backend]) for backend in selected_backends)
        lines.append("| " + " | ".join(cells) + " |")

    lines.extend(
        [
            "",
            "## Estimated production-workload warm evaluation",
            "",
            (
                f"Median of {arguments.batches} calibrated batches. In each batch, "
                "an optimized-backend estimate is "
                "`f_MHV t_MHV + (1-f_MHV) t_BG`; the table reports the median of "
                "those mixtures. Here `f_MHV` is the global canonical "
                "leading-colour proxy fraction. A default-BG backend uses its "
                "representative nonzero BG timing directly."
            ),
            "",
            header,
            divider,
        ]
    )
    for result in results:
        timings = {backend: float("nan") for backend in selected_backends}
        for backend in result.runs:
            batches = result.runs[backend].weighted_batches
            if batches:
                timings[backend] = statistics.median(batches)
        lines.append("| " + " | ".join(_backend_timing_cells(result, selected_backends, timings)) + " |")

    lines.extend(
        [
            "",
            "## Warm MHV/anti-MHV path",
            "",
            (
                f"Median of {arguments.batches} batches over "
                f"{arguments.mhv_samples} deterministic proxy-weighted sample(s)."
            ),
            "",
            header,
            divider,
        ]
    )
    for result in results:
        timings = {backend: float("nan") for backend in selected_backends}
        for backend, run in result.runs.items():
            timings[backend] = _median_batch_value(run.mhv_batches)
        lines.append(
            "| "
            + " | ".join(_backend_timing_cells(result, selected_backends, timings))
            + " |"
        )

    lines.extend(
        [
            "",
            "## Warm Berends-Giele path",
            "",
            (
                f"Median of {arguments.batches} batches. Optimized backends use "
                "the selected general-helicity samples. Default-BG backends use "
                "one representative nonzero row because their operation count "
                "does not depend on the nonzero helicity; at N=4 and N=5 that "
                "representative is necessarily MHV."
            ),
            "",
            header,
            divider,
        ]
    )
    for result in results:
        timings = {backend: float("nan") for backend in selected_backends}
        for backend in result.runs:
            timings[backend] = _median_batch_value(result.runs[backend].bg_batches)
        lines.append(
            "| "
            + " | ".join(_backend_timing_cells(result, selected_backends, timings))
            + " |"
        )

    lines.extend(
        [
            "",
            "## Numerical agreement",
            "",
            "| Total gluons | Compared backends | Max relative spread (nonzero) | Max absolute value (zero sector) |",
            "|---:|:---|---:|---:|",
        ]
    )
    for result in results:
        compared = (
            ", ".join(
                BACKEND_LABELS.get(name, name) for name in result.compared_backends
            )
            if result.compared_backends
            else "N/A"
        )
        max_relative = (
            f"{result.maximum_relative_difference:.3e}"
            if math.isfinite(result.maximum_relative_difference)
            else "N/A"
        )
        max_zero = (
            f"{result.maximum_zero_value:.3e}"
            if math.isfinite(result.maximum_zero_value)
            else "N/A"
        )
        lines.append(
            f"| {result.total_gluons} | {compared} | {max_relative} | {max_zero} |"
        )

    lines.extend(
        [
            "",
            "## Leading-colour helicity proxy",
            "",
            (
                "The fraction is computed from the exhaustive canonical-order "
                "proxy before sampling: `sum_p,h in MHV |A_p,h|^2 / "
                "sum_p,h |A_p,h|^2`. The global ratio represents points drawn "
                "according to the proxy helicity-summed result. Because this is "
                "one colour order rather than the full leading-colour sum, the "
                "combined timing is explicitly an approximation."
            ),
            "",
            "| Total gluons | Global MHV fraction | Point-to-point range |",
            "|---:|---:|---:|",
        ]
    )
    for result in results:
        fractions = list(result.proxy_point_mhv_fractions.values())
        lines.append(
            f"| {result.total_gluons} | {result.proxy_mhv_fraction:.6f} | "
            f"{min(fractions):.6f} – {max(fractions):.6f} |"
        )

    lines.extend(
        [
            "",
            "## Selected helicity samples",
            "",
            (
                "Original event/configuration indices, crossed-helicity path, and "
                "canonical leading-colour proxy weight in the selected pool. "
                "Sample 1 is drawn from the complete nonzero proxy distribution. "
                "Optimized backends use the full pool; default-BG backends use "
                "samples 1 and 2. Repeated draws are retained."
            ),
            "",
            "| Gluons | Sample | Source point | Source helicity | Path | Proxy weight |",
            "|---:|---:|---:|---:|:---|---:|",
        ]
    )
    for result in results:
        for sample_index, sample in enumerate(result.samples, start=1):
            lines.append(
                f"| {result.total_gluons} | {sample_index} | "
                f"{sample.source_event} | {sample.source_configuration} | "
                f"{sample.path} | {sample.proxy_weight:.12e} |"
            )

    lines.extend(
        [
            "",
            "## Skipped backends",
            "",
            "| Total gluons | Backend | Reason |",
            "|---:|---|---|",
        ]
    )
    skipped_any = False
    for result in results:
        for backend, reason in result.skipped_backends.items():
            skipped_any = True
            lines.append(
                f"| {result.total_gluons} | {BACKEND_LABELS.get(backend, backend)} | "
                f"{_markdown_table_text(reason)} |"
            )
    if not skipped_any:
        lines.append("| - | - | no backends were skipped |")

    lines.extend(
        [
            "",
            "## Reproducibility",
            "",
            "| Total gluons | Generator seed | Combined event SHA-256 |",
            "|---:|---:|:---|",
        ]
    )
    for result in results:
        lines.append(
            f"| {result.total_gluons} | {arguments.seed + result.total_gluons} | "
            f"`{result.event_hash}` |"
        )

    lines.extend(
        [
            "",
            f"Compiler: `{compiler_name}`",
            "",
            f"Flags for selected backends: `{arguments.fflags}`",
            "",
            (
                f"Backend address-space cap: `{arguments.max_memory_gib:.3g} GiB`; "
                "peak resident memory is measured with GNU time"
            ),
            f"Backend invocation timeout: `{arguments.backend_timeout:.3g} s`",
            (
                "Separate initialization preflight: "
                f"`{'disabled' if arguments.skip_initialization_preflight else 'enabled'}`"
            ),
            "Minimum calibration target per sampled timing batch: "
            f"`{arguments.target_seconds:.3g} s`",
            (
                "Requested expensive samples per multiplicity: "
                f"`{arguments.mhv_samples}` MHV/anti-MHV and "
                f"`{arguments.non_mhv_samples}` general-helicity; one additional "
                "analytic-zero row is used only for validation."
            ),
            "",
            (
                "The selected analytically vanishing row must be below "
                f"{arguments.zero_tolerance:.3e}; selected nonzero rows are "
                "compared relatively. Zero sectors have zero proxy weight and are "
                "excluded from the calibrated timing cells and production "
                "estimate."
            ),
        ]
    )
    if arguments.table is not None:
        lines.extend(
            [
                f"Wigner table: `{arguments.table}`",
                (
                    f"Wigner prefix depth: `{arguments.table_depth}` "
                    f"(supports through `{2 * arguments.table_depth + 1}` total gluons)"
                ),
                f"Wigner table SHA-256: `{table_hash}`",
            ]
        )

    return "\n".join(lines) + "\n"


def main() -> int:
    try:
        arguments = parse_arguments()
        validate_arguments(arguments)
        compiler = shlex.split(arguments.fc)
        flags = shlex.split(arguments.fflags)
        if not compiler:
            raise BenchmarkError("--fc must not be empty")
        arguments.build_dir.mkdir(parents=True, exist_ok=True)
        compiler_name = compiler_version(compiler).strip()
        table_hash = (
            hashlib.sha256(arguments.table.read_bytes()).hexdigest()
            if arguments.table is not None
            else "n/a"
        )
        rambo_executable, proxy_executable, executables = build_executables(
            arguments, compiler, compiler_name, flags
        )
        results = collect_results(
            arguments,
            rambo_executable,
            proxy_executable,
            executables,
            compiler_name,
            table_hash,
        )
        report = render_report(
            arguments=arguments,
            results=results,
            compiler_name=compiler_name,
            table_hash=table_hash,
        )
        if arguments.output is not None:
            arguments.output.parent.mkdir(parents=True, exist_ok=True)
            arguments.output.write_text(report, encoding="utf-8")
            print(f"Wrote {arguments.output}", file=sys.stderr)
        print(report, end="")
        return 0
    except (BenchmarkError, OSError, ValueError) as error:
        print(f"benchmark: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
