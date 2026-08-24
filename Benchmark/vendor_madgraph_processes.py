#!/usr/bin/env python3
"""Vendor lightweight fixed-helicity all-gluon MadGraph standalones."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile


BENCHMARK_DIR = Path(__file__).resolve().parent
DEFAULT_DESTINATION = BENCHMARK_DIR / "MadGraph5"
SHARED_DHELAS_SOURCES = (
    "aloha_functions.f",
    "VVV1P0_1.f",
    "VVV1_0.f",
    "VVVV1P0_1.f",
    "VVVV1_0.f",
    "VVVV3P0_1.f",
    "VVVV3_0.f",
    "VVVV4P0_1.f",
    "VVVV4_0.f",
)


def git_object(repository: Path, revision: str) -> str | None:
    """Return a Git object name when generation uses a checkout."""

    result = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "--verify", revision],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if result.returncode != 0 or not re.fullmatch(r"[0-9a-f]{40}\n?", result.stdout):
        return None
    return result.stdout.strip()


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mg5",
        type=Path,
        required=True,
        help="path to the MG5_aMC executable used for generation",
    )
    parser.add_argument("--min-gluons", type=int, default=4)
    parser.add_argument("--max-gluons", type=int, default=6)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--timeout", type=float, default=3600.0)
    return parser.parse_args()


def process_command(total_gluons: int, output: Path) -> str:
    final_state = " ".join("g" for _ in range(total_gluons - 2))
    return "\n".join(
        (
            "set automatic_html_opening False",
            "import model sm",
            f"generate g g > {final_state} QED=0",
            f"output standalone {output} -f",
            "quit",
            "",
        )
    )


def copy_file(source: Path, destination: Path) -> None:
    """Copy generated text with stable whitespace and one final newline."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    text = source.read_text(encoding="utf-8")
    normalized = "\n".join(line.rstrip() for line in text.splitlines()) + "\n"
    destination.write_text(normalized, encoding="utf-8")


def sanitize_process_card(card: Path, process: str) -> None:
    """Replace transient output paths and strip terminal colour escapes."""

    text = card.read_text(encoding="utf-8")
    text = re.sub(r"\x1b\[[0-9;]*m", "", text)
    text = re.sub(
        r"(?m)^output standalone .* -f$",
        f"output standalone standalone_{process.replace(' ', '_')} -f",
        text,
    )
    card.write_text(text, encoding="utf-8")


def main() -> int:
    arguments = parse_arguments()
    mg5 = arguments.mg5.expanduser().resolve()
    destination = arguments.destination.expanduser().resolve()
    if not mg5.is_file():
        raise SystemExit(f"MG5_aMC executable not found: {mg5}")
    if arguments.min_gluons < 4 or arguments.max_gluons < arguments.min_gluons:
        raise SystemExit("invalid total-gluon range")
    mg5_root = mg5.parent.parent
    version_file = mg5_root / "VERSION"
    version_text = version_file.read_text(encoding="utf-8")
    copy_file(version_file, destination / "VERSION")
    generator_metadata = {
        "reported_version": re.search(r"version\s*=\s*([^\n]+)", version_text)
        .group(1)
        .strip(),
        "source_git_commit": git_object(mg5_root, "HEAD"),
        "sm_model_git_tree": git_object(mg5_root, "HEAD:models/sm"),
        "generation_constraint": "QED=0",
    }
    (destination / "GENERATOR.json").write_text(
        json.dumps(generator_metadata, indent=2) + "\n", encoding="utf-8"
    )

    with tempfile.TemporaryDirectory(prefix="ampligluon-mg5-") as temporary:
        temporary_root = Path(temporary)
        for total_gluons in range(arguments.min_gluons, arguments.max_gluons + 1):
            standalone = temporary_root / f"standalone_N{total_gluons}"
            command_file = temporary_root / f"generate_N{total_gluons}.mg5"
            command_file.write_text(
                process_command(total_gluons, standalone), encoding="utf-8"
            )
            process = subprocess.run(
                [str(mg5), str(command_file)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
                timeout=arguments.timeout,
            )
            if process.returncode != 0:
                raise SystemExit(
                    f"MadGraph generation failed at N={total_gluons}:\n{process.stdout}"
                )
            subprocesses = sorted(
                path
                for path in (standalone / "SubProcesses").glob("P*")
                if path.is_dir() and (path / "matrix.f").is_file()
            )
            if len(subprocesses) != 1:
                raise SystemExit(
                    f"expected one generated subprocess at N={total_gluons}, "
                    f"found {len(subprocesses)}"
                )
            generated = subprocesses[0]
            process_dir = destination / "Processes" / f"N{total_gluons}"
            for filename in ("matrix.f", "nexternal.inc", "ngraphs.inc"):
                copy_file(generated / filename, process_dir / filename)
            copy_file(
                standalone / "Cards" / "proc_card_mg5.dat",
                process_dir / "proc_card_mg5.dat",
            )
            sanitize_process_card(
                process_dir / "proc_card_mg5.dat",
                f"gg_to_{total_gluons - 2}g",
            )
            matrix = process_dir / "matrix.f"
            matrix_text = matrix.read_text(encoding="utf-8")
            diagrams_match = re.search(r"Process has (\d+) diagrams", process.stdout)
            matrix_graphs_match = re.search(r"PARAMETER \(NGRAPHS=(\d+)\)", matrix_text)
            colours_match = re.search(
                r"PARAMETER \(NWAVEFUNCS=\d+, NCOLOR=(\d+)\)", matrix_text
            )
            metadata = {
                "total_gluons": total_gluons,
                "process": f"g g > {' '.join('g' for _ in range(total_gluons - 2))}",
                "madgraph_version": re.search(r"version\s*=\s*([^\n]+)", version_text)
                .group(1)
                .strip(),
                "feynman_diagrams": (
                    int(diagrams_match.group(1)) if diagrams_match else None
                ),
                "generated_matrix_graphs": (
                    int(matrix_graphs_match.group(1)) if matrix_graphs_match else None
                ),
                "colour_flows": (
                    int(colours_match.group(1)) if colours_match else None
                ),
                "matrix_sha256": hashlib.sha256(matrix.read_bytes()).hexdigest(),
            }
            (process_dir / "process.json").write_text(
                json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
            )

            source = standalone / "Source"
            shared_files = [
                (
                    standalone / "SubProcesses" / "coupl.inc",
                    destination / "Source" / "coupl.inc",
                ),
                *(
                    (
                        source / "DHELAS" / filename,
                        destination / "Source" / "DHELAS" / filename,
                    )
                    for filename in SHARED_DHELAS_SOURCES
                ),
            ]
            for shared_source, shared_destination in shared_files:
                if not shared_source.is_file():
                    continue
                if (
                    shared_destination.is_file()
                    and shared_destination.read_text(encoding="utf-8")
                    != "\n".join(
                        line.rstrip()
                        for line in shared_source.read_text(
                            encoding="utf-8"
                        ).splitlines()
                    )
                    + "\n"
                ):
                    raise SystemExit(
                        f"inconsistent shared MadGraph source: {shared_source.name}"
                    )
                copy_file(shared_source, shared_destination)
            print(
                f"Vendored N={total_gluons}: {metadata['feynman_diagrams']} "
                f"diagrams, {matrix.stat().st_size / 1024.0:.1f} KiB matrix.f"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
