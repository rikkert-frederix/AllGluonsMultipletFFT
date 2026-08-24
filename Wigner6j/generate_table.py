#!/usr/bin/env python3
"""Generate the SU(3) adjacent-adjoint swap table used by Eq. (8)."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import tempfile
import time

try:
    import numpy as np
except ImportError as error:  # pragma: no cover - exercised only on bad setups
    raise SystemExit(
        "NumPy is required. Run this program with a Python 3.10+ environment "
        "containing NumPy (see requirements.txt)."
    ) from error

from su3wigner import (
    DimensionOnlySwapTableBuilder,
    DirectSpechtSwapTableBuilder,
    FundamentalSplitSwapTableBuilder,
    RecursiveReductionSwapTableBuilder,
    SwapTableBuilder,
    braid_residuals,
    read_table,
    write_table,
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate normalized SU(3) recoupling matrices for swapping two "
            "adjacent adjoints after every gluonic prefix representation."
        )
    )
    parser.add_argument(
        "--method",
        choices=(
            "explicit-index",
            "fundamental-split",
            "recursive-reduction",
            "dimension-only",
            "direct-specht",
        ),
        default="explicit-index",
        help=(
            "coefficient construction: trace full adjoint-index chains, or "
            "split each adjoint into a traceless 3 x 3-bar pair and evaluate "
            "one highest-weight state, or recursively reduce two-adjoint "
            "recouplings to elementary fundamental ones, or use analytic "
            "fusion labels, dimensions, Casimirs, and Young-content swap "
            "matrices without constructing CG tensors, or use the direct "
            "six-box skew-Specht representation "
            "(default: explicit-index)"
        ),
    )
    parser.add_argument(
        "--cross-check",
        action="store_true",
        help=(
            "also build with every other full-table method"
        ),
    )
    parser.add_argument(
        "--max-prefix-gluons",
        type=int,
        default=2,
        help=(
            "include every left/prefix irrep reachable from the singlet with "
            "at most this many adjoints (default: 2)"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="output table (default: data/su3_adjoint_swap_prefix_N.tbl)",
    )
    parser.add_argument(
        "--zero-tolerance",
        type=float,
        default=5.0e-14,
        help="omit values whose absolute magnitude is at most this number",
    )
    parser.add_argument(
        "--braid-check-depth",
        type=int,
        default=1,
        help=(
            "check the three-adjoint braid recursion for prefix irreps with "
            "minimum depth at most this value; use -1 to skip (default: 1)"
        ),
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if arguments.max_prefix_gluons < 0:
        raise SystemExit("--max-prefix-gluons must be non-negative")
    if arguments.braid_check_depth < -1:
        raise SystemExit("--braid-check-depth must be -1 or non-negative")
    if (
        not np.isfinite(arguments.zero_tolerance)
        or arguments.zero_tolerance < 0.0
    ):
        raise SystemExit("--zero-tolerance must be finite and non-negative")

    output = arguments.output
    if output is None:
        output = (
            Path(__file__).resolve().parent
            / "data"
            / f"su3_adjoint_swap_prefix_{arguments.max_prefix_gluons}.tbl"
        )
    output.parent.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()

    def make_builder(method: str, couplings=None):
        if method == "dimension-only":
            return DimensionOnlySwapTableBuilder()
        if method == "explicit-index":
            return SwapTableBuilder(couplings)
        if method == "fundamental-split":
            return FundamentalSplitSwapTableBuilder(couplings)
        if method == "recursive-reduction":
            return RecursiveReductionSwapTableBuilder(couplings)
        if method == "direct-specht":
            return DirectSpechtSwapTableBuilder()
        raise AssertionError(f"unhandled coefficient method {method}")

    builder = make_builder(arguments.method)
    reachable = builder.reachable_irreps(arguments.max_prefix_gluons)
    blocks = builder.build(arguments.max_prefix_gluons)

    cross_check_errors: dict[str, float] = {}
    if arguments.cross_check:
        methods = (
            "explicit-index",
            "fundamental-split",
            "recursive-reduction",
            "dimension-only",
            "direct-specht",
        )
        shared_couplings = (
            builder.couplings
            if arguments.method in (
                "explicit-index",
                "fundamental-split",
                "recursive-reduction",
            )
            else None
        )
        for method in methods:
            if method == arguments.method:
                continue
            independent = make_builder(method, shared_couplings)
            if (
                method not in ("dimension-only", "direct-specht")
                and shared_couplings is None
            ):
                shared_couplings = independent.couplings
            comparison = independent.build(arguments.max_prefix_gluons)
            if len(comparison) != len(blocks):
                raise ArithmeticError(
                    f"{method} changed the number of blocks"
                )
            cross_check_error = 0.0
            for primary, copy in zip(blocks, comparison, strict=True):
                if (
                    primary.left != copy.left
                    or primary.right != copy.right
                    or primary.paths != copy.paths
                ):
                    raise ArithmeticError(
                        f"{method} changed a block key or path"
                    )
                cross_check_error = max(
                    cross_check_error,
                    float(np.max(np.abs(primary.matrix - copy.matrix))),
                )
            if cross_check_error > 5.0e-11:
                raise ArithmeticError(
                    f"{method} residual {cross_check_error:.3e}"
                )
            cross_check_errors[method] = cross_check_error

    worst_braid = 0.0
    checked_prefixes = 0
    if arguments.braid_check_depth >= 0:
        labels = tuple(
            label
            for label, depth in reachable.items()
            if depth <= arguments.braid_check_depth
        )
        for label in labels:
            residuals = braid_residuals(builder, label)
            if any(not np.isfinite(value) for value in residuals.values()):
                raise ArithmeticError(
                    f"non-finite braid recursion residual for prefix {label}"
                )
            checked_prefixes += 1
            worst_braid = max(worst_braid, max(residuals.values(), default=0.0))
        if worst_braid > 5.0e-8:
            raise ArithmeticError(f"braid recursion residual {worst_braid:.3e}")

    # Build and validate beside the destination so the final replacement is
    # atomic.  In particular, a failed round-trip must not corrupt an existing
    # production table.
    with tempfile.TemporaryDirectory(
        dir=output.parent, prefix=f".{output.name}."
    ) as temporary_directory:
        temporary_output = Path(temporary_directory) / output.name
        write_table(
            temporary_output,
            blocks,
            max_prefix_gluons=arguments.max_prefix_gluons,
            zero_tolerance=arguments.zero_tolerance,
        )

        # A write/read comparison catches indexing mistakes in the
        # Fortran-facing serialization separately from the group-theory checks.
        restored = read_table(temporary_output)
        if len(restored) != len(blocks):
            raise ArithmeticError("table round-trip changed the number of blocks")
        round_trip_error = 0.0
        for original, copy in zip(blocks, restored, strict=True):
            if (
                original.left != copy.left
                or original.right != copy.right
                or original.paths != copy.paths
            ):
                raise ArithmeticError("table round-trip changed a block key or path")
            round_trip_error = max(
                round_trip_error,
                float(np.max(np.abs(original.matrix - copy.matrix))),
            )
            builder.validate_block(copy, tolerance=5.0e-10)
        if round_trip_error > max(arguments.zero_tolerance, 1.0e-13):
            raise ArithmeticError(
                f"table round-trip residual {round_trip_error:.3e}"
            )

        temporary_output.replace(output)

    nonzero = sum(np.count_nonzero(block.matrix) for block in blocks)
    elapsed = time.perf_counter() - started
    print(f"wrote {output}")
    print(f"method={arguments.method}")
    print(
        f"{len(reachable)} prefix irreps, {len(blocks)} blocks, "
        f"{sum(block.size for block in blocks)} paths, {nonzero} nonzero values"
    )
    validation_basis = {
        "explicit-index": "full adjoint indices",
        "fundamental-split": "split fundamental indices",
        "recursive-reduction": "lower-complexity fundamental recouplings",
        "dimension-only": "analytic labels, dimensions, and Young-content swaps",
        "direct-specht": "direct six-box skew-Specht projectors",
    }[arguments.method]
    print(
        f"checked every block on {validation_basis}; "
        f"braid prefixes={checked_prefixes}, "
        f"worst braid residual={worst_braid:.3e}, elapsed={elapsed:.2f}s"
    )
    for method, error in cross_check_errors.items():
        print(f"cross-check[{method}] residual={error:.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
