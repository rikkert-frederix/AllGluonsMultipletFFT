"""Independent identities for the recursively generated swap coefficients."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .recoupling import Path, SwapBlock, SwapTableBuilder
from .representations import Label


@dataclass(frozen=True, order=True, slots=True)
class ThreeAdjointPath:
    """``left --8-> first --8-> second --8-> right``."""

    first: Label
    first_multiplicity: int
    second: Label
    second_multiplicity: int
    third_multiplicity: int


def _path_index(block: SwapBlock, middle: Label, left_mult: int, right_mult: int) -> int:
    matches = [
        index
        for index, path in enumerate(block.paths)
        if path.middle == middle
        and path.left_multiplicity == left_mult
        and path.right_multiplicity == right_mult
    ]
    if len(matches) != 1:
        raise ArithmeticError(
            f"expected one path through {middle}[{left_mult},{right_mult}], "
            f"found {len(matches)}"
        )
    return matches[0]


def braid_residuals(
    builder: SwapTableBuilder, left: Label
) -> dict[Label, float]:
    """Check ``s1 s2 s1 = s2 s1 s2`` on three adjacent adjoints.

    This is the component form of the associativity/permutation recursion.
    It links blocks with prefix ``left`` to blocks whose prefixes occur in
    ``left tensor 8`` and is independent of any closed-graph normalization.
    """

    paths_by_right: dict[Label, list[ThreeAdjointPath]] = {}
    for first_edge in builder.couplings.decompose(left).couplings:
        for second_edge in builder.couplings.decompose(first_edge.target).couplings:
            for third_edge in builder.couplings.decompose(second_edge.target).couplings:
                path = ThreeAdjointPath(
                    first=first_edge.target,
                    first_multiplicity=first_edge.multiplicity,
                    second=second_edge.target,
                    second_multiplicity=second_edge.multiplicity,
                    third_multiplicity=third_edge.multiplicity,
                )
                paths_by_right.setdefault(third_edge.target, []).append(path)

    block_cache: dict[Label, dict[Label, SwapBlock]] = {}

    def block(prefix: Label, suffix: Label) -> SwapBlock:
        if prefix not in block_cache:
            block_cache[prefix] = {
                candidate.right: candidate
                for candidate in builder.blocks_for_left(prefix)
            }
        return block_cache[prefix][suffix]

    residuals: dict[Label, float] = {}
    for right, unsorted_paths in paths_by_right.items():
        paths = tuple(sorted(unsorted_paths))
        lookup = {path: index for index, path in enumerate(paths)}
        size = len(paths)
        swap_1 = np.zeros((size, size), dtype=float)
        swap_2 = np.zeros((size, size), dtype=float)
        for input_index, input_path in enumerate(paths):
            first_block = block(left, input_path.second)
            first_input = _path_index(
                first_block,
                input_path.first,
                input_path.first_multiplicity,
                input_path.second_multiplicity,
            )
            for output_local, output_path in enumerate(first_block.paths):
                output = ThreeAdjointPath(
                    first=output_path.middle,
                    first_multiplicity=output_path.left_multiplicity,
                    second=input_path.second,
                    second_multiplicity=output_path.right_multiplicity,
                    third_multiplicity=input_path.third_multiplicity,
                )
                swap_1[lookup[output], input_index] = first_block.matrix[
                    output_local, first_input
                ]

            second_block = block(input_path.first, right)
            second_input = _path_index(
                second_block,
                input_path.second,
                input_path.second_multiplicity,
                input_path.third_multiplicity,
            )
            for output_local, output_path in enumerate(second_block.paths):
                output = ThreeAdjointPath(
                    first=input_path.first,
                    first_multiplicity=input_path.first_multiplicity,
                    second=output_path.middle,
                    second_multiplicity=output_path.left_multiplicity,
                    third_multiplicity=output_path.right_multiplicity,
                )
                swap_2[lookup[output], input_index] = second_block.matrix[
                    output_local, second_input
                ]

        left_side = swap_1 @ swap_2 @ swap_1
        right_side = swap_2 @ swap_1 @ swap_2
        braid = float(np.max(np.abs(left_side - right_side)))
        involution_1 = float(np.max(np.abs(swap_1 @ swap_1 - np.eye(size))))
        involution_2 = float(np.max(np.abs(swap_2 @ swap_2 - np.eye(size))))
        residuals[right] = max(braid, involution_1, involution_2)
    return residuals


def assert_braid_relations(
    builder: SwapTableBuilder,
    left_labels: tuple[Label, ...],
    *,
    tolerance: float = 5.0e-8,
) -> None:
    for left in left_labels:
        residuals = braid_residuals(builder, left)
        for right, residual in residuals.items():
            if not np.isfinite(residual) or residual > tolerance:
                raise ArithmeticError(
                    f"braid recursion failed for {left} -> {right}: {residual:.3e}"
                )
