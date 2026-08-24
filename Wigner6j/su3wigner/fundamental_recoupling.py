"""Adjacent-adjoint recoupling through fundamental line splitting.

This backend is deliberately independent of the full-index contraction in
``recoupling.py``.  It resolves each adjoint as the traceless subspace of
``3 tensor 3-bar``, swaps the two resulting index pairs, and evaluates the
invariant endomorphism on one highest-weight state of the final irrep.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from .coupling import AdjointCouplings, Coupling
from .fundamental import adjoint_split_isometry
from .recoupling import Path, SwapBlock, SwapTableBuilder
from .representations import Label


class FundamentalSplitSwapTableBuilder(SwapTableBuilder):
    """Build the same swap blocks by resolving both gluons into ``3 x 3-bar``."""

    def __init__(self, couplings: AdjointCouplings | None = None) -> None:
        super().__init__(couplings)
        adjoint = self.couplings.irreps.get((1, 1))
        self.splitting = adjoint_split_isometry(adjoint)
        self._pair_splitting: np.ndarray | None = None

    @property
    def pair_splitting(self) -> np.ndarray:
        """Two-adjoint split isometry, retained as a lazy public attribute."""

        if self._pair_splitting is None:
            self._pair_splitting = np.einsum(
                "aij,bkl->abijkl", self.splitting, self.splitting
            ).reshape(self.splitting.shape[0] ** 2, 3**4)
        return self._pair_splitting

    @lru_cache(maxsize=None)
    def blocks_for_left(self, left: Label) -> tuple[SwapBlock, ...]:
        path_data = self._path_couplings_for_left(left)
        blocks: list[SwapBlock] = []
        for right in sorted(path_data):
            ordered = sorted(path_data[right], key=lambda item: item[0])
            paths = tuple(item[0] for item in ordered)
            left_dimension = self.couplings.irreps.get(left).dim
            adjoint_dimension = self.splitting.shape[0]
            original = np.empty(
                (
                    len(ordered),
                    left_dimension,
                    adjoint_dimension,
                    adjoint_dimension,
                ),
                dtype=float,
            )
            for index, (_, first, second) in enumerate(ordered):
                original[index] = self._compose_highest(first, second)
            exchanged = np.swapaxes(original, 2, 3)
            original_flat = original.reshape(len(ordered), -1)
            exchanged_flat = exchanged.reshape(len(ordered), -1)
            matrix = exchanged_flat @ original_flat.T
            self._clean_matrix(matrix)
            block = SwapBlock(left=left, right=right, paths=paths, matrix=matrix)
            self.validate_block(block)
            self._validate_split_equation(
                block,
                original_flat,
                exchanged_flat,
            )
            blocks.append(block)
        return tuple(blocks)

    def _compose_highest(
        self, first: Coupling, second: Coupling
    ) -> np.ndarray:
        """Return one normalized final-state column before the fixed split.

        Splitting both adjoints applies the same isometry ``S tensor S`` to
        every path.  It cancels exactly from both the overlap and reconstruction
        equations, including after pair exchange.  Keeping the compact 8-by-8
        column therefore evaluates the identical fundamental-split coefficient
        without materializing its 9-by-9 image.
        """

        left = self.couplings.irreps.get(first.source)
        right = self.couplings.irreps.get(second.target)
        adjoint_dimension = self.splitting.shape[0]
        highest_candidates = [
            index for index, weight in enumerate(right.weights) if weight == right.label
        ]
        if len(highest_candidates) != 1:
            raise ArithmeticError(
                f"{right.label} has {len(highest_candidates)} highest states"
            )
        highest = highest_candidates[0]
        # Contract only the exact-weight entries that can contribute to this
        # highest-state column.  Materializing the complete first vertex here
        # would defeat packed CG storage for every path through a large middle
        # representation.
        result = np.zeros(
            (left.dim * adjoint_dimension, adjoint_dimension), dtype=float
        )
        second_rows, second_values = second._column_entries(highest)
        first_columns: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        for second_row, second_value in zip(
            second_rows, second_values, strict=True
        ):
            middle_index, second_adjoint = divmod(
                int(second_row), adjoint_dimension
            )
            first_column = first_columns.get(middle_index)
            if first_column is None:
                first_column = first._column_entries(middle_index)
                first_columns[middle_index] = first_column
            first_rows, first_values = first_column
            result[first_rows, second_adjoint] += first_values * second_value
        return result.reshape(left.dim, adjoint_dimension, adjoint_dimension)

    @staticmethod
    def _validate_split_equation(
        block: SwapBlock,
        original_flat: np.ndarray,
        exchanged_flat: np.ndarray,
    ) -> None:
        if any(
            not np.all(np.isfinite(array))
            for array in (block.matrix, original_flat, exchanged_flat)
        ):
            raise ArithmeticError(
                "fundamental-split swap equation contains a non-finite value"
            )
        identity = np.eye(block.size)
        original_gram = original_flat @ original_flat.T
        exchanged_gram = exchanged_flat @ exchanged_flat.T
        reconstructed = block.matrix.T @ exchanged_flat
        if any(
            not np.all(np.isfinite(array))
            for array in (original_gram, exchanged_gram, reconstructed)
        ):
            raise ArithmeticError(
                "fundamental-split swap residual contains a non-finite value"
            )
        basis_error = max(
            float(np.max(np.abs(original_gram - identity))),
            float(np.max(np.abs(exchanged_gram - identity))),
        )
        equation_error = float(np.max(np.abs(original_flat - reconstructed)))
        if max(basis_error, equation_error) > 2.0e-8 * max(1, block.size):
            raise ArithmeticError(
                "fundamental-split swap equation failed for "
                f"{block.left}->{block.right}: basis={basis_error:.3e}, "
                f"equation={equation_error:.3e}"
            )
