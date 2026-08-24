from __future__ import annotations

import unittest

import numpy as np

from su3wigner import (
    AdjointCouplings,
    FundamentalSplitSwapTableBuilder,
    Path,
    SwapBlock,
    SwapTableBuilder,
    adjoint_split_isometry,
)


class FundamentalRecouplingTests(unittest.TestCase):
    def test_adjoint_resolution_is_the_traceless_projector(self) -> None:
        couplings = AdjointCouplings()
        splitting = adjoint_split_isometry(couplings.irreps.get((1, 1)))
        flattened = splitting.reshape(8, 9)
        singlet = np.eye(3).reshape(9) / np.sqrt(3.0)
        expected = np.eye(9) - np.outer(singlet, singlet)
        np.testing.assert_allclose(flattened @ flattened.T, np.eye(8), atol=2.0e-14)
        np.testing.assert_allclose(flattened.T @ flattened, expected, atol=2.0e-14)

    def test_cutoff_two_matches_full_adjoint_index_contraction(self) -> None:
        couplings = AdjointCouplings()
        explicit = SwapTableBuilder(couplings).build(2)
        split = FundamentalSplitSwapTableBuilder(couplings).build(2)
        self.assertEqual(len(explicit), len(split))
        for expected, actual in zip(explicit, split, strict=True):
            self.assertEqual(expected.left, actual.left)
            self.assertEqual(expected.right, actual.right)
            self.assertEqual(expected.paths, actual.paths)
            np.testing.assert_allclose(actual.matrix, expected.matrix, atol=3.0e-14)

    def test_split_equation_validation_rejects_nonfinite_data(self) -> None:
        block = SwapBlock(
            left=(0, 0),
            right=(0, 0),
            paths=(
                Path(
                    middle=(1, 1),
                    left_multiplicity=0,
                    right_multiplicity=0,
                ),
            ),
            matrix=np.ones((1, 1), dtype=float),
        )
        finite = np.ones((1, 1), dtype=float)
        nonfinite = np.full((1, 1), np.nan)
        for name, matrix, original, exchanged in (
            ("matrix", nonfinite, finite, finite),
            ("original", finite, nonfinite, finite),
            ("exchanged", finite, finite, nonfinite),
        ):
            with self.subTest(name=name):
                candidate = SwapBlock(
                    left=block.left,
                    right=block.right,
                    paths=block.paths,
                    matrix=matrix,
                )
                with self.assertRaisesRegex(ArithmeticError, "non-finite"):
                    FundamentalSplitSwapTableBuilder._validate_split_equation(
                        candidate,
                        original,
                        exchanged,
                    )


if __name__ == "__main__":
    unittest.main()
