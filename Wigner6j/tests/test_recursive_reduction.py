from __future__ import annotations

import unittest

import numpy as np

from su3wigner import (
    AdjointCouplings,
    RecursiveReductionSwapTableBuilder,
    SwapTableBuilder,
)
from su3wigner.recursive_reduction import (
    FUNDAMENTAL,
    FundamentalCouplings,
    FundamentalRecouplings,
)


class RecursiveReductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.adjoint_couplings = AdjointCouplings()
        cls.direct = SwapTableBuilder(cls.adjoint_couplings)
        cls.recursive = RecursiveReductionSwapTableBuilder(
            cls.adjoint_couplings
        )

    def test_terminal_three_quark_recouplings(self) -> None:
        """The terminal cases follow directly from S3 permutation symmetry."""

        elementary = FundamentalRecouplings(
            FundamentalCouplings(self.adjoint_couplings.irreps)
        )
        singlet = elementary.block(
            (1, 0), FUNDAMENTAL, FUNDAMENTAL, (0, 0)
        )
        decuplet = elementary.block(
            (1, 0), FUNDAMENTAL, FUNDAMENTAL, (3, 0)
        )
        octets = elementary.block(
            (1, 0), FUNDAMENTAL, FUNDAMENTAL, (1, 1)
        )
        np.testing.assert_array_equal(singlet.matrix, [[-1.0]])
        np.testing.assert_array_equal(decuplet.matrix, [[1.0]])
        np.testing.assert_allclose(
            octets.matrix,
            [[0.5, np.sqrt(3.0) / 2.0], [np.sqrt(3.0) / 2.0, -0.5]],
            atol=2.0e-14,
        )

    def test_four_terminal_moves_equal_literal_pair_swap(self) -> None:
        recouplings = self.recursive.fundamental_recouplings
        singlet_targets = [
            block.right for block in self.recursive.blocks_for_left((0, 0))
        ]
        residuals = [
            recouplings.pair_swap_direct_residual((0, 0), target)
            for target in singlet_targets
        ]
        residuals.append(
            recouplings.pair_swap_direct_residual((1, 1), (2, 2))
        )
        self.assertLess(max(residuals), 8.0e-14)

    def test_repeated_adjoint_vertex_gauge_is_preserved(self) -> None:
        reductions = self.recursive.vertex_reductions
        symmetric = reductions.reduce((1, 1), (1, 1), 0)
        antisymmetric = reductions.reduce((1, 1), (1, 1), 1)
        self.assertEqual(symmetric.intermediates, antisymmetric.intermediates)
        np.testing.assert_allclose(
            np.dot(symmetric.coefficients, antisymmetric.coefficients),
            0.0,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            np.dot(symmetric.coefficients, symmetric.coefficients),
            1.0,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            np.dot(antisymmetric.coefficients, antisymmetric.coefficients),
            1.0,
            atol=2.0e-14,
        )

    def test_terminal_validation_rejects_nonfinite_data(self) -> None:
        finite = np.ones((1, 1), dtype=float)
        nonfinite = np.full((1, 1), np.nan)
        for name, matrix, original, exchanged in (
            ("matrix", nonfinite, finite, finite),
            ("original", finite, nonfinite, finite),
            ("exchanged", finite, finite, nonfinite),
        ):
            with self.subTest(name=name):
                with self.assertRaisesRegex(ArithmeticError, "non-finite"):
                    FundamentalRecouplings._validate_terminal(
                        matrix,
                        original,
                        exchanged,
                        (0, 0),
                        (0, 0),
                    )

    def test_cutoff_two_matches_full_adjoint_contraction(self) -> None:
        expected = self.direct.build(2)
        actual = self.recursive.build(2)
        self.assertEqual(len(expected), len(actual))
        for direct, reduced in zip(expected, actual, strict=True):
            self.assertEqual(direct.left, reduced.left)
            self.assertEqual(direct.right, reduced.right)
            self.assertEqual(direct.paths, reduced.paths)
            np.testing.assert_allclose(
                reduced.matrix, direct.matrix, atol=3.0e-14
            )


if __name__ == "__main__":
    unittest.main()
