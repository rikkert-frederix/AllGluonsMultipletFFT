from __future__ import annotations

from collections import Counter
from pathlib import Path
from unittest import mock
import unittest

import numpy as np

from su3wigner import (
    DimensionOnlySwapTableBuilder,
    SwapTableBuilder,
    braid_residuals,
    read_table,
)
from su3wigner.coupling import AdjointCouplings
from su3wigner.recursive_reduction import FundamentalCouplings
from su3wigner.representations import IrrepFactory, dimension
from su3wigner.dimension_only import (
    ADJOINT,
    ANTIFUNDAMENTAL,
    FUNDAMENTAL,
    _FUSION_PATH_CACHE_SIZE,
    DimensionOnlyRecouplings,
)


class DimensionOnlyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.builder = DimensionOnlySwapTableBuilder()
        cls.direct_blocks = SwapTableBuilder().build(2)

    def test_backend_never_constructs_a_cg_or_representation_object(self) -> None:
        forbidden = AssertionError("dimension-only backend entered a CG backend")
        with (
            mock.patch.object(IrrepFactory, "__init__", side_effect=forbidden),
            mock.patch.object(AdjointCouplings, "__init__", side_effect=forbidden),
            mock.patch.object(FundamentalCouplings, "__init__", side_effect=forbidden),
        ):
            guarded = DimensionOnlySwapTableBuilder()
            blocks = guarded.build(2)
        self.assertEqual(len(blocks), len(self.direct_blocks))
        self.assertFalse(
            hasattr(guarded.couplings.decompose((1, 1)).couplings[0], "embedding")
        )

    def test_elementary_fusion_rules_are_exact(self) -> None:
        fusion = self.builder.elementary_fusion
        self.assertEqual(
            fusion.decompose((2, 1), FUNDAMENTAL).targets,
            ((1, 2), (2, 0), (3, 1)),
        )
        self.assertEqual(
            fusion.decompose((2, 1), ANTIFUNDAMENTAL).targets,
            ((1, 1), (2, 2), (3, 0)),
        )

    def test_transient_swap_matrices_are_not_retained(self) -> None:
        recouplings = DimensionOnlyRecouplings()
        first = recouplings.adjoint_pair_swap((2, 1), (2, 1))
        second = recouplings.adjoint_pair_swap((2, 1), (2, 1))

        self.assertIsNot(first[1], second[1])
        self.assertFalse(hasattr(recouplings.adjoint_pair_swap, "cache_info"))
        self.assertFalse(hasattr(recouplings.line_across_adjoint, "cache_info"))
        self.assertFalse(hasattr(recouplings.adjacent_swap, "cache_info"))
        self.assertLessEqual(
            recouplings.fusion_paths.cache_info().currsize,
            _FUSION_PATH_CACHE_SIZE,
        )
        for name, cached in (
            ("adjoint decompositions", self.builder.couplings.decompose),
            ("elementary decompositions", recouplings.fusion.decompose),
            ("elementary frames", recouplings.frame),
            ("elementary blocks", recouplings.block),
            ("vertex reductions", self.builder.vertex_reductions.reduce),
        ):
            with self.subTest(cache=name):
                info = cached.cache_info()
                self.assertIsNotNone(info.maxsize)
                self.assertLessEqual(info.currsize, info.maxsize)

    def test_sparse_pair_application_matches_dense_matrix(self) -> None:
        recouplings = DimensionOnlyRecouplings()
        for source, target in (
            ((1, 1), (1, 1)),
            ((2, 1), (2, 1)),
            ((2, 1), (3, 2)),
        ):
            with self.subTest(source=source, target=target):
                paths, dense = recouplings.adjoint_pair_swap(source, target)
                values = np.arange(len(paths) * 3, dtype=float).reshape(len(paths), 3)
                sparse = recouplings.apply_adjoint_pair_swap(source, target, values)
                np.testing.assert_allclose(
                    sparse,
                    dense @ values,
                    rtol=0.0,
                    atol=2.0e-14,
                )

    def test_structural_adjoint_rule_accounts_for_every_dimension(self) -> None:
        for p in range(5):
            for q in range(5):
                source = (p, q)
                expected: list[tuple[int, int]] = [(p + 1, q + 1)]
                if q:
                    expected.append((p + 2, q - 1))
                if p:
                    expected.append((p - 1, q + 2))
                expected.extend([source] * (int(p > 0) + int(q > 0)))
                if q >= 2:
                    expected.append((p + 1, q - 2))
                if p >= 2:
                    expected.append((p - 2, q + 1))
                if p and q:
                    expected.append((p - 1, q - 1))
                decomposition = self.builder.couplings.decompose(source)
                self.assertEqual(
                    Counter(edge.target for edge in decomposition.couplings),
                    Counter(expected),
                )
                self.assertEqual(
                    sum(dimension(edge.target) for edge in decomposition.couplings),
                    8 * dimension(source),
                )

    def test_adjoint_seed_parities_are_structural(self) -> None:
        observed = {
            (edge.target, edge.multiplicity): edge.exchange_parity
            for edge in self.builder.couplings.decompose(ADJOINT).couplings
        }
        self.assertEqual(
            observed,
            {
                ((0, 0), 0): 1,
                ((0, 3), 0): -1,
                ((1, 1), 0): 1,
                ((1, 1), 1): -1,
                ((2, 2), 0): 1,
                ((3, 0), 0): -1,
            },
        )

    def test_young_seminormal_terminal_examples(self) -> None:
        recouplings = DimensionOnlyRecouplings()
        singlet = recouplings.block((1, 0), FUNDAMENTAL, FUNDAMENTAL, (0, 0))
        decuplet = recouplings.block((1, 0), FUNDAMENTAL, FUNDAMENTAL, (3, 0))
        octets = recouplings.block((1, 0), FUNDAMENTAL, FUNDAMENTAL, (1, 1))
        np.testing.assert_array_equal(singlet.matrix, [[-1.0]])
        np.testing.assert_array_equal(decuplet.matrix, [[1.0]])
        np.testing.assert_allclose(
            octets.matrix,
            [[0.5, np.sqrt(3.0) / 2.0], [np.sqrt(3.0) / 2.0, -0.5]],
            atol=2.0e-15,
        )

        # Conjugation, including the induced reordering of path labels.
        anti_octets = recouplings.block(
            (0, 1), ANTIFUNDAMENTAL, ANTIFUNDAMENTAL, (1, 1)
        )
        np.testing.assert_allclose(
            anti_octets.matrix,
            [[-0.5, np.sqrt(3.0) / 2.0], [np.sqrt(3.0) / 2.0, 0.5]],
            atol=2.0e-15,
        )

    def test_dimension_casimir_frames_fix_vertex_gauges(self) -> None:
        recouplings = self.builder.elementary_recouplings
        frame = recouplings.frame((2, 1), FUNDAMENTAL)
        np.testing.assert_allclose(frame.matrix.T @ frame.matrix, np.eye(3), atol=2e-15)

        generic_zero = self.builder.vertex_reductions.reduce((2, 1), (2, 1), 0)
        generic_one = self.builder.vertex_reductions.reduce((2, 1), (2, 1), 1)
        np.testing.assert_allclose(generic_zero.coefficients, frame.action, atol=2e-15)
        np.testing.assert_allclose(
            generic_one.coefficients,
            np.cross(frame.singlet, frame.action),
            atol=2e-15,
        )

        octet_frame = recouplings.frame(ADJOINT, FUNDAMENTAL)
        d_vertex = self.builder.vertex_reductions.reduce(ADJOINT, ADJOINT, 0)
        f_vertex = self.builder.vertex_reductions.reduce(ADJOINT, ADJOINT, 1)
        np.testing.assert_allclose(
            f_vertex.coefficients, -octet_frame.action, atol=2e-15
        )
        np.testing.assert_allclose(
            d_vertex.coefficients,
            np.cross(-octet_frame.action, octet_frame.singlet),
            atol=2e-15,
        )
        np.testing.assert_allclose(
            np.dot(d_vertex.coefficients, f_vertex.coefficients), 0.0, atol=2e-15
        )
        unique = self.builder.vertex_reductions.reduce((2, 1), (3, 2), 0)
        np.testing.assert_array_equal(unique.coefficients, [1.0])

    def test_all_blocks_through_cutoff_two_match_direct_contraction(self) -> None:
        actual = self.builder.build(2)
        self.assertEqual(len(actual), len(self.direct_blocks))
        for dimension_only, direct in zip(actual, self.direct_blocks, strict=True):
            self.assertEqual(dimension_only.left, direct.left)
            self.assertEqual(dimension_only.right, direct.right)
            self.assertEqual(dimension_only.paths, direct.paths)
            np.testing.assert_allclose(
                dimension_only.matrix,
                direct.matrix,
                rtol=0.0,
                atol=3.0e-14,
            )

    def test_all_checked_tables_through_cutoff_six(self) -> None:
        data = Path(__file__).resolve().parents[1] / "data"
        for cutoff in range(2, 7):
            with self.subTest(cutoff=cutoff):
                checked = read_table(data / f"su3_adjoint_swap_prefix_{cutoff}.tbl")
                actual = self.builder.build(cutoff)
                self.assertEqual(len(actual), len(checked))
                for dimension_only, expected in zip(actual, checked, strict=True):
                    self.assertEqual(dimension_only.left, expected.left)
                    self.assertEqual(dimension_only.right, expected.right)
                    self.assertEqual(dimension_only.paths, expected.paths)
                    np.testing.assert_allclose(
                        dimension_only.matrix,
                        expected.matrix,
                        rtol=0.0,
                        atol=4.0e-14,
                    )

    def test_all_prefix_sixteen_braid_relations(self) -> None:
        for left in self.builder.reachable_irreps(16):
            residuals = braid_residuals(self.builder, left)
            self.assertLess(max(residuals.values()), 4.0e-14)


if __name__ == "__main__":
    unittest.main()
