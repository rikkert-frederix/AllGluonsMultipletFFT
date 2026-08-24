from __future__ import annotations

from collections import Counter
from dataclasses import replace
from pathlib import Path
import time
import unittest
from unittest import mock

import numpy as np

from su3wigner import DirectSpechtSwapOracle, DirectSpechtSwapTableBuilder, read_table
from su3wigner.coupling import AdjointCouplings
from su3wigner.direct_specht import (
    _path_vectors,
    _projector_kernel,
    _SkewRepresentation,
)
from su3wigner.recursive_reduction import FundamentalCouplings
from su3wigner.representations import IrrepFactory


class DirectSpechtOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.oracle = DirectSpechtSwapOracle()
        filename = (
            Path(__file__).resolve().parents[1]
            / "data"
            / "su3_adjoint_swap_prefix_2.tbl"
        )
        cls.reference = {
            (block.left, block.right): block for block in read_table(filename)
        }

    def assert_matches_reference(
        self, left: tuple[int, int], right: tuple[int, int]
    ) -> None:
        actual = self.oracle.block(left, right)
        expected = self.reference[(left, right)]
        self.assertEqual(actual.paths, expected.paths)
        np.testing.assert_allclose(
            actual.matrix, expected.matrix, rtol=0.0, atol=4.0e-14
        )

    def test_complete_seed_parities(self) -> None:
        blocks = {block.right: block for block in self.oracle.blocks_for_left((0, 0))}
        self.assertEqual(set(blocks), {(0, 0), (0, 3), (1, 1), (2, 2), (3, 0)})
        np.testing.assert_array_equal(blocks[(0, 0)].matrix, [[1.0]])
        np.testing.assert_array_equal(blocks[(2, 2)].matrix, [[1.0]])
        np.testing.assert_array_equal(blocks[(3, 0)].matrix, [[-1.0]])
        np.testing.assert_array_equal(blocks[(0, 3)].matrix, [[-1.0]])
        np.testing.assert_array_equal(blocks[(1, 1)].matrix, [[1.0, 0.0], [0.0, -1.0]])

    def test_multiplicity_free_radical_block(self) -> None:
        self.assert_matches_reference((1, 1), (1, 4))
        diagnostics = self.oracle.diagnostics((1, 1), (1, 4))
        self.assertEqual(diagnostics.final_partition, (5, 4, 0))
        self.assertEqual(diagnostics.skew_dimension, 14)
        self.assertEqual(diagnostics.product_rank, 2)
        self.assertEqual(diagnostics.prefix_ranks, (((0, 3), 1), ((2, 2), 1)))
        expected = np.array([[0.5, np.sqrt(3.0) / 2.0], [np.sqrt(3.0) / 2.0, -0.5]])
        np.testing.assert_allclose(
            self.oracle.block((1, 1), (1, 4)).matrix,
            expected,
            rtol=0.0,
            atol=3.0e-15,
        )

    def test_generic_action_and_complement_block(self) -> None:
        self.assert_matches_reference((2, 2), (0, 3))
        diagnostics = self.oracle.diagnostics((2, 2), (0, 3))
        self.assertEqual(diagnostics.final_partition, (5, 5, 2))
        self.assertEqual(diagnostics.skew_dimension, 45)
        self.assertEqual(diagnostics.product_rank, 5)
        self.assertEqual(
            diagnostics.prefix_ranks,
            (((0, 3), 1), ((1, 1), 1), ((1, 4), 1), ((2, 2), 2)),
        )

    def test_skew_generators_projectors_and_literal_swap(self) -> None:
        representation = _SkewRepresentation((4, 2, 0), (5, 5, 2))
        identity = np.eye(representation.dimension)
        for generator in representation.generators:
            np.testing.assert_allclose(generator, generator.T, rtol=0.0, atol=2e-15)
            np.testing.assert_allclose(
                generator @ generator, identity, rtol=0.0, atol=3e-15
            )
        for index in range(4):
            left = (
                representation.generators[index]
                @ representation.generators[index + 1]
                @ representation.generators[index]
            )
            right = (
                representation.generators[index + 1]
                @ representation.generators[index]
                @ representation.generators[index + 1]
            )
            np.testing.assert_allclose(left, right, rtol=0.0, atol=4e-15)

        first = representation.adjoint(1)
        second = representation.adjoint(4)
        for projector in (first, second):
            np.testing.assert_allclose(projector, projector.T, rtol=0.0, atol=2e-15)
            np.testing.assert_allclose(
                projector @ projector, projector, rtol=0.0, atol=3e-15
            )
        np.testing.assert_allclose(first @ second, second @ first, rtol=0.0, atol=3e-15)

        swap = representation.block_swap()
        np.testing.assert_allclose(swap, swap.T, rtol=0.0, atol=5e-15)
        np.testing.assert_allclose(swap @ swap, identity, rtol=0.0, atol=6e-15)

    def test_prefix_masks_match_dense_projectors(self) -> None:
        representation = _SkewRepresentation((4, 2, 0), (5, 5, 2))
        partition = (5, 3, 1)
        mask = representation.prefix_mask(partition)
        np.testing.assert_array_equal(
            representation.prefix_projector(partition), np.diag(mask)
        )
        self.assertFalse(mask.flags.writeable)

    def test_oracle_never_calls_recursive_dimension_or_cg_backends(self) -> None:
        import su3wigner.direct_specht as module

        source = Path(module.__file__).read_text(encoding="utf-8")
        self.assertNotIn("dimension_only", source)
        self.assertNotIn("recursive_reduction", source)
        self.assertNotIn("AdjointCouplings", source)
        self.assertNotIn("FundamentalCouplings", source)

        forbidden = AssertionError("direct Specht oracle entered a CG backend")
        with (
            mock.patch.object(IrrepFactory, "__init__", side_effect=forbidden),
            mock.patch.object(AdjointCouplings, "__init__", side_effect=forbidden),
            mock.patch.object(FundamentalCouplings, "__init__", side_effect=forbidden),
        ):
            guarded = DirectSpechtSwapOracle()
            guarded.blocks_for_left((0, 0))
            guarded.block((1, 1), (1, 4))
            guarded.block((2, 2), (0, 3))

    def test_formerly_unsupported_rank_four_block_matches(self) -> None:
        self.assert_matches_reference((1, 1), (1, 1))

    def test_non_singlet_block_family_is_complete(self) -> None:
        actual = self.oracle.blocks_for_left((1, 1))
        expected = [block for block in self.reference.values() if block.left == (1, 1)]
        self.assertEqual([(block.left, block.right) for block in actual],
                         [(block.left, block.right) for block in expected])

    def test_unexpected_middle_pattern_is_rejected(self) -> None:
        kernel = _projector_kernel((2, 2), (0, 3))
        paths = list(kernel.paths)
        repeated = next(
            index
            for index, path in enumerate(paths)
            if path.middle == (2, 2) and path.left_multiplicity == 1
        )
        paths[repeated] = replace(paths[repeated], right_multiplicity=1)
        with self.assertRaisesRegex(ArithmeticError, "unexpected multiplicity pattern"):
            _path_vectors(replace(kernel, paths=tuple(paths)))

    def test_prefix_six_raw_projectors_cover_every_block(self) -> None:
        filename = (
            Path(__file__).resolve().parents[1]
            / "data"
            / "su3_adjoint_swap_prefix_6.tbl"
        )
        blocks = read_table(filename)
        dimensions: list[int] = []
        for block in blocks:
            diagnostics = self.oracle.diagnostics(block.left, block.right)
            self.assertEqual(diagnostics.product_rank, block.size)
            eigenvalues = np.linalg.eigvalsh(block.matrix)
            self.assertEqual(
                diagnostics.swap_signature,
                (
                    int(np.count_nonzero(eigenvalues > 0.5)),
                    int(np.count_nonzero(eigenvalues < -0.5)),
                ),
            )
            self.assertEqual(
                diagnostics.prefix_ranks,
                tuple(sorted(Counter(path.middle for path in block.paths).items())),
            )
            dimensions.append(diagnostics.skew_dimension)
        self.assertEqual(len(dimensions), 354)
        self.assertEqual(max(dimensions), 90)
        self.assertTrue(
            all(
                self.oracle.supported(block.left, block.right)
                for block in blocks
            )
        )

    def test_full_prefix_six_signed_table(self) -> None:
        filename = (
            Path(__file__).resolve().parents[1]
            / "data"
            / "su3_adjoint_swap_prefix_6.tbl"
        )
        expected = read_table(filename)
        actual = DirectSpechtSwapTableBuilder().build(6)
        self.assertEqual(len(actual), len(expected))
        for left, right in zip(actual, expected, strict=True):
            self.assertEqual((left.left, left.right, left.paths),
                             (right.left, right.right, right.paths))
            np.testing.assert_allclose(
                left.matrix, right.matrix, rtol=0.0, atol=4.0e-14
            )

    def test_support_caches_are_bounded_and_tableau_recursion_stays_local(
        self,
    ) -> None:
        import su3wigner.direct_specht as module

        expected_bounds = (
            (module._adjoint_edges, module._EDGE_CACHE_SIZE),
            (module._paths, module._PATH_CACHE_SIZE),
            (module._skew_tableaux, module._TABLEAU_CACHE_SIZE),
            (module._local_edge_frame, module._LOCAL_FRAME_CACHE_SIZE),
            (self.oracle.block, module._SIGNED_BLOCK_CACHE_SIZE),
            (self.oracle.diagnostics, module._DIAGNOSTIC_CACHE_SIZE),
        )
        for cached, expected in expected_bounds:
            information = cached.cache_info()
            self.assertEqual(information.maxsize, expected)
            self.assertLessEqual(information.currsize, expected)

        module._skew_tableaux.cache_clear()
        self.assertEqual(module.skew_dimension((2, 2), (2, 2)), 90)
        information = module._skew_tableaux.cache_info()
        self.assertGreater(information.hits, 0)
        self.assertLess(information.currsize, module._TABLEAU_CACHE_SIZE)

    def test_audited_blocks_are_small_and_fast(self) -> None:
        fresh = DirectSpechtSwapOracle()
        started = time.perf_counter()
        fresh.blocks_for_left((0, 0))
        fresh.block((1, 1), (1, 4))
        fresh.block((2, 2), (0, 3))
        elapsed = time.perf_counter() - started
        # This is deliberately generous enough for unoptimized CI machines;
        # the local reference run is below one tenth of a second.
        self.assertLess(elapsed, 5.0)


if __name__ == "__main__":
    unittest.main()
