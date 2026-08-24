from __future__ import annotations

from unittest import mock
import unittest

import numpy as np

from su3wigner import AdjointCouplings, SwapTableBuilder
from su3wigner import recoupling


class ExplicitRecouplingTests(unittest.TestCase):
    def test_dense_decode_lru_is_bounded_and_private(self) -> None:
        couplings = AdjointCouplings()
        vertices = couplings.decompose((1, 1)).couplings[:2]
        dense_sizes = [
            coupling._packed_embedding.shape[0]
            * coupling._packed_embedding.shape[1]
            * coupling._packed_embedding.dtype.itemsize
            for coupling in vertices
            if coupling._packed_embedding is not None
        ]
        self.assertEqual(len(dense_sizes), 2)

        with mock.patch.object(
            recoupling, "_DENSE_DECODE_CACHE_BYTES", max(dense_sizes)
        ):
            builder = SwapTableBuilder(couplings)
        first = builder._decoded_coupling_tensor(vertices[0])
        self.assertIs(first, builder._decoded_coupling_tensor(vertices[0]))
        second = builder._decoded_coupling_tensor(vertices[1])

        self.assertIsNotNone(first)
        self.assertIsNotNone(second)
        self.assertNotIn(id(vertices[0]), builder._dense_decode_cache)
        self.assertIn(id(vertices[1]), builder._dense_decode_cache)
        self.assertLessEqual(
            builder._dense_decode_bytes,
            builder._dense_decode_budget,
        )
        self.assertTrue(
            all(coupling._embedding_cache is None for coupling in vertices)
        )

    def test_index_equation_residual_is_relative_to_chain_norm(self) -> None:
        gram = np.eye(2)
        right_dimension = 648
        relative_error = 7.1e-9
        residual_squared = np.longdouble(
            relative_error**2 * right_dimension * len(gram)
        )

        basis_error, equation_error = SwapTableBuilder._index_residuals(
            gram,
            residual_squared,
            right_dimension,
        )

        self.assertEqual(basis_error, 0.0)
        self.assertAlmostEqual(equation_error, relative_error)
        self.assertLess(equation_error, 2.0e-8)

    def test_bounded_slabs_match_single_slab(self) -> None:
        dense_couplings = AdjointCouplings()
        dense_builder = SwapTableBuilder(dense_couplings)
        dense_path_data = dense_builder._path_couplings_for_left((1, 1))
        for entries in dense_path_data.values():
            for _path, first, second in entries:
                first.embedding
                second.embedding
        with mock.patch.object(recoupling, "_CHAIN_SLAB_BYTES", 2**30):
            expected = dense_builder.blocks_for_left((1, 1))

        packed_couplings = AdjointCouplings()
        with mock.patch.object(recoupling, "_DENSE_DECODE_CACHE_BYTES", 0):
            packed_builder = SwapTableBuilder(packed_couplings)
        packed_path_data = packed_builder._path_couplings_for_left((1, 1))
        packed_vertices = [
            coupling
            for entries in packed_path_data.values()
            for _path, first, second in entries
            for coupling in (first, second)
        ]
        self.assertTrue(packed_vertices)
        self.assertTrue(
            all(coupling._embedding_cache is None for coupling in packed_vertices)
        )
        with mock.patch.object(recoupling, "_CHAIN_SLAB_BYTES", 4096):
            actual = packed_builder.blocks_for_left((1, 1))

        self.assertEqual(packed_builder._dense_decode_budget, 0)
        self.assertFalse(packed_builder._dense_decode_cache)
        self.assertTrue(
            all(coupling._embedding_cache is None for coupling in packed_vertices)
        )

        self.assertEqual(len(actual), len(expected))
        for bounded, single in zip(actual, expected, strict=True):
            self.assertEqual(bounded.left, single.left)
            self.assertEqual(bounded.right, single.right)
            self.assertEqual(bounded.paths, single.paths)
            np.testing.assert_allclose(
                bounded.matrix,
                single.matrix,
                rtol=0.0,
                atol=5.0e-15,
            )


if __name__ == "__main__":
    unittest.main()
