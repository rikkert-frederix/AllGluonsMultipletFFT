from __future__ import annotations

from dataclasses import FrozenInstanceError
import unittest

import numpy as np

from su3wigner.packed_embedding import PackedEmbedding


class PackedEmbeddingTests(unittest.TestCase):
    def setUp(self) -> None:
        # Deliberately interleave equal weights in both row and column order.
        self.product_weights = ((1, 0), (0, 1), (1, 0), (2, -1), (0, 1), (1, 0))
        self.target_weights = ((0, 1), (1, 0), (2, -1), (1, 0), (3, 3))
        self.dense = np.zeros((6, 5), dtype=np.float64)
        for row, row_weight in enumerate(self.product_weights):
            for column, column_weight in enumerate(self.target_weights):
                if row_weight == column_weight:
                    self.dense[row, column] = 10.0 * row + column + 0.25
        self.packed = PackedEmbedding(
            self.dense,
            self.product_weights,
            self.target_weights,
            particle_dimension=2,
        )

    def test_round_trip_and_noncontiguous_sector_access(self) -> None:
        np.testing.assert_array_equal(self.packed.to_dense(), self.dense)
        np.testing.assert_array_equal(
            self.packed.source_slab(1, 3), self.dense.reshape(3, 2, 5)[1:3]
        )
        np.testing.assert_array_equal(
            self.packed.target_columns(1, 4),
            self.dense.reshape(3, 2, 5)[:, :, 1:4],
        )
        np.testing.assert_array_equal(
            self.packed.target_column(3), self.dense[:, 3].reshape(3, 2)
        )

        rows, values = self.packed.column_entries(3)
        np.testing.assert_array_equal(rows, np.array([0, 2, 5]))
        np.testing.assert_array_equal(values, self.dense[rows, 3])
        self.assertFalse(rows.flags.writeable)
        self.assertFalse(values.flags.writeable)

        # A target weight absent from the product has an exactly empty column.
        empty_rows, empty_values = self.packed.column_entries(4)
        self.assertEqual(empty_rows.size, 0)
        self.assertEqual(empty_values.size, 0)
        np.testing.assert_array_equal(self.packed.target_column(4), 0.0)

    def test_rejects_any_nonzero_off_sector_without_a_tolerance(self) -> None:
        invalid = self.dense.copy()
        invalid[1, 1] = np.nextafter(0.0, 1.0)
        with self.assertRaisesRegex(ValueError, r"off-sector entry at \(1, 1\)"):
            PackedEmbedding(
                invalid,
                self.product_weights,
                self.target_weights,
                particle_dimension=2,
            )

    def test_storage_is_one_read_only_buffer_with_byte_accounting(self) -> None:
        expected_values = sum(
            sum(weight == row_weight for row_weight in self.product_weights)
            * sum(weight == column_weight for column_weight in self.target_weights)
            for weight in set(self.product_weights) & set(self.target_weights)
        )
        self.assertEqual(self.packed.value_count, expected_values)
        self.assertEqual(self.packed.packed_values.shape, (expected_values,))
        self.assertEqual(
            self.packed.nbytes,
            self.packed.data_nbytes + self.packed.descriptor_nbytes,
        )
        self.assertEqual(
            self.packed.data_nbytes, expected_values * self.dense.dtype.itemsize
        )
        self.assertLess(self.packed.data_nbytes, self.dense.nbytes)
        self.assertFalse(self.packed.packed_values.flags.writeable)
        with self.assertRaises(ValueError):
            self.packed.packed_values[0] = 0.0
        with self.assertRaises(FrozenInstanceError):
            self.packed._target_dimension = 99

    def test_validation_and_access_bounds(self) -> None:
        with self.assertRaisesRegex(ValueError, "not divisible"):
            PackedEmbedding(
                self.dense,
                self.product_weights,
                self.target_weights,
                particle_dimension=4,
            )
        with self.assertRaisesRegex(ValueError, "non-finite"):
            invalid = self.dense.copy()
            invalid[0, 1] = np.inf
            PackedEmbedding(
                invalid,
                self.product_weights,
                self.target_weights,
                particle_dimension=2,
            )
        with self.assertRaises(IndexError):
            self.packed.source_slab(-1, 1)
        with self.assertRaises(IndexError):
            self.packed.target_columns(2, 6)
        with self.assertRaises(IndexError):
            self.packed.target_column(5)


if __name__ == "__main__":
    unittest.main()
