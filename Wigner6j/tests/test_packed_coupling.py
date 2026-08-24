from __future__ import annotations

from dataclasses import asdict, fields, replace
import unittest

import numpy as np

from su3wigner import AdjointCouplings, Coupling


class PackedCouplingTests(unittest.TestCase):
    def test_public_dataclass_contract_is_preserved(self) -> None:
        dense = np.array([1.0])
        coupling = Coupling((0, 0), (1, 0), 3, dense, -1)

        self.assertEqual(
            [item.name for item in fields(Coupling)],
            [
                "source",
                "target",
                "multiplicity",
                "embedding",
                "exchange_parity",
            ],
        )
        self.assertEqual(
            Coupling.__match_args__,
            (
                "source",
                "target",
                "multiplicity",
                "embedding",
                "exchange_parity",
            ),
        )
        match coupling:
            case Coupling(source, target, multiplicity, embedding, parity):
                matched = (source, target, multiplicity, embedding, parity)
            case _:
                self.fail("public positional Coupling pattern did not match")
        self.assertEqual(matched[:3], ((0, 0), (1, 0), 3))
        self.assertIs(matched[3], dense)
        self.assertEqual(matched[4], -1)

        serialized = asdict(coupling)
        self.assertEqual(list(serialized), [item.name for item in fields(Coupling)])
        np.testing.assert_array_equal(serialized["embedding"], dense)
        self.assertFalse(any(name.startswith("_") for name in serialized))

        changed = replace(coupling, multiplicity=4)
        self.assertEqual(changed.multiplicity, 4)
        self.assertIs(changed.embedding, dense)
        self.assertNotEqual(
            coupling,
            Coupling((0, 0), (1, 0), 3, np.array([2.0]), -1),
        )
        with self.assertRaises(TypeError):
            hash(coupling)

    def test_public_dense_constructor_and_property_are_preserved(self) -> None:
        dense = np.arange(16.0).reshape(8, 2)
        coupling = Coupling(
            source=(0, 0),
            target=(1, 0),
            multiplicity=3,
            embedding=dense,
            exchange_parity=-1,
        )
        self.assertIs(coupling.embedding, dense)
        self.assertEqual(coupling.source, (0, 0))
        self.assertEqual(coupling.target, (1, 0))
        self.assertEqual(coupling.multiplicity, 3)
        self.assertEqual(coupling.exchange_parity, -1)

    def test_real_coupling_is_packed_and_accessors_round_trip(self) -> None:
        couplings = AdjointCouplings()
        coupling = couplings.decompose((2, 2)).for_target((2, 2))[0]
        packed = coupling._packed_embedding
        self.assertIsNotNone(packed)
        assert packed is not None
        self.assertIsNone(coupling._embedding_cache)

        expected = packed.to_dense()
        self.assertLess(packed.nbytes, expected.nbytes // 10)
        source_dimension = packed.source_dimension
        target_dimension = packed.target_dimension
        tensor = expected.reshape(source_dimension, 8, target_dimension)
        np.testing.assert_array_equal(coupling._source_slab(2, 7), tensor[2:7])
        np.testing.assert_array_equal(
            coupling._target_columns(3, 9), tensor[:, :, 3:9]
        )
        np.testing.assert_array_equal(coupling._target_column(4), tensor[:, :, 4])
        rows, values = coupling._column_entries(4)
        np.testing.assert_array_equal(values, expected[rows, 4])
        self.assertIsNone(coupling._embedding_cache)

        column_with_sector_zeros = next(
            column
            for column in range(target_dimension)
            if (
                (stored := packed.column_entries(column)[1]).size
                and np.any(stored == 0)
                and np.any(stored != 0)
            )
        )
        expected_rows = np.flatnonzero(expected[:, column_with_sector_zeros])
        packed_rows, packed_values = coupling._column_entries(
            column_with_sector_zeros
        )
        np.testing.assert_array_equal(packed_rows, expected_rows)
        np.testing.assert_array_equal(
            packed_values, expected[expected_rows, column_with_sector_zeros]
        )

        dense = coupling.embedding
        self.assertIs(coupling.embedding, dense)
        np.testing.assert_array_equal(dense, expected)
        dense_rows, dense_values = coupling._column_entries(
            column_with_sector_zeros
        )
        np.testing.assert_array_equal(dense_rows, packed_rows)
        np.testing.assert_array_equal(dense_values, packed_values)

        # Once callers request the mutable public array, private accessors must
        # observe that same cached ndarray rather than stale packed storage.
        dense[0, 4] = 17.0
        self.assertEqual(coupling._target_column(4)[0, 0], 17.0)

    def test_decomposition_is_complete_ordered_and_packed(self) -> None:
        couplings = AdjointCouplings()
        source = couplings.irreps.get((2, 2))
        decomposition = couplings.decompose(source.label)

        self.assertEqual(
            [(item.target, item.multiplicity) for item in decomposition.couplings],
            [
                ((4, 1), 0),
                ((3, 3), 0),
                ((3, 0), 0),
                ((2, 2), 0),
                ((2, 2), 1),
                ((1, 4), 0),
                ((1, 1), 0),
                ((0, 3), 0),
            ],
        )
        self.assertEqual(
            sum(
                couplings.irreps.get(item.target).dim
                for item in decomposition.couplings
            ),
            source.dim * 8,
        )
        self.assertTrue(
            all(item._packed_embedding is not None for item in decomposition.couplings)
        )
        self.assertTrue(
            all(item._embedding_cache is None for item in decomposition.couplings)
        )


if __name__ == "__main__":
    unittest.main()
