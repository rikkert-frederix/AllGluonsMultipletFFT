from __future__ import annotations

from dataclasses import asdict, fields, replace
import unittest

import numpy as np

from su3wigner import RecursiveReductionSwapTableBuilder
from su3wigner.recursive_reduction import (
    ADJOINT_SPLIT_SEQUENCE,
    ANTIFUNDAMENTAL,
    FUNDAMENTAL,
    ElementaryCoupling,
    FundamentalCouplings,
)


class PackedElementaryCouplingTests(unittest.TestCase):
    def test_public_dataclass_contract_is_preserved(self) -> None:
        dense = np.array([1.0])
        coupling = ElementaryCoupling((0, 0), FUNDAMENTAL, FUNDAMENTAL, dense)

        self.assertEqual(
            [item.name for item in fields(ElementaryCoupling)],
            ["source", "particle", "target", "embedding"],
        )
        self.assertEqual(
            ElementaryCoupling.__match_args__,
            ("source", "particle", "target", "embedding"),
        )
        match coupling:
            case ElementaryCoupling(source, particle, target, embedding):
                matched = (source, particle, target, embedding)
            case _:
                self.fail("public positional ElementaryCoupling pattern did not match")
        self.assertEqual(matched[:3], ((0, 0), FUNDAMENTAL, FUNDAMENTAL))
        self.assertIs(matched[3], dense)

        serialized = asdict(coupling)
        self.assertEqual(
            list(serialized), [item.name for item in fields(ElementaryCoupling)]
        )
        np.testing.assert_array_equal(serialized["embedding"], dense)
        self.assertFalse(any(name.startswith("_") for name in serialized))

        changed = replace(coupling, target=(2, 0))
        self.assertEqual(changed.target, (2, 0))
        self.assertIs(changed.embedding, dense)
        self.assertNotEqual(
            coupling,
            ElementaryCoupling(
                (0, 0), FUNDAMENTAL, FUNDAMENTAL, np.array([2.0])
            ),
        )
        with self.assertRaises(TypeError):
            hash(coupling)

    def test_public_dense_constructor_and_property_are_preserved(self) -> None:
        dense = np.arange(9.0).reshape(3, 3)
        coupling = ElementaryCoupling(
            source=(0, 0),
            particle=FUNDAMENTAL,
            target=FUNDAMENTAL,
            embedding=dense,
        )
        self.assertIs(coupling.embedding, dense)
        self.assertEqual(coupling.source, (0, 0))
        self.assertEqual(coupling.particle, FUNDAMENTAL)
        self.assertEqual(coupling.target, FUNDAMENTAL)

    def test_real_elementary_coupling_round_trips_without_promotion(self) -> None:
        couplings = FundamentalCouplings()
        decomposition = couplings.decompose((3, 2), FUNDAMENTAL)
        packed_bytes = sum(
            coupling._packed_embedding.nbytes
            for coupling in decomposition.couplings
            if coupling._packed_embedding is not None
        )
        dense_bytes = sum(
            packed.shape[0] * packed.shape[1] * packed.dtype.itemsize
            for coupling in decomposition.couplings
            if (packed := coupling._packed_embedding) is not None
        )
        self.assertLess(packed_bytes, dense_bytes // 10)

        coupling = max(
            decomposition.couplings,
            key=lambda item: item._target_dimension,
        )
        packed = coupling._packed_embedding
        self.assertIsNotNone(packed)
        assert packed is not None
        self.assertIsNone(coupling._embedding_cache)
        expected = packed.to_dense()
        tensor = expected.reshape(packed.source_dimension, 3, packed.target_dimension)
        np.testing.assert_array_equal(coupling._source_slab(2, 7), tensor[2:7])
        np.testing.assert_array_equal(
            coupling._target_columns(3, 9), tensor[:, :, 3:9]
        )
        np.testing.assert_array_equal(coupling._target_column(4), tensor[:, :, 4])
        rows, values = coupling._column_entries(4)
        np.testing.assert_array_equal(values, expected[rows, 4])
        self.assertIsNone(coupling._embedding_cache)

        dense = coupling.embedding
        self.assertIs(coupling.embedding, dense)
        np.testing.assert_array_equal(dense, expected)

    def test_recursive_consumers_leave_elementary_cache_compact(self) -> None:
        builder = RecursiveReductionSwapTableBuilder()
        couplings = builder.fundamental_couplings
        recouplings = builder.fundamental_recouplings
        retained: list[ElementaryCoupling] = list(
            couplings.decompose(FUNDAMENTAL, ANTIFUNDAMENTAL).couplings
        )

        terminal_options = recouplings._two_step_options(
            (2, 2), ANTIFUNDAMENTAL, FUNDAMENTAL, (2, 2)
        )
        retained.extend(
            coupling
            for _middle, first, second in terminal_options
            for coupling in (first, second)
        )
        recouplings.block(
            (2, 2), ANTIFUNDAMENTAL, FUNDAMENTAL, (2, 2)
        )

        for first in couplings.decompose((2, 2), FUNDAMENTAL).couplings:
            retained.append(first)
            retained.extend(
                couplings.decompose(
                    first.target, ANTIFUNDAMENTAL
                ).for_target((2, 2))
            )
        builder.vertex_reductions.reduce((2, 2), (2, 2), 0)

        paths = recouplings.fusion_paths(
            (0, 0), ADJOINT_SPLIT_SEQUENCE, (0, 0)
        )
        self.assertTrue(paths)
        audit_path = paths[0]
        for edge, particle in enumerate(ADJOINT_SPLIT_SEQUENCE):
            options = couplings.decompose(
                audit_path[edge], particle
            ).for_target(audit_path[edge + 1])
            self.assertEqual(len(options), 1)
            retained.append(options[0])
        recouplings.path_tensor(audit_path, ADJOINT_SPLIT_SEQUENCE)

        unique = {id(coupling): coupling for coupling in retained}.values()
        for coupling in unique:
            with self.subTest(
                source=coupling.source,
                particle=coupling.particle,
                target=coupling.target,
            ):
                self.assertIsNotNone(coupling._packed_embedding)
                self.assertIsNone(coupling._embedding_cache)
        self.assertEqual(couplings.decompose.cache_parameters()["maxsize"], 256)


if __name__ == "__main__":
    unittest.main()
