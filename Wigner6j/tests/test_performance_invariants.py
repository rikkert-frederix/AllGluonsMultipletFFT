from __future__ import annotations

from dataclasses import asdict, fields, replace
from inspect import signature
import unittest

import numpy as np

from su3wigner import Irrep, IrrepFactory
from su3wigner.coupling import AdjointCouplings, _TensorProductAction
from su3wigner.packed_embedding import PackedEmbedding
from su3wigner.representations import _harmonic_basis


class PerformanceInvariantTests(unittest.TestCase):
    def test_irrep_public_dataclass_contract_is_preserved(self) -> None:
        generators = tuple(
            tuple(np.zeros((1, 1)) for _column in range(3)) for _row in range(3)
        )
        states = ((0, 0, 0, 0, 0, 0),)
        embedding = np.array([[1.0]])
        irrep = Irrep((0, 0), generators, ((0, 0),), (), states, embedding)

        public_fields = [
            "label",
            "generators",
            "weights",
            "lowering_recipe",
            "ambient_states",
            "ambient_embedding",
        ]
        self.assertEqual([item.name for item in fields(Irrep)], public_fields)
        self.assertEqual(Irrep.__match_args__, tuple(public_fields))
        match irrep:
            case Irrep(
                label, matched_generators, weights, recipe, matched_states, dense
            ):
                matched = (
                    label,
                    matched_generators,
                    weights,
                    recipe,
                    matched_states,
                    dense,
                )
            case _:
                self.fail("public positional Irrep pattern did not match")
        self.assertEqual(matched[0], (0, 0))
        self.assertIs(matched[1], generators)
        self.assertEqual(matched[2:5], (((0, 0),), (), states))
        self.assertIs(matched[5], embedding)

        serialized = asdict(irrep)
        self.assertEqual(list(serialized), public_fields)
        np.testing.assert_array_equal(serialized["ambient_embedding"], embedding)
        self.assertFalse(any(name.startswith("_") for name in serialized))

        changed = replace(irrep, label=(1, 0))
        self.assertEqual(changed.label, (1, 0))
        self.assertIs(changed.ambient_states, states)
        self.assertIs(changed.ambient_embedding, embedding)
        different = Irrep(
            (0, 0),
            generators,
            ((0, 0),),
            (),
            states,
            np.array([[2.0]]),
        )
        self.assertNotEqual(irrep, different)
        with self.assertRaises(TypeError):
            hash(irrep)

    def test_nonadjoint_ambient_embedding_is_lazy_and_public(self) -> None:
        irrep = IrrepFactory().get((2, 1))
        self.assertIsNone(irrep._ambient_embedding_cache)

        embedding = irrep.ambient_embedding
        self.assertIs(irrep.ambient_embedding, embedding)
        self.assertEqual(embedding.shape, (len(irrep.ambient_states), irrep.dim))
        np.testing.assert_allclose(
            embedding.T @ embedding,
            np.eye(irrep.dim),
            atol=2.0e-13,
        )

    def test_factory_ambient_transform_is_packed_and_bitwise_compatible(self) -> None:
        irrep = IrrepFactory().get((3, 3))
        packed = irrep._ambient_transform
        self.assertIsInstance(packed, PackedEmbedding)
        assert isinstance(packed, PackedEmbedding)
        dense_transform = packed.to_dense()
        self.assertLess(packed.nbytes, dense_transform.nbytes // 8)
        self.assertIsNone(irrep._ambient_embedding_cache)

        harmonic, _weights = _harmonic_basis(*irrep.label)
        expected_embedding = harmonic @ dense_transform
        np.testing.assert_array_equal(irrep.ambient_embedding, expected_embedding)

        # Irrep is public, so retain its original dense constructor exactly.
        self.assertEqual(
            list(signature(Irrep).parameters),
            [
                "label",
                "generators",
                "weights",
                "lowering_recipe",
                "ambient_states",
                "ambient_embedding",
            ],
        )
        ambient_states = irrep.ambient_states
        dense_irrep = Irrep(
            label=irrep.label,
            generators=irrep.generators,
            weights=irrep.weights,
            lowering_recipe=irrep.lowering_recipe,
            ambient_states=ambient_states,
            ambient_embedding=expected_embedding,
        )
        self.assertIs(dense_irrep.ambient_states, ambient_states)
        self.assertIs(dense_irrep.ambient_embedding, expected_embedding)
        np.testing.assert_array_equal(
            irrep.ambient_embedding, dense_irrep.ambient_embedding
        )
        self.assertIs(irrep.ambient_embedding, irrep.ambient_embedding)

    def test_multiply_connected_weight_replay_remains_intertwining(self) -> None:
        couplings = AdjointCouplings()
        coupling = couplings.decompose((4, 4)).for_target((3, 6))[0]
        source = couplings.irreps.get((4, 4))
        adjoint = couplings.irreps.get((1, 1))
        target = couplings.irreps.get((3, 6))

        worst = 0.0
        for row, col in ((0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)):
            product = _TensorProductAction(
                source.generator(row, col), adjoint.generator(row, col)
            )
            residual = (
                product @ coupling.embedding
                - coupling.embedding @ target.generator(row, col)
            )
            worst = max(worst, float(np.max(np.abs(residual))))
        self.assertLess(worst, 2.0e-12)


if __name__ == "__main__":
    unittest.main()
