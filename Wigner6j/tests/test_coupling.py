from __future__ import annotations

from collections import Counter
import unittest

import numpy as np

from su3wigner import AdjointCouplings
from su3wigner.coupling import (
    _action_embedding,
    _adjoint_product_embedding,
    _product_action,
    _product_generator,
)


class CouplingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.couplings = AdjointCouplings()

    def test_singlet_times_adjoint(self) -> None:
        decomposition = self.couplings.decompose((0, 0))
        self.assertEqual(
            [(item.target, item.multiplicity) for item in decomposition.couplings],
            [((1, 1), 0)],
        )

    def test_implicit_product_action_matches_dense_kronecker_sum(self) -> None:
        source = self.couplings.irreps.get((2, 1))
        adjoint = self.couplings.irreps.get((1, 1))
        action = _product_action(source, adjoint, 1, 0)
        dense = _product_generator(source, adjoint, 1, 0)
        random = np.random.default_rng(1729)
        vectors = random.normal(size=(source.dim * adjoint.dim, 4))
        np.testing.assert_allclose(action @ vectors, dense @ vectors, atol=2.0e-14)

        inputs = (0, 3, 7, 11)
        outputs = (1, 5, 9, 13, 17)
        np.testing.assert_array_equal(
            action.restricted(outputs, inputs), dense[np.ix_(outputs, inputs)]
        )

    def test_adjoint_times_adjoint_and_parities(self) -> None:
        decomposition = self.couplings.decompose((1, 1))
        observed = {
            (item.target, item.multiplicity): item.exchange_parity
            for item in decomposition.couplings
        }
        expected = {
            ((0, 0), 0): 1,
            ((1, 1), 0): 1,
            ((1, 1), 1): -1,
            ((3, 0), 0): -1,
            ((0, 3), 0): -1,
            ((2, 2), 0): 1,
        }
        self.assertEqual(observed, expected)

    def test_generic_double_copy(self) -> None:
        decomposition = self.couplings.decompose((2, 2))
        copies = decomposition.for_target((2, 2))
        self.assertEqual([item.multiplicity for item in copies], [0, 1])
        source = self.couplings.irreps.get((2, 2))
        adjoint = self.couplings.irreps.get((1, 1))
        np.testing.assert_allclose(
            copies[0].embedding, _action_embedding(source, adjoint), atol=2.0e-13
        )

    def test_large_double_copy_replays_target_lowering_recipe(self) -> None:
        decomposition = self.couplings.decompose((5, 5))
        copies = decomposition.for_target((5, 5))
        self.assertEqual([item.multiplicity for item in copies], [0, 1])

    def test_octet_vertex_phases_are_product_anchored(self) -> None:
        adjoint = self.couplings.irreps.get((1, 1))
        copies = self.couplings.decompose((1, 1)).for_target((1, 1))
        np.testing.assert_allclose(
            copies[0].embedding,
            _adjoint_product_embedding(adjoint, 1),
            atol=2.0e-13,
        )
        np.testing.assert_allclose(
            copies[1].embedding,
            _adjoint_product_embedding(adjoint, -1),
            atol=2.0e-13,
        )

    def test_discovered_decomposition_through_small_generic_reps(self) -> None:
        # This product rule follows independently by writing 8=3 x 3bar - 1.
        # It is deliberately a test only; the generator discovers highest
        # weights from raising-operator kernels and does not call this formula.
        for p in range(4):
            for q in range(4):
                expected: list[tuple[int, int]] = [(p + 1, q + 1)]
                if q >= 1:
                    expected.append((p + 2, q - 1))
                if p >= 1:
                    expected.append((p - 1, q + 2))
                expected.extend([(p, q)] * ((p >= 1) + (q >= 1)))
                if q >= 2:
                    expected.append((p + 1, q - 2))
                if p >= 2:
                    expected.append((p - 2, q + 1))
                if p >= 1 and q >= 1:
                    expected.append((p - 1, q - 1))
                observed = [
                    item.target for item in self.couplings.decompose((p, q)).couplings
                ]
                self.assertEqual(Counter(observed), Counter(expected), (p, q))


if __name__ == "__main__":
    unittest.main()
