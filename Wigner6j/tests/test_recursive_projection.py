from __future__ import annotations

import unittest

import numpy as np

from su3wigner import RecursiveReductionSwapTableBuilder
from su3wigner.coupling import _product_action
from su3wigner.recursive_reduction import (
    ANTIFUNDAMENTAL,
    FUNDAMENTAL,
    FundamentalCouplings,
)


class RecursiveProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.builder = RecursiveReductionSwapTableBuilder()

    def full_elementary_chain(self, first, second) -> np.ndarray:
        irreps = self.builder.fundamental_couplings.irreps
        source_dimension = irreps.get(first.source).dim
        first_dimension = irreps.get(first.particle).dim
        middle_dimension = irreps.get(first.target).dim
        second_dimension = irreps.get(second.particle).dim
        target_dimension = irreps.get(second.target).dim
        first_tensor = first.embedding.reshape(
            source_dimension, first_dimension, middle_dimension
        )
        second_tensor = second.embedding.reshape(
            middle_dimension, second_dimension * target_dimension
        )
        return (
            first_tensor.reshape(source_dimension * first_dimension, middle_dimension)
            @ second_tensor
        ).reshape(
            source_dimension,
            first_dimension,
            second_dimension,
            target_dimension,
        )

    def test_terminal_highest_projection_matches_full_trace(self) -> None:
        recouplings = self.builder.fundamental_recouplings
        source = FUNDAMENTAL
        target = (1, 1)
        options = recouplings._two_step_options(
            source, FUNDAMENTAL, FUNDAMENTAL, target
        )
        full = np.stack(
            [self.full_elementary_chain(first, second) for _, first, second in options]
        )
        exchanged = np.swapaxes(full, 2, 3)
        target_dimension = self.builder.couplings.irreps.get(target).dim
        expected = (
            exchanged.reshape(len(options), -1)
            @ full.reshape(len(options), -1).T
        ) / target_dimension

        actual = recouplings.block(
            source, FUNDAMENTAL, FUNDAMENTAL, target
        )
        np.testing.assert_allclose(actual.matrix, expected, rtol=0.0, atol=2.0e-14)

    def test_vertex_highest_projection_matches_full_trace(self) -> None:
        source = target = (1, 1)
        reductions = self.builder.vertex_reductions
        options = []
        for first in self.builder.fundamental_couplings.decompose(
            source, FUNDAMENTAL
        ).couplings:
            for second in self.builder.fundamental_couplings.decompose(
                first.target, ANTIFUNDAMENTAL
            ).for_target(target):
                options.append((first.target, first, second))
        options.sort(key=lambda item: item[0])
        chains = np.stack(
            [self.full_elementary_chain(first, second) for _, first, second in options]
        )
        chains_flat = chains.reshape(len(options), -1)
        source_dimension = self.builder.couplings.irreps.get(source).dim
        target_dimension = self.builder.couplings.irreps.get(target).dim
        splitting = reductions.splitting.reshape(8, 9)

        for adjoint_vertex in self.builder.couplings.decompose(source).for_target(
            target
        ):
            with self.subTest(multiplicity=adjoint_vertex.multiplicity):
                adjoint_tensor = adjoint_vertex.embedding.reshape(
                    source_dimension, 8, target_dimension
                )
                embedded = (
                    adjoint_tensor.transpose(0, 2, 1).reshape(-1, 8) @ splitting
                ).reshape(source_dimension, target_dimension, 3, 3).transpose(
                    0, 2, 3, 1
                )
                expected = (chains_flat @ embedded.reshape(-1)) / target_dimension
                actual = reductions.reduce(
                    source, target, adjoint_vertex.multiplicity
                )
                self.assertEqual(
                    actual.intermediates, tuple(item[0] for item in options)
                )
                np.testing.assert_allclose(
                    actual.coefficients, expected, rtol=0.0, atol=2.0e-14
                )

    def test_vertex_reduction_does_not_materialize_packed_adjoint(self) -> None:
        builder = RecursiveReductionSwapTableBuilder()
        source = target = (2, 2)
        vertex = builder.couplings.decompose(source).for_target(target)[0]
        self.assertIsNotNone(vertex._packed_embedding)
        self.assertIsNone(vertex._embedding_cache)

        builder.vertex_reductions.reduce(source, target, vertex.multiplicity)

        self.assertIsNone(vertex._embedding_cache)

    def test_deep_fundamental_replay_obeys_reference_lowerings(self) -> None:
        couplings = FundamentalCouplings()
        source_label = (5, 5)
        source = couplings.irreps.get(source_label)
        particle = couplings.irreps.get(FUNDAMENTAL)
        decomposition = couplings.decompose(source_label, FUNDAMENTAL)

        worst = 0.0
        for coupling in decomposition.couplings:
            target = couplings.irreps.get(coupling.target)
            for row, column in ((1, 0), (2, 1)):
                residual = (
                    _product_action(source, particle, row, column)
                    @ coupling.embedding
                    - coupling.embedding @ target.generator(row, column)
                )
                worst = max(worst, float(np.max(np.abs(residual))))
        self.assertLess(worst, 1.0e-12)


if __name__ == "__main__":
    unittest.main()
