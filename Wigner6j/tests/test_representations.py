from __future__ import annotations

import unittest

import numpy as np

from su3wigner import IrrepFactory, dimension


class RepresentationTests(unittest.TestCase):
    def test_dimension_and_weights(self) -> None:
        factory = IrrepFactory()
        for label in ((0, 0), (1, 0), (0, 1), (1, 1), (3, 0), (2, 2)):
            irrep = factory.get(label)
            self.assertEqual(irrep.dim, dimension(label))
            self.assertEqual(irrep.weights[0], label)

    def test_simple_root_commutator(self) -> None:
        irrep = IrrepFactory().get((2, 2))
        e12 = irrep.generator(0, 1)
        e21 = irrep.generator(1, 0)
        h1 = irrep.generator(0, 0) - irrep.generator(1, 1)
        np.testing.assert_allclose(e12 @ e21 - e21 @ e12, h1, atol=2.0e-12)


if __name__ == "__main__":
    unittest.main()
