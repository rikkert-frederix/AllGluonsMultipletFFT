from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from su3wigner import (
    SwapBlock,
    SwapTableBuilder,
    braid_residuals,
    read_table,
    write_table,
)


class RecouplingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.builder = SwapTableBuilder()

    def test_eq8_seed_signs(self) -> None:
        blocks = {block.right: block for block in self.builder.blocks_for_left((0, 0))}
        np.testing.assert_array_equal(blocks[(0, 0)].matrix, [[1.0]])
        np.testing.assert_array_equal(blocks[(2, 2)].matrix, [[1.0]])
        np.testing.assert_array_equal(blocks[(3, 0)].matrix, [[-1.0]])
        np.testing.assert_array_equal(blocks[(0, 3)].matrix, [[-1.0]])
        np.testing.assert_array_equal(
            blocks[(1, 1)].matrix, [[1.0, 0.0], [0.0, -1.0]]
        )

    def test_nontrivial_block_is_symmetric_involution(self) -> None:
        blocks = {block.right: block for block in self.builder.blocks_for_left((1, 1))}
        block = blocks[(2, 2)]
        np.testing.assert_allclose(block.matrix, block.matrix.T, atol=2.0e-14)
        np.testing.assert_allclose(
            block.matrix @ block.matrix, np.eye(block.size), atol=2.0e-14
        )

    def test_nonfinite_block_is_rejected(self) -> None:
        original = self.builder.blocks_for_left((0, 0))[0]
        for value in (np.nan, np.inf, -np.inf):
            with self.subTest(value=value):
                block = SwapBlock(
                    left=original.left,
                    right=original.right,
                    paths=original.paths,
                    matrix=np.array([[value]]),
                )
                with self.assertRaisesRegex(ArithmeticError, "non-finite"):
                    self.builder.validate_block(block)

    def test_low_nontrivial_radical_block(self) -> None:
        blocks = {block.right: block for block in self.builder.blocks_for_left((1, 1))}
        expected = np.array(
            [[0.5, np.sqrt(3.0) / 2.0], [np.sqrt(3.0) / 2.0, -0.5]]
        )
        np.testing.assert_allclose(blocks[(1, 4)].matrix, expected, atol=2.0e-14)

    def test_braid_recursion(self) -> None:
        for left in ((0, 0), (1, 1)):
            residuals = braid_residuals(self.builder, left)
            self.assertLess(max(residuals.values()), 2.0e-12)

    def test_table_round_trip(self) -> None:
        blocks = self.builder.build(1)
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "swap.tbl"
            write_table(filename, blocks, max_prefix_gluons=1)
            restored = read_table(filename)
        self.assertEqual(len(blocks), len(restored))
        for original, copy in zip(blocks, restored, strict=True):
            self.assertEqual(original.left, copy.left)
            self.assertEqual(original.right, copy.right)
            self.assertEqual(original.paths, copy.paths)
            np.testing.assert_array_equal(original.matrix, copy.matrix)


if __name__ == "__main__":
    unittest.main()
