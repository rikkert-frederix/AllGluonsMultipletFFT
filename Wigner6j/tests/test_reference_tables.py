from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np

from su3wigner import SwapTableBuilder, read_table


class ReferenceTableTests(unittest.TestCase):
    def test_checked_in_tables_are_valid(self) -> None:
        data_directory = Path(__file__).resolve().parents[1] / "data"

        for prefix in range(2, 7):
            filename = data_directory / f"su3_adjoint_swap_prefix_{prefix}.tbl"
            with self.subTest(prefix=prefix):
                blocks = read_table(filename)
                self.assertTrue(blocks)
                self.assertEqual(
                    len({(block.left, block.right) for block in blocks}),
                    len(blocks),
                )

                for block in blocks:
                    self.assertGreater(block.size, 0)
                    self.assertEqual(block.matrix.shape, (block.size, block.size))
                    self.assertEqual(len(set(block.paths)), block.size)
                    self.assertTrue(np.isfinite(block.matrix).all())
                    for label in (block.left, block.right):
                        self.assertTrue(all(value >= 0 for value in label))
                    for path in block.paths:
                        self.assertTrue(all(value >= 0 for value in path.middle))
                        self.assertGreaterEqual(path.left_multiplicity, 0)
                        self.assertGreaterEqual(path.right_multiplicity, 0)
                        self.assertIn(path.left_exchange_parity, (None, -1, 1))
                        self.assertIn(path.right_exchange_parity, (None, -1, 1))
                    SwapTableBuilder.validate_block(block, tolerance=5.0e-10)


if __name__ == "__main__":
    unittest.main()
