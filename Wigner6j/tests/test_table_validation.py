from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from su3wigner import Path as FusionPath
from su3wigner import SwapBlock, read_table, write_table


VALID_TABLE = """\
SU3_ADJOINT_SWAP_TABLE_V1
NC 3
ADJOINT_PQ 1 1
MAX_PREFIX_GLUONS 0
NBLOCKS 1
NPATHS 1
NVALUES 1
BEGIN_BLOCKS
0 0 0 0 0 1
END_BLOCKS
BEGIN_PATHS
0 0 1 1 0 0 0 1
END_PATHS
BEGIN_VALUES
0 0 0 1.0
END_VALUES
END_TABLE
"""


class TableValidationTests(unittest.TestCase):
    def read_text(self, contents: str) -> tuple[SwapBlock, ...]:
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "table.tbl"
            filename.write_text(contents, encoding="ascii")
            return read_table(filename)

    def assert_invalid(self, contents: str) -> None:
        with self.assertRaises(ValueError):
            self.read_text(contents)

    def test_minimal_and_checked_in_tables_are_accepted(self) -> None:
        blocks = self.read_text(VALID_TABLE)
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks[0].paths[0].right_exchange_parity, 1)
        np.testing.assert_array_equal(blocks[0].matrix, [[1.0]])

        data_directory = Path(__file__).parents[1] / "data"
        for filename in sorted(data_directory.glob("*.tbl")):
            with self.subTest(filename=filename.name):
                self.assertGreater(len(read_table(filename)), 0)

    def test_required_metadata_and_declared_counts_are_validated(self) -> None:
        required_rows = (
            "NC 3",
            "ADJOINT_PQ 1 1",
            "MAX_PREFIX_GLUONS 0",
            "NBLOCKS 1",
            "NPATHS 1",
            "NVALUES 1",
        )
        for row in required_rows:
            with self.subTest(missing=row):
                self.assert_invalid(VALID_TABLE.replace(f"{row}\n", "", 1))

        replacements = {
            "wrong NC": ("NC 3", "NC 2"),
            "wrong adjoint": ("ADJOINT_PQ 1 1", "ADJOINT_PQ 1 0"),
            "negative cutoff": ("MAX_PREFIX_GLUONS 0", "MAX_PREFIX_GLUONS -1"),
            "block count": ("NBLOCKS 1", "NBLOCKS 2"),
            "path count": ("NPATHS 1", "NPATHS 2"),
            "value count": ("NVALUES 1", "NVALUES 2"),
            "negative count": ("NVALUES 1", "NVALUES -1"),
            "unknown metadata": ("NC 3", "FUTURE_FIELD 3"),
            "duplicate metadata": ("NC 3", "NC 3\nNC 3"),
        }
        for name, (old, new) in replacements.items():
            with self.subTest(name=name):
                self.assert_invalid(VALID_TABLE.replace(old, new, 1))

    def test_sections_must_be_ordered_complete_and_terminated(self) -> None:
        variants = {
            "missing block end": VALID_TABLE.replace("END_BLOCKS\n", "", 1),
            "missing paths begin": VALID_TABLE.replace("BEGIN_PATHS\n", "", 1),
            "missing paths end": VALID_TABLE.replace("END_PATHS\n", "", 1),
            "missing values begin": VALID_TABLE.replace("BEGIN_VALUES\n", "", 1),
            "missing values end": VALID_TABLE.replace("END_VALUES\n", "", 1),
            "missing table end": VALID_TABLE.replace("END_TABLE\n", "", 1),
            "trailing data": VALID_TABLE + "EXTRA\n",
            "misordered marker": VALID_TABLE.replace(
                "END_BLOCKS\nBEGIN_PATHS", "BEGIN_PATHS\nEND_BLOCKS", 1
            ),
        }
        for name, contents in variants.items():
            with self.subTest(name=name):
                self.assert_invalid(contents)

    def test_block_ids_labels_and_sizes_are_validated(self) -> None:
        block_row = "0 0 0 0 0 1"
        duplicate = VALID_TABLE.replace("NBLOCKS 1", "NBLOCKS 2", 1).replace(
            block_row, f"{block_row}\n{block_row}", 1
        )
        variants = {
            "duplicate ID": duplicate,
            "noncontiguous ID": VALID_TABLE.replace(block_row, "1 0 0 0 0 1", 1),
            "negative left label": VALID_TABLE.replace(
                block_row, "0 -1 0 0 0 1", 1
            ),
            "negative right label": VALID_TABLE.replace(
                block_row, "0 0 0 0 -1 1", 1
            ),
            "zero size": VALID_TABLE.replace(block_row, "0 0 0 0 0 0", 1),
        }
        for name, contents in variants.items():
            with self.subTest(name=name):
                self.assert_invalid(contents)

    def test_path_indices_and_fields_are_validated(self) -> None:
        path_row = "0 0 1 1 0 0 0 1"
        duplicate = (
            VALID_TABLE.replace("NPATHS 1", "NPATHS 2", 1)
            .replace("0 0 0 0 0 1", "0 0 0 0 0 2", 1)
            .replace(path_row, f"{path_row}\n{path_row}", 1)
        )
        variants = {
            "duplicate index": duplicate,
            "unknown block": VALID_TABLE.replace(path_row, "1 0 1 1 0 0 0 1", 1),
            "out-of-range index": VALID_TABLE.replace(
                path_row, "0 1 1 1 0 0 0 1", 1
            ),
            "negative middle label": VALID_TABLE.replace(
                path_row, "0 0 -1 1 0 0 0 1", 1
            ),
            "negative multiplicity": VALID_TABLE.replace(
                path_row, "0 0 1 1 -1 0 0 1", 1
            ),
            "invalid left parity": VALID_TABLE.replace(
                path_row, "0 0 1 1 0 0 2 1", 1
            ),
            "invalid right parity": VALID_TABLE.replace(
                path_row, "0 0 1 1 0 0 0 -2", 1
            ),
        }
        for name, contents in variants.items():
            with self.subTest(name=name):
                self.assert_invalid(contents)

    def test_value_indices_and_coefficients_are_validated(self) -> None:
        value_row = "0 0 0 1.0"
        duplicate = VALID_TABLE.replace("NVALUES 1", "NVALUES 2", 1).replace(
            value_row, f"{value_row}\n{value_row}", 1
        )
        variants = {
            "duplicate index": duplicate,
            "unknown block": VALID_TABLE.replace(value_row, "1 0 0 1.0", 1),
            "bad output index": VALID_TABLE.replace(value_row, "0 1 0 1.0", 1),
            "bad input index": VALID_TABLE.replace(value_row, "0 0 -1 1.0", 1),
            "NaN": VALID_TABLE.replace(value_row, "0 0 0 nan", 1),
            "positive infinity": VALID_TABLE.replace(value_row, "0 0 0 inf", 1),
            "negative infinity": VALID_TABLE.replace(value_row, "0 0 0 -inf", 1),
        }
        for name, contents in variants.items():
            with self.subTest(name=name):
                self.assert_invalid(contents)

    def test_writer_rejects_nonfinite_coefficients_before_writing(self) -> None:
        path = FusionPath(
            middle=(1, 1),
            left_multiplicity=0,
            right_multiplicity=0,
        )
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "table.tbl"
            for coefficient in (np.nan, np.inf, -np.inf):
                with self.subTest(coefficient=coefficient):
                    filename.write_text("unchanged", encoding="ascii")
                    block = SwapBlock(
                        left=(0, 0),
                        right=(0, 0),
                        paths=(path,),
                        matrix=np.array([[coefficient]]),
                    )
                    with self.assertRaises(ValueError):
                        write_table(filename, (block,), max_prefix_gluons=0)
                    self.assertEqual(filename.read_text(encoding="ascii"), "unchanged")


if __name__ == "__main__":
    unittest.main()
