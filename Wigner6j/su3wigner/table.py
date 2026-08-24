"""A deliberately simple, round-trippable text format for Fortran input."""

from __future__ import annotations

from pathlib import Path as FilePath
from typing import Iterable

import numpy as np

from .recoupling import Path, SwapBlock

MAGIC = "SU3_ADJOINT_SWAP_TABLE_V1"

_METADATA_WIDTHS = {
    "NC": 1,
    "ADJOINT_PQ": 2,
    "MAX_PREFIX_GLUONS": 1,
    "NBLOCKS": 1,
    "NPATHS": 1,
    "NVALUES": 1,
}
_SECTION_MARKERS = {
    "BEGIN_BLOCKS",
    "END_BLOCKS",
    "BEGIN_PATHS",
    "END_PATHS",
    "BEGIN_VALUES",
    "END_VALUES",
    "END_TABLE",
}


def _parity_value(value: int | None) -> int:
    return 0 if value is None else value


def write_table(
    filename: str | FilePath,
    blocks: Iterable[SwapBlock],
    *,
    max_prefix_gluons: int,
    zero_tolerance: float = 0.0,
) -> None:
    """Write swap blocks as whitespace-delimited ASCII.

    Values use 17 digits after the decimal point (18 significant decimal
    digits), sufficient to recover the generating IEEE binary64 number
    exactly.  A zero tolerance may be used to omit numerical zeros, but matrix
    dimensions and paths remain explicit.
    """

    block_tuple = tuple(blocks)
    for block_index, block in enumerate(block_tuple):
        if not np.all(np.isfinite(block.matrix)):
            raise ValueError(f"non-finite coefficient in block {block_index}")
    path_count = sum(block.size for block in block_tuple)
    values = [
        (block_index, row, col, float(value))
        for block_index, block in enumerate(block_tuple)
        for row in range(block.size)
        for col, value in enumerate(block.matrix[row])
        if abs(value) > zero_tolerance
    ]
    lines = [
        MAGIC,
        "# Eq8 convention: B(g1,g2)[in] = sum_out W[out,in] B(g2,g1)[out]",
        "# Irrep convention: SU(3) Dynkin labels p q; adjoint=(1,1)",
        "# Vertex convention: isometric real CG maps; multiplicities are zero-based",
        "# Self vertex mult=0: normalized action v tensor X -> rho(X)v (except 8x8)",
        "# 8x8->8: mult=0 Jordan/symmetric, mult=1 Lie-bracket/antisymmetric",
        "NC 3",
        "ADJOINT_PQ 1 1",
        f"MAX_PREFIX_GLUONS {max_prefix_gluons}",
        f"NBLOCKS {len(block_tuple)}",
        f"NPATHS {path_count}",
        f"NVALUES {len(values)}",
        "BEGIN_BLOCKS",
        "# block left_p left_q right_p right_q size",
    ]
    lines.extend(
        f"{index} {block.left[0]} {block.left[1]} "
        f"{block.right[0]} {block.right[1]} {block.size}"
        for index, block in enumerate(block_tuple)
    )
    lines.extend(
        (
            "END_BLOCKS",
            "BEGIN_PATHS",
            "# block path middle_p middle_q left_mult right_mult left_parity right_parity",
        )
    )
    for block_index, block in enumerate(block_tuple):
        lines.extend(
            f"{block_index} {path_index} {path.middle[0]} {path.middle[1]} "
            f"{path.left_multiplicity} {path.right_multiplicity} "
            f"{_parity_value(path.left_exchange_parity)} "
            f"{_parity_value(path.right_exchange_parity)}"
            for path_index, path in enumerate(block.paths)
        )
    lines.extend(
        (
            "END_PATHS",
            "BEGIN_VALUES",
            "# block out_path in_path coefficient",
        )
    )
    lines.extend(
        f"{block} {row} {col} {value:.17e}" for block, row, col, value in values
    )
    lines.extend(("END_VALUES", "END_TABLE", ""))
    FilePath(filename).write_text("\n".join(lines), encoding="ascii")


def _next_data_line(
    lines: list[str], index: int
) -> tuple[int, tuple[int, str] | None]:
    while index < len(lines):
        line_number = index + 1
        line = lines[index].strip()
        index += 1
        if line and not line.startswith("#"):
            return index, (line_number, line)
    return index, None


def _read_section(
    lines: list[str], index: int, end_marker: str
) -> tuple[int, list[tuple[int, str]]]:
    rows: list[tuple[int, str]] = []
    while True:
        index, item = _next_data_line(lines, index)
        if item is None:
            raise ValueError(f"missing {end_marker}")
        line_number, line = item
        if line == end_marker:
            return index, rows
        if line in _SECTION_MARKERS:
            raise ValueError(
                f"unexpected {line} at line {line_number}; expected {end_marker}"
            )
        rows.append(item)


def _expect_marker(lines: list[str], index: int, marker: str) -> int:
    index, item = _next_data_line(lines, index)
    if item is None:
        raise ValueError(f"missing {marker}")
    line_number, line = item
    if line != marker:
        raise ValueError(f"expected {marker} at line {line_number}, found {line!r}")
    return index


def _parse_int_row(
    item: tuple[int, str], width: int, description: str
) -> tuple[int, ...]:
    line_number, line = item
    fields = line.split()
    if len(fields) != width:
        raise ValueError(
            f"invalid {description} row at line {line_number}: "
            f"expected {width} fields"
        )
    try:
        return tuple(map(int, fields))
    except ValueError as error:
        raise ValueError(
            f"invalid integer in {description} row at line {line_number}"
        ) from error


def read_table(filename: str | FilePath) -> tuple[SwapBlock, ...]:
    """Read a table written by :func:`write_table`."""

    lines = FilePath(filename).read_text(encoding="ascii").splitlines()
    if not lines or lines[0] != MAGIC:
        raise ValueError(f"not a {MAGIC} file")

    metadata_values: dict[str, tuple[int, ...]] = {}
    index = 1
    while True:
        index, item = _next_data_line(lines, index)
        if item is None:
            raise ValueError("missing BEGIN_BLOCKS")
        line_number, line = item
        if line == "BEGIN_BLOCKS":
            break
        if line in _SECTION_MARKERS:
            raise ValueError(f"unexpected {line} at line {line_number}")
        fields = line.split()
        name = fields[0]
        if name not in _METADATA_WIDTHS:
            raise ValueError(f"unknown metadata {name!r} at line {line_number}")
        if name in metadata_values:
            raise ValueError(f"duplicate metadata {name} at line {line_number}")
        values = _parse_int_row(
            (line_number, " ".join(fields[1:])),
            _METADATA_WIDTHS[name],
            name,
        )
        metadata_values[name] = values

    missing_metadata = set(_METADATA_WIDTHS) - set(metadata_values)
    if missing_metadata:
        missing = ", ".join(sorted(missing_metadata))
        raise ValueError(f"missing required metadata: {missing}")
    if metadata_values["NC"] != (3,):
        raise ValueError("NC must be 3")
    if metadata_values["ADJOINT_PQ"] != (1, 1):
        raise ValueError("ADJOINT_PQ must be 1 1")

    max_prefix_gluons = metadata_values["MAX_PREFIX_GLUONS"][0]
    declared_blocks = metadata_values["NBLOCKS"][0]
    declared_paths = metadata_values["NPATHS"][0]
    declared_values = metadata_values["NVALUES"][0]
    if max_prefix_gluons < 0:
        raise ValueError("MAX_PREFIX_GLUONS must be non-negative")
    if min(declared_blocks, declared_paths, declared_values) < 0:
        raise ValueError("NBLOCKS, NPATHS, and NVALUES must be non-negative")

    index, block_rows = _read_section(lines, index, "END_BLOCKS")
    index = _expect_marker(lines, index, "BEGIN_PATHS")
    index, path_rows = _read_section(lines, index, "END_PATHS")
    index = _expect_marker(lines, index, "BEGIN_VALUES")
    index, value_rows = _read_section(lines, index, "END_VALUES")
    index = _expect_marker(lines, index, "END_TABLE")
    index, trailing = _next_data_line(lines, index)
    if trailing is not None:
        line_number, line = trailing
        raise ValueError(
            f"unexpected data after END_TABLE at line {line_number}: {line!r}"
        )

    if len(block_rows) != declared_blocks:
        raise ValueError(
            f"NBLOCKS declares {declared_blocks}, found {len(block_rows)} rows"
        )
    if len(path_rows) != declared_paths:
        raise ValueError(
            f"NPATHS declares {declared_paths}, found {len(path_rows)} rows"
        )
    if len(value_rows) != declared_values:
        raise ValueError(
            f"NVALUES declares {declared_values}, found {len(value_rows)} rows"
        )

    metadata: dict[int, tuple[tuple[int, int], tuple[int, int], int]] = {}
    for item in block_rows:
        block, lp, lq, rp, rq, size = _parse_int_row(item, 6, "block")
        if block < 0 or block >= declared_blocks:
            raise ValueError(f"block ID {block} is outside 0..{declared_blocks - 1}")
        if block in metadata:
            raise ValueError(f"duplicate block ID {block}")
        if min(lp, lq, rp, rq) < 0:
            raise ValueError(f"negative irrep label in block {block}")
        if size < 1:
            raise ValueError(f"block {block} must have positive size")
        metadata[block] = ((lp, lq), (rp, rq), size)
    if set(metadata) != set(range(declared_blocks)):
        raise ValueError("block IDs must be contiguous starting at zero")

    expected_paths = sum(data[2] for data in metadata.values())
    if declared_paths != expected_paths:
        raise ValueError(
            f"NPATHS declares {declared_paths}, block sizes require {expected_paths}"
        )
    paths: dict[int, list[Path | None]] = {
        block: [None] * data[2] for block, data in metadata.items()
    }
    for item in path_rows:
        block, path_index, mp, mq, left_mult, right_mult, left_parity, right_parity = (
            _parse_int_row(item, 8, "path")
        )
        if block not in metadata:
            raise ValueError(f"path references unknown block {block}")
        if path_index < 0 or path_index >= metadata[block][2]:
            raise ValueError(f"path index {path_index} is outside block {block}")
        if paths[block][path_index] is not None:
            raise ValueError(f"duplicate path index {path_index} in block {block}")
        if min(mp, mq) < 0:
            raise ValueError(f"negative middle irrep label in block {block}")
        if min(left_mult, right_mult) < 0:
            raise ValueError(
                f"negative multiplicity in block {block} path {path_index}"
            )
        if left_parity not in (-1, 0, 1) or right_parity not in (-1, 0, 1):
            raise ValueError(f"invalid parity in block {block} path {path_index}")
        paths[block][path_index] = Path(
            middle=(mp, mq),
            left_multiplicity=left_mult,
            right_multiplicity=right_mult,
            left_exchange_parity=None if left_parity == 0 else left_parity,
            right_exchange_parity=None if right_parity == 0 else right_parity,
        )
    matrices = {
        block: np.zeros((data[2], data[2]), dtype=float)
        for block, data in metadata.items()
    }
    seen_values: set[tuple[int, int, int]] = set()
    for line_number, row in value_rows:
        fields = row.split()
        if len(fields) != 4:
            raise ValueError(
                f"invalid value row at line {line_number}: expected 4 fields"
            )
        try:
            block, out_path, in_path = map(int, fields[:3])
            value = float(fields[3])
        except ValueError as error:
            raise ValueError(f"invalid value row at line {line_number}") from error
        if block not in metadata:
            raise ValueError(f"value references unknown block {block}")
        block_size = metadata[block][2]
        if not (0 <= out_path < block_size and 0 <= in_path < block_size):
            raise ValueError(f"matrix index outside block {block}")
        key = (block, out_path, in_path)
        if key in seen_values:
            raise ValueError(
                f"duplicate value index ({out_path}, {in_path}) in block {block}"
            )
        if not np.isfinite(value):
            raise ValueError(f"non-finite coefficient at line {line_number}")
        seen_values.add(key)
        matrices[block][out_path, in_path] = value

    result: list[SwapBlock] = []
    for block in sorted(metadata):
        left, right, _size = metadata[block]
        if any(path is None for path in paths[block]):
            raise ValueError(f"missing path in block {block}")
        result.append(
            SwapBlock(
                left=left,
                right=right,
                paths=tuple(path for path in paths[block] if path is not None),
                matrix=matrices[block],
            )
        )
    return tuple(result)
