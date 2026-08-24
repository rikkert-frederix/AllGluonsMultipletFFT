"""Weight-sector storage for sparse Clebsch--Gordan embeddings.

The embeddings used by the recoupling backends are dense only within exact
weight sectors.  :class:`PackedEmbedding` stores those rectangular blocks in
one contiguous value buffer and keeps compact integer descriptors for their
rows and columns.  It is intentionally an internal building block and is not
re-exported from :mod:`su3wigner`.
"""

from __future__ import annotations

from dataclasses import dataclass
from operator import index
from typing import Iterable, Sequence

import numpy as np


Weight = tuple[int, int]


def _normalize_weights(
    weights: Iterable[Sequence[int]], expected: int, name: str
) -> tuple[Weight, ...]:
    """Validate and freeze a sequence of two-component integral weights."""

    supplied = tuple(weights)
    if len(supplied) != expected:
        raise ValueError(f"{name} has length {len(supplied)}, expected {expected}")

    result: list[Weight] = []
    for position, weight in enumerate(supplied):
        try:
            first, second = weight
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"{name}[{position}] must have exactly two components"
            ) from error
        try:
            result.append((index(first), index(second)))
        except TypeError as error:
            raise TypeError(
                f"{name}[{position}] components must be integers"
            ) from error
    return tuple(result)


def _compact_index_dtype(maximum: int) -> np.dtype[np.signedinteger]:
    """Use 32-bit descriptors unless an unusually large array requires 64."""

    if maximum <= np.iinfo(np.int32).max:
        return np.dtype(np.int32)
    return np.dtype(np.int64)


def _readonly(array: np.ndarray) -> np.ndarray:
    result = np.ascontiguousarray(array)
    result.setflags(write=False)
    return result


def _half_open_range(
    start: int, stop: int, limit: int, description: str
) -> tuple[int, int]:
    try:
        normalized_start = index(start)
        normalized_stop = index(stop)
    except TypeError as error:
        raise TypeError(f"{description} bounds must be integers") from error
    if not 0 <= normalized_start <= normalized_stop <= limit:
        raise IndexError(
            f"invalid {description} [{normalized_start}, {normalized_stop}) "
            f"for extent {limit}"
        )
    return normalized_start, normalized_stop


@dataclass(frozen=True, slots=True, init=False, eq=False, repr=False)
class PackedEmbedding:
    """Immutable exact-weight-sector representation of a matrix.

    Product rows are assumed to be ordered as ``source x particle``.  Slab
    access therefore returns the three-dimensional view expected by the
    recoupling contractions, while :meth:`to_dense` reconstructs the original
    two-dimensional matrix.
    """

    _product_dimension: int
    _target_dimension: int
    _particle_dimension: int
    _values: np.ndarray
    _value_offsets: np.ndarray
    _row_offsets: np.ndarray
    _column_offsets: np.ndarray
    _row_indices: np.ndarray
    _column_indices: np.ndarray
    _column_sector: np.ndarray
    _column_local: np.ndarray

    def __init__(
        self,
        embedding: np.ndarray,
        product_weights: Iterable[Sequence[int]],
        target_weights: Iterable[Sequence[int]],
        *,
        particle_dimension: int,
    ) -> None:
        """Pack ``embedding`` after verifying every off-sector entry is zero."""

        try:
            particle_dimension = index(particle_dimension)
        except TypeError as error:
            raise TypeError("particle_dimension must be an integer") from error
        if particle_dimension <= 0:
            raise ValueError("particle_dimension must be positive")

        dense = np.asarray(embedding)
        if dense.ndim != 2:
            raise ValueError("embedding must be a two-dimensional array")
        if not np.issubdtype(dense.dtype, np.number):
            raise TypeError("embedding must have a numeric dtype")
        product_dimension, target_dimension = dense.shape
        if product_dimension % particle_dimension:
            raise ValueError(
                f"product dimension {product_dimension} is not divisible by "
                f"particle dimension {particle_dimension}"
            )
        if not np.all(np.isfinite(dense)):
            raise ValueError("embedding contains a non-finite entry")

        normalized_product_weights = _normalize_weights(
            product_weights, product_dimension, "product_weights"
        )
        normalized_target_weights = _normalize_weights(
            target_weights, target_dimension, "target_weights"
        )

        rows_by_weight: dict[Weight, list[int]] = {}
        for row, weight in enumerate(normalized_product_weights):
            rows_by_weight.setdefault(weight, []).append(row)
        columns_by_weight: dict[Weight, list[int]] = {}
        for column, weight in enumerate(normalized_target_weights):
            columns_by_weight.setdefault(weight, []).append(column)

        # This is deliberately an exact check, not a tolerance.  Packing must
        # never hide a numerically small violation of weight preservation.
        for column, weight in enumerate(normalized_target_weights):
            invalid = dense[:, column] != 0
            allowed = rows_by_weight.get(weight)
            if allowed:
                invalid[allowed] = False
            if np.any(invalid):
                row = int(np.flatnonzero(invalid)[0])
                raise ValueError(
                    "embedding has a nonzero off-sector entry at "
                    f"({row}, {column}): product weight "
                    f"{normalized_product_weights[row]} != target weight {weight}"
                )

        sector_weights = tuple(
            sorted(rows_by_weight.keys() & columns_by_weight.keys())
        )
        row_blocks = [rows_by_weight[weight] for weight in sector_weights]
        column_blocks = [columns_by_weight[weight] for weight in sector_weights]

        value_offsets = [0]
        row_offsets = [0]
        column_offsets = [0]
        for rows, columns in zip(row_blocks, column_blocks, strict=True):
            value_offsets.append(value_offsets[-1] + len(rows) * len(columns))
            row_offsets.append(row_offsets[-1] + len(rows))
            column_offsets.append(column_offsets[-1] + len(columns))

        descriptor_maximum = max(
            product_dimension,
            target_dimension,
            len(sector_weights),
            value_offsets[-1],
        )
        descriptor_dtype = _compact_index_dtype(descriptor_maximum)
        value_offsets_array = np.asarray(value_offsets, dtype=descriptor_dtype)
        row_offsets_array = np.asarray(row_offsets, dtype=descriptor_dtype)
        column_offsets_array = np.asarray(column_offsets, dtype=descriptor_dtype)
        row_indices = np.asarray(
            [row for rows in row_blocks for row in rows], dtype=descriptor_dtype
        )
        column_indices = np.asarray(
            [column for columns in column_blocks for column in columns],
            dtype=descriptor_dtype,
        )

        values = np.empty(value_offsets[-1], dtype=dense.dtype)
        column_sector = np.full(target_dimension, -1, dtype=descriptor_dtype)
        column_local = np.full(target_dimension, -1, dtype=descriptor_dtype)
        for sector, (rows, columns) in enumerate(
            zip(row_blocks, column_blocks, strict=True)
        ):
            value_start = value_offsets[sector]
            value_stop = value_offsets[sector + 1]
            values[value_start:value_stop] = dense[np.ix_(rows, columns)].reshape(-1)
            for local, column in enumerate(columns):
                column_sector[column] = sector
                column_local[column] = local

        object.__setattr__(self, "_product_dimension", product_dimension)
        object.__setattr__(self, "_target_dimension", target_dimension)
        object.__setattr__(self, "_particle_dimension", particle_dimension)
        object.__setattr__(self, "_values", _readonly(values))
        object.__setattr__(self, "_value_offsets", _readonly(value_offsets_array))
        object.__setattr__(self, "_row_offsets", _readonly(row_offsets_array))
        object.__setattr__(self, "_column_offsets", _readonly(column_offsets_array))
        object.__setattr__(self, "_row_indices", _readonly(row_indices))
        object.__setattr__(self, "_column_indices", _readonly(column_indices))
        object.__setattr__(self, "_column_sector", _readonly(column_sector))
        object.__setattr__(self, "_column_local", _readonly(column_local))

    @classmethod
    def from_dense(
        cls,
        embedding: np.ndarray,
        product_weights: Iterable[Sequence[int]],
        target_weights: Iterable[Sequence[int]],
        *,
        particle_dimension: int,
    ) -> PackedEmbedding:
        """Named constructor mirroring the operation performed by ``__init__``."""

        return cls(
            embedding,
            product_weights,
            target_weights,
            particle_dimension=particle_dimension,
        )

    @property
    def shape(self) -> tuple[int, int]:
        return self._product_dimension, self._target_dimension

    @property
    def dtype(self) -> np.dtype:
        return self._values.dtype

    @property
    def source_dimension(self) -> int:
        return self._product_dimension // self._particle_dimension

    @property
    def particle_dimension(self) -> int:
        return self._particle_dimension

    @property
    def target_dimension(self) -> int:
        return self._target_dimension

    @property
    def sector_count(self) -> int:
        return self._value_offsets.size - 1

    @property
    def value_count(self) -> int:
        return self._values.size

    @property
    def packed_values(self) -> np.ndarray:
        """Read-only view of the sole coefficient buffer."""

        result = self._values.view()
        result.setflags(write=False)
        return result

    @property
    def data_nbytes(self) -> int:
        return self._values.nbytes

    @property
    def descriptor_nbytes(self) -> int:
        return sum(
            array.nbytes
            for array in (
                self._value_offsets,
                self._row_offsets,
                self._column_offsets,
                self._row_indices,
                self._column_indices,
                self._column_sector,
                self._column_local,
            )
        )

    @property
    def nbytes(self) -> int:
        """NumPy payload bytes, excluding small Python object headers."""

        return self.data_nbytes + self.descriptor_nbytes

    def __repr__(self) -> str:
        return (
            f"PackedEmbedding(shape={self.shape}, dtype={self.dtype}, "
            f"sectors={self.sector_count}, values={self.value_count})"
        )

    def _sector_parts(
        self, sector: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        row_start = int(self._row_offsets[sector])
        row_stop = int(self._row_offsets[sector + 1])
        column_start = int(self._column_offsets[sector])
        column_stop = int(self._column_offsets[sector + 1])
        value_start = int(self._value_offsets[sector])
        value_stop = int(self._value_offsets[sector + 1])
        rows = self._row_indices[row_start:row_stop]
        columns = self._column_indices[column_start:column_stop]
        block = self._values[value_start:value_stop].reshape(rows.size, columns.size)
        return rows, columns, block

    def to_dense(self) -> np.ndarray:
        """Reconstruct a new dense matrix in the original basis order."""

        result = np.zeros(self.shape, dtype=self.dtype)
        for sector in range(self.sector_count):
            rows, columns, block = self._sector_parts(sector)
            result[np.ix_(rows, columns)] = block
        return result

    def source_slab(self, start: int, stop: int) -> np.ndarray:
        """Return source rows ``[start, stop)`` as ``source x particle x target``."""

        start, stop = _half_open_range(
            start, stop, self.source_dimension, "source slab"
        )
        flat_start = start * self.particle_dimension
        flat_stop = stop * self.particle_dimension
        result = np.zeros(
            (stop - start, self.particle_dimension, self.target_dimension),
            dtype=self.dtype,
        )
        flat_result = result.reshape(-1, self.target_dimension)
        for sector in range(self.sector_count):
            rows, columns, block = self._sector_parts(sector)
            local_start = int(np.searchsorted(rows, flat_start, side="left"))
            local_stop = int(np.searchsorted(rows, flat_stop, side="left"))
            if local_start == local_stop:
                continue
            slab_rows = rows[local_start:local_stop] - flat_start
            flat_result[np.ix_(slab_rows, columns)] = block[local_start:local_stop]
        return result

    def target_columns(self, start: int, stop: int) -> np.ndarray:
        """Return target columns as ``source x particle x selected-target``."""

        start, stop = _half_open_range(
            start, stop, self.target_dimension, "target columns"
        )
        result = np.zeros(
            (self.source_dimension, self.particle_dimension, stop - start),
            dtype=self.dtype,
        )
        flat_result = result.reshape(self._product_dimension, stop - start)
        for sector in range(self.sector_count):
            rows, columns, block = self._sector_parts(sector)
            local_start = int(np.searchsorted(columns, start, side="left"))
            local_stop = int(np.searchsorted(columns, stop, side="left"))
            if local_start == local_stop:
                continue
            selected_columns = columns[local_start:local_stop] - start
            flat_result[np.ix_(rows, selected_columns)] = block[
                :, local_start:local_stop
            ]
        return result

    def column_entries(self, column: int) -> tuple[np.ndarray, np.ndarray]:
        """Return flat product-row indices and stored values for one column."""

        try:
            column = index(column)
        except TypeError as error:
            raise TypeError("target column must be an integer") from error
        if not 0 <= column < self.target_dimension:
            raise IndexError(
                f"target column {column} is outside extent {self.target_dimension}"
            )
        sector = int(self._column_sector[column])
        if sector < 0:
            return self._row_indices[:0], self._values[:0]
        rows, _columns, block = self._sector_parts(sector)
        local = int(self._column_local[column])
        return rows, block[:, local]

    def target_column(self, column: int) -> np.ndarray:
        """Return one column as a dense ``source x particle`` array."""

        rows, values = self.column_entries(column)
        result = np.zeros(self._product_dimension, dtype=self.dtype)
        result[rows] = values
        return result.reshape(self.source_dimension, self.particle_dimension)
