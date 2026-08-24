"""Reduce adjacent-adjoint recoupling to elementary quark recouplings.

The construction in this module is deliberately different from both direct
backends.  Each adjoint vertex is first expanded in the complete set of
``3``--``3-bar`` fusion paths.  Swapping two adjoints then becomes the
permutation

    (3, 3-bar, 3, 3-bar) -> (3, 3-bar, 3, 3-bar)

which exchanges the two ordered pairs.  That permutation is evaluated as
four adjacent, elementary two-line recouplings.  Thus a two-adjoint
coefficient is a sum of products of one-adjoint/fundamental coefficients,
and those in turn are products of terminal two-fundamental coefficients.

No formula from an external convention is used.  The terminal coefficients
are normalized overlaps of multiplicity-free fundamental coupling trees,
and the adjoint-vertex expansion follows from orthogonal completeness.  The
already fixed adjoint vertices supply only the unavoidable phase and
multiplicity gauge required by the public table convention.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from operator import index
from typing import Sequence

import numpy as np

from .coupling import (
    AdjointCouplings,
    _TensorProductAction,
    _product_action,
    _raising_annihilation,
)
from .linalg import canonical_nullspace, descendant_basis
from .packed_embedding import PackedEmbedding
from .recoupling import Path, SwapBlock, SwapTableBuilder
from .representations import Irrep, IrrepFactory, Label


FUNDAMENTAL: Label = (1, 0)
ANTIFUNDAMENTAL: Label = (0, 1)
ADJOINT_SPLIT_SEQUENCE: tuple[Label, ...] = (
    FUNDAMENTAL,
    ANTIFUNDAMENTAL,
    FUNDAMENTAL,
    ANTIFUNDAMENTAL,
)


def _require_finite(description: str, *arrays: np.ndarray) -> None:
    """Reject non-finite intermediate data before residual comparisons."""

    if any(not np.all(np.isfinite(array)) for array in arrays):
        raise ArithmeticError(f"{description} contains a non-finite value")


class _ElementaryCouplingStorage:
    """Private slots kept out of the public dataclass field contract."""

    __slots__ = ("_packed_embedding", "_embedding_cache", "_particle_dimension")
    _packed_embedding: PackedEmbedding | None
    _embedding_cache: np.ndarray | None
    _particle_dimension: int


@dataclass(frozen=True, slots=True)
class ElementaryCoupling(_ElementaryCouplingStorage):
    """An isometry ``target -> source tensor particle``."""

    source: Label
    particle: Label
    target: Label
    embedding: np.ndarray

    def __post_init__(self) -> None:
        """Initialize private storage for a public dense construction."""

        object.__setattr__(self, "_packed_embedding", None)
        object.__setattr__(
            self, "_embedding_cache", object.__getattribute__(self, "embedding")
        )
        object.__setattr__(self, "_particle_dimension", 3)

    def __getattribute__(self, name: str):
        """Materialize the public dataclass field on its first actual access."""

        if name == "embedding":
            try:
                cached = object.__getattribute__(self, "_embedding_cache")
            except AttributeError:
                # The generated dataclass initializer sets ``embedding`` before
                # ``__post_init__`` has initialized the private cache slots.
                return object.__getattribute__(self, "embedding")
            if cached is None:
                packed = object.__getattribute__(self, "_packed_embedding")
                if packed is None:
                    return object.__getattribute__(self, "embedding")
                cached = packed.to_dense()
                object.__setattr__(self, "embedding", cached)
                object.__setattr__(self, "_embedding_cache", cached)
            return cached
        return object.__getattribute__(self, name)

    def __getstate__(self) -> tuple[Label, Label, Label, np.ndarray]:
        """Serialize exactly the original dense public dataclass state."""

        return self.source, self.particle, self.target, self.embedding

    def __setstate__(
        self, state: tuple[Label, Label, Label, np.ndarray]
    ) -> None:
        source, particle, target, embedding = state
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "particle", particle)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "embedding", embedding)
        object.__setattr__(self, "_packed_embedding", None)
        object.__setattr__(self, "_embedding_cache", embedding)
        object.__setattr__(self, "_particle_dimension", 3)

    @classmethod
    def _from_packed(
        cls,
        source: Label,
        particle: Label,
        target: Label,
        packed: PackedEmbedding,
    ) -> ElementaryCoupling:
        result = cls.__new__(cls)
        object.__setattr__(result, "source", source)
        object.__setattr__(result, "particle", particle)
        object.__setattr__(result, "target", target)
        object.__setattr__(result, "embedding", None)
        object.__setattr__(result, "_packed_embedding", packed)
        object.__setattr__(result, "_embedding_cache", None)
        object.__setattr__(result, "_particle_dimension", packed.particle_dimension)
        return result

    @property
    def _source_dimension(self) -> int:
        packed = self._packed_embedding
        if packed is not None and self._embedding_cache is None:
            return packed.source_dimension
        return self.embedding.shape[0] // self._particle_dimension

    @property
    def _target_dimension(self) -> int:
        packed = self._packed_embedding
        if packed is not None and self._embedding_cache is None:
            return packed.target_dimension
        return self.embedding.shape[1]

    def _pack(
        self,
        product_weights: Sequence[Label],
        target_weights: Sequence[Label],
    ) -> ElementaryCoupling:
        packed = PackedEmbedding.from_dense(
            self.embedding,
            product_weights,
            target_weights,
            particle_dimension=self._particle_dimension,
        )
        return self._from_packed(
            self.source,
            self.particle,
            self.target,
            packed,
        )

    def _source_slab(self, start: int, stop: int) -> np.ndarray:
        """Return ``source[start:stop] x particle x target`` coefficients."""

        packed = self._packed_embedding
        if packed is not None and self._embedding_cache is None:
            return packed.source_slab(start, stop)
        source_dimension = self.embedding.shape[0] // self._particle_dimension
        start = index(start)
        stop = index(stop)
        if not 0 <= start <= stop <= source_dimension:
            raise IndexError(
                f"invalid source slab [{start}, {stop}) for extent "
                f"{source_dimension}"
            )
        return self.embedding.reshape(
            source_dimension, self._particle_dimension, self.embedding.shape[1]
        )[start:stop].copy()

    def _target_columns(self, start: int, stop: int) -> np.ndarray:
        """Return selected columns as ``source x particle x width``."""

        packed = self._packed_embedding
        if packed is not None and self._embedding_cache is None:
            return packed.target_columns(start, stop)
        start = index(start)
        stop = index(stop)
        target_dimension = self.embedding.shape[1]
        if not 0 <= start <= stop <= target_dimension:
            raise IndexError(
                f"invalid target columns [{start}, {stop}) for extent "
                f"{target_dimension}"
            )
        source_dimension = self.embedding.shape[0] // self._particle_dimension
        return self.embedding.reshape(
            source_dimension, self._particle_dimension, target_dimension
        )[:, :, start:stop].copy()

    def _target_column(self, column: int) -> np.ndarray:
        """Return one column as ``source x particle`` coefficients."""

        packed = self._packed_embedding
        if packed is not None and self._embedding_cache is None:
            return packed.target_column(column)
        column = index(column)
        if not 0 <= column < self.embedding.shape[1]:
            raise IndexError(
                f"target column {column} is outside extent "
                f"{self.embedding.shape[1]}"
            )
        source_dimension = self.embedding.shape[0] // self._particle_dimension
        return self.embedding[:, column].reshape(
            source_dimension, self._particle_dimension
        ).copy()

    def _column_entries(self, column: int) -> tuple[np.ndarray, np.ndarray]:
        """Return stored flat product rows and values for one target column."""

        packed = self._packed_embedding
        if packed is not None and self._embedding_cache is None:
            return packed.column_entries(column)
        column = index(column)
        if not 0 <= column < self.embedding.shape[1]:
            raise IndexError(
                f"target column {column} is outside extent "
                f"{self.embedding.shape[1]}"
            )
        values = self.embedding[:, column]
        rows = np.flatnonzero(values).astype(np.intp, copy=False)
        rows.setflags(write=False)
        selected = np.ascontiguousarray(values[rows])
        selected.setflags(write=False)
        return rows, selected


_ELEMENTARY_SOURCE_SLAB_BYTES = 8 * 1024 * 1024


def _compose_elementary_highest(
    first: ElementaryCoupling,
    second: ElementaryCoupling,
    highest: int,
) -> np.ndarray:
    """Compose two elementary vertices on one target column in source slabs."""

    source_dimension = first._source_dimension
    first_particle_dimension = first._particle_dimension
    middle_dimension = first._target_dimension
    second_particle_dimension = second._particle_dimension
    second_column = second._target_column(highest)
    if second_column.shape[0] != middle_dimension:
        raise ArithmeticError("elementary coupling dimensions do not compose")
    result = np.empty(
        (
            source_dimension,
            first_particle_dimension,
            second_particle_dimension,
        ),
        dtype=float,
    )
    bytes_per_source = (
        first_particle_dimension * middle_dimension * np.dtype(float).itemsize
    )
    rows_per_slab = max(1, _ELEMENTARY_SOURCE_SLAB_BYTES // bytes_per_source)
    for start in range(0, source_dimension, rows_per_slab):
        stop = min(source_dimension, start + rows_per_slab)
        first_slab = first._source_slab(start, stop)
        result[start:stop] = (
            first_slab.reshape(
                (stop - start) * first_particle_dimension, middle_dimension
            )
            @ second_column
        ).reshape(
            stop - start,
            first_particle_dimension,
            second_particle_dimension,
        )
    return result


@dataclass(frozen=True, slots=True)
class ElementaryDecomposition:
    source: Label
    particle: Label
    couplings: tuple[ElementaryCoupling, ...]

    @property
    def targets(self) -> tuple[Label, ...]:
        return tuple(coupling.target for coupling in self.couplings)

    def for_target(self, target: Label) -> tuple[ElementaryCoupling, ...]:
        return tuple(c for c in self.couplings if c.target == target)


class FundamentalCouplings:
    """Multiplicity-free couplings to ``3`` or ``3-bar`` from first principles."""

    def __init__(self, irreps: IrrepFactory | None = None) -> None:
        self.irreps = irreps or IrrepFactory()

    # Packing makes a complete depth-seven working set cheaper than a handful
    # of the former dense elementary maps.  Retain it so evicted decompositions
    # do not force expensive irrep reconstruction, while keeping arbitrary
    # exploratory use bounded.
    @lru_cache(maxsize=256)
    def decompose(
        self, source_label: Label, particle_label: Label
    ) -> ElementaryDecomposition:
        if particle_label not in (FUNDAMENTAL, ANTIFUNDAMENTAL):
            raise ValueError(
                "elementary reduction supports only 3 and 3-bar particles"
            )
        source = self.irreps.get(source_label)
        particle = self.irreps.get(particle_label)
        product_dimension = source.dim * particle.dim
        weights = tuple(
            (
                source_weight[0] + particle_weight[0],
                source_weight[1] + particle_weight[1],
            )
            for source_weight in source.weights
            for particle_weight in particle.weights
        )
        sectors: dict[Label, list[int]] = {}
        for index, weight in enumerate(weights):
            sectors.setdefault(weight, []).append(index)

        product_actions = tuple(
            tuple(
                _product_action(source, particle, row, col)
                for col in range(3)
            )
            for row in range(3)
        )
        raising_1 = product_actions[0][1]
        raising_2 = product_actions[1][2]
        lowering_1 = product_actions[1][0]
        lowering_2 = product_actions[2][1]

        highest_by_target: dict[Label, np.ndarray] = {}
        for weight in sorted(sectors, reverse=True):
            if weight[0] < 0 or weight[1] < 0:
                continue
            indices = sectors[weight]
            annihilation = _raising_annihilation(
                sectors, weight, indices, raising_1, raising_2
            )
            kernel = canonical_nullspace(annihilation)
            if kernel.shape[1] == 0:
                continue
            if kernel.shape[1] != 1:
                raise ArithmeticError(
                    f"{source_label} x {particle_label} has unexpected "
                    f"highest-weight multiplicity {kernel.shape[1]} at {weight}"
                )
            highest_by_target[weight] = kernel[:, 0]

        couplings: list[ElementaryCoupling] = []
        accounted_dimension = 0
        for target_label in sorted(highest_by_target, reverse=True):
            target = self.irreps.get(target_label)
            highest = np.zeros(product_dimension, dtype=float)
            highest[sectors[target_label]] = highest_by_target[target_label]
            descendants, weights_found, recipe = descendant_basis(
                lowering_1,
                lowering_2,
                highest,
                target_label,
                target.dim,
                recipe=target.lowering_recipe,
                coordinate_weights=weights,
                reference_lowerings=(
                    target.generator(1, 0),
                    target.generator(2, 1),
                ),
            )
            if weights_found != target.weights or recipe != target.lowering_recipe:
                raise ArithmeticError(
                    f"basis recursion mismatch for {source_label} x "
                    f"{particle_label} -> {target_label}"
                )
            self._validate_intertwiner(
                source,
                particle,
                target,
                descendants,
                product_actions,
            )
            couplings.append(
                ElementaryCoupling(
                    source=source_label,
                    particle=particle_label,
                    target=target_label,
                    embedding=descendants,
                )._pack(weights, target.weights)
            )
            accounted_dimension += target.dim

        if accounted_dimension != product_dimension:
            raise ArithmeticError(
                f"decomposition of {source_label} x {particle_label} accounts "
                f"for {accounted_dimension}/{product_dimension} dimensions"
            )
        return ElementaryDecomposition(
            source_label,
            particle_label,
            tuple(sorted(couplings, key=lambda coupling: coupling.target)),
        )

    @staticmethod
    def _validate_intertwiner(
        source: Irrep,
        particle: Irrep,
        target: Irrep,
        embedding: np.ndarray,
        product_actions: tuple[tuple[_TensorProductAction, ...], ...],
    ) -> None:
        _require_finite("elementary CG embedding", embedding)
        isometry_error = float(
            np.max(np.abs(embedding.T @ embedding - np.eye(target.dim)))
        )
        if not np.isfinite(isometry_error):
            raise ArithmeticError("elementary CG isometry residual is non-finite")
        worst = 0.0
        for row in range(3):
            for col in range(3):
                if row == col:
                    continue
                product = product_actions[row][col]
                residual = product @ embedding - embedding @ target.generator(
                    row, col
                )
                _require_finite("elementary CG intertwining residual", residual)
                worst = max(worst, float(np.max(np.abs(residual))))
        for first, second in ((0, 1), (1, 2)):
            product = _TensorProductAction(
                source.generator(first, first)
                - source.generator(second, second),
                particle.generator(first, first)
                - particle.generator(second, second),
            )
            target_cartan = target.generator(first, first) - target.generator(
                second, second
            )
            residual = product @ embedding - embedding @ target_cartan
            _require_finite("elementary CG Cartan residual", residual)
            worst = max(worst, float(np.max(np.abs(residual))))
        scale = max(1.0, float(target.dim))
        if isometry_error > 3.0e-9 * scale:
            raise ArithmeticError(
                f"elementary CG isometry residual {isometry_error:.3e}"
            )
        if worst > 3.0e-8 * scale:
            raise ArithmeticError(
                f"elementary CG intertwining residual {worst:.3e}"
            )


@dataclass(frozen=True, slots=True)
class ElementarySwapBlock:
    """Terminal two-fundamental recoupling coefficient matrix."""

    source: Label
    target: Label
    first_particle: Label
    second_particle: Label
    input_middles: tuple[Label, ...]
    output_middles: tuple[Label, ...]
    matrix: np.ndarray


@dataclass(frozen=True, slots=True)
class PathSwap:
    """One adjacent exchange acting on a complete left-associated path basis."""

    input_sequence: tuple[Label, ...]
    output_sequence: tuple[Label, ...]
    input_paths: tuple[tuple[Label, ...], ...]
    output_paths: tuple[tuple[Label, ...], ...]
    matrix: np.ndarray


class FundamentalRecouplings:
    """Terminal recouplings and their recursive action on longer paths."""

    def __init__(self, couplings: FundamentalCouplings) -> None:
        self.couplings = couplings

    @lru_cache(maxsize=None)
    def block(
        self,
        source: Label,
        first_particle: Label,
        second_particle: Label,
        target: Label,
    ) -> ElementarySwapBlock:
        if first_particle > second_particle:
            reverse = self.block(
                source, second_particle, first_particle, target
            )
            return ElementarySwapBlock(
                source=source,
                target=target,
                first_particle=first_particle,
                second_particle=second_particle,
                input_middles=reverse.output_middles,
                output_middles=reverse.input_middles,
                matrix=reverse.matrix.T,
            )

        inputs = self._two_step_options(
            source, first_particle, second_particle, target
        )
        outputs = self._two_step_options(
            source, second_particle, first_particle, target
        )
        if not inputs or not outputs or len(inputs) != len(outputs):
            raise ArithmeticError(
                "elementary recoupling path mismatch for "
                f"{source} x {first_particle} x {second_particle} -> {target}"
            )

        source_dimension = self.couplings.irreps.get(source).dim
        target_irrep = self.couplings.irreps.get(target)
        highest_candidates = [
            index
            for index, weight in enumerate(target_irrep.weights)
            if weight == target
        ]
        if len(highest_candidates) != 1:
            raise ArithmeticError(
                f"{target} has {len(highest_candidates)} highest states"
            )
        highest = highest_candidates[0]
        first_dimension = self.couplings.irreps.get(first_particle).dim
        second_dimension = self.couplings.irreps.get(second_particle).dim
        # The two trees are validated intertwiners from the same irreducible
        # target, so their overlap is a scalar by Schur's lemma.  Evaluate and
        # reconstruct it on the normalized highest-weight column instead of
        # carrying every descendant and amplifying long-lowering roundoff.
        original = np.empty(
            (
                len(inputs),
                source_dimension,
                first_dimension,
                second_dimension,
            ),
            dtype=float,
        )
        for index, (_, first, second) in enumerate(inputs):
            original[index] = self._compose_highest(first, second, highest)
        if first_particle == second_particle:
            exchanged = np.swapaxes(original, 2, 3)
        else:
            exchanged = np.empty_like(original)
            for index, (_, first, second) in enumerate(outputs):
                exchanged[index] = np.swapaxes(
                    self._compose_highest(first, second, highest), 1, 2
                )
        original_flat = original.reshape(len(inputs), -1)
        exchanged_flat = exchanged.reshape(len(outputs), -1)
        matrix = exchanged_flat @ original_flat.T
        SwapTableBuilder._clean_matrix(matrix)
        self._validate_terminal(
            matrix,
            original_flat,
            exchanged_flat,
            source,
            target,
        )
        return ElementarySwapBlock(
            source=source,
            target=target,
            first_particle=first_particle,
            second_particle=second_particle,
            input_middles=tuple(item[0] for item in inputs),
            output_middles=tuple(item[0] for item in outputs),
            matrix=matrix,
        )

    def _two_step_options(
        self,
        source: Label,
        first_particle: Label,
        second_particle: Label,
        target: Label,
    ) -> tuple[
        tuple[Label, ElementaryCoupling, ElementaryCoupling], ...
    ]:
        options: list[tuple[Label, ElementaryCoupling, ElementaryCoupling]] = []
        for first in self.couplings.decompose(source, first_particle).couplings:
            for second in self.couplings.decompose(
                first.target, second_particle
            ).for_target(target):
                options.append((first.target, first, second))
        return tuple(sorted(options, key=lambda item: item[0]))

    def _compose_highest(
        self,
        first: ElementaryCoupling,
        second: ElementaryCoupling,
        highest: int,
    ) -> np.ndarray:
        return _compose_elementary_highest(first, second, highest)

    @staticmethod
    def _validate_terminal(
        matrix: np.ndarray,
        original_flat: np.ndarray,
        exchanged_flat: np.ndarray,
        source: Label,
        target: Label,
    ) -> None:
        _require_finite(
            "terminal recoupling data", matrix, original_flat, exchanged_flat
        )
        size = matrix.shape[0]
        identity = np.eye(size)
        original_gram = original_flat @ original_flat.T
        exchanged_gram = exchanged_flat @ exchanged_flat.T
        reconstructed = matrix.T @ exchanged_flat
        column_gram = matrix.T @ matrix
        row_gram = matrix @ matrix.T
        _require_finite(
            "terminal recoupling residual",
            original_gram,
            exchanged_gram,
            reconstructed,
            column_gram,
            row_gram,
        )
        worst = max(
            float(np.max(np.abs(original_gram - identity))),
            float(np.max(np.abs(exchanged_gram - identity))),
            float(np.max(np.abs(column_gram - identity))),
            float(np.max(np.abs(row_gram - identity))),
            float(np.max(np.abs(original_flat - reconstructed))),
        )
        if worst > 3.0e-8 * max(1, size):
            raise ArithmeticError(
                f"terminal recoupling failed for {source}->{target}: "
                f"residual {worst:.3e}"
            )

    @lru_cache(maxsize=None)
    def _fusion_paths_by_target(
        self,
        source: Label,
        sequence: tuple[Label, ...],
    ) -> dict[Label, tuple[tuple[Label, ...], ...]]:
        paths: tuple[tuple[Label, ...], ...] = ((source,),)
        for particle in sequence:
            extended: list[tuple[Label, ...]] = []
            for path in paths:
                for coupling in self.couplings.decompose(
                    path[-1], particle
                ).couplings:
                    extended.append(path + (coupling.target,))
            paths = tuple(extended)
        grouped: dict[Label, list[tuple[Label, ...]]] = {}
        for path in paths:
            grouped.setdefault(path[-1], []).append(path)
        return {
            target: tuple(sorted(target_paths))
            for target, target_paths in grouped.items()
        }

    def fusion_paths(
        self,
        source: Label,
        sequence: tuple[Label, ...],
        target: Label,
    ) -> tuple[tuple[Label, ...], ...]:
        return self._fusion_paths_by_target(source, sequence).get(target, ())

    def adjacent_swap(
        self,
        source: Label,
        target: Label,
        sequence: tuple[Label, ...],
        position: int,
    ) -> PathSwap:
        if position < 0 or position + 1 >= len(sequence):
            raise IndexError("adjacent-swap position is outside the particle sequence")
        output_sequence = list(sequence)
        output_sequence[position], output_sequence[position + 1] = (
            output_sequence[position + 1],
            output_sequence[position],
        )
        output_sequence_tuple = tuple(output_sequence)
        input_paths = self.fusion_paths(source, sequence, target)
        output_paths = self.fusion_paths(source, output_sequence_tuple, target)
        if len(input_paths) != len(output_paths):
            raise ArithmeticError("adjacent exchange changed fusion-space dimension")
        output_lookup = {path: index for index, path in enumerate(output_paths)}
        matrix = np.zeros((len(output_paths), len(input_paths)), dtype=float)

        for input_index, path in enumerate(input_paths):
            local = self.block(
                path[position],
                sequence[position],
                sequence[position + 1],
                path[position + 2],
            )
            local_input = local.input_middles.index(path[position + 1])
            for local_output, middle in enumerate(local.output_middles):
                changed = list(path)
                changed[position + 1] = middle
                output_index = output_lookup[tuple(changed)]
                matrix[output_index, input_index] = local.matrix[
                    local_output, local_input
                ]

        _require_finite("path-level elementary exchange", matrix)
        identity = np.eye(len(input_paths))
        input_gram = matrix.T @ matrix
        output_gram = matrix @ matrix.T
        _require_finite(
            "path-level elementary exchange residual", input_gram, output_gram
        )
        worst = max(
            float(np.max(np.abs(input_gram - identity))),
            float(np.max(np.abs(output_gram - identity))),
        )
        if worst > 3.0e-8 * max(1, len(input_paths)):
            raise ArithmeticError(
                f"path-level elementary exchange residual {worst:.3e}"
            )
        return PathSwap(
            input_sequence=sequence,
            output_sequence=output_sequence_tuple,
            input_paths=input_paths,
            output_paths=output_paths,
            matrix=matrix,
        )

    def path_tensor(
        self,
        path: tuple[Label, ...],
        sequence: tuple[Label, ...],
    ) -> np.ndarray:
        """Construct one elementary fusion tree for an independent audit."""

        if len(path) != len(sequence) + 1:
            raise ValueError("fusion path and particle sequence lengths disagree")
        first_options = self.couplings.decompose(path[0], sequence[0]).for_target(
            path[1]
        )
        if len(first_options) != 1:
            raise ArithmeticError("elementary path edge is not multiplicity-free")
        source_dimension = self.couplings.irreps.get(path[0]).dim
        tensor = first_options[0]._source_slab(0, source_dimension)
        for edge in range(1, len(sequence)):
            options = self.couplings.decompose(
                path[edge], sequence[edge]
            ).for_target(path[edge + 1])
            if len(options) != 1:
                raise ArithmeticError("elementary path edge is not multiplicity-free")
            next_tensor = options[0]._source_slab(
                0, self.couplings.irreps.get(path[edge]).dim
            )
            tensor = np.tensordot(tensor, next_tensor, axes=([-1], [0]))
        return tensor

    def pair_swap_direct_residual(self, source: Label, target: Label) -> float:
        """Compare the recursive four-move product with a literal index swap."""

        paths, recursive = self.adjoint_pair_swap(source, target)
        chains = np.stack(
            [self.path_tensor(path, ADJOINT_SPLIT_SEQUENCE) for path in paths]
        )
        # Output trees see (q2,b2,q1,b1).  Reorder their explicit axes back
        # to the input order (q1,b1,q2,b2) before taking the overlap.
        exchanged = np.transpose(chains, (0, 1, 4, 5, 2, 3, 6))
        target_dimension = self.couplings.irreps.get(target).dim
        direct = np.einsum(
            "orabcdt,irabcdt->oi", exchanged, chains, optimize=True
        ) / target_dimension
        reconstruction = np.einsum(
            "oi,orabcdt->irabcdt", recursive, exchanged, optimize=True
        )
        return max(
            float(np.max(np.abs(recursive - direct))),
            float(np.max(np.abs(chains - reconstruction))),
        )

    def line_across_adjoint(
        self,
        source: Label,
        target: Label,
        sequence: tuple[Label, ...],
        pair_position: int,
    ) -> PathSwap:
        """Move one elementary line left across one split adjoint pair."""

        if sequence[pair_position : pair_position + 2] != (
            FUNDAMENTAL,
            ANTIFUNDAMENTAL,
        ):
            raise ValueError("selected positions do not contain a split adjoint")
        if pair_position + 2 >= len(sequence):
            raise ValueError("there is no elementary line to move across the pair")
        input_sequence = sequence
        input_paths = self.fusion_paths(source, input_sequence, target)
        total = np.eye(len(input_paths))
        current_sequence = input_sequence
        current_paths = input_paths
        # (q,b,x) -> (q,x,b) -> (x,q,b).
        for position in (pair_position + 1, pair_position):
            step = self.adjacent_swap(
                source, target, current_sequence, position
            )
            if step.input_paths != current_paths:
                raise ArithmeticError("one-adjoint path ordering changed")
            total = step.matrix @ total
            current_sequence = step.output_sequence
            current_paths = step.output_paths
        return PathSwap(
            input_sequence=input_sequence,
            output_sequence=current_sequence,
            input_paths=input_paths,
            output_paths=current_paths,
            matrix=total,
        )

    def adjoint_pair_swap(
        self, source: Label, target: Label
    ) -> tuple[tuple[tuple[Label, ...], ...], np.ndarray]:
        """Exchange two ``(3,3-bar)`` pairs using four terminal moves."""

        sequence = ADJOINT_SPLIT_SEQUENCE
        paths = self.fusion_paths(source, sequence, target)
        # First move q2 across (q1,b1), then move b2 across the same pair.
        # Each cached move is a one-adjoint recoupling made from exactly two
        # terminal two-fundamental recouplings.
        first = self.line_across_adjoint(source, target, sequence, 0)
        second = self.line_across_adjoint(
            source, target, first.output_sequence, 1
        )
        total = second.matrix @ first.matrix
        current_sequence = second.output_sequence
        current_paths = second.output_paths
        if current_sequence != sequence or current_paths != paths:
            raise ArithmeticError("adjoint pair permutation ended in the wrong basis")
        SwapTableBuilder._clean_matrix(total)
        _require_finite("four-line pair permutation", total)
        identity = np.eye(len(paths))
        squared = total @ total
        _require_finite("four-line pair permutation residual", squared)
        worst = max(
            float(np.max(np.abs(total - total.T))),
            float(np.max(np.abs(squared - identity))),
        )
        if worst > 5.0e-8 * max(1, len(paths)):
            raise ArithmeticError(
                f"four-line pair permutation residual {worst:.3e}"
            )
        return paths, total


@dataclass(frozen=True, slots=True)
class VertexReduction:
    """Expansion of one fixed adjoint vertex in ``3``--``3-bar`` paths."""

    intermediates: tuple[Label, ...]
    coefficients: np.ndarray


class AdjointVertexReductions:
    """Derive the convention-preserving adjoint vertex expansion."""

    def __init__(
        self,
        adjoint_couplings: AdjointCouplings,
        fundamental_couplings: FundamentalCouplings,
    ) -> None:
        self.adjoint_couplings = adjoint_couplings
        self.fundamental_couplings = fundamental_couplings
        split_vertices = fundamental_couplings.decompose(
            FUNDAMENTAL, ANTIFUNDAMENTAL
        )
        adjoint_vertices = split_vertices.for_target((1, 1))
        singlet_vertices = split_vertices.for_target((0, 0))
        if len(adjoint_vertices) != 1 or len(singlet_vertices) != 1:
            raise ArithmeticError("3 x 3-bar did not resolve as 1 plus 8")
        # The elementary coupling is an isometry 8 -> 3 tensor 3-bar in
        # the same abstract particle bases used by every terminal move.
        self.splitting = np.transpose(
            adjoint_vertices[0]._source_slab(0, 3), (2, 0, 1)
        )
        flattened = self.splitting.reshape(8, 9)
        singlet = singlet_vertices[0]._target_column(0).reshape(-1)
        _require_finite("recursive 3 x 3-bar resolution", flattened, singlet)
        resolution_error = max(
            float(np.max(np.abs(flattened @ flattened.T - np.eye(8)))),
            float(
                np.max(
                    np.abs(
                        flattened.T @ flattened
                        + np.outer(singlet, singlet)
                        - np.eye(9)
                    )
                )
            ),
        )
        if resolution_error > 2.0e-12:
            raise ArithmeticError(
                "recursive 3 x 3-bar resolution failed: "
                f"{resolution_error:.3e}"
            )

    def reduce(
        self, source: Label, target: Label, multiplicity: int
    ) -> VertexReduction:
        matches = [
            reduction
            for candidate_multiplicity, reduction in self._reduce_all(source, target)
            if candidate_multiplicity == multiplicity
        ]
        if len(matches) != 1:
            raise ArithmeticError(
                f"expected one adjoint vertex {source}->{target}[{multiplicity}]"
            )
        return matches[0]

    @lru_cache(maxsize=None)
    def _reduce_all(
        self, source: Label, target: Label
    ) -> tuple[tuple[int, VertexReduction], ...]:
        adjoint_vertices = self.adjoint_couplings.decompose(source).for_target(
            target
        )
        if not adjoint_vertices:
            raise ArithmeticError(f"no adjoint vertex {source}->{target}")
        options: list[
            tuple[Label, ElementaryCoupling, ElementaryCoupling]
        ] = []
        for first in self.fundamental_couplings.decompose(
            source, FUNDAMENTAL
        ).couplings:
            for second in self.fundamental_couplings.decompose(
                first.target, ANTIFUNDAMENTAL
            ).for_target(target):
                options.append((first.target, first, second))
        options.sort(key=lambda item: item[0])
        if not options:
            raise ArithmeticError(
                f"no fundamental expansion paths for adjoint vertex {source}->{target}"
            )

        source_dimension = self.adjoint_couplings.irreps.get(source).dim
        target_irrep = self.adjoint_couplings.irreps.get(target)
        # Each path and adjoint vertex is an already-validated intertwiner from
        # the same irreducible target.  Schur's lemma therefore makes their
        # overlap a scalar on the target space, so one normalized highest-weight
        # column determines it exactly.  This avoids amplifying roundoff from
        # long lowering chains in otherwise equivalent descendant columns.
        highest_candidates = [
            index
            for index, weight in enumerate(target_irrep.weights)
            if weight == target
        ]
        if len(highest_candidates) != 1:
            raise ArithmeticError(
                f"{target} has {len(highest_candidates)} highest states"
            )
        highest = highest_candidates[0]
        chains = np.empty(
            (len(options), source_dimension, 3, 3),
            dtype=float,
        )
        for index, (_, first, second) in enumerate(options):
            chains[index] = self._compose_fundamental_highest(
                first, second, highest
            )
        _require_finite("adjoint vertex path expansion", chains)
        chains_flat = chains.reshape(len(options), -1)
        adjoint_dimension = self.splitting.shape[0]
        path_gram = chains_flat @ chains_flat.T
        _require_finite("adjoint vertex path Gram matrix", path_gram)
        path_gram_error = float(
            np.max(np.abs(path_gram - np.eye(len(options))))
        )
        reductions: list[tuple[int, VertexReduction]] = []
        for adjoint_vertex in adjoint_vertices:
            adjoint_column = adjoint_vertex._target_column(highest)
            embedded = (
                adjoint_column
                @ self.splitting.reshape(adjoint_dimension, 9)
            ).reshape(source_dimension, 3, 3)
            embedded_flat = embedded.reshape(-1)
            coefficients = chains_flat @ embedded_flat
            reconstructed = coefficients @ chains_flat
            _require_finite(
                "adjoint vertex reduction",
                adjoint_column,
                embedded,
                coefficients,
                reconstructed,
            )
            coefficient_norm = float(np.dot(coefficients, coefficients))
            if not np.isfinite(coefficient_norm):
                raise ArithmeticError(
                    "adjoint vertex reduction norm is non-finite"
                )
            worst = max(
                path_gram_error,
                abs(coefficient_norm - 1.0),
                float(np.max(np.abs(embedded_flat - reconstructed))),
            )
            if worst > 3.0e-8 * max(1, len(options)):
                raise ArithmeticError(
                    f"adjoint vertex reduction failed for {source}->{target}"
                    f"[{adjoint_vertex.multiplicity}]: residual {worst:.3e}"
                )
            coefficients[np.abs(coefficients) < 5.0e-14] = 0.0
            reductions.append(
                (
                    adjoint_vertex.multiplicity,
                    VertexReduction(
                        intermediates=tuple(item[0] for item in options),
                        coefficients=coefficients,
                    ),
                )
            )
        return tuple(reductions)

    def _compose_fundamental_highest(
        self,
        first: ElementaryCoupling,
        second: ElementaryCoupling,
        highest: int,
    ) -> np.ndarray:
        return _compose_elementary_highest(first, second, highest)


class RecursiveReductionSwapTableBuilder(SwapTableBuilder):
    """Build adjoint swap blocks entirely from lower-complexity recouplings."""

    def __init__(self, couplings: AdjointCouplings | None = None) -> None:
        super().__init__(couplings)
        self.fundamental_couplings = FundamentalCouplings(self.couplings.irreps)
        self.fundamental_recouplings = FundamentalRecouplings(
            self.fundamental_couplings
        )
        self.vertex_reductions = AdjointVertexReductions(
            self.couplings, self.fundamental_couplings
        )

    @lru_cache(maxsize=None)
    def blocks_for_left(self, left: Label) -> tuple[SwapBlock, ...]:
        path_data = self._path_couplings_for_left(left)
        blocks: list[SwapBlock] = []
        for right in sorted(path_data):
            ordered = sorted(path_data[right], key=lambda item: item[0])
            paths = tuple(item[0] for item in ordered)
            fundamental_paths, pair_swap = (
                self.fundamental_recouplings.adjoint_pair_swap(left, right)
            )
            lookup = {
                path: index for index, path in enumerate(fundamental_paths)
            }
            expansion = np.zeros((len(fundamental_paths), len(paths)), dtype=float)
            for column, (_, first, second) in enumerate(ordered):
                first_reduction = self.vertex_reductions.reduce(
                    first.source, first.target, first.multiplicity
                )
                second_reduction = self.vertex_reductions.reduce(
                    second.source, second.target, second.multiplicity
                )
                for first_middle, first_coefficient in zip(
                    first_reduction.intermediates,
                    first_reduction.coefficients,
                    strict=True,
                ):
                    for second_middle, second_coefficient in zip(
                        second_reduction.intermediates,
                        second_reduction.coefficients,
                        strict=True,
                    ):
                        fundamental_path = (
                            left,
                            first_middle,
                            first.target,
                            second_middle,
                            right,
                        )
                        expansion[lookup[fundamental_path], column] += (
                            first_coefficient * second_coefficient
                        )

            _require_finite("recursive adjoint expansion", expansion, pair_swap)
            expansion_gram = expansion.T @ expansion
            _require_finite(
                "recursive adjoint expansion Gram matrix", expansion_gram
            )
            expansion_error = float(
                np.max(np.abs(expansion_gram - np.eye(len(paths))))
            )
            transformed = pair_swap @ expansion
            matrix = expansion.T @ transformed
            self._clean_matrix(matrix)
            _require_finite("recursive adjoint projection", transformed, matrix)
            reconstructed = expansion @ matrix
            _require_finite("recursive adjoint reconstruction", reconstructed)
            closure_error = float(
                np.max(np.abs(transformed - reconstructed))
            )
            if max(expansion_error, closure_error) > 5.0e-8 * max(1, len(paths)):
                raise ArithmeticError(
                    f"recursive adjoint reduction failed for {left}->{right}: "
                    f"basis={expansion_error:.3e}, closure={closure_error:.3e}"
                )
            block = SwapBlock(left=left, right=right, paths=paths, matrix=matrix)
            self.validate_block(block)
            blocks.append(block)
        return tuple(blocks)
