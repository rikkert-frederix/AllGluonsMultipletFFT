"""Clebsch--Gordan intertwiners for ``(p,q) tensor (1,1)``.

The decomposition and every coefficient are discovered from the Lie-algebra
action.  In particular, the familiar adjoint-product rule is *checked*, not
used as input: highest-weight vectors are the simultaneous kernel of the two
raising operators in each product-weight sector.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from operator import index
from typing import Sequence

import numpy as np

from .fundamental import adjoint_split_isometry
from .linalg import (
    canonical_nullspace,
    canonical_subspace_basis,
    descendant_basis,
)
from .packed_embedding import PackedEmbedding
from .representations import Irrep, IrrepFactory, Label


class _CouplingStorage:
    """Private slots kept out of the public dataclass field contract."""

    __slots__ = ("_packed_embedding", "_embedding_cache", "_particle_dimension")
    _packed_embedding: PackedEmbedding | None
    _embedding_cache: np.ndarray | None
    _particle_dimension: int


@dataclass(frozen=True, slots=True)
class Coupling(_CouplingStorage):
    """An isometric embedding ``target -> source tensor adjoint``."""

    source: Label
    target: Label
    multiplicity: int
    embedding: np.ndarray
    exchange_parity: int | None = None

    def __post_init__(self) -> None:
        """Initialize private storage for a public dense construction."""

        object.__setattr__(self, "_packed_embedding", None)
        object.__setattr__(
            self, "_embedding_cache", object.__getattribute__(self, "embedding")
        )
        object.__setattr__(self, "_particle_dimension", 8)

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

    def __getstate__(self) -> tuple[Label, Label, int, np.ndarray, int | None]:
        """Serialize exactly the original dense public dataclass state."""

        return (
            self.source,
            self.target,
            self.multiplicity,
            self.embedding,
            self.exchange_parity,
        )

    def __setstate__(
        self, state: tuple[Label, Label, int, np.ndarray, int | None]
    ) -> None:
        source, target, multiplicity, embedding, exchange_parity = state
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "multiplicity", multiplicity)
        object.__setattr__(self, "embedding", embedding)
        object.__setattr__(self, "exchange_parity", exchange_parity)
        object.__setattr__(self, "_packed_embedding", None)
        object.__setattr__(self, "_embedding_cache", embedding)
        object.__setattr__(self, "_particle_dimension", 8)

    @classmethod
    def _from_packed(
        cls,
        source: Label,
        target: Label,
        multiplicity: int,
        exchange_parity: int | None,
        packed: PackedEmbedding,
    ) -> Coupling:
        result = cls.__new__(cls)
        object.__setattr__(result, "source", source)
        object.__setattr__(result, "target", target)
        object.__setattr__(result, "multiplicity", multiplicity)
        object.__setattr__(result, "embedding", None)
        object.__setattr__(result, "exchange_parity", exchange_parity)
        object.__setattr__(result, "_packed_embedding", packed)
        object.__setattr__(result, "_embedding_cache", None)
        object.__setattr__(result, "_particle_dimension", packed.particle_dimension)
        return result

    def _pack(
        self,
        product_weights: Sequence[Label],
        target_weights: Sequence[Label],
    ) -> Coupling:
        packed = PackedEmbedding.from_dense(
            self.embedding,
            product_weights,
            target_weights,
            particle_dimension=self._particle_dimension,
        )
        return self._from_packed(
            self.source,
            self.target,
            self.multiplicity,
            self.exchange_parity,
            packed,
        )

    def _source_slab(self, start: int, stop: int) -> np.ndarray:
        """Return ``source[start:stop] x adjoint x target`` coefficients."""

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
        """Return selected columns as ``source x adjoint x width``."""

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
        """Return one column as ``source x adjoint`` coefficients."""

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
        """Return nonzero flat product rows and values for one column."""

        packed = self._packed_embedding
        if packed is not None and self._embedding_cache is None:
            rows, values = packed.column_entries(column)
            nonzero = values != 0
            selected_rows = np.ascontiguousarray(rows[nonzero], dtype=np.intp)
            selected_values = np.ascontiguousarray(values[nonzero])
            selected_rows.setflags(write=False)
            selected_values.setflags(write=False)
            return selected_rows, selected_values
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


@dataclass(frozen=True, slots=True)
class ProductDecomposition:
    """Complete decomposition of one irrep tensored with the adjoint."""

    source: Label
    couplings: tuple[Coupling, ...]

    @property
    def targets(self) -> tuple[Label, ...]:
        return tuple(sorted({coupling.target for coupling in self.couplings}))

    def for_target(self, target: Label) -> tuple[Coupling, ...]:
        return tuple(c for c in self.couplings if c.target == target)


@dataclass(frozen=True, slots=True)
class _TensorProductAction:
    """Kronecker-sum action without materializing its dense matrix."""

    left: np.ndarray
    right: np.ndarray

    @property
    def left_dimension(self) -> int:
        return self.left.shape[0]

    @property
    def right_dimension(self) -> int:
        return self.right.shape[0]

    def __matmul__(self, value: np.ndarray) -> np.ndarray:
        array = np.asarray(value)
        left_dimension = self.left_dimension
        right_dimension = self.right_dimension
        product_dimension = left_dimension * right_dimension
        if array.shape[0] != product_dimension:
            raise ValueError(
                f"tensor-product action expected leading dimension "
                f"{product_dimension}, got {array.shape[0]}"
            )
        if array.ndim == 1:
            tensor = array.reshape(left_dimension, right_dimension)
            return (self.left @ tensor + tensor @ self.right.T).reshape(-1)
        if array.ndim != 2:
            raise ValueError("tensor-product action accepts vectors or matrices")
        columns = array.shape[1]
        tensor = array.reshape(left_dimension, right_dimension, columns)
        left_result = (
            self.left @ tensor.reshape(left_dimension, right_dimension * columns)
        ).reshape(left_dimension, right_dimension, columns)
        right_result = self.right @ tensor
        return (left_result + right_result).reshape(product_dimension, columns)

    def restricted(
        self, output_indices: Sequence[int], input_indices: Sequence[int]
    ) -> np.ndarray:
        """Return a weight-sector block of the represented operator."""

        outputs = np.asarray(output_indices, dtype=int)
        inputs = np.asarray(input_indices, dtype=int)
        if outputs.size == 0:
            return np.zeros((0, inputs.size), dtype=float)
        right_dimension = self.right_dimension
        output_left, output_right = np.divmod(outputs, right_dimension)
        input_left, input_right = np.divmod(inputs, right_dimension)
        result = self.left[np.ix_(output_left, input_left)] * (
            output_right[:, None] == input_right[None, :]
        )
        result += self.right[np.ix_(output_right, input_right)] * (
            output_left[:, None] == input_left[None, :]
        )
        return result

    def dense(self) -> np.ndarray:
        return np.kron(self.left, np.eye(self.right_dimension)) + np.kron(
            np.eye(self.left_dimension), self.right
        )


def _product_action(
    left: Irrep, right: Irrep, row: int, col: int
) -> _TensorProductAction:
    return _TensorProductAction(
        left.generator(row, col), right.generator(row, col)
    )


def _product_generator(left: Irrep, right: Irrep, row: int, col: int) -> np.ndarray:
    """Materialize a product generator for diagnostics outside hot paths."""

    return _product_action(left, right, row, col).dense()


def _raising_annihilation(
    sectors: dict[Label, list[int]],
    weight: Label,
    input_indices: Sequence[int],
    raising_1: _TensorProductAction,
    raising_2: _TensorProductAction,
) -> np.ndarray:
    """Stack only the nonzero weight-sector rows of both simple raisings."""

    target_1 = sectors.get((weight[0] + 2, weight[1] - 1), ())
    target_2 = sectors.get((weight[0] - 1, weight[1] + 2), ())
    return np.vstack(
        (
            raising_1.restricted(target_1, input_indices),
            raising_2.restricted(target_2, input_indices),
        )
    )


def _swap_on_weight_sector(
    indices: list[int], dimension: int
) -> np.ndarray:
    """Permutation of equal tensor factors, restricted to one weight sector."""

    position = {full_index: local for local, full_index in enumerate(indices)}
    result = np.zeros((len(indices), len(indices)), dtype=float)
    for source_local, source_full in enumerate(indices):
        left, right = divmod(source_full, dimension)
        target_full = right * dimension + left
        result[position[target_full], source_local] = 1.0
    return result


def _adjoint_matrix_basis(adjoint: Irrep) -> tuple[np.ndarray, ...]:
    """Interpret the harmonic ``(1,1)`` tensors as traceless 3x3 matrices."""

    return tuple(adjoint_split_isometry(adjoint))


def _normalize_embedding(embedding: np.ndarray) -> np.ndarray:
    target_dimension = embedding.shape[1]
    gram = embedding.T @ embedding
    normalization = float(np.trace(gram) / target_dimension)
    residual = float(np.max(np.abs(gram - normalization * np.eye(target_dimension))))
    if normalization <= 0.0 or residual > 2.0e-9 * max(1, target_dimension):
        raise ArithmeticError(
            f"canonical vertex has non-Schur norm: n={normalization:.3e}, "
            f"residual={residual:.3e}"
        )
    return embedding / np.sqrt(normalization)


def _action_embedding(
    source: Irrep,
    adjoint: Irrep,
    matrices: tuple[np.ndarray, ...] | None = None,
) -> np.ndarray:
    """Normalized adjoint of ``v tensor X -> rho_source(X) v``."""

    if matrices is None:
        matrices = _adjoint_matrix_basis(adjoint)
    projection = np.zeros((source.dim, source.dim * adjoint.dim), dtype=float)
    for adjoint_index, matrix in enumerate(matrices):
        action = sum(
            (
                matrix[row, col] * source.generator(row, col)
                for row in range(3)
                for col in range(3)
            ),
            start=np.zeros((source.dim, source.dim), dtype=float),
        )
        projection[:, adjoint_index :: adjoint.dim] = action
    return _normalize_embedding(projection.T)


def _adjoint_product_embedding(
    adjoint: Irrep,
    parity: int,
    matrices: tuple[np.ndarray, ...] | None = None,
) -> np.ndarray:
    """Anchor 8s to the Jordan product and 8a to the Lie bracket."""

    if matrices is None:
        matrices = _adjoint_matrix_basis(adjoint)
    projection = np.zeros((adjoint.dim, adjoint.dim * adjoint.dim), dtype=float)
    identity = np.eye(3)
    for left_index, left in enumerate(matrices):
        for right_index, right in enumerate(matrices):
            if parity == -1:
                product = left @ right - right @ left
            elif parity == 1:
                product = left @ right + right @ left
                product -= (2.0 / 3.0) * np.trace(left @ right) * identity
            else:
                raise ValueError("adjoint exchange parity must be +1 or -1")
            input_index = left_index * adjoint.dim + right_index
            for output_index, output in enumerate(matrices):
                projection[output_index, input_index] = np.sum(output * product)
    return _normalize_embedding(projection.T)


class AdjointCouplings:
    """Construct and cache all ``R tensor 8`` coupling vertices."""

    def __init__(self, irreps: IrrepFactory | None = None) -> None:
        self.irreps = irreps or IrrepFactory()
        self._matrix_basis: tuple[np.ndarray, ...] | None = None

    def _adjoint_matrix_basis(self, adjoint: Irrep) -> tuple[np.ndarray, ...]:
        if self._matrix_basis is None:
            self._matrix_basis = _adjoint_matrix_basis(adjoint)
        return self._matrix_basis

    @lru_cache(maxsize=None)
    def decompose(self, source_label: Label) -> ProductDecomposition:
        source = self.irreps.get(source_label)
        adjoint = self.irreps.get((1, 1))
        product_dimension = source.dim * adjoint.dim
        weights = tuple(
            (
                source_weight[0] + adjoint_weight[0],
                source_weight[1] + adjoint_weight[1],
            )
            for source_weight in source.weights
            for adjoint_weight in adjoint.weights
        )
        sectors: dict[Label, list[int]] = {}
        for index, weight in enumerate(weights):
            sectors.setdefault(weight, []).append(index)

        product_actions = tuple(
            tuple(
                _product_action(source, adjoint, row, col)
                for col in range(3)
            )
            for row in range(3)
        )
        raising_1 = product_actions[0][1]
        raising_2 = product_actions[1][2]

        highest_by_target: dict[Label, list[tuple[np.ndarray, int | None]]] = {}
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

            vectors_with_parity: list[tuple[np.ndarray, int | None]] = []
            if source.label == adjoint.label:
                swap = _swap_on_weight_sector(indices, source.dim)
                compressed_swap = kernel.T @ swap @ kernel
                eigenvalues, eigenvectors = np.linalg.eigh(compressed_swap)
                # Symmetric precedes antisymmetric.  This fixes the two octet
                # vertices as d-like (multiplicity 0) and f-like (1).
                for parity in (1, -1):
                    selected = np.flatnonzero(
                        np.abs(eigenvalues - parity) < 2.0e-8
                    )
                    if not selected.size:
                        continue
                    subspace = kernel @ eigenvectors[:, selected]
                    fixed = canonical_subspace_basis(subspace)
                    vectors_with_parity.extend(
                        (fixed[:, column], parity)
                        for column in range(fixed.shape[1])
                    )
            else:
                vectors_with_parity.extend(
                    (kernel[:, column], None)
                    for column in range(kernel.shape[1])
                )
            highest_by_target[weight] = vectors_with_parity

        packed_couplings: list[Coupling] = []
        accounted_dimension = 0
        for target_label in sorted(highest_by_target, reverse=True):
            target = self.irreps.get(target_label)
            target_couplings = self._decompose_target(
                source,
                adjoint,
                target,
                highest_by_target[target_label],
                sectors,
                weights,
                product_actions,
            )
            packed_couplings.extend(target_couplings)
            accounted_dimension += len(target_couplings) * target.dim

        if accounted_dimension != product_dimension:
            decomposition = ", ".join(
                f"{coupling.target}[{coupling.multiplicity}]"
                for coupling in packed_couplings
            )
            raise ArithmeticError(
                f"decomposition of {source_label} x (1,1) accounts for "
                f"{accounted_dimension}/{product_dimension}: {decomposition}"
            )
        return ProductDecomposition(source_label, tuple(packed_couplings))

    def _decompose_target(
        self,
        source: Irrep,
        adjoint: Irrep,
        target: Irrep,
        highest_vectors: Sequence[tuple[np.ndarray, int | None]],
        sectors: dict[Label, list[int]],
        product_weights: Sequence[Label],
        product_actions: tuple[tuple[_TensorProductAction, ...], ...],
    ) -> tuple[Coupling, ...]:
        """Construct, validate, and pack all copies of one target irrep.

        Dense embeddings are needed together only when anchoring and checking
        equivalent copies.  Keeping that lifetime inside this method releases
        one target group before :meth:`decompose` starts the next one.
        """

        lowering_1 = product_actions[1][0]
        lowering_2 = product_actions[2][1]
        dense_couplings: list[Coupling] = []
        for multiplicity, (local_highest, parity) in enumerate(highest_vectors):
            highest = np.zeros(len(product_weights), dtype=float)
            highest[sectors[target.label]] = local_highest
            descendants, descendant_weights, descendant_recipe = descendant_basis(
                lowering_1,
                lowering_2,
                highest,
                target.label,
                target.dim,
                recipe=target.lowering_recipe,
                coordinate_weights=product_weights,
                reference_lowerings=(
                    target.generator(1, 0),
                    target.generator(2, 1),
                ),
            )
            if descendant_weights != target.weights:
                raise ArithmeticError(
                    f"basis recursion mismatch for {source.label} x 8 -> "
                    f"{target.label}"
                )
            if descendant_recipe != target.lowering_recipe:
                raise ArithmeticError(
                    f"lowering recipe mismatch for {source.label} x 8 -> "
                    f"{target.label}"
                )
            dense_couplings.append(
                Coupling(
                    source=source.label,
                    target=target.label,
                    multiplicity=multiplicity,
                    embedding=descendants,
                    exchange_parity=parity,
                )
            )

        dense_couplings = self._anchor_canonical_vertices(
            source, adjoint, dense_couplings
        )
        for coupling in dense_couplings:
            self._validate_intertwiner(
                source,
                adjoint,
                target,
                coupling.embedding,
                product_actions,
            )
        self._validate_mutual_orthogonality(dense_couplings)
        return tuple(
            coupling._pack(product_weights, target.weights)
            for coupling in dense_couplings
        )

    def _anchor_canonical_vertices(
        self, source: Irrep, adjoint: Irrep, couplings: list[Coupling]
    ) -> list[Coupling]:
        """Tie otherwise arbitrary multiplicity phases to invariant products."""

        result = list(couplings)
        self_indices = [
            index
            for index, coupling in enumerate(result)
            if coupling.target == source.label
        ]
        if not self_indices:
            return result

        if source.label == adjoint.label:
            # Preserve the documented order 8s, 8a, but orient them relative
            # to matrix multiplication rather than an arbitrary coordinate.
            for index in self_indices:
                old = result[index]
                if old.exchange_parity is None:
                    raise ArithmeticError("8 x 8 -> 8 vertex lacks exchange parity")
                anchor = _adjoint_product_embedding(
                    adjoint,
                    old.exchange_parity,
                    self._adjoint_matrix_basis(adjoint),
                )
                overlap = float(np.sum(old.embedding * anchor) / adjoint.dim)
                if abs(abs(overlap) - 1.0) > 3.0e-8:
                    raise ArithmeticError(
                        f"8 x 8 vertex/product overlap is {overlap:.12g}"
                    )
                result[index] = Coupling(
                    source=old.source,
                    target=old.target,
                    multiplicity=old.multiplicity,
                    embedding=anchor,
                    exchange_parity=old.exchange_parity,
                )
            return result

        # The representation action selects a canonical copy of R in R x 8.
        # If there is a second copy, keep the orthogonal complement in the
        # already-fixed lexicographic gauge and orient its highest component
        # positively.
        action = _action_embedding(
            source, adjoint, self._adjoint_matrix_basis(adjoint)
        )
        old_copies = [result[index] for index in self_indices]
        target_dimension = source.dim
        coefficients = np.array(
            [
                np.sum(coupling.embedding * action) / target_dimension
                for coupling in old_copies
            ],
            dtype=float,
        )
        if abs(float(np.dot(coefficients, coefficients)) - 1.0) > 4.0e-8:
            raise ArithmeticError(
                f"generator-action vertex is incomplete in {source.label} x 8"
            )
        anchored = [action]
        if len(old_copies) == 2:
            complement = (
                -coefficients[1] * old_copies[0].embedding
                + coefficients[0] * old_copies[1].embedding
            )
            complement = _normalize_embedding(complement)
            pivot = np.flatnonzero(np.abs(complement[:, 0]) > 2.0e-11)
            if pivot.size and complement[pivot[0], 0] < 0.0:
                complement = -complement
            anchored.append(complement)
        elif len(old_copies) != 1:
            raise ArithmeticError(
                f"unexpected self multiplicity {len(old_copies)} for {source.label} x 8"
            )

        for multiplicity, (index, embedding) in enumerate(
            zip(self_indices, anchored, strict=True)
        ):
            old = result[index]
            result[index] = Coupling(
                source=old.source,
                target=old.target,
                multiplicity=multiplicity,
                embedding=embedding,
                exchange_parity=old.exchange_parity,
            )
        return result

    @staticmethod
    def _validate_intertwiner(
        source: Irrep,
        adjoint: Irrep,
        target: Irrep,
        embedding: np.ndarray,
        product_actions: tuple[tuple[_TensorProductAction, ...], ...],
    ) -> None:
        isometry_error = float(
            np.max(np.abs(embedding.T @ embedding - np.eye(target.dim)))
        )
        worst = 0.0
        # Off-diagonal E_ij and the two traceless Cartan generators span
        # su(3).  Individual E_ii must not be compared: tensor realizations
        # whose U(3) highest weights differ by a determinant column restrict
        # to the same SU(3) irrep (for example 3 x 8 contains 6-bar).
        for row in range(3):
            for col in range(3):
                if row == col:
                    continue
                product = product_actions[row][col]
                residual = product @ embedding - embedding @ target.generator(row, col)
                worst = max(worst, float(np.max(np.abs(residual))))
        for first, second in ((0, 1), (1, 2)):
            product = _TensorProductAction(
                source.generator(first, first)
                - source.generator(second, second),
                adjoint.generator(first, first)
                - adjoint.generator(second, second),
            )
            target_cartan = target.generator(first, first) - target.generator(
                second, second
            )
            residual = product @ embedding - embedding @ target_cartan
            worst = max(worst, float(np.max(np.abs(residual))))
        scale = max(1.0, float(target.dim))
        if isometry_error > 3.0e-9 * scale:
            raise ArithmeticError(f"CG isometry residual {isometry_error:.3e}")
        if worst > 3.0e-8 * scale:
            raise ArithmeticError(f"CG intertwining residual {worst:.3e}")

    @staticmethod
    def _validate_mutual_orthogonality(couplings: list[Coupling]) -> None:
        for left_index, left in enumerate(couplings):
            for right in couplings[left_index + 1 :]:
                if left.target != right.target:
                    continue
                overlap = float(np.max(np.abs(left.embedding.T @ right.embedding)))
                if overlap > 3.0e-8 * max(1, left.embedding.shape[1]):
                    raise ArithmeticError(
                        "equivalent CG copies are not orthogonal: "
                        f"{left.source} -> {left.target}, residual {overlap:.3e}"
                    )
