"""Unitary SU(3) irreps built directly as traceless symmetric tensors.

No tabulated Clebsch--Gordan data enter here.  The irrep ``(p,q)`` is the
kernel of contraction in

    Sym^p(3) tensor Sym^q(3-bar).

Normalized occupation-number states make the tensor-product scalar product
Euclidean.  This gives all square roots and signs of the Lie algebra action
without a convention imported from an external table.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import sqrt

import numpy as np

from .linalg import LoweringRecipe, canonical_nullspace, descendant_basis
from .packed_embedding import PackedEmbedding

Label = tuple[int, int]
Weight = tuple[int, int]
Occupation = tuple[int, int, int, int, int, int]


def dimension(label: Label) -> int:
    """Dimension of the SU(3) irrep with Dynkin labels ``(p,q)``."""

    p, q = label
    if p < 0 or q < 0:
        raise ValueError("Dynkin labels must be non-negative")
    return (p + 1) * (q + 1) * (p + q + 2) // 2


@lru_cache(maxsize=None)
def _compositions3(total: int) -> tuple[tuple[int, int, int], ...]:
    # Highest fundamental weights occur early in this ordering.  Nothing in
    # the mathematics depends on it, but it makes the gauge easy to inspect.
    return tuple(
        (first, second, total - first - second)
        for first in range(total, -1, -1)
        for second in range(total - first, -1, -1)
    )


@lru_cache(maxsize=None)
def _ambient_basis(p: int, q: int) -> tuple[Occupation, ...]:
    return tuple(
        upper + lower
        for upper in _compositions3(p)
        for lower in reversed(_compositions3(q))
    )


def _weight(state: Occupation) -> Weight:
    u0, u1, u2, l0, l1, l2 = state
    return ((u0 - l0) - (u1 - l1), (u1 - l1) - (u2 - l2))


def _ambient_generator_action(
    basis: tuple[Occupation, ...],
    lookup: dict[Occupation, int],
    row: int,
    col: int,
    vectors: np.ndarray,
) -> np.ndarray:
    """Apply ``E[row,col]`` to several ambient vectors at once."""

    size = len(basis)
    result = np.zeros((size, vectors.shape[1]), dtype=float)
    for source, state in enumerate(basis):
        upper = list(state[:3])
        lower = list(state[3:])
        if row == col:
            result[source] = (upper[row] - lower[row]) * vectors[source]
            continue
        if upper[col]:
            target = upper.copy()
            coefficient = sqrt((target[row] + 1) * target[col])
            target[row] += 1
            target[col] -= 1
            new_state = tuple(target) + tuple(lower)
            result[lookup[new_state]] += coefficient * vectors[source]
        if lower[row]:
            target = lower.copy()
            coefficient = -sqrt((target[col] + 1) * target[row])
            target[col] += 1
            target[row] -= 1
            new_state = tuple(upper) + tuple(target)
            result[lookup[new_state]] += coefficient * vectors[source]
    return result


def _harmonic_basis(p: int, q: int) -> tuple[np.ndarray, tuple[Weight, ...]]:
    """Canonical orthonormal basis of the contraction kernel."""

    source_basis = _ambient_basis(p, q)
    source_by_weight: dict[Weight, list[int]] = {}
    source_weights: list[Weight] = []
    source_local = np.empty(len(source_basis), dtype=int)
    for index, state in enumerate(source_basis):
        weight = _weight(state)
        source_weights.append(weight)
        indices = source_by_weight.setdefault(weight, [])
        source_local[index] = len(indices)
        indices.append(index)

    if p == 0 or q == 0:
        target_basis: tuple[Occupation, ...] = ()
        target_by_weight: dict[Weight, list[int]] = {}
        target_local = np.empty(0, dtype=int)
    else:
        target_basis = _ambient_basis(p - 1, q - 1)
        target_by_weight = {}
        target_local = np.empty(len(target_basis), dtype=int)
        for index, state in enumerate(target_basis):
            indices = target_by_weight.setdefault(_weight(state), [])
            target_local[index] = len(indices)
            indices.append(index)

    target_lookup = {state: index for index, state in enumerate(target_basis)}
    contraction_by_weight = {
        weight: np.zeros(
            (len(target_by_weight.get(weight, ())), len(source_indices)),
            dtype=float,
        )
        for weight, source_indices in source_by_weight.items()
    }
    if target_basis:
        for source, state in enumerate(source_basis):
            upper = list(state[:3])
            lower = list(state[3:])
            block = contraction_by_weight[source_weights[source]]
            for color in range(3):
                if upper[color] and lower[color]:
                    new_upper = upper.copy()
                    new_lower = lower.copy()
                    coefficient = sqrt(new_upper[color] * new_lower[color])
                    new_upper[color] -= 1
                    new_lower[color] -= 1
                    target = tuple(new_upper) + tuple(new_lower)
                    target_index = target_lookup[target]
                    block[
                        target_local[target_index], source_local[source]
                    ] += coefficient

    expected = dimension((p, q))
    harmonic = np.zeros((len(source_basis), expected), dtype=float)
    column = 0
    weights: list[Weight] = []
    residual_squared = 0.0
    # Descending weights put the highest-weight sector first and make the
    # resulting convention human-readable.
    for weight in sorted(source_by_weight, reverse=True):
        source_indices = source_by_weight[weight]
        block = contraction_by_weight[weight]
        local_kernel = canonical_nullspace(block)
        local_columns = local_kernel.shape[1]
        harmonic[source_indices, column : column + local_columns] = local_kernel
        column += local_columns
        weights.extend([weight] * local_columns)
        residual_squared += float(np.linalg.norm(block @ local_kernel) ** 2)

    if column != expected:
        raise ArithmeticError(
            f"ker(contraction) for ({p},{q}) has dimension {column}, "
            f"expected {expected}"
        )
    residual = sqrt(residual_squared)
    if residual > 2.0e-9:
        raise ArithmeticError(f"tracelessness residual is {residual:.3e}")
    return harmonic, tuple(weights)


class _IrrepStorage:
    """Private slots kept out of the public dataclass field contract."""

    __slots__ = (
        "_ambient_states_cache",
        "_ambient_transform",
        "_ambient_embedding_cache",
    )
    _ambient_states_cache: tuple[Occupation, ...] | None
    _ambient_transform: PackedEmbedding | None
    _ambient_embedding_cache: np.ndarray | None


@dataclass(frozen=True, slots=True)
class Irrep(_IrrepStorage):
    """A concrete, real, orthonormal realization of an SU(3) irrep."""

    label: Label
    generators: tuple[tuple[np.ndarray, ...], ...]
    weights: tuple[Weight, ...]
    lowering_recipe: LoweringRecipe
    ambient_states: tuple[Occupation, ...]
    ambient_embedding: np.ndarray

    def __post_init__(self) -> None:
        """Initialize private storage for a public dense construction."""

        object.__setattr__(
            self,
            "_ambient_states_cache",
            object.__getattribute__(self, "ambient_states"),
        )
        object.__setattr__(self, "_ambient_transform", None)
        object.__setattr__(
            self,
            "_ambient_embedding_cache",
            object.__getattribute__(self, "ambient_embedding"),
        )

    def __getattribute__(self, name: str):
        """Materialize lazy public dataclass fields on their first access."""

        if name == "ambient_states":
            try:
                cached = object.__getattribute__(self, "_ambient_states_cache")
            except AttributeError:
                return object.__getattribute__(self, "ambient_states")
            if cached is None:
                transform = object.__getattribute__(self, "_ambient_transform")
                if transform is None:
                    return object.__getattribute__(self, "ambient_states")
                label = object.__getattribute__(self, "label")
                cached = _ambient_basis(label[0], label[1])
                object.__setattr__(self, "ambient_states", cached)
                object.__setattr__(self, "_ambient_states_cache", cached)
            return cached
        if name == "ambient_embedding":
            try:
                cached = object.__getattribute__(self, "_ambient_embedding_cache")
            except AttributeError:
                return object.__getattribute__(self, "ambient_embedding")
            if cached is None:
                transform = object.__getattribute__(self, "_ambient_transform")
                if transform is None:
                    return object.__getattribute__(self, "ambient_embedding")
                label = object.__getattribute__(self, "label")
                harmonic, _weights = _harmonic_basis(label[0], label[1])
                cached = harmonic @ transform.to_dense()
                object.__setattr__(self, "ambient_embedding", cached)
                object.__setattr__(self, "_ambient_embedding_cache", cached)
            return cached
        return object.__getattribute__(self, name)

    def __getstate__(
        self,
    ) -> tuple[
        Label,
        tuple[tuple[np.ndarray, ...], ...],
        tuple[Weight, ...],
        LoweringRecipe,
        tuple[Occupation, ...],
        np.ndarray,
    ]:
        """Serialize exactly the original dense public dataclass state."""

        return (
            self.label,
            self.generators,
            self.weights,
            self.lowering_recipe,
            self.ambient_states,
            self.ambient_embedding,
        )

    def __setstate__(
        self,
        state: tuple[
            Label,
            tuple[tuple[np.ndarray, ...], ...],
            tuple[Weight, ...],
            LoweringRecipe,
            tuple[Occupation, ...],
            np.ndarray,
        ],
    ) -> None:
        (
            label,
            generators,
            weights,
            lowering_recipe,
            ambient_states,
            ambient_embedding,
        ) = state
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "generators", generators)
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "lowering_recipe", lowering_recipe)
        object.__setattr__(self, "ambient_states", ambient_states)
        object.__setattr__(self, "ambient_embedding", ambient_embedding)
        object.__setattr__(self, "_ambient_states_cache", ambient_states)
        object.__setattr__(self, "_ambient_transform", None)
        object.__setattr__(self, "_ambient_embedding_cache", ambient_embedding)

    @classmethod
    def _from_packed(
        cls,
        label: Label,
        generators: tuple[tuple[np.ndarray, ...], ...],
        weights: tuple[Weight, ...],
        lowering_recipe: LoweringRecipe,
        ambient_transform: PackedEmbedding,
        ambient_embedding: np.ndarray | None,
    ) -> Irrep:
        result = cls.__new__(cls)
        object.__setattr__(result, "label", label)
        object.__setattr__(result, "generators", generators)
        object.__setattr__(result, "weights", weights)
        object.__setattr__(result, "lowering_recipe", lowering_recipe)
        object.__setattr__(result, "ambient_states", None)
        object.__setattr__(result, "ambient_embedding", ambient_embedding)
        object.__setattr__(result, "_ambient_states_cache", None)
        object.__setattr__(result, "_ambient_transform", ambient_transform)
        object.__setattr__(result, "_ambient_embedding_cache", ambient_embedding)
        return result

    @property
    def p(self) -> int:
        return self.label[0]

    @property
    def q(self) -> int:
        return self.label[1]

    @property
    def dim(self) -> int:
        return len(self.weights)

    def generator(self, row: int, col: int) -> np.ndarray:
        return self.generators[row][col]


class IrrepFactory:
    """Cache and validate concrete ``(p,q)`` representations."""

    @lru_cache(maxsize=None)
    def get(self, label: Label) -> Irrep:
        p, q = label
        expected = dimension(label)
        ambient = _ambient_basis(p, q)
        ambient_lookup = {state: index for index, state in enumerate(ambient)}
        harmonic, harmonic_weights = _harmonic_basis(p, q)
        harmonic_generators = tuple(
            tuple(
                harmonic.T
                @ _ambient_generator_action(
                    ambient,
                    ambient_lookup,
                    row,
                    col,
                    harmonic,
                )
                for col in range(3)
            )
            for row in range(3)
        )

        highest_candidates = [
            index for index, weight in enumerate(harmonic_weights) if weight == label
        ]
        if len(highest_candidates) != 1:
            raise ArithmeticError(
                f"({p},{q}) has {len(highest_candidates)} canonical highest states"
            )
        highest = np.zeros(expected, dtype=float)
        highest[highest_candidates[0]] = 1.0
        descendants, weights, lowering_recipe = descendant_basis(
            harmonic_generators[1][0],
            harmonic_generators[2][1],
            highest,
            label,
            expected,
            coordinate_weights=harmonic_weights,
        )
        generators = tuple(
            tuple(
                descendants.T @ harmonic_generators[row][col] @ descendants
                for col in range(3)
            )
            for row in range(3)
        )
        ambient_transform = PackedEmbedding.from_dense(
            descendants,
            harmonic_weights,
            weights,
            particle_dimension=1,
        )
        result = Irrep._from_packed(
            label=label,
            generators=generators,
            weights=weights,
            lowering_recipe=lowering_recipe,
            ambient_transform=ambient_transform,
            ambient_embedding=(
                harmonic @ descendants if label == (1, 1) else None
            ),
        )
        self._validate(result)
        return result

    @staticmethod
    def _validate(irrep: Irrep) -> None:
        scale = max(1.0, float(irrep.dim))
        worst_adjoint = 0.0
        worst_commutator = 0.0
        indexed_generators = tuple(
            (row, col, irrep.generator(row, col))
            for row in range(3)
            for col in range(3)
        )
        for i in range(3):
            for j in range(3):
                worst_adjoint = max(
                    worst_adjoint,
                    float(
                        np.max(
                            np.abs(
                                irrep.generator(i, j).T - irrep.generator(j, i)
                            )
                        )
                    ),
                )
        # A commutator is antisymmetric in its two arguments, so checking
        # each unordered generator pair once covers exactly the same defining
        # relations without repeating all reverse-order products.
        for left_index, (i, j, left) in enumerate(indexed_generators):
            for k, ell, right in indexed_generators[left_index + 1 :]:
                residual = left @ right - right @ left
                if j == k:
                    residual -= irrep.generator(i, ell)
                if i == ell:
                    residual += irrep.generator(k, j)
                worst_commutator = max(
                    worst_commutator, float(np.max(np.abs(residual)))
                )
        if worst_adjoint > 2.0e-9 * scale:
            raise ArithmeticError(f"E_ij adjoint residual {worst_adjoint:.3e}")
        if worst_commutator > 2.0e-8 * scale:
            raise ArithmeticError(
                f"SU(3) commutator residual {worst_commutator:.3e}"
            )
