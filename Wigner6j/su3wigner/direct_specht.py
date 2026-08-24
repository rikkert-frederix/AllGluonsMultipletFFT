"""Direct, Clebsch--Gordan-free adjacent-adjoint swap construction.

This module supplies an independent check on the tensor and recursively
reduced backends.  It realizes an SU(3) irrep as a Hermitian Young projector
inside a tensor power of the fundamental.  The two external adjoints are two
copies of the ``(2, 1)`` Young projector.  Their exchange is evaluated once,
as the literal permutation of the two three-index blocks.

Only standard *skew* tableaux for the six newly added boxes are required.
The very large Specht space belonging to the final partition is never built.
No representation matrices for SU(3), Clebsch--Gordan embeddings, highest
weight kernels, elementary-line recouplings, or lower 6j symbols are used.

The repository vertex gauge is fixed locally in the same permutation algebra.
Three-box edge frames use the Jucys--Murphy action direction, its oriented
complement, and the special even/odd convention for ``8 x 8 -> 8``.  Tensoring
two such frames in the skew-tableau branching basis therefore fixes both the
magnitudes and signs of every block without consulting another backend.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .recoupling import Path, SwapBlock, SwapTableBuilder
from .representations import Label


ADJOINT: Label = (1, 1)
_EDGE_CACHE_SIZE = 512
_PATH_CACHE_SIZE = 256
_TABLEAU_CACHE_SIZE = 256
_SIGNED_BLOCK_CACHE_SIZE = 16
_DIAGNOSTIC_CACHE_SIZE = 256
_LOCAL_FRAME_CACHE_SIZE = 512


@dataclass(frozen=True, slots=True)
class SpechtDiagnostics:
    """Small-space dimensions used to construct one oracle block."""

    left_partition: tuple[int, int, int]
    final_partition: tuple[int, int, int]
    skew_dimension: int
    product_rank: int
    swap_signature: tuple[int, int]
    prefix_ranks: tuple[tuple[Label, int], ...]


@dataclass(frozen=True, slots=True)
class _Edge:
    source: Label
    target: Label
    multiplicity: int
    exchange_parity: int | None


@lru_cache(maxsize=_EDGE_CACHE_SIZE)
def _adjoint_edges(source: Label) -> tuple[_Edge, ...]:
    """The label rule for ``source tensor 8``, with no coupling tensors."""

    p, q = source
    targets: list[Label] = [(p + 1, q + 1)]
    if q:
        targets.append((p + 2, q - 1))
    if p:
        targets.append((p - 1, q + 2))
    targets.extend([source] * (int(p > 0) + int(q > 0)))
    if q >= 2:
        targets.append((p + 1, q - 2))
    if p >= 2:
        targets.append((p - 2, q + 1))
    if p and q:
        targets.append((p - 1, q - 1))

    multiplicities: dict[Label, int] = {}
    edges: list[_Edge] = []
    for target in targets:
        multiplicity = multiplicities.get(target, 0)
        multiplicities[target] = multiplicity + 1
        parity: int | None = None
        if source == ADJOINT:
            if target in ((0, 0), (2, 2)):
                parity = 1
            elif target in ((3, 0), (0, 3)):
                parity = -1
            elif target == ADJOINT:
                parity = 1 if multiplicity == 0 else -1
        edges.append(_Edge(source, target, multiplicity, parity))
    return tuple(sorted(edges, key=lambda edge: (edge.target, edge.multiplicity)))


@lru_cache(maxsize=_PATH_CACHE_SIZE)
def _paths(left: Label, right: Label) -> tuple[Path, ...]:
    result: list[Path] = []
    for first in _adjoint_edges(left):
        for second in _adjoint_edges(first.target):
            if second.target != right:
                continue
            result.append(
                Path(
                    middle=first.target,
                    left_multiplicity=first.multiplicity,
                    right_multiplicity=second.multiplicity,
                    left_exchange_parity=first.exchange_parity,
                    right_exchange_parity=second.exchange_parity,
                )
            )
    return tuple(sorted(result))


def _minimal_partition(label: Label) -> tuple[int, int, int]:
    p, q = label
    return p + q, q, 0


def _lifted_partition(label: Label, degree: int) -> tuple[int, int, int]:
    p, q = label
    remainder = degree - p - 2 * q
    if remainder < 0 or remainder % 3:
        raise ValueError(f"{label} does not occur at polynomial degree {degree}")
    determinant_columns = remainder // 3
    return (
        p + q + determinant_columns,
        q + determinant_columns,
        determinant_columns,
    )


@lru_cache(maxsize=_TABLEAU_CACHE_SIZE)
def _skew_tableaux(
    inner: tuple[int, int, int], outer: tuple[int, int, int]
) -> tuple[tuple[int, ...], ...]:
    """Rows of labels in standard tableaux of ``outer / inner``.

    Removal of the largest label reproduces the ordering of the repository's
    Young-orthogonal implementation while never enumerating the full Specht
    basis.
    """

    if any(inner[row] > outer[row] for row in range(3)):
        return ()
    if inner == outer:
        return ((),)
    result: list[tuple[int, ...]] = []
    for row in range(3):
        if outer[row] <= inner[row]:
            continue
        if row < 2 and outer[row] <= outer[row + 1]:
            continue
        reduced = list(outer)
        reduced[row] -= 1
        for tableau in _skew_tableaux(inner, tuple(reduced)):
            result.append(tableau + (row,))
    return tuple(result)


def skew_dimension(left: Label, right: Label) -> int:
    """Dimension of the exact six-box space for a fixed outer block."""

    inner = _minimal_partition(left)
    degree = sum(inner) + 6
    outer = _lifted_partition(right, degree)
    return len(_skew_tableaux(inner, outer))


class _SkewRepresentation:
    """Young's seminormal representation restricted to six added boxes."""

    def __init__(self, inner: tuple[int, int, int], outer: tuple[int, int, int]):
        self.inner = inner
        self.outer = outer
        self.tableaux = _skew_tableaux(inner, outer)
        if not self.tableaux:
            raise ValueError(f"empty skew shape {outer}/{inner}")
        self.box_count = sum(outer) - sum(inner)
        if self.box_count not in (3, 6):
            raise ValueError(
                "direct adjoint construction requires three or six added boxes"
            )
        self.dimension = len(self.tableaux)
        self._index = {tableau: index for index, tableau in enumerate(self.tableaux)}
        self.contents = tuple(self._contents(tableau) for tableau in self.tableaux)
        self._prefix_partitions = tuple(
            self._prefix_partition(tableau) for tableau in self.tableaux
        )
        self.identity = np.eye(self.dimension)
        self.generators = tuple(
            self._generator(index) for index in range(self.box_count - 1)
        )
        self._transpositions: dict[tuple[int, int], np.ndarray] = {}
        self._determinants: dict[int, np.ndarray] = {}
        self._adjoints: dict[int, np.ndarray] = {}
        self._prefix_masks: dict[tuple[int, int, int], np.ndarray] = {}
        self._block_swap: np.ndarray | None = None

    def _contents(self, tableau: tuple[int, ...]) -> tuple[int, ...]:
        shape = list(self.inner)
        result: list[int] = []
        for row in tableau:
            result.append(shape[row] - row)
            shape[row] += 1
        if tuple(shape) != self.outer:
            raise ArithmeticError("skew tableau does not end in the outer shape")
        return tuple(result)

    def _prefix_partition(self, tableau: tuple[int, ...]) -> tuple[int, int, int]:
        shape = list(self.inner)
        for row in tableau[:3]:
            shape[row] += 1
        return tuple(shape)

    def _generator(self, index: int) -> np.ndarray:
        matrix = np.zeros((self.dimension, self.dimension), dtype=float)
        for column, tableau in enumerate(self.tableaux):
            axial_distance = (
                self.contents[column][index + 1] - self.contents[column][index]
            )
            if axial_distance == 0:
                raise ArithmeticError("equal contents in a standard skew tableau")
            matrix[column, column] = 1.0 / axial_distance
            if abs(axial_distance) == 1:
                continue
            partner = list(tableau)
            partner[index], partner[index + 1] = (
                partner[index + 1],
                partner[index],
            )
            partner_index = self._index.get(tuple(partner))
            if partner_index is None:
                raise ArithmeticError("missing seminormal tableau partner")
            matrix[partner_index, column] = np.sqrt(
                1.0 - 1.0 / axial_distance**2
            )
        return matrix

    def _transposition_matrix(self, first: int, second: int) -> np.ndarray:
        if first > second:
            first, second = second, first
        if first == second:
            return self.identity
        if not 1 <= first < second <= self.box_count:
            raise IndexError(
                f"transposition labels must lie between 1 and {self.box_count}"
            )
        if second == first + 1:
            return self.generators[first - 1]
        key = first, second
        cached = self._transpositions.get(key)
        if cached is not None:
            return cached
        result = self.identity.copy()
        generator_indices = list(range(first - 1, second - 1))
        generator_indices.extend(range(second - 3, first - 2, -1))
        for generator in generator_indices:
            result = result @ self.generators[generator]
        self._transpositions[key] = result
        return result

    def transposition(self, first: int, second: int) -> np.ndarray:
        """Representation of a transposition of 1-based added labels."""

        return self._transposition_matrix(first, second).copy()

    def _determinant_matrix(self, start: int) -> np.ndarray:
        cached = self._determinants.get(start)
        if cached is not None:
            return cached
        if not 1 <= start <= self.box_count - 2:
            raise IndexError("three-line block lies outside the skew representation")
        first = (
            self.identity - self._transposition_matrix(start, start + 1)
        ) / 2.0
        second = (
            self.identity
            - self._transposition_matrix(start, start + 2)
            - self._transposition_matrix(start + 1, start + 2)
        ) / 3.0
        result = first @ second
        result = (result + result.T) / 2.0
        self._determinants[start] = result
        return result

    def determinant(self, start: int) -> np.ndarray:
        """Three-line antisymmetrizer beginning at 1-based ``start``."""

        return self._determinant_matrix(start).copy()

    def _adjoint_matrix(self, start: int) -> np.ndarray:
        cached = self._adjoints.get(start)
        if cached is not None:
            return cached
        anti_pair = (
            self.identity - self._transposition_matrix(start + 1, start + 2)
        ) / 2.0
        result = anti_pair - self._determinant_matrix(start)
        result = (result + result.T) / 2.0
        self._adjoints[start] = result
        return result

    def adjoint(self, start: int) -> np.ndarray:
        """The traceless part of ``3 tensor 3-bar`` on three lines."""

        return self._adjoint_matrix(start).copy()

    def prefix_mask(self, partition: tuple[int, int, int]) -> np.ndarray:
        """Boolean mask for the irrep carried by the first three new boxes."""

        cached = self._prefix_masks.get(partition)
        if cached is None:
            cached = np.asarray(
                [item == partition for item in self._prefix_partitions],
                dtype=bool,
            )
            cached.setflags(write=False)
            self._prefix_masks[partition] = cached
        return cached

    def prefix_projector(self, partition: tuple[int, int, int]) -> np.ndarray:
        return np.diag(self.prefix_mask(partition).astype(float))

    def block_swap(self) -> np.ndarray:
        # [1,2,3,4,5,6] -> [4,5,6,1,2,3].
        if self.box_count != 6:
            raise ValueError("block swap requires exactly six added boxes")
        if self._block_swap is not None:
            return self._block_swap.copy()
        operations = (3, 2, 1, 4, 3, 2, 5, 4, 3)
        result = self.identity.copy()
        for operation in operations:
            result = self.generators[operation - 1] @ result
        self._block_swap = result
        return result.copy()


def _matrix_rank(matrix: np.ndarray, tolerance: float = 2.0e-10) -> int:
    eigenvalues = np.linalg.eigvalsh((matrix + matrix.T) / 2.0)
    return int(np.count_nonzero(eigenvalues > tolerance))


def _prefix_restriction(
    representation: _SkewRepresentation,
    product: np.ndarray,
    partition: tuple[int, int, int],
) -> np.ndarray:
    """Restrict a product projector without constructing a dense diagonal."""

    result = product * representation.prefix_mask(partition)[None, :]
    return (result + result.T) / 2.0


def _canonical_ray(
    projector: np.ndarray,
    *,
    orthogonal_to: tuple[np.ndarray, ...] = (),
    tolerance: float = 2.0e-11,
) -> np.ndarray:
    """Project coordinate vectors in order, matching the canonical gauge rule."""

    for coordinate in range(projector.shape[0]):
        vector = projector[:, coordinate].copy()
        for _ in range(2):
            for old in orthogonal_to:
                vector -= old * np.dot(old, vector)
        norm = float(np.linalg.norm(vector))
        if norm <= tolerance:
            continue
        vector /= norm
        pivot = np.flatnonzero(np.abs(vector) > tolerance)
        if pivot.size and vector[pivot[0]] < 0.0:
            vector = -vector
        return vector
    raise ArithmeticError("could not extract a ray from a nonzero path projector")


@dataclass(frozen=True, slots=True)
class _ProjectorKernel:
    """Gauge-independent six-box data for one outer block."""

    paths: tuple[Path, ...]
    inner: tuple[int, int, int]
    outer: tuple[int, int, int]
    degree: int
    representation: _SkewRepresentation
    product: np.ndarray
    swap: np.ndarray
    swap_signature: tuple[int, int]


def _projector_residual(projector: np.ndarray) -> float:
    symmetry = float(np.max(np.abs(projector - projector.T)))
    idempotence = float(np.max(np.abs(projector @ projector - projector)))
    return max(symmetry, idempotence)


def _projector_kernel(left: Label, right: Label) -> _ProjectorKernel:
    """Construct and validate the gauge-independent part of one block."""

    paths = _paths(left, right)
    if not paths:
        raise ValueError(f"no two-adjoint paths from {left} to {right}")
    inner = _minimal_partition(left)
    degree = sum(inner) + 6
    outer = _lifted_partition(right, degree)
    representation = _SkewRepresentation(inner, outer)
    product = (
        representation._adjoint_matrix(1)
        @ representation._adjoint_matrix(4)
    )
    product = (product + product.T) / 2.0
    projector_error = _projector_residual(product)
    if projector_error > 3.0e-12:
        raise ArithmeticError(
            f"adjoint product is not a projector: residual {projector_error:.3e}"
        )
    observed_rank = _matrix_rank(product)
    if observed_rank != len(paths):
        raise ArithmeticError(
            f"product projector rank {observed_rank}, expected {len(paths)}"
        )

    swap = representation.block_swap()
    swap_error = max(
        float(np.max(np.abs(swap - swap.T))),
        float(np.max(np.abs(swap @ swap - representation.identity))),
    )
    if swap_error > 3.0e-12:
        raise ArithmeticError(
            f"literal block swap is not an involution: residual {swap_error:.3e}"
        )
    swap_product = swap @ product
    product_swap = product @ swap
    commute_error = float(np.max(np.abs(swap_product - product_swap)))
    if commute_error > 2.0e-12:
        raise ArithmeticError(
            "adjoint product projector does not commute with swap: "
            f"residual {commute_error:.3e}"
        )
    restricted_norm_squared = float(np.sum(product_swap * product_swap))
    if abs(restricted_norm_squared - len(paths)) > 3.0e-11 * len(paths):
        raise ArithmeticError(
            "restricted swap has the wrong Frobenius norm: "
            f"{restricted_norm_squared:.12g}, expected {len(paths)}"
        )
    signed_trace = float(np.trace(product_swap))
    rounded_trace = round(signed_trace)
    if abs(signed_trace - rounded_trace) > 3.0e-11 * len(paths):
        raise ArithmeticError(
            f"restricted swap trace {signed_trace:.12g} is not integral"
        )
    if (len(paths) + rounded_trace) % 2:
        raise ArithmeticError("restricted swap trace has the wrong parity")
    positive = (len(paths) + rounded_trace) // 2
    negative = len(paths) - positive
    return _ProjectorKernel(
        paths=paths,
        inner=inner,
        outer=outer,
        degree=degree,
        representation=representation,
        product=product,
        swap=swap,
        swap_signature=(positive, negative),
    )


def _valid_label(label: object) -> bool:
    return (
        isinstance(label, tuple)
        and len(label) == 2
        and all(isinstance(item, int) and item >= 0 for item in label)
    )


def _partition_label(partition: tuple[int, int, int]) -> Label:
    return partition[0] - partition[1], partition[1] - partition[2]


def _repository_edge_phase(source: Label, target: Label) -> float:
    """Pieri/determinant phase relating lifted Young and repository gauges.

    This is the explicit three-row form of the lexicographic phase convention
    in DERIVATION.md Sections 6.4--6.5.  Only gluonic labels (triality zero)
    reach this routine.
    """

    delta = target[0] - source[0], target[1] - source[1]
    negative = (
        (source[0] == 0 and target != source)
        or (
            source[0] == 1
            and source[1] >= 4
            and delta == (-1, -1)
        )
        or (
            source[0] == 2
            and source[1] >= 2
            and delta == (-2, 1)
        )
        or (
            source[1] == 1
            and source[0] >= 4
            and delta == (-1, -1)
        )
    )
    return -1.0 if negative else 1.0


@lru_cache(maxsize=_LOCAL_FRAME_CACHE_SIZE)
def _local_edge_frame(
    source: Label,
    target: Label,
    source_partition: tuple[int, int, int],
    target_partition: tuple[int, int, int],
) -> tuple[tuple[np.ndarray, ...], tuple[tuple[int, ...], ...]]:
    """Return convention-fixed rays for one three-box adjoint edge."""

    representation = _SkewRepresentation(source_partition, target_partition)
    adjoint = representation._adjoint_matrix(1)
    expected = len([edge for edge in _adjoint_edges(source) if edge.target == target])
    observed = _matrix_rank(adjoint)
    if observed != expected or observed not in (1, 2):
        raise ArithmeticError(
            f"local edge {source}->{target} has rank {observed}, expected {expected}"
        )
    phase = _repository_edge_phase(source, target)
    if observed == 1:
        vectors = (phase * _canonical_ray(adjoint),)
    else:
        if source != target:
            raise ArithmeticError("only a self edge may have multiplicity two")
        # J_{m+1} is diagonal in the seminormal basis.  Sandwiching it
        # between determinant and adjoint projectors selects the action copy.
        jucys_murphy = np.asarray(
            [contents[0] for contents in representation.contents]
        )
        determinant_ray = _canonical_ray(representation._determinant_matrix(1))
        action = adjoint @ (jucys_murphy * determinant_ray)
        action_norm = float(np.linalg.norm(action))
        if action_norm <= 2.0e-11:
            raise ArithmeticError(f"vanishing action ray for self edge {source}")
        action /= action_norm
        complement = _canonical_ray(adjoint, orthogonal_to=(action,))
        if source == ADJOINT:
            # The repository uses the even Jordan ray first and the bracket
            # with left-then-right operand order second.
            vectors = (phase * complement, -phase * action)
        else:
            vectors = (phase * action, phase * complement)
    for vector in vectors:
        vector.setflags(write=False)
    return vectors, representation.tableaux


def _middle_groups(
    paths: tuple[Path, ...],
) -> tuple[tuple[Label, tuple[tuple[int, Path], ...]], ...]:
    grouped: dict[Label, list[tuple[int, Path]]] = {}
    for index, path in enumerate(paths):
        grouped.setdefault(path.middle, []).append((index, path))
    return tuple((middle, tuple(items)) for middle, items in sorted(grouped.items()))


def _path_vectors(kernel: _ProjectorKernel) -> np.ndarray:
    """Assemble the path-ordered six-box basis from two local edge frames."""

    allowed_patterns = {
        ((0, 0),),
        ((0, 0), (1, 0)),
        ((0, 0), (0, 1)),
        ((0, 0), (0, 1), (1, 0), (1, 1)),
    }
    for middle, items in _middle_groups(kernel.paths):
        pattern = tuple(
            (path.left_multiplicity, path.right_multiplicity)
            for _index, path in items
        )
        if pattern not in allowed_patterns:
            raise ArithmeticError(
                f"unexpected multiplicity pattern through {middle}: {pattern}"
            )
        partition = _lifted_partition(middle, kernel.degree - 3)
        sector = _prefix_restriction(
            kernel.representation, kernel.product, partition
        )
        if _matrix_rank(sector) != len(items):
            raise ArithmeticError(
                f"middle sector {middle} has the wrong path rank"
            )

    left = _partition_label(kernel.inner)
    right = _partition_label(kernel.outer)
    columns: list[np.ndarray] = []
    for path in kernel.paths:
        middle_partition = _lifted_partition(path.middle, kernel.degree - 3)
        left_frame, left_tableaux = _local_edge_frame(
            left, path.middle, kernel.inner, middle_partition
        )
        right_frame, right_tableaux = _local_edge_frame(
            path.middle, right, middle_partition, kernel.outer
        )
        left_ray = left_frame[path.left_multiplicity]
        right_ray = right_frame[path.right_multiplicity]
        left_lookup = {tableau: index for index, tableau in enumerate(left_tableaux)}
        right_lookup = {
            tableau: index for index, tableau in enumerate(right_tableaux)
        }
        vector = np.zeros(kernel.representation.dimension, dtype=float)
        for index, tableau in enumerate(kernel.representation.tableaux):
            first = left_lookup.get(tableau[:3])
            second = right_lookup.get(tableau[3:])
            if first is not None and second is not None:
                vector[index] = left_ray[first] * right_ray[second]
        columns.append(vector)
    return np.column_stack(columns)


class DirectSpechtSwapOracle:
    """Block-level direct-Specht oracle in the repository vertex gauge."""

    @staticmethod
    def supported(left: Label, right: Label) -> bool:
        """Return whether a nonempty two-adjoint path space exists."""

        if not _valid_label(left) or not _valid_label(right):
            return False
        return bool(_paths(left, right))

    @lru_cache(maxsize=_SIGNED_BLOCK_CACHE_SIZE)
    def block(self, left: Label, right: Label) -> SwapBlock:
        if not self.supported(left, right):
            raise ValueError(f"no two-adjoint paths from {left} to {right}")
        kernel = _projector_kernel(left, right)
        vectors = _path_vectors(kernel)
        expected_rank = len(kernel.paths)
        gram_error = float(np.max(np.abs(vectors.T @ vectors - np.eye(expected_rank))))
        closure_error = float(np.max(np.abs(vectors @ vectors.T - kernel.product)))
        if max(gram_error, closure_error) > 3.0e-11:
            raise ArithmeticError(
                f"path basis failed: gram={gram_error:.3e}, "
                f"closure={closure_error:.3e}"
            )
        matrix = vectors.T @ kernel.swap @ vectors
        SwapTableBuilder._clean_matrix(matrix)
        block = SwapBlock(left=left, right=right, paths=kernel.paths, matrix=matrix)
        SwapTableBuilder.validate_block(block, tolerance=3.0e-10)
        return block

    @lru_cache(maxsize=_DIAGNOSTIC_CACHE_SIZE)
    def diagnostics(self, left: Label, right: Label) -> SpechtDiagnostics:
        """Return gauge-independent diagnostics for any valid outer block."""

        kernel = _projector_kernel(left, right)
        middle_counts: dict[Label, int] = {}
        for path in kernel.paths:
            middle_counts[path.middle] = middle_counts.get(path.middle, 0) + 1
        prefix_ranks: list[tuple[Label, int]] = []
        for middle, expected_rank in sorted(middle_counts.items()):
            partition = _lifted_partition(middle, kernel.degree - 3)
            projector = _prefix_restriction(
                kernel.representation, kernel.product, partition
            )
            projector_error = _projector_residual(projector)
            if projector_error > 3.0e-12:
                raise ArithmeticError(
                    f"prefix {middle} is not a projector: residual "
                    f"{projector_error:.3e}"
                )
            observed_rank = _matrix_rank(projector)
            if observed_rank != expected_rank:
                raise ArithmeticError(
                    f"prefix {middle} rank {observed_rank}, expected "
                    f"{expected_rank}"
                )
            prefix_ranks.append((middle, observed_rank))
        return SpechtDiagnostics(
            left_partition=kernel.inner,
            final_partition=kernel.outer,
            skew_dimension=kernel.representation.dimension,
            product_rank=len(kernel.paths),
            swap_signature=kernel.swap_signature,
            prefix_ranks=tuple(prefix_ranks),
        )

    def blocks_for_left(self, left: Label) -> tuple[SwapBlock, ...]:
        if not _valid_label(left):
            raise ValueError(f"invalid SU(3) label {left}")
        rights = sorted(
            {
                second.target
                for first in _adjoint_edges(left)
                for second in _adjoint_edges(first.target)
            }
        )
        return tuple(self.block(left, right) for right in rights)


@dataclass(frozen=True, slots=True)
class _DirectDecomposition:
    couplings: tuple[_Edge, ...]


class _DirectFusionRules:
    """Small structural protocol consumed by braid validation."""

    @staticmethod
    def decompose(source: Label) -> _DirectDecomposition:
        return _DirectDecomposition(_adjoint_edges(source))


class DirectSpechtSwapTableBuilder(SwapTableBuilder):
    """Build complete swap tables with the direct skew-Specht construction."""

    def __init__(self) -> None:
        # This backend deliberately has no carrier-space coupling object.
        self.couplings = _DirectFusionRules()
        self.oracle = DirectSpechtSwapOracle()

    def reachable_irreps(self, max_prefix_gluons: int) -> dict[Label, int]:
        if max_prefix_gluons < 0:
            raise ValueError("max_prefix_gluons must be non-negative")
        depths: dict[Label, int] = {(0, 0): 0}
        frontier = {(0, 0)}
        for depth in range(1, max_prefix_gluons + 1):
            next_frontier: set[Label] = set()
            for source in frontier:
                for edge in _adjoint_edges(source):
                    if edge.target not in depths:
                        depths[edge.target] = depth
                        next_frontier.add(edge.target)
            frontier = next_frontier
            if not frontier:
                break
        return depths

    @lru_cache(maxsize=None)
    def blocks_for_left(self, left: Label) -> tuple[SwapBlock, ...]:
        return self.oracle.blocks_for_left(left)
