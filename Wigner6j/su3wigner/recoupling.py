"""Eq. (8) as an index-level change of normalized coupling-tree basis."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Iterator
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .coupling import AdjointCouplings, Coupling
from .representations import Label


# ``chains`` has one path and two adjoint axes in addition to the outer irrep
# indices.  Keeping a complete high-cutoff block can therefore consume several
# hundred MiB.  Slabs target this much chain payload for all paths together;
# decoded coupling operands and overlap/validation temporaries are separate, so
# this constant is deliberately not a bound on total process memory.
_CHAIN_SLAB_BYTES = 32 * 1024 * 1024
_DENSE_DECODE_CACHE_BYTES = 64 * 1024 * 1024


@dataclass(frozen=True, order=True, slots=True)
class Path:
    """One two-adjoint fusion path ``left -> middle -> right``."""

    middle: Label
    left_multiplicity: int
    right_multiplicity: int
    left_exchange_parity: int | None = None
    right_exchange_parity: int | None = None


@dataclass(frozen=True, slots=True)
class SwapBlock:
    """Local adjacent-gluon swap matrix for fixed outer irreps."""

    left: Label
    right: Label
    paths: tuple[Path, ...]
    matrix: np.ndarray

    @property
    def size(self) -> int:
        return len(self.paths)


class SwapTableBuilder:
    """Build normalized adjoint swap coefficients and check Eq. (8)."""

    def __init__(self, couplings: AdjointCouplings | None = None) -> None:
        self.couplings = couplings or AdjointCouplings()
        self._dense_decode_budget = _DENSE_DECODE_CACHE_BYTES
        self._dense_decode_bytes = 0
        self._dense_decode_cache: OrderedDict[
            int, tuple[Coupling, np.ndarray]
        ] = OrderedDict()

    def reachable_irreps(self, max_prefix_gluons: int) -> dict[Label, int]:
        """Return gluonic irreps and their minimum adjoint-fusion depth."""

        if max_prefix_gluons < 0:
            raise ValueError("max_prefix_gluons must be non-negative")
        minimum_depth: dict[Label, int] = {(0, 0): 0}
        frontier = {(0, 0)}
        for depth in range(1, max_prefix_gluons + 1):
            new_frontier: set[Label] = set()
            for label in frontier:
                for target in self.couplings.decompose(label).targets:
                    if target not in minimum_depth:
                        minimum_depth[target] = depth
                        new_frontier.add(target)
            frontier = new_frontier
            if not frontier:
                break
        return minimum_depth

    @lru_cache(maxsize=None)
    def build(self, max_prefix_gluons: int) -> tuple[SwapBlock, ...]:
        left_labels = self.reachable_irreps(max_prefix_gluons)
        blocks: list[SwapBlock] = []
        for left in sorted(left_labels, key=lambda x: (left_labels[x], x)):
            blocks.extend(self.blocks_for_left(left))
        return tuple(blocks)

    @lru_cache(maxsize=None)
    def blocks_for_left(self, left: Label) -> tuple[SwapBlock, ...]:
        path_data = self._path_couplings_for_left(left)

        blocks: list[SwapBlock] = []
        for right in sorted(path_data):
            ordered = sorted(path_data[right], key=lambda item: item[0])
            paths = tuple(item[0] for item in ordered)
            left_dimension = self.couplings.irreps.get(left).dim
            right_dimension = self.couplings.irreps.get(right).dim
            adjoint_dimension = self.couplings.irreps.get((1, 1)).dim
            retain_single_slab = (
                len(ordered)
                * left_dimension
                * adjoint_dimension**2
                * right_dimension
                * np.dtype(float).itemsize
                <= _CHAIN_SLAB_BYTES
            )
            matrix_accumulator = np.zeros(
                (len(ordered), len(ordered)), dtype=np.longdouble
            )
            gram_accumulator = np.zeros_like(matrix_accumulator)
            for chains in self._chain_slabs(
                ordered,
                left_dimension,
                right_dimension,
                adjoint_dimension,
            ):
                chains_flat = chains.reshape(len(ordered), -1)
                swapped_flat = np.swapaxes(chains, 2, 3).reshape(
                    len(ordered), -1
                )
                matrix_accumulator += (
                    swapped_flat @ chains_flat.T
                ).astype(np.longdouble)
                gram_accumulator += (
                    chains_flat @ chains_flat.T
                ).astype(np.longdouble)
            if not retain_single_slab:
                del chains, chains_flat, swapped_flat
            matrix = np.asarray(
                matrix_accumulator / np.longdouble(right_dimension),
                dtype=float,
            )
            gram = np.asarray(
                gram_accumulator / np.longdouble(right_dimension),
                dtype=float,
            )
            self._clean_matrix(matrix)
            block = SwapBlock(left=left, right=right, paths=paths, matrix=matrix)
            self.validate_block(block)

            equation_residual_squared = np.longdouble(0.0)
            if retain_single_slab:
                reconstructed = matrix.T @ swapped_flat
                residual = chains_flat - reconstructed
                equation_residual_squared += np.longdouble(
                    residual.ravel() @ residual.ravel()
                )
            else:
                for chains in self._chain_slabs(
                    ordered,
                    left_dimension,
                    right_dimension,
                    adjoint_dimension,
                ):
                    chains_flat = chains.reshape(len(ordered), -1)
                    swapped_flat = np.swapaxes(chains, 2, 3).reshape(
                        len(ordered), -1
                    )
                    reconstructed = matrix.T @ swapped_flat
                    residual = chains_flat - reconstructed
                    equation_residual_squared += np.longdouble(
                        residual.ravel() @ residual.ravel()
                    )
            del chains, chains_flat, swapped_flat, reconstructed, residual
            basis_error, equation_error = self._index_residuals(
                gram,
                equation_residual_squared,
                right_dimension,
            )
            if (
                not np.isfinite(basis_error)
                or basis_error > 2.0e-8 * max(1, block.size)
                or not np.isfinite(equation_error)
                or equation_error > 2.0e-8
            ):
                raise ArithmeticError(
                    f"index form of Eq. (8) failed for {block.left}->{block.right}: "
                    f"basis {basis_error:.3e}, relative equation "
                    f"{equation_error:.3e}"
                )
            blocks.append(block)
        return tuple(blocks)

    def _chain_slabs(
        self,
        ordered: list[tuple[Path, Coupling, Coupling]],
        left_dimension: int,
        right_dimension: int,
        adjoint_dimension: int,
    ) -> Iterator[np.ndarray]:
        """Yield all path tensors in bounded outer-index slabs."""

        bytes_per_outer_pair = (
            len(ordered) * adjoint_dimension**2 * np.dtype(float).itemsize
        )
        outer_pairs = max(1, _CHAIN_SLAB_BYTES // bytes_per_outer_pair)
        if outer_pairs >= left_dimension:
            left_chunk = left_dimension
            right_chunk = max(
                1, min(right_dimension, outer_pairs // left_dimension)
            )
        else:
            left_chunk = outer_pairs
            right_chunk = 1

        for left_start in range(0, left_dimension, left_chunk):
            left_stop = min(left_start + left_chunk, left_dimension)
            for right_start in range(0, right_dimension, right_chunk):
                right_stop = min(right_start + right_chunk, right_dimension)
                chains = np.empty(
                    (
                        len(ordered),
                        left_stop - left_start,
                        adjoint_dimension,
                        adjoint_dimension,
                        right_stop - right_start,
                    ),
                    dtype=float,
                )
                for index, (_path, first, second) in enumerate(ordered):
                    chains[index] = self._compose_chain_slab(
                        first,
                        second,
                        left_start,
                        left_stop,
                        right_start,
                        right_stop,
                    )
                yield chains

    @staticmethod
    def _index_residuals(
        gram: np.ndarray,
        equation_residual_squared: np.longdouble,
        right_dimension: int,
    ) -> tuple[float, float]:
        """Return basis and relative Frobenius errors for the index equation."""

        identity = np.eye(gram.shape[0])
        basis_error = float(np.max(np.abs(gram - identity)))
        basis_squared = np.longdouble(right_dimension) * np.sum(
            np.diag(gram), dtype=np.longdouble
        )
        if (
            not np.isfinite(equation_residual_squared)
            or equation_residual_squared < 0.0
            or not np.isfinite(basis_squared)
            or basis_squared <= 0.0
        ):
            equation_error = float("inf")
        else:
            equation_error = float(
                np.sqrt(equation_residual_squared / basis_squared)
            )
        return basis_error, equation_error

    def _path_couplings_for_left(
        self, left: Label
    ) -> dict[Label, list[tuple[Path, Coupling, Coupling]]]:
        """Enumerate path labels together with their two vertex tensors."""

        path_data: dict[Label, list[tuple[Path, Coupling, Coupling]]] = {}
        for first in self.couplings.decompose(left).couplings:
            for second in self.couplings.decompose(first.target).couplings:
                path = Path(
                    middle=first.target,
                    left_multiplicity=first.multiplicity,
                    right_multiplicity=second.multiplicity,
                    left_exchange_parity=first.exchange_parity,
                    right_exchange_parity=second.exchange_parity,
                )
                path_data.setdefault(second.target, []).append((path, first, second))
        return path_data

    @staticmethod
    def _clean_matrix(matrix: np.ndarray) -> None:
        """Remove roundoff at the three exact values used most often."""

        matrix[np.abs(matrix) < 5.0e-14] = 0.0
        matrix[np.abs(matrix - 1.0) < 5.0e-14] = 1.0
        matrix[np.abs(matrix + 1.0) < 5.0e-14] = -1.0

    def chain_tensor(self, left: Label, right: Label, path: Path) -> np.ndarray:
        """Return ``B[r,g1,g2,t]`` for direct checking of Eq. (8)."""

        first_options = self.couplings.decompose(left).for_target(path.middle)
        first = next(
            coupling
            for coupling in first_options
            if coupling.multiplicity == path.left_multiplicity
        )
        second_options = self.couplings.decompose(path.middle).for_target(right)
        second = next(
            coupling
            for coupling in second_options
            if coupling.multiplicity == path.right_multiplicity
        )
        return self._compose_chain(first, second)

    def _compose_chain(self, first: Coupling, second: Coupling) -> np.ndarray:
        left_dimension = self.couplings.irreps.get(first.source).dim
        right_dimension = self.couplings.irreps.get(second.target).dim
        return self._compose_chain_slab(
            first, second, 0, left_dimension, 0, right_dimension
        )

    def _compose_chain_slab(
        self,
        first: Coupling,
        second: Coupling,
        left_start: int,
        left_stop: int,
        right_start: int,
        right_stop: int,
    ) -> np.ndarray:
        middle_dimension = self.couplings.irreps.get(first.target).dim
        adjoint_dimension = self.couplings.irreps.get((1, 1)).dim
        first_dense = self._decoded_coupling_tensor(first)
        if first_dense is None:
            first_tensor = first._source_slab(left_start, left_stop)
        else:
            first_tensor = first_dense[left_start:left_stop]
        second_dense = self._decoded_coupling_tensor(second)
        if second_dense is None:
            second_tensor = second._target_columns(right_start, right_stop)
        else:
            second_tensor = second_dense[:, :, right_start:right_stop]
        return (
            first_tensor.reshape(-1, middle_dimension)
            @ second_tensor.reshape(
                middle_dimension, -1
            )
        ).reshape(
            left_stop - left_start,
            adjoint_dimension,
            adjoint_dimension,
            right_stop - right_start,
        )

    def _decoded_coupling_tensor(self, coupling: Coupling) -> np.ndarray | None:
        """Return an LRU-cached dense tensor without populating Coupling caches."""

        key = id(coupling)
        public_dense = coupling._embedding_cache
        if public_dense is not None:
            stale = self._dense_decode_cache.pop(key, None)
            if stale is not None:
                self._dense_decode_bytes -= stale[1].nbytes
            source_dimension = (
                public_dense.shape[0] // coupling._particle_dimension
            )
            return public_dense.reshape(
                source_dimension,
                coupling._particle_dimension,
                public_dense.shape[1],
            )

        cached = self._dense_decode_cache.get(key)
        if cached is not None and cached[0] is coupling:
            self._dense_decode_cache.move_to_end(key)
            return cached[1]
        if cached is not None:
            self._dense_decode_bytes -= cached[1].nbytes
            del self._dense_decode_cache[key]

        packed = coupling._packed_embedding
        if packed is None:
            raise AssertionError("coupling has neither dense nor packed storage")
        dense_bytes = (
            packed.shape[0] * packed.shape[1] * packed.dtype.itemsize
        )
        if dense_bytes > self._dense_decode_budget:
            return None

        while (
            self._dense_decode_cache
            and self._dense_decode_bytes + dense_bytes
            > self._dense_decode_budget
        ):
            _old_key, (_old_coupling, old_dense) = (
                self._dense_decode_cache.popitem(last=False)
            )
            self._dense_decode_bytes -= old_dense.nbytes
        dense = packed.to_dense().reshape(
            packed.source_dimension,
            packed.particle_dimension,
            packed.target_dimension,
        )
        dense.setflags(write=False)
        self._dense_decode_cache[key] = (coupling, dense)
        self._dense_decode_bytes += dense.nbytes
        return dense

    @staticmethod
    def validate_block(
        block: SwapBlock,
        *,
        chains: np.ndarray | None = None,
        swapped_flat: np.ndarray | None = None,
        tolerance: float = 2.0e-8,
    ) -> None:
        matrix = block.matrix
        expected_shape = (block.size, block.size)
        if matrix.shape != expected_shape:
            raise ArithmeticError(
                f"invalid swap block {block.left}->{block.right}: matrix shape "
                f"{matrix.shape}, expected {expected_shape}"
            )
        if not np.all(np.isfinite(matrix)):
            raise ArithmeticError(
                f"invalid swap block {block.left}->{block.right}: "
                "matrix contains a non-finite value"
            )
        identity = np.eye(block.size)
        symmetry = float(np.max(np.abs(matrix - matrix.T)))
        orthogonality = float(np.max(np.abs(matrix.T @ matrix - identity)))
        involution = float(np.max(np.abs(matrix @ matrix - identity)))
        worst = max(symmetry, orthogonality, involution)
        if worst > tolerance * max(1, block.size):
            raise ArithmeticError(
                f"invalid swap block {block.left}->{block.right}: residual {worst:.3e}"
            )
        if chains is not None:
            right_dimension = chains.shape[-1]
            chains_flat = chains.reshape(block.size, -1)
            gram = (chains_flat @ chains_flat.T) / right_dimension
            if swapped_flat is None:
                swapped_flat = np.swapaxes(chains, 2, 3).reshape(block.size, -1)
            reconstructed = matrix.T @ swapped_flat
            residual = chains_flat - reconstructed
            equation_residual_squared = np.longdouble(
                residual.ravel() @ residual.ravel()
            )
            basis_error, equation_error = SwapTableBuilder._index_residuals(
                gram,
                equation_residual_squared,
                right_dimension,
            )
            if (
                not np.isfinite(basis_error)
                or basis_error > tolerance * max(1, block.size)
                or not np.isfinite(equation_error)
                or equation_error > tolerance
            ):
                raise ArithmeticError(
                    f"index form of Eq. (8) failed for {block.left}->{block.right}: "
                    f"basis {basis_error:.3e}, relative equation "
                    f"{equation_error:.3e}"
                )
