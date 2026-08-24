"""Dimension-only SU(3) adjacent-adjoint recoupling.

This module constructs the same :class:`~su3wigner.recoupling.SwapBlock`
objects as the tensor backends, but it never constructs a representation,
generator, Clebsch--Gordan embedding, or nullspace.  Its only scalar inputs
are the closed Weyl dimension and quadratic-Casimir formulae.  The remaining
ingredients are finite label rules:

* the ``(p,q)`` fusion rules with ``3``, ``3-bar``, and ``8``;
* Young's seminormal adjacent-transposition matrix for two equal elementary
  lines;
* dimension/Casimir formulae for the two mixed elementary frames; and
* the symmetric/antisymmetric convention for the two ``8 x 8 -> 8``
  vertices.

Resolving an adjoint label as ``3,3-bar`` is therefore only a combinatorial
resolution of fusion paths.  It is not a Clebsch--Gordan decomposition.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import sqrt

import numpy as np

from .recoupling import Path, SwapBlock, SwapTableBuilder
from .representations import Label, dimension


FUNDAMENTAL: Label = (1, 0)
ANTIFUNDAMENTAL: Label = (0, 1)
ADJOINT: Label = (1, 1)
ADJOINT_SPLIT_SEQUENCE: tuple[Label, ...] = (
    FUNDAMENTAL,
    ANTIFUNDAMENTAL,
    FUNDAMENTAL,
    ANTIFUNDAMENTAL,
)
_STRUCTURAL_CACHE_SIZE = 4096
_ELEMENTARY_BLOCK_CACHE_SIZE = 4096
_FRAME_CACHE_SIZE = 1024
_VERTEX_REDUCTION_CACHE_SIZE = 4096
_FUSION_PATH_CACHE_SIZE = 256


def _conjugate(label: Label) -> Label:
    return label[1], label[0]


def _casimir(label: Label) -> float:
    """Quadratic Casimir in the normalization ``C_F=4/3``."""

    p, q = label
    return (p * p + p * q + q * q + 3 * p + 3 * q) / 3.0


@dataclass(frozen=True, slots=True)
class StructuralCoupling:
    """One label-level edge in ``source tensor 8``.

    Deliberately absent is an ``embedding`` field: this backend has no vector
    spaces on which such an object could act.
    """

    source: Label
    target: Label
    multiplicity: int
    exchange_parity: int | None = None


@dataclass(frozen=True, slots=True)
class StructuralDecomposition:
    """The analytic representation-ring decomposition of ``source x 8``."""

    source: Label
    couplings: tuple[StructuralCoupling, ...]

    @property
    def targets(self) -> tuple[Label, ...]:
        return tuple(sorted({coupling.target for coupling in self.couplings}))

    def for_target(self, target: Label) -> tuple[StructuralCoupling, ...]:
        return tuple(c for c in self.couplings if c.target == target)


class AdjointFusionRules:
    """Analytic ``(p,q) tensor (1,1)`` fusion edges only."""

    @lru_cache(maxsize=_STRUCTURAL_CACHE_SIZE)
    def decompose(self, source: Label) -> StructuralDecomposition:
        p, q = source
        if p < 0 or q < 0:
            raise ValueError("Dynkin labels must be non-negative")

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
        edges: list[StructuralCoupling] = []
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
            edges.append(StructuralCoupling(source, target, multiplicity, parity))
        accounted_dimension = sum(dimension(edge.target) for edge in edges)
        expected_dimension = 8 * dimension(source)
        if accounted_dimension != expected_dimension:
            raise ArithmeticError(
                f"structural decomposition of {source} x 8 accounts for "
                f"{accounted_dimension}/{expected_dimension} dimensions"
            )
        return StructuralDecomposition(
            source,
            tuple(sorted(edges, key=lambda edge: (edge.target, edge.multiplicity))),
        )


@dataclass(frozen=True, slots=True)
class ElementaryDecomposition:
    """Multiplicity-free label edges for tensoring by ``3`` or ``3-bar``."""

    source: Label
    particle: Label
    targets: tuple[Label, ...]


class ElementaryFusionRules:
    """Pieri rules for an elementary fundamental or antifundamental line."""

    @lru_cache(maxsize=_STRUCTURAL_CACHE_SIZE)
    def decompose(self, source: Label, particle: Label) -> ElementaryDecomposition:
        p, q = source
        if p < 0 or q < 0:
            raise ValueError("Dynkin labels must be non-negative")
        if particle == FUNDAMENTAL:
            targets: list[Label] = [(p + 1, q)]
            if p:
                targets.append((p - 1, q + 1))
            if q:
                targets.append((p, q - 1))
        elif particle == ANTIFUNDAMENTAL:
            # Conjugate of the preceding rule.
            targets = [(p, q + 1)]
            if q:
                targets.append((p + 1, q - 1))
            if p:
                targets.append((p - 1, q))
        else:
            raise ValueError("elementary fusion supports only 3 and 3-bar")
        return ElementaryDecomposition(source, particle, tuple(sorted(targets)))


@dataclass(frozen=True, slots=True)
class ElementarySwapBlock:
    """An adjacent swap in a two-line elementary fusion space."""

    source: Label
    target: Label
    first_particle: Label
    second_particle: Label
    input_middles: tuple[Label, ...]
    output_middles: tuple[Label, ...]
    matrix: np.ndarray


@dataclass(frozen=True, slots=True)
class ElementaryFrame:
    """Singlet/adjoint frame for ``R -> middle -> R`` mixed paths."""

    source: Label
    first_particle: Label
    intermediates: tuple[Label, ...]
    singlet: np.ndarray
    action: np.ndarray | None
    complement: np.ndarray | None
    matrix: np.ndarray


@dataclass(frozen=True, slots=True)
class PathSwap:
    """One adjacent exchange on complete left-associated label paths."""

    input_sequence: tuple[Label, ...]
    output_sequence: tuple[Label, ...]
    input_paths: tuple[tuple[Label, ...], ...]
    output_paths: tuple[tuple[Label, ...], ...]
    matrix: np.ndarray


class DimensionOnlyRecouplings:
    """Elementary recouplings obtained without constructing state vectors."""

    def __init__(self, fusion: ElementaryFusionRules | None = None) -> None:
        self.fusion = fusion or ElementaryFusionRules()

    def _two_step_options(
        self,
        source: Label,
        first_particle: Label,
        second_particle: Label,
        target: Label,
    ) -> tuple[Label, ...]:
        return tuple(
            sorted(
                middle
                for middle in self.fusion.decompose(source, first_particle).targets
                if target in self.fusion.decompose(middle, second_particle).targets
            )
        )

    @staticmethod
    def _fundamental_row(source: Label, target: Label) -> int:
        p, q = source
        branches = {
            (p + 1, q): 0,
            (p - 1, q + 1): 1,
            (p, q - 1): 2,
        }
        try:
            row = branches[target]
        except KeyError as error:
            raise ArithmeticError(f"{source} -> {target} is not a 3 edge") from error
        if row == 1 and not p:
            raise ArithmeticError(f"forbidden second-row edge {source} -> {target}")
        if row == 2 and not q:
            raise ArithmeticError(f"forbidden third-row edge {source} -> {target}")
        return row

    @classmethod
    def _box_contents(
        cls, source: Label, first_middle: Label, target: Label
    ) -> tuple[int, int]:
        """Contents of the two added boxes, retaining determinant columns."""

        p, q = source
        partition = [p + q, q, 0]
        first_row = cls._fundamental_row(source, first_middle)
        first_content = partition[first_row] - first_row
        partition[first_row] += 1
        second_row = cls._fundamental_row(first_middle, target)
        second_content = partition[second_row] - second_row
        partition[second_row] += 1
        if any(partition[index] < partition[index + 1] for index in range(2)):
            raise ArithmeticError("fundamental path did not produce a partition")
        return first_content, second_content

    @classmethod
    def _ff_matrix(
        cls, source: Label, target: Label, middles: tuple[Label, ...]
    ) -> np.ndarray:
        if len(middles) not in (1, 2):
            raise ArithmeticError(
                f"two-fundamental space has unexpected size {len(middles)}"
            )
        diagonal = []
        for middle in middles:
            first_content, second_content = cls._box_contents(source, middle, target)
            axial_distance = second_content - first_content
            if not axial_distance:
                raise ArithmeticError("two added boxes have equal content")
            diagonal.append(1.0 / axial_distance)
        matrix = np.diag(diagonal)
        if len(middles) == 2:
            off_diagonal = sqrt(max(0.0, 1.0 - diagonal[0] ** 2))
            matrix[0, 1] = off_diagonal
            matrix[1, 0] = off_diagonal
        return matrix

    @staticmethod
    def _branch_phase(source: Label, middle: Label, particle: Label) -> float:
        p, q = source
        if particle == FUNDAMENTAL:
            return -1.0 if middle == (p - 1, q + 1) else 1.0
        if particle == ANTIFUNDAMENTAL:
            return -1.0 if middle == (p + 1, q - 1) else 1.0
        raise ValueError("branch phase is defined only for 3 and 3-bar")

    @lru_cache(maxsize=_FRAME_CACHE_SIZE)
    def frame(self, source: Label, first_particle: Label) -> ElementaryFrame:
        if first_particle == FUNDAMENTAL:
            second_particle = ANTIFUNDAMENTAL
        elif first_particle == ANTIFUNDAMENTAL:
            second_particle = FUNDAMENTAL
        else:
            raise ValueError("mixed frame starts with 3 or 3-bar")
        middles = self._two_step_options(
            source, first_particle, second_particle, source
        )
        if not middles:
            raise ArithmeticError(f"no mixed return paths for {source}")

        source_dimension = dimension(source)
        phases = np.array(
            [self._branch_phase(source, middle, first_particle) for middle in middles]
        )
        middle_dimensions = np.array(
            [dimension(middle) for middle in middles], dtype=float
        )
        singlet = phases * np.sqrt(middle_dimensions / (3.0 * source_dimension))
        singlet /= np.linalg.norm(singlet)

        columns = [singlet]
        action: np.ndarray | None = None
        complement: np.ndarray | None = None
        if len(middles) > 1:
            action = (
                phases
                * np.sqrt(middle_dimensions)
                * np.array(
                    [
                        (_casimir(middle) - _casimir(source) - 4.0 / 3.0) / 2.0
                        for middle in middles
                    ]
                )
            )
            action /= np.linalg.norm(action)
            if abs(float(np.dot(singlet, action))) > 3.0e-14:
                raise ArithmeticError("dimension/Casimir frame is not orthogonal")
            columns.append(action)
        if len(middles) == 3:
            assert action is not None
            complement = np.cross(singlet, action)
            complement /= np.linalg.norm(complement)
            columns.append(complement)
        if len(middles) > 3:
            raise ArithmeticError("mixed elementary frame has dimension above three")
        matrix = np.column_stack(columns)
        residual = float(np.max(np.abs(matrix.T @ matrix - np.eye(len(middles)))))
        if residual > 3.0e-14:
            raise ArithmeticError(f"elementary frame residual {residual:.3e}")
        return ElementaryFrame(
            source,
            first_particle,
            middles,
            singlet,
            action,
            complement,
            matrix,
        )

    @lru_cache(maxsize=_ELEMENTARY_BLOCK_CACHE_SIZE)
    def block(
        self,
        source: Label,
        first_particle: Label,
        second_particle: Label,
        target: Label,
    ) -> ElementarySwapBlock:
        inputs = self._two_step_options(source, first_particle, second_particle, target)
        outputs = self._two_step_options(
            source, second_particle, first_particle, target
        )
        if not inputs or len(inputs) != len(outputs):
            raise ArithmeticError(
                "elementary recoupling path mismatch for "
                f"{source} x {first_particle} x {second_particle} -> {target}"
            )

        if first_particle == second_particle == FUNDAMENTAL:
            matrix = self._ff_matrix(source, target, inputs)
        elif first_particle == second_particle == ANTIFUNDAMENTAL:
            conjugate_block = self.block(
                _conjugate(source),
                FUNDAMENTAL,
                FUNDAMENTAL,
                _conjugate(target),
            )
            conjugate_lookup = {
                middle: index
                for index, middle in enumerate(conjugate_block.input_middles)
            }
            indices = [conjugate_lookup[_conjugate(middle)] for middle in inputs]
            matrix = conjugate_block.matrix[np.ix_(indices, indices)].copy()
        elif first_particle == FUNDAMENTAL and second_particle == ANTIFUNDAMENTAL:
            if target != source:
                if len(inputs) != 1:
                    raise ArithmeticError("non-self mixed fusion is not unique")
                matrix = np.ones((1, 1), dtype=float)
            else:
                input_frame = self.frame(source, FUNDAMENTAL)
                output_frame = self.frame(source, ANTIFUNDAMENTAL)
                if (
                    input_frame.intermediates != inputs
                    or output_frame.intermediates != outputs
                ):
                    raise ArithmeticError("mixed-frame path ordering mismatch")
                eigenvalues = np.ones(len(inputs), dtype=float)
                if len(inputs) > 1:
                    eigenvalues[1] = -1.0
                matrix = (
                    output_frame.matrix @ np.diag(eigenvalues) @ input_frame.matrix.T
                )
        elif first_particle == ANTIFUNDAMENTAL and second_particle == FUNDAMENTAL:
            reverse = self.block(source, FUNDAMENTAL, ANTIFUNDAMENTAL, target)
            if reverse.input_middles != outputs or reverse.output_middles != inputs:
                raise ArithmeticError("reverse mixed-path ordering mismatch")
            matrix = reverse.matrix.T.copy()
        else:
            raise ValueError("elementary swaps support only 3 and 3-bar")

        SwapTableBuilder._clean_matrix(matrix)
        identity = np.eye(len(inputs))
        worst = max(
            float(np.max(np.abs(matrix.T @ matrix - identity))),
            float(np.max(np.abs(matrix @ matrix.T - identity))),
        )
        if worst > 5.0e-13:
            raise ArithmeticError(f"elementary swap residual {worst:.3e}")
        return ElementarySwapBlock(
            source,
            target,
            first_particle,
            second_particle,
            inputs,
            outputs,
            matrix,
        )

    @lru_cache(maxsize=_FUSION_PATH_CACHE_SIZE)
    def fusion_paths(
        self,
        source: Label,
        sequence: tuple[Label, ...],
        target: Label,
    ) -> tuple[tuple[Label, ...], ...]:
        paths: tuple[tuple[Label, ...], ...] = ((source,),)
        for particle in sequence:
            paths = tuple(
                path + (next_label,)
                for path in paths
                for next_label in self.fusion.decompose(path[-1], particle).targets
            )
        return tuple(sorted(path for path in paths if path[-1] == target))

    def adjacent_swap(
        self,
        source: Label,
        target: Label,
        sequence: tuple[Label, ...],
        position: int,
    ) -> PathSwap:
        if position < 0 or position + 1 >= len(sequence):
            raise IndexError("adjacent-swap position is outside the sequence")
        output_sequence_list = list(sequence)
        output_sequence_list[position], output_sequence_list[position + 1] = (
            output_sequence_list[position + 1],
            output_sequence_list[position],
        )
        output_sequence = tuple(output_sequence_list)
        input_paths = self.fusion_paths(source, sequence, target)
        output_paths = self.fusion_paths(source, output_sequence, target)
        if len(input_paths) != len(output_paths):
            raise ArithmeticError("adjacent swap changed the path-space dimension")
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
        identity = np.eye(len(input_paths))
        worst = max(
            float(np.max(np.abs(matrix.T @ matrix - identity))),
            float(np.max(np.abs(matrix @ matrix.T - identity))),
        )
        if worst > 8.0e-13 * max(1, len(input_paths)):
            raise ArithmeticError(f"path-swap residual {worst:.3e}")
        return PathSwap(
            sequence,
            output_sequence,
            input_paths,
            output_paths,
            matrix,
        )

    def _apply_adjacent_swap(
        self,
        sequence: tuple[Label, ...],
        position: int,
        input_paths: tuple[tuple[Label, ...], ...],
        values: np.ndarray,
    ) -> tuple[tuple[Label, ...], tuple[tuple[Label, ...], ...], np.ndarray]:
        """Apply one sparse local move without materializing its square matrix."""

        if position < 0 or position + 1 >= len(sequence):
            raise IndexError("adjacent-swap position is outside the sequence")
        if values.shape[0] != len(input_paths):
            raise ValueError("adjacent-swap input basis does not match its values")
        output_sequence_list = list(sequence)
        output_sequence_list[position], output_sequence_list[position + 1] = (
            output_sequence_list[position + 1],
            output_sequence_list[position],
        )
        output_sequence = tuple(output_sequence_list)

        # Every elementary block lists the complete local output basis.  Build
        # the globally ordered output basis from those sparse transitions so
        # that applying a move does not enumerate the complete fusion tree a
        # second time.
        transitions: list[tuple[int, tuple[Label, ...], float]] = []
        output_path_set: set[tuple[Label, ...]] = set()
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
                output_path = tuple(changed)
                transitions.append(
                    (
                        input_index,
                        output_path,
                        float(local.matrix[local_output, local_input]),
                    )
                )
                output_path_set.add(output_path)
        output_paths = tuple(sorted(output_path_set))
        if len(output_paths) != len(input_paths):
            raise ArithmeticError("adjacent swap changed the path-space dimension")
        output_lookup = {path: index for index, path in enumerate(output_paths)}
        output = np.zeros((len(output_paths), values.shape[1]), dtype=float)
        for input_index, output_path, coefficient in transitions:
            output[output_lookup[output_path]] += coefficient * values[input_index]
        return output_sequence, output_paths, output

    def apply_adjoint_pair_swap(
        self, source: Label, target: Label, values: np.ndarray
    ) -> np.ndarray:
        """Apply the four sparse elementary moves to one or more columns."""

        array = np.asarray(values, dtype=float)
        if array.ndim != 2:
            raise ValueError("pair-swap values must be a two-dimensional array")
        sequence = ADJOINT_SPLIT_SEQUENCE
        paths = self.fusion_paths(source, sequence, target)
        current_sequence = sequence
        current_paths = paths
        current = array
        # Move the first (3,3-bar) pair right, then the second pair left.
        for position in (1, 0, 2, 1):
            current_sequence, current_paths, current = self._apply_adjacent_swap(
                current_sequence,
                position,
                current_paths,
                current,
            )
        if current_sequence != sequence or current_paths != paths:
            raise ArithmeticError("sparse pair permutation ended in the wrong basis")
        return current

    def line_across_adjoint(
        self,
        source: Label,
        target: Label,
        sequence: tuple[Label, ...],
        pair_position: int,
    ) -> PathSwap:
        if sequence[pair_position : pair_position + 2] != (
            FUNDAMENTAL,
            ANTIFUNDAMENTAL,
        ):
            raise ValueError("selected positions are not a split adjoint")
        if pair_position + 2 >= len(sequence):
            raise ValueError("there is no elementary line to move")
        input_paths = self.fusion_paths(source, sequence, target)
        total = np.eye(len(input_paths))
        current_sequence = sequence
        current_paths = input_paths
        for position in (pair_position + 1, pair_position):
            step = self.adjacent_swap(source, target, current_sequence, position)
            if step.input_paths != current_paths:
                raise ArithmeticError("one-adjoint path ordering changed")
            total = step.matrix @ total
            current_sequence = step.output_sequence
            current_paths = step.output_paths
        return PathSwap(
            sequence,
            current_sequence,
            input_paths,
            current_paths,
            total,
        )

    def adjoint_pair_swap(
        self, source: Label, target: Label
    ) -> tuple[tuple[tuple[Label, ...], ...], np.ndarray]:
        sequence = ADJOINT_SPLIT_SEQUENCE
        paths = self.fusion_paths(source, sequence, target)
        first = self.line_across_adjoint(source, target, sequence, 0)
        second = self.line_across_adjoint(source, target, first.output_sequence, 1)
        total = second.matrix @ first.matrix
        if second.output_sequence != sequence or second.output_paths != paths:
            raise ArithmeticError("pair permutation ended in the wrong path basis")
        SwapTableBuilder._clean_matrix(total)
        identity = np.eye(len(paths))
        worst = max(
            float(np.max(np.abs(total - total.T))),
            float(np.max(np.abs(total @ total - identity))),
        )
        if worst > 2.0e-12 * max(1, len(paths)):
            raise ArithmeticError(f"elementary pair-swap residual {worst:.3e}")
        return paths, total


@dataclass(frozen=True, slots=True)
class VertexReduction:
    """Dimension-only expansion of a fixed adjoint edge in ``3,3-bar`` paths."""

    intermediates: tuple[Label, ...]
    coefficients: np.ndarray


class DimensionOnlyVertexReductions:
    """Analytic convention frames for a single adjoint vertex."""

    def __init__(self, recouplings: DimensionOnlyRecouplings) -> None:
        self.recouplings = recouplings

    @lru_cache(maxsize=_VERTEX_REDUCTION_CACHE_SIZE)
    def reduce(
        self, source: Label, target: Label, multiplicity: int
    ) -> VertexReduction:
        middles = self.recouplings._two_step_options(
            source, FUNDAMENTAL, ANTIFUNDAMENTAL, target
        )
        if target != source:
            if multiplicity != 0 or len(middles) != 1:
                raise ArithmeticError("non-self adjoint edge is not unique")
            coefficients = np.ones(1, dtype=float)
        else:
            frame = self.recouplings.frame(source, FUNDAMENTAL)
            if frame.intermediates != middles or frame.action is None:
                raise ArithmeticError("self-adjoint frame is incomplete")
            if source == ADJOINT:
                if frame.complement is None:
                    raise ArithmeticError("8 x 8 frame lacks its third direction")
                if multiplicity == 0:
                    # d = cross(f,g), with f=-action.
                    coefficients = np.cross(-frame.action, frame.singlet)
                elif multiplicity == 1:
                    coefficients = -frame.action.copy()
                else:
                    raise ArithmeticError("8 x 8 -> 8 has two copies")
            elif multiplicity == 0:
                # The ordinary representation-action vertex.
                coefficients = frame.action.copy()
            elif multiplicity == 1 and frame.complement is not None:
                coefficients = np.cross(frame.singlet, frame.action)
            else:
                raise ArithmeticError(
                    f"invalid self-edge multiplicity {multiplicity} for {source}"
                )
        SwapTableBuilder._clean_matrix(coefficients)
        norm_error = abs(float(np.dot(coefficients, coefficients)) - 1.0)
        if norm_error > 5.0e-13:
            raise ArithmeticError(f"vertex-frame norm residual {norm_error:.3e}")
        return VertexReduction(middles, coefficients)


class DimensionOnlySwapTableBuilder(SwapTableBuilder):
    """Build adjacent-adjoint swap tables from labels and dimensions alone."""

    def __init__(self) -> None:
        # Do not call SwapTableBuilder.__init__: its default constructs the
        # full Clebsch--Gordan backend.  ``couplings`` intentionally supplies
        # only the small structural protocol used by clients and validators.
        self.couplings = AdjointFusionRules()
        self.elementary_fusion = ElementaryFusionRules()
        self.elementary_recouplings = DimensionOnlyRecouplings(self.elementary_fusion)
        self.vertex_reductions = DimensionOnlyVertexReductions(
            self.elementary_recouplings
        )

    def chain_tensor(self, left: Label, right: Label, path: Path) -> np.ndarray:
        """Reject the carrier-space audit that this backend intentionally lacks."""

        del left, right, path
        raise NotImplementedError(
            "dimension-only construction has no carrier-space chain tensor"
        )

    def _path_edges_for_left(self, left: Label) -> dict[
        Label,
        list[tuple[Path, StructuralCoupling, StructuralCoupling]],
    ]:
        data: dict[
            Label,
            list[tuple[Path, StructuralCoupling, StructuralCoupling]],
        ] = {}
        for first in self.couplings.decompose(left).couplings:
            for second in self.couplings.decompose(first.target).couplings:
                path = Path(
                    middle=first.target,
                    left_multiplicity=first.multiplicity,
                    right_multiplicity=second.multiplicity,
                    left_exchange_parity=first.exchange_parity,
                    right_exchange_parity=second.exchange_parity,
                )
                data.setdefault(second.target, []).append((path, first, second))
        return data

    @lru_cache(maxsize=None)
    def blocks_for_left(self, left: Label) -> tuple[SwapBlock, ...]:
        path_data = self._path_edges_for_left(left)
        blocks: list[SwapBlock] = []
        for right in sorted(path_data):
            ordered = sorted(path_data[right], key=lambda item: item[0])
            paths = tuple(item[0] for item in ordered)
            elementary_paths = self.elementary_recouplings.fusion_paths(
                left, ADJOINT_SPLIT_SEQUENCE, right
            )
            lookup = {path: index for index, path in enumerate(elementary_paths)}
            expansion = np.zeros((len(elementary_paths), len(paths)), dtype=float)
            for column, (_, first, second) in enumerate(ordered):
                first_frame = self.vertex_reductions.reduce(
                    first.source, first.target, first.multiplicity
                )
                second_frame = self.vertex_reductions.reduce(
                    second.source, second.target, second.multiplicity
                )
                for first_middle, first_coefficient in zip(
                    first_frame.intermediates,
                    first_frame.coefficients,
                    strict=True,
                ):
                    for second_middle, second_coefficient in zip(
                        second_frame.intermediates,
                        second_frame.coefficients,
                        strict=True,
                    ):
                        elementary_path = (
                            left,
                            first_middle,
                            first.target,
                            second_middle,
                            right,
                        )
                        expansion[lookup[elementary_path], column] += (
                            first_coefficient * second_coefficient
                        )

            identity = np.eye(len(paths))
            basis_error = float(np.max(np.abs(expansion.T @ expansion - identity)))
            transformed = self.elementary_recouplings.apply_adjoint_pair_swap(
                left, right, expansion
            )
            matrix = expansion.T @ transformed
            self._clean_matrix(matrix)
            closure_error = float(np.max(np.abs(transformed - expansion @ matrix)))
            if max(basis_error, closure_error) > 4.0e-12 * max(1, len(paths)):
                raise ArithmeticError(
                    f"dimension-only reduction failed for {left}->{right}: "
                    f"basis={basis_error:.3e}, closure={closure_error:.3e}"
                )
            block = SwapBlock(left, right, paths, matrix)
            self.validate_block(block, tolerance=4.0e-12)
            blocks.append(block)
        return tuple(blocks)
