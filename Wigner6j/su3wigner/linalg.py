"""Small deterministic linear-algebra helpers used to fix basis gauges."""

from __future__ import annotations

import numpy as np

LoweringRecipe = tuple[tuple[int, int], ...]


def numerical_rank(singular_values: np.ndarray, rows: int, cols: int) -> int:
    """Return the standard scale-aware SVD rank."""

    if singular_values.size == 0:
        return 0
    tolerance = max(rows, cols) * np.finfo(float).eps * singular_values[0]
    return int(np.count_nonzero(singular_values > tolerance))


def canonical_subspace_basis(
    spanning_basis: np.ndarray,
    dimension: int | None = None,
    *,
    tolerance: float = 2.0e-11,
) -> np.ndarray:
    """Choose a deterministic orthonormal basis for a supplied subspace.

    ``spanning_basis`` may contain any orthonormal basis in its columns (for
    example, vectors returned by an SVD).  The orthogonal projector is
    independent of that arbitrary basis.  Projecting coordinate vectors in
    lexicographic order and applying modified Gram--Schmidt therefore fixes a
    reproducible real gauge.  The first significant component of every new
    vector is made positive.
    """

    rows = spanning_basis.shape[0]
    if dimension is None:
        dimension = spanning_basis.shape[1]
    if dimension == 0:
        return np.zeros((rows, 0), dtype=float)

    projector = spanning_basis @ spanning_basis.T
    vectors: list[np.ndarray] = []
    for coordinate in range(rows):
        vector = projector[:, coordinate].copy()
        # Two passes are cheap here and suppress loss of orthogonality in
        # nearly parallel projected coordinate vectors.
        for _ in range(2):
            for old in vectors:
                vector -= old * np.dot(old, vector)
        norm = np.linalg.norm(vector)
        if norm <= tolerance:
            continue
        vector /= norm
        pivot = np.flatnonzero(np.abs(vector) > tolerance)
        if pivot.size and vector[pivot[0]] < 0.0:
            vector = -vector
        vectors.append(vector)
        if len(vectors) == dimension:
            break

    if len(vectors) != dimension:
        raise ArithmeticError(
            f"could construct only {len(vectors)} of {dimension} basis vectors"
        )
    return np.column_stack(vectors)


def canonical_nullspace(
    matrix: np.ndarray, *, tolerance: float = 2.0e-11
) -> np.ndarray:
    """Return a deterministic real orthonormal basis for ``ker(matrix)``."""

    rows, cols = matrix.shape
    if cols == 0:
        return np.zeros((0, 0), dtype=float)
    if rows == 0:
        return np.eye(cols, dtype=float)
    if not np.any(matrix):
        return np.eye(cols, dtype=float)
    if cols == 1:
        return np.zeros((1, 0), dtype=float)
    _u, singular_values, vh = np.linalg.svd(matrix, full_matrices=True)
    rank = numerical_rank(singular_values, rows, cols)
    raw = vh[rank:, :].T
    return canonical_subspace_basis(raw, cols - rank, tolerance=tolerance)


def descendant_basis(
    lowering_1: np.ndarray,
    lowering_2: np.ndarray,
    highest: np.ndarray,
    highest_weight: tuple[int, int],
    expected_dimension: int,
    *,
    tolerance: float = 2.0e-10,
    recipe: LoweringRecipe | None = None,
    coordinate_weights: tuple[tuple[int, int], ...] | None = None,
    reference_lowerings: tuple[np.ndarray, np.ndarray] | None = None,
) -> tuple[np.ndarray, tuple[tuple[int, int], ...], LoweringRecipe]:
    """Generate an irrep basis recursively from its highest-weight state.

    The two simple-root lowering operators are always attempted in the order
    ``E_21``, ``E_32`` (zero-based matrices ``E[1,0]``, ``E[2,1]``).
    Modified Gram--Schmidt removes descendants already present at the target
    weight.  In discovery mode the accepted parent/lowering-operation pairs
    are returned as a recipe.  Replaying that recipe in an equivalent
    realization fixes the unique unitary intertwiner, including its phase,
    without letting roundoff change which nominally dependent candidate is
    accepted first.  When the reference representation's lowering matrices
    are supplied, replay uses their already-fixed expansion coefficients
    rather than independently renormalizing every descendant.  This prevents
    small nonlinear normalization errors from accumulating along long paths.
    """

    vector = np.asarray(highest, dtype=float).copy()
    coordinate_indices: dict[tuple[int, int], np.ndarray] | None = None
    if coordinate_weights is not None:
        if len(coordinate_weights) != vector.size:
            raise ValueError(
                "coordinate weight count does not match the vector dimension"
            )
        grouped: dict[tuple[int, int], list[int]] = {}
        for index, weight in enumerate(coordinate_weights):
            grouped.setdefault(weight, []).append(index)
        coordinate_indices = {
            weight: np.asarray(indices, dtype=int)
            for weight, indices in grouped.items()
        }

    vector /= np.linalg.norm(vector)
    basis = np.empty((vector.size, expected_dimension), dtype=float)
    basis[:, 0] = vector
    vectors = [basis[:, 0]]
    weights = [highest_weight]
    vectors_by_weight: dict[tuple[int, int], list[np.ndarray]] = {
        highest_weight: [vectors[0]]
    }
    indices_by_weight: dict[tuple[int, int], list[int]] = {highest_weight: [0]}
    operations = (
        (lowering_1, (-2, 1)),
        (lowering_2, (1, -2)),
    )
    accepted_steps: list[tuple[int, int]] = []

    def lowered_candidate(
        parent: int,
        lowering: np.ndarray,
        target_weight: tuple[int, int],
    ) -> np.ndarray:
        if coordinate_indices is None:
            return lowering @ vectors[parent]
        source_indices = coordinate_indices.get(weights[parent])
        target_indices = coordinate_indices.get(target_weight)
        candidate = np.zeros_like(vectors[parent])
        if source_indices is None or target_indices is None:
            return candidate
        restricted = getattr(lowering, "restricted", None)
        if restricted is None:
            block = lowering[np.ix_(target_indices, source_indices)]
        else:
            block = restricted(target_indices, source_indices)
        candidate[target_indices] = block @ vectors[parent][source_indices]
        return candidate

    def normalize_candidate(
        candidate: np.ndarray, target_weight: tuple[int, int]
    ) -> tuple[np.ndarray, float]:
        indices: slice | np.ndarray
        if coordinate_indices is None:
            indices = slice(None)
        else:
            indices = coordinate_indices.get(
                target_weight, np.empty(0, dtype=int)
            )
        for _ in range(2):
            for old in vectors_by_weight.get(target_weight, ()):
                candidate[indices] -= old[indices] * np.dot(
                    old[indices], candidate[indices]
                )
        norm = float(np.linalg.norm(candidate[indices]))
        if norm > tolerance:
            candidate[indices] /= norm
        return candidate, norm

    def replay_weight_blocks(
        replay_recipe: LoweringRecipe,
        reference: tuple[np.ndarray, np.ndarray],
    ) -> tuple[np.ndarray, tuple[tuple[int, int], ...], LoweringRecipe]:
        """Replay all incoming lowering relations one weight space at a time.

        For a target weight space with embedding ``X``, every already-built
        parent supplies an equation ``X B = A``.  Here ``A`` is the product
        lowering applied to the parent embedding and ``B`` is the matching
        block of the reference lowering matrix.  The orthogonal Procrustes
        factor of ``A B.T`` is the isometry that best satisfies all incoming
        equations simultaneously.  It avoids selecting and repeatedly
        renormalizing one path through a multiply connected weight diagram.
        """

        if coordinate_indices is None:
            raise AssertionError("weight-block replay needs coordinate weights")
        replay_weights = [highest_weight]
        for step, (parent, operation_index) in enumerate(replay_recipe, start=1):
            if not 0 <= parent < step:
                raise ValueError(
                    f"lowering recipe parent {parent} is unavailable at step {step}"
                )
            if operation_index not in (0, 1):
                raise ValueError(
                    "lowering recipe operation must be 0 or 1, got "
                    f"{operation_index}"
                )
            shift = operations[operation_index][1]
            replay_weights.append(
                (
                    replay_weights[parent][0] + shift[0],
                    replay_weights[parent][1] + shift[1],
                )
            )

        target_indices_by_weight: dict[tuple[int, int], list[int]] = {}
        first_index: dict[tuple[int, int], int] = {}
        for index, weight in enumerate(replay_weights):
            target_indices_by_weight.setdefault(weight, []).append(index)
            first_index.setdefault(weight, index)

        replayed = np.zeros((vector.size, expected_dimension), dtype=float)
        replayed[:, 0] = vector
        ordered_weights = sorted(
            target_indices_by_weight,
            key=lambda weight: (
                highest_weight[0]
                + highest_weight[1]
                - weight[0]
                - weight[1],
                first_index[weight],
            ),
        )
        for target_weight in ordered_weights:
            if target_weight == highest_weight:
                continue
            output_coordinates = coordinate_indices.get(target_weight)
            if output_coordinates is None:
                raise ArithmeticError(
                    f"product has no coordinate sector at weight {target_weight}"
                )
            output_targets = np.asarray(
                target_indices_by_weight[target_weight], dtype=int
            )
            product_blocks: list[np.ndarray] = []
            reference_blocks: list[np.ndarray] = []
            for operation_index, (lowering, shift) in enumerate(operations):
                parent_weight = (
                    target_weight[0] - shift[0],
                    target_weight[1] - shift[1],
                )
                parent_targets_list = target_indices_by_weight.get(parent_weight)
                parent_coordinates = coordinate_indices.get(parent_weight)
                if parent_targets_list is None or parent_coordinates is None:
                    continue
                parent_targets = np.asarray(parent_targets_list, dtype=int)
                lowering_restricted = getattr(lowering, "restricted", None)
                if lowering_restricted is None:
                    product_lowering = lowering[
                        np.ix_(output_coordinates, parent_coordinates)
                    ]
                else:
                    product_lowering = lowering_restricted(
                        output_coordinates, parent_coordinates
                    )
                product_blocks.append(
                    product_lowering
                    @ replayed[np.ix_(parent_coordinates, parent_targets)]
                )
                reference_blocks.append(
                    reference[operation_index][
                        np.ix_(output_targets, parent_targets)
                    ]
                )
            if not product_blocks:
                raise ArithmeticError(
                    f"weight {target_weight} has no constructed parent sector"
                )
            product_data = np.concatenate(product_blocks, axis=1)
            reference_data = np.concatenate(reference_blocks, axis=1)
            correlation = product_data @ reference_data.T
            left, singular_values, right_transpose = np.linalg.svd(
                correlation, full_matrices=False
            )
            rank = numerical_rank(
                singular_values, correlation.shape[0], correlation.shape[1]
            )
            if rank != output_targets.size:
                raise ArithmeticError(
                    f"lowering relations span rank {rank}/{output_targets.size} "
                    f"at weight {target_weight}"
                )
            replayed[np.ix_(output_coordinates, output_targets)] = (
                left @ right_transpose
            )
        return replayed, tuple(replay_weights), replay_recipe

    if (
        recipe is not None
        and reference_lowerings is not None
        and coordinate_indices is not None
    ):
        if len(recipe) != expected_dimension - 1:
            raise ValueError(
                f"lowering recipe has {len(recipe)} steps, expected "
                f"{expected_dimension - 1}"
            )
        return replay_weight_blocks(recipe, reference_lowerings)

    if recipe is not None:
        if len(recipe) != expected_dimension - 1:
            raise ValueError(
                f"lowering recipe has {len(recipe)} steps, expected "
                f"{expected_dimension - 1}"
            )
        for parent, operation_index in recipe:
            if not 0 <= parent < len(vectors):
                raise ValueError(
                    f"lowering recipe parent {parent} is unavailable at step "
                    f"{len(vectors)}"
                )
            if operation_index not in (0, 1):
                raise ValueError(
                    f"lowering recipe operation must be 0 or 1, got "
                    f"{operation_index}"
                )
            lowering, shift = operations[operation_index]
            target_weight = (
                weights[parent][0] + shift[0],
                weights[parent][1] + shift[1],
            )
            candidate = lowered_candidate(parent, lowering, target_weight)
            if reference_lowerings is None:
                candidate, norm = normalize_candidate(candidate, target_weight)
                if norm <= tolerance:
                    raise ArithmeticError(
                        "lowering recipe produced a dependent descendant at "
                        f"step {len(vectors)}: parent={parent}, "
                        f"operation={operation_index}, norm={norm:.3e}"
                    )
            else:
                reference = reference_lowerings[operation_index]
                new_index = len(vectors)
                if reference.shape != (expected_dimension, expected_dimension):
                    raise ValueError(
                        "reference lowering matrix has the wrong shape"
                    )
                for old_index in indices_by_weight.get(target_weight, ()):
                    candidate -= (
                        vectors[old_index] * reference[old_index, parent]
                    )
                # Suppress only floating-point leakage along existing columns;
                # their exact coefficients were removed immediately above.
                for _ in range(2):
                    for old_index in indices_by_weight.get(target_weight, ()):
                        old = vectors[old_index]
                        candidate -= old * np.dot(old, candidate)
                pivot = float(reference[new_index, parent])
                if abs(pivot) <= tolerance:
                    raise ArithmeticError(
                        "reference lowering coefficient vanished at "
                        f"step {new_index}: parent={parent}, "
                        f"operation={operation_index}, coefficient={pivot:.3e}"
                    )
                candidate /= pivot
                norm = float(np.linalg.norm(candidate))
                if norm <= tolerance:
                    raise ArithmeticError(
                        "reference lowering replay produced a zero descendant at "
                        f"step {new_index}"
                    )
            basis[:, len(vectors)] = candidate
            candidate = basis[:, len(vectors)]
            vectors.append(candidate)
            weights.append(target_weight)
            vectors_by_weight.setdefault(target_weight, []).append(candidate)
            indices_by_weight.setdefault(target_weight, []).append(
                len(vectors) - 1
            )
        return basis, tuple(weights), recipe

    cursor = 0
    while cursor < len(vectors) and len(vectors) < expected_dimension:
        source_weight = weights[cursor]
        for operation_index, (lowering, shift) in enumerate(operations):
            target_weight = (
                source_weight[0] + shift[0],
                source_weight[1] + shift[1],
            )
            candidate = lowered_candidate(cursor, lowering, target_weight)
            candidate, norm = normalize_candidate(candidate, target_weight)
            if norm <= tolerance:
                continue
            # Do not flip this sign: it is inherited from the lowering
            # operator and is precisely the relative phase convention.
            basis[:, len(vectors)] = candidate
            candidate = basis[:, len(vectors)]
            vectors.append(candidate)
            weights.append(target_weight)
            vectors_by_weight.setdefault(target_weight, []).append(candidate)
            indices_by_weight.setdefault(target_weight, []).append(
                len(vectors) - 1
            )
            accepted_steps.append((cursor, operation_index))
            if len(vectors) == expected_dimension:
                break
        cursor += 1

    if len(vectors) != expected_dimension:
        raise ArithmeticError(
            "lowering recursion generated "
            f"{len(vectors)} states, expected {expected_dimension}"
        )
    return basis, tuple(weights), tuple(accepted_steps)
