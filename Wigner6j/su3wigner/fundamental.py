"""Resolve the adjoint into a traceless fundamental--antifundamental pair.

The only input is the concrete adjoint module constructed in
``representations.py``.  In particular, no Fierz coefficient or tabulated
birdtrack normalization is used: both completeness identities are checked
against the orthogonal decomposition of the nine-dimensional matrix space
into the singlet and its traceless complement.
"""

from __future__ import annotations

import numpy as np

from .representations import Irrep


def adjoint_split_isometry(adjoint: Irrep) -> np.ndarray:
    """Return ``S[a,i,j]`` embedding ``8`` in ``3 tensor 3-bar``.

    The adjoint basis used everywhere else in the package is recovered from
    the traceless matrices ``S[a]``.  Rows of the flattened result are
    orthonormal, while their outer products resolve the projector onto the
    traceless part of the full matrix space.
    """

    if adjoint.label != (1, 1):
        raise ValueError("fundamental splitting requested for a non-adjoint irrep")

    splitting = np.zeros((adjoint.dim, 3, 3), dtype=float)
    for adjoint_index in range(adjoint.dim):
        for coefficient, state in zip(
            adjoint.ambient_embedding[:, adjoint_index],
            adjoint.ambient_states,
            strict=True,
        ):
            if abs(coefficient) < 1.0e-15:
                continue
            upper = state[:3].index(1)
            lower = state[3:].index(1)
            splitting[adjoint_index, upper, lower] += coefficient

    flattened = splitting.reshape(adjoint.dim, 9)
    singlet = np.eye(3, dtype=float).reshape(9) / np.sqrt(3.0)
    traceless_projector = np.eye(9) - np.outer(singlet, singlet)
    isometry_error = float(
        np.max(np.abs(flattened @ flattened.T - np.eye(adjoint.dim)))
    )
    completeness_error = float(
        np.max(np.abs(flattened.T @ flattened - traceless_projector))
    )
    equivariance_error = 0.0
    for row in range(3):
        for col in range(3):
            elementary = np.zeros((3, 3), dtype=float)
            elementary[row, col] = 1.0
            for source in range(adjoint.dim):
                matrix_action = (
                    elementary @ splitting[source]
                    - splitting[source] @ elementary
                )
                adjoint_action = np.einsum(
                    "b,bij->ij",
                    adjoint.generator(row, col)[:, source],
                    splitting,
                )
                equivariance_error = max(
                    equivariance_error,
                    float(np.max(np.abs(matrix_action - adjoint_action))),
                )
    if max(isometry_error, completeness_error, equivariance_error) > 2.0e-12:
        raise ArithmeticError(
            "adjoint/fundamental resolution failed: "
            f"isometry={isometry_error:.3e}, "
            f"completeness={completeness_error:.3e}, "
            f"equivariance={equivariance_error:.3e}"
        )
    return splitting
