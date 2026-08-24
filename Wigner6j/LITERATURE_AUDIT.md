# Literature audit for the CG-free construction

The relevant representation-only route is developed in two papers:

- J. Alcock-Zeilinger, S. Keppeler, S. Plätzer and M. Sjödahl,
  [*Wigner 6j symbols for SU(N): Symbols with at least two quark-lines*](https://arxiv.org/abs/2209.15013),
  [J. Math. Phys. 64 (2023) 023504](https://doi.org/10.1063/5.0131538).
- S. Keppeler, S. Plätzer and M. Sjödahl,
  [*Wigner 6j symbols with gluon lines: completing the set of 6j symbols required for color decomposition*](https://arxiv.org/abs/2312.16688),
  [JHEP 05 (2024) 051](https://doi.org/10.1007/JHEP05(2024)051).

The first paper reduces the magnitudes of recouplings with two opposing
fundamental lines to representation dimensions.  The second splits adjoint
lines into fundamental--antifundamental pairs, obtains repeated adjoint
vertices by Gram--Schmidt from a small dimension-only Gram matrix, and
reduces the required gluon symbols to the fundamental-line symbols.  This is
enough to eliminate carrier-space Clebsch--Gordan tensors.  Fusion *labels*
are still required: a recoupling coefficient has no meaning without saying
which irreducible summands and multiplicity channels it connects.

The code does not transcribe the papers' displayed gluon formulae.  It derives
the same fixed-depth reduction in the repository's isometric-vertex gauge and
checks the completed matrices against the explicit index construction.  This
avoids silently importing a different closed-vertex normalization, a vertex
orientation, or a repeated-channel basis.

## Normalization conversion

The first paper normalizes each nonzero closed three-vertex graph to one.
Its symbols are therefore not the isometric coupling-tree matrices stored by
this repository.  If `a` and `b` are the input and output middle irreps and
$S_{a,b}$ is a unit-3j symbol in the paper's convention, then

$$
 W_{b,a}=\sqrt{d_b d_a}\,S_{a,b}.
$$

Comparing a displayed paper symbol directly with a table entry would import
spurious square roots of dimensions.  With this conversion, the magnitude in
the paper's Eq. (32) is exactly the off-diagonal magnitude of the seminormal
matrix below.

## A confirmed sign error in the quark-line paper

The magnitude theorem in Eq. (32) of the published paper (Theorem 1; Eq. (32)
in the arXiv version as well) is consistent with orthogonal recoupling.  The
linear sign identity Eq. (26d) and the sign argument built from it in
Appendix B.3 are not.

An explicit SU(3) counterexample starts from the Young diagram
`alpha=[3]`.  Adding boxes in rows 1 and 2 gives

| diagram | dimension |
|---|---:|
| `alpha=[3]` | 10 |
| `M_1=[4]` | 15 |
| `M_2=[3,1]` | 15 |
| `M^{11}=[5]` | 21 |
| `M^{22}=[3,2]` | 15 |
| `M^{12}=[4,1]` | 24 |

Equation (32) gives

$$
 d_1 S_{1,1}^{12}=\pm\sqrt{1-
 \frac{d_1d_2}{d_\alpha d_{12}}}=\pm\frac14,
 \qquad
 d_2 S_{2,2}^{12}=-d_1 S_{1,1}^{12}.
$$

Thus the two nonzero symbols must be `+/- 1/60` with opposite signs.
Equation (26d) for row 1 instead requires

$$
 1=\frac{21}{15}+24S_{1,1}^{12},
 \qquad S_{1,1}^{12}=-\frac1{60},
$$

whereas the same equation for row 2 requires

$$
 1=\frac{15}{15}+24S_{2,2}^{12},
 \qquad S_{2,2}^{12}=0.
$$

This is a direct contradiction.  The algebraic failure is also visible in
Appendix B: Eqs. (B17)--(B21) define `chi_ij=chi_ji` and a positive quantity
`A_ij` proportional to `d_ij/d_i`, then assert `A_ji=-A_ij`.  That
antisymmetry does not follow from the definition.  Equation (B22) therefore
builds a skew matrix from a non-skew quantity, and the inference in (B25)
does not repair the inconsistency.

The likely failed step is the removal of the barred/conjugated vertex just
before Eq. (26d).  Whatever its graphical diagnosis, Eq. (26d) and
Appendix B.3 are not used here.

## Independent sign prescription

Signs for two identical fundamental lines are fixed instead in Young's
orthogonal (seminormal) basis.  If the two added boxes have contents `c_1`
and `c_2`, set `r=c_2-c_1`.  On the two possible Young paths the normalized
adjacent swap is

$$
 \begin{pmatrix}
  1/r & \sqrt{1-1/r^2}\\
  \sqrt{1-1/r^2} & -1/r
 \end{pmatrix}.
$$

When only one path is admissible, `r=+/-1` and this reduces to its symmetric
or antisymmetric exchange parity.  The off-diagonal magnitude is the
dimension formula of Eq. (32) after converting the paper's unit-3j symbols
to normalized coupling trees.  Conjugation gives the antifundamental case.

This fixes signs locally, without the inconsistent global sign recursion.
The dimension-only backend then reconstructs every requested adjoint swap
matrix and compares it with the explicit-index backend, so the later gluon
paper is tested as a structural method rather than trusted as a source of
normalizations or signs.

The axial-content rule was also compared directly with all 56 nontrivial
two-fundamental terminal blocks for $0\leq p,q\leq4$; the worst binary64
difference was $1.33\times10^{-15}$.

## What is and is not taken from the gluon paper

The 2024 paper's fixed-depth chain is

$$
 \text{two-gluon symbols (PDF Eq. 34)}
 \longrightarrow
 \text{gluon--quark symbols (PDF Eq. 33)}
 \longrightarrow
 \text{two-quark symbols (2023 Eq. 32)}.
$$

For a repeated adjoint channel $\alpha\otimes A\to\alpha$, its PDF Eq. (27)
uses the small Gram matrix

$$
 s_{jk}=\frac1{N^2-1}
 \left(\frac{\delta_{jk}}{d_{\lambda_j}}
       -\frac1{N d_\alpha}\right).
$$

An ordered Cholesky or Gram--Schmidt factor fixes a multiplicity frame.  This
shows that carrier-space CG tensors are unnecessary, but that frame is not
automatically the repository's generator-action/complement frame.  In
particular, the paper's two triple-adjoint vectors must be rotated to the
repository order `mult0=d`, `mult1=f`; barred-vertex and reversal phases must
also be mapped explicitly.

The mixed $3/\bar3$ dimension--Casimir frame used by
`su3wigner/dimension_only.py` is an independent derivation in the repository
gauge.  It should not be attributed to a numbered formula in either paper.
The implementation uses the papers for the structural fixed-depth idea, then
checks all normalizations, rotations, and signs against literal permutation
operators.

## Audit of the gluon paper

One hard typographical error occurs in Appendix A of arXiv:2312.16688:

$$
 F\otimes\bar F=A\otimes\bullet
$$

must use a direct sum,

$$
 F\otimes\bar F=A\oplus\bullet.
$$

The next equality uses the direct sum correctly, so this typo does not
propagate.  No numerical or algebraic contradiction was found in its PDF
Eqs. (22), (27)--(29), (33)--(34), or (A8).  The experimental arXiv HTML
currently numbers the two reduction formulae differently from the PDF, so
PDF equation numbers are used here.

The remaining hazards are convention differences, not demonstrated paper
errors: unit-3j versus isometric-tree normalization, row-ordered
Gram--Schmidt versus `d/f` or action/complement multiplicity bases, vertex
reversal phases, and conjugate-arrow orientation.  Every one of these can
change signs or dimension factors while leaving a diagram topologically
similar.
