# Colour and recursion derivation

This note fixes the conventions used by the Fortran implementation.  The
recursion idea comes from the Multiplet Recursion work, but the coefficient,
orientation, and normalization rules below were derived directly from the
isometric Clebsch--Gordan tensors used to generate the local `Wigner6j` table.
They were then checked by explicit adjoint-index contraction before being
translated to Fortran.

## 1. Open multiplet paths

For a current containing `k` gluons, start with the singlet and successively
fuse `k` adjoints.  Only paths ending in an adjoint can connect to the rest of
a pure-gluon tree.  A path records every intermediate Dynkin label `(p,q)` and
every outer-multiplicity label.

Each three-vertex in `Wigner6j` is an isometric embedding.  Consequently, for
a fixed open path, the map from the final adjoint to the tensor product of the
`k` external adjoints has orthonormal columns.  Distinct paths are orthogonal.
Taking the trace over the open adjoint gives norm squared `dim(8)=8` for each
closed colour tensor.

The path counts ending in an adjoint are

```text
k = 1, 2, 3, 4,   5,   6,    7,     8
    1, 2, 8, 32, 145, 702, 3598, 19280.
```

These counts are discovered from the transitions present in the table, not
stored as special cases.

### Wigner-table cutoff used by the consumer

Close a `k`-gluon current with its open adjoint and label the representation
after `j` external adjoints by `R_j`.  The resulting colour chain has `k+1`
adjoints and singlets at both ends.  Reading the suffix from the other end
gives a fusion path to the conjugate of `R_j`.  Since the adjoint is
self-conjugate, `R_j` and its conjugate have the same minimum adjoint-fusion
depth `d`, and hence

```text
d(R_j) <= min(j, k+1-j).
```

The largest value of this bound over the chain is `floor((k+1)/2)`.  Every
outer representation used by the current one-sided path catalog and local
swap lookup is one of these cuts.  The table generator's
`MAX_PREFIX_GLUONS` selects labels by minimum fusion depth and, for each
selected left outer label, includes every possible result after the two
gluons being swapped.  The consumer therefore enforces the sufficient prefix
cutoff `floor((k+1)/2)` for a `k`-gluon current.  For the top current of an
amplitude with `N=k+1` total gluons this becomes

```text
required prefix depth = floor(N/2).
```

Equivalently, the current implementation accepts a prefix-`p` table for
amplitudes through `2*p+1` total gluons.  This is a conservative contract for
this path and swap construction, rather than a minimality statement about all
possible recoupling algorithms.  It concerns the Wigner data only and does
not change the number of complete multiplet paths.

## 2. Local adjacent swaps

For fixed irreps immediately outside two neighbouring gluons, the table gives
the complete real change-of-basis matrix `W` with the convention

```text
B(g1,g2)[input] = sum_output W(output,input) B(g2,g1)[output].
```

The Fortran loader validates every local matrix as a symmetric orthogonal
involution.  It then embeds each local block in the complete open-path space,
producing one sparse matrix for every possible adjacent position.  The test
suite independently checks that these full-path matrices are symmetric
involutions too.

## 3. Joining two currents

Let the left and right child paths both end in an adjoint.  The colour part of
a gluonic binary merge is the antisymmetric `8 x 8 -> 8` vertex, which is
multiplicity 1 in the table.

Draw the left chain toward the central vertex and reflect the right chain away
from it.  Reading the resulting closed chain from the singlet gives a seed
path consisting of:

1. the complete left path;
2. the central adjoint with multiplicity 1; and
3. the conjugated intermediate irreps of the right path in reverse order,
   with its multiplicities reversed.

Reflecting a child is convention-sensitive when complex irreps or repeated
self channels occur.  In the real, isometric vertex convention fixed by the
table generator, that sign can be obtained locally.  Write a path as vertices

```text
R_(j-1) x 8 -> R_j  with outer multiplicity mu_j.
```

Reflection replaces each vertex by its charge-conjugate vertex.  The canonical
vertex gauges of `Wigner6j` give the following factors:

- a non-self transition has factor `+1`;
- `8 x 8 -> 8` has factor `(-1)^mu`, since multiplicities zero and one are
  respectively the Jordan-product and Lie-bracket vertices; and
- every other self transition `R x 8 -> R` has factor `(-1)^(mu+1)`.
  Multiplicity zero is the generator-action vertex; its minus sign follows
  from the contragredient identity
  `rho_conjugate(X) = -transpose(rho_R(X))`.  The canonically oriented
  orthogonal copy has the opposite reflection eigenvalue.

Multiplying the local factors gives

```text
rho(path) = (-1)^(
    sum_j mu_j
    + number of self transitions R_(j-1) = R_j with R_j /= 8).
```

The reflected path itself is still constructed and looked up explicitly.  If
the original open path has length `L`, then for `j=2,...,L` its labels are

```text
R'_j  = conjugate(R_(L+1-j)),
mu'_j = mu_(L+2-j),
```

with the unique initial `1 -> 8` vertex at `j=1`.  Initialization checks that
this map exists, squares to the identity, and assigns the same sign to both
members of every reflected pair.  It now costs `O(P*L)` operations and
`O(P+L)` memory for `P` paths.  A regression retains the former complete-swap
construction only through tractable open length six and compares every sign;
the two constructions agree there path by path.  Fixed positive/negative sign
counts and order-sensitive checksums cover the complete length-eight space
without recreating the old quadratic workspace there.

If this computed reflection sign is `rho`, the seed coefficient is `-rho`.
The minus sign is the parity of the one antisymmetric central vertex reversed
by the reflection.  A direct program which contracted every adjoint index for
all merges through four external gluons agreed with this construction to a
maximum absolute residual of `2.4e-15`.

Before sorting, the external order of a merge is

```text
left legs in ascending order, parent marker, right legs in descending order.
```

Repeated adjacent swaps transform this into canonical ascending order.  All
unordered partitions of the same external subset share intermediate orderings,
so initialization interns them into one directed acyclic graph.  Runtime
sorting is therefore a sparse pass through the graph, not a separate bubble
sort for every diagram.

The number of inversions supplies a natural layering of this graph: every
sorting edge goes from layer `k` to layer `k-1`.  Let `V_k` be the number of
interned nodes in layer `k`, `W = max_k V_k`, and `P` the number of colour
paths.  Runtime keeps two reusable coefficient buffers, each with `P*W`
entries, rather than one `P`-entry current for every graph node.  Seed currents
are inserted into the appropriate destination layer before its incoming
sorting edges, so merging and floating-point accumulation follow the original
all-node evaluation order.  Only the active slots of each layer are cleared;
over a complete pass the clearing work remains proportional to
`P*sum_k V_k`, while persistent sorting storage is reduced from
`O(P*sum_k V_k)` to `O(P*W)`.

## 4. Colour normalization

Resolve an adjoint state as an orthonormal traceless `3 x 3` matrix.  Acting
with the ordinary matrix commutator on two such states gives an adjoint map.
Direct contraction of that map with its transpose yields `6` times the
identity.  The multiplicity-1 table vertex is isometric, so the ordinary colour
bracket is `sqrt(6)` times that normalized vertex.

Every binary QCD merge therefore contributes

```text
strong_coupling * sqrt(6)
```

in addition to the normalized seed and Wigner coefficients.  This result is
obtained from the table's own adjoint realization and is not a hard-coded
translation to a different generator convention.

## 5. Kinematic recursion

All momenta are converted internally to the all-outgoing convention: the two
physical incoming momenta are negated and the final-state momenta are not.
External polarization vectors and the colour-ordered Lorentz kernels use the
same phase and metric conventions as AmpliCol.

Two current species suffice:

- a four-component gluon current;
- a six-component antisymmetric auxiliary-tensor current.

The auxiliary field decomposes the four-gluon contact term into two binary
vertices.  It has an identity, nonpropagating connection.  Thus every runtime
operation combines exactly two child currents.  For one unordered partition,
the possible combinations are

```text
gluon + gluon  -> gluon
gluon + gluon  -> tensor
tensor + gluon -> gluon
gluon + tensor -> gluon.
```

Both the colour bracket and each corresponding kinematic rule reverse sign
when the two sides are exchanged.  Their product is symmetric, which proves
that one orientation of each unordered bipartition is sufficient.
AmpliGluonMultiplet chooses the orientation whose lowest-ranked external leg
lies on the left.

Currents are evaluated by increasing subset cardinality.  Proper gluon
currents receive the massless propagator; the current containing all but the
closing external gluon does not.  Tensor currents never receive a propagator.
The top tensor current cannot close the amplitude and is not constructed.

## 6. Closure and squared matrix element

The current containing the first `N-1` external gluons is contracted with the
polarization of gluon `N`.  As established in Section 1, its closed path tensor
has norm squared 8.  Multiplying every raw path coefficient by `sqrt(8)` gives
amplitudes in a unit-normalized colour basis.

The full-colour result is consequently the ordinary Euclidean norm of the
complex basis-amplitude vector.  No colour Gram matrix remains at runtime.
Initial-colour averaging divides this sum by `8^2`; helicities remain fixed.

## 7. MHV radiation inside the diagonal multiplet basis

The all-outgoing MHV and anti-MHV sectors have exactly two helicities of one
sign. Their Parke--Taylor insertion identity can be applied without changing
colour basis. Let a closed source basis contain `m+1` gluons, represented here
by an open path of length `m`. Define `E_m` by copying a source path and
appending the normalized antisymmetric transition

```text
8 x 8 -> 8,  outer multiplicity = 1.
```

If `S_j^(m)` is the already validated adjacent-swap operator on the source
space, radiation from old leg `i` is the matrix-free map

```text
R_i = (target swaps back) E_m (source swaps i to the closure)
    = S_i^(m+1) ... S_(m-1)^(m+1)
      E_m
      S_(m-1)^(m) ... S_i^(m).
```

Both sides of this equation are coordinates in normalized multiplet path
bases. No trace, DDM, adjoint-index, or non-diagonal colour representation is
introduced. Since the swaps and `E_m` are isometries, each `R_i` is an
isometry. The implementation additionally checks the complete colour-charge
identity

```text
sum_i R_i = 0
```

column by column through target path length seven.

Choose one minority-helicity leg `r`, the other minority-helicity leg, and one
majority-helicity leg as the three-gluon seed. For each remaining
majority-helicity leg `q`, the angle-bracket MHV update used by the code is

```text
x_(m+1) = sum_(i /= r) <i r>/(<i q><q r>) R_i x_m.
```

The anti-MHV update replaces angle brackets by square brackets. This is an
exact inverse-soft form of the Parke--Taylor identity, not a soft-limit
approximation. After all insertions, ordinary full-path swaps restore the
canonical external ordering expected by the public basis vector.

The unit seed is the unique length-two path ending through the antisymmetric
adjoint copy. The public normalization is

```text
sqrt(8) * (strong_coupling*sqrt(6))^(N-2)
```

times the formal three-point helicity amplitude and the common phase relating
the analytic spinor chart to `external_massless_vector`. A single serial,
planar Berends--Giele amplitude determines that phase as
`A_BG(1,...,N)/A_PT(1,...,N)`. It supplies no colour vector or colour
contraction and costs polynomial rather than factorial work. Keeping the full
complex ratio makes crossing, polarization, and `i` conventions automatic.

If any external pair is too close to a collinear configuration, or a spinor
denominator is too small for a stable inverse-soft evaluation, the evaluator
uses the general multiplet Berends--Giele recursion below. The pair guard is
needed because the single planar phase calibration cannot see every
nonadjacent collinear channel. The construction follows the
multiplet-radiation idea of
[Du, Sjodahl and Thoren, arXiv:1503.00530](https://arxiv.org/abs/1503.00530),
but derives the maps from this repository's sequential path basis rather than
importing that paper's balanced-basis matrices.

## 8. NMHV path-basis prototype

The same radiation maps also give an exact colour-dressed CSW construction
for NMHV and anti-NMHV amplitudes without leaving the diagonal multiplet
basis. For each admissible partition, both off-shell MHV vertices are built
as multiplet vectors by the recursion above. Contracting their internal
adjoint amounts to a direct path stitch: retain the left path, append the
conjugate-reversed right path with its reflection phase, and use the existing
Wigner swaps to restore canonical external order. The channel normalization
is

```text
V_left * V_right / (sqrt(8) * P_channel^2).
```

No ordered-colour vector or colour Gram matrix enters this construction. An
independent prototype reproduced the complete generic multiplet coefficient
vector for every NMHV and anti-NMHV helicity placement at six and seven
gluons, with worst relative residuals `4.1e-14` and `2.7e-13` respectively.

The direct prototype is intentionally not a production fast path: sorting a
full top-level vector independently for all
`3*(2^(N-3)-1)` channels made it six to eight times slower than the existing
generic recursion through nine gluons. A competitive implementation must
memoize repeated off-shell MHV subset vectors, precompute the signed stitch
maps, inject all channel seeds into shared ordering nodes, and run one
top-level layered Wigner DAG. This is the basis-preserving NMHV direction;
replacing it by DDM partial amplitudes would change the public colour basis
and is not an implementation option here.

## 9. Independent checks

The convention was checked at three levels:

1. explicit SU(3) Clebsch--Gordan index contractions versus the seed and swap
   construction;
2. structural Fortran tests of every path count, reflection, and seed through
   open length eight, and dense sparse-swap properties through length seven,
   including direct reflection signs against the former complete-swap
   construction through open length six plus fixed sentinels for all 19,280
   length-eight paths; and
3. full fixed-helicity matrix elements versus AmpliCol's independent
   colour-flow amplitudes and colour Gram matrix for total multiplicities
   four, five, six, seven, and eight.

The MHV path additionally has an exhaustive column-by-column test of colour
conservation and isometry for the radiation operators through target path
length seven. Independent event fuzzing compared 7,705,350 complex multiplet
coefficients, covering all MHV and anti-MHV minority-leg placements through
nine gluons and positive, negative, and zero couplings. A permanent
near-collinear event verifies that the numerical guard enters the generic
multiplet fallback before inverse-soft cancellation becomes inaccurate.

The last comparison exercises interference between cubic diagrams and the
auxiliary representation of the four-gluon interaction, so it also fixes their
relative phase rather than only an overall normalization.
