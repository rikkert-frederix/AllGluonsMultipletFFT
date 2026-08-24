# Direct Specht, from zero group theory to the swap table

This document explains the `direct-specht` method used by the Wigner-table
generator. It assumes that you know ordinary matrices and vectors, but it does
not assume that you know group theory, particle physics, Young diagrams, or
Wigner symbols.

The short version is:

> Direct Specht replaces a very large calculation involving color tensors by
> a small calculation about moving six numbered boxes. It keeps exactly the
> box arrangements that behave like two gluons, literally swaps the first
> three boxes with the last three, and records that swap as a matrix.

The implementation is in
[`su3wigner/direct_specht.py`](su3wigner/direct_specht.py).

## 1. What problem is being solved?

In this project, a many-gluon color state is built one gluon at a time. There
can be several valid ways to reach the same final state. Those ways are called
**fusion paths**.

Suppose two neighboring gluons are called `A` and `B`. A path might first add
`A` and then `B`:

```text
left state --A--> middle state --B--> right state
```

After exchanging the two gluons, the same physical state is described by
paths that add `B` first and `A` second. The coefficients relating the two path
descriptions form a small square matrix:

```text
old path coefficients  --swap matrix-->  new path coefficients
```

That is the matrix stored in each `SwapBlock`. The complete table contains one
block for every valid pair of outer labels `(left, right)`.

This is a change-of-basis problem. It is not merely a yes/no question about
which states are allowed. The signs and square roots in the matrix matter.

## 2. Why are there several methods?

The most direct physics calculation constructs explicit color tensors and
contracts their indices. That is conceptually reliable, but those tensors grow
quickly.

`direct-specht` takes a different route. It uses the fact that the same
information can be encoded by permutations: instructions that exchange
numbered slots. This produces a small, independent calculation with no
explicit SU(3) carrier-space matrices or Clebsch--Gordan tensors.

Having independent methods is valuable. If a tensor contraction and a
box-permutation calculation produce the same thousands of matrices, including
all signs, it is unlikely that they share the same implementation mistake.

## 3. The first mental model: boxes on shelves

A **Young diagram** is a collection of left-aligned boxes arranged in rows.
For this SU(3) calculation, at most three rows are needed. We store a diagram
as its row lengths:

```text
(4, 2, 0) means

[] [] [] []
[] []
```

An SU(3) label is written `(p, q)`. Its smallest three-row diagram is

```text
(p + q, q, 0).
```

For example:

```text
label (2, 1)  ->  diagram (3, 1, 0)

[] [] []
[]
```

Adding a full column of three boxes does not change the SU(3) state relevant
here. Such a column represents a determinant factor, which is trivial for
SU(3). Therefore the same label may also be represented by

```text
(p + q + c, q + c, c)
```

for a nonnegative integer `c`. The code calls this a **lifted partition**.

You can think of the extra three-box columns as harmless packing material. We
sometimes need them so that the total number of boxes comes out correctly.

## 4. Why does one gluon become three boxes?

The color representation carried by a gluon is the SU(3) adjoint, commonly
called the `8` because it has eight components. In Young-diagram language it
has shape `(2, 1)`, which contains three boxes:

```text
[] []
[]
```

The method studies the exchange of two gluons, so it adds two groups of three
boxes:

```text
first gluon:   boxes 1, 2, 3
second gluon:  boxes 4, 5, 6
```

This is the origin of the “six-box space” mentioned in code and diagnostics.

## 5. Standard tableaux: legal histories of adding boxes

Write the numbers `1` through `6` in the six new boxes. A **standard tableau**
is a numbering that increases from left to right in every row and from top to
bottom in every column.

Instead of storing the picture, the code stores the row chosen for each new
box. For example,

```text
(0, 1, 0, 2, 1, 0)
```

means:

```text
box 1 went in row 0
box 2 went in row 1
box 3 went in row 0
box 4 went in row 2
box 5 went in row 1
box 6 went in row 0
```

Only histories that remain valid Young diagrams are retained.

A **skew tableau** describes only the new boxes between a starting diagram and
an ending diagram. This is crucial for performance: the method does not build
the enormous space associated with all old boxes. It tracks only the six boxes
that participate in the local swap.

For every block in the checked-in depth-6 table, this skew space has dimension
at most 90, even though the full physical representation can be much larger.

## 6. Turning box exchanges into matrices

Every adjacent exchange, such as swapping box `2` with box `3`, is represented
by a real matrix. The matrix acts on the list of standard skew tableaux.

The diagonal part depends on the **content** of a box:

```text
content = column number - row number.
```

If two consecutive boxes have content difference `r`, their two-dimensional
exchange block has the familiar orthogonal form

```text
[ 1/r                 sqrt(1 - 1/r^2) ]
[ sqrt(1 - 1/r^2)    -1/r              ]
```

In one-dimensional cases the answer is simply `+1` or `-1`.

These are Young's seminormal, or Young-orthogonal, permutation matrices. The
implementation constructs the five adjacent exchanges `(1 2)` through
`(5 6)`. Longer permutations are products of these small matrices.

The code checks the expected permutation rules numerically:

```text
swap_i squared = identity
swap_i swap_(i+1) swap_i = swap_(i+1) swap_i swap_(i+1)
```

The second identity says that two different sequences of neighboring moves
which perform the same rearrangement must agree.

## 7. Selecting the part that really is a gluon

Three ordinary boxes contain more possibilities than the gluon `8`. We need a
filter that keeps the adjoint-shaped part and rejects the determinant part.

In linear algebra, such a filter is a **projector**: a matrix `P` satisfying

```text
P @ P = P.
```

The determinant projector antisymmetrizes all three boxes. In plain language,
it keeps the part that changes sign whenever any two of the three boxes are
exchanged.

The adjoint projector keeps the antisymmetric pair in the last two slots but
subtracts the fully antisymmetric determinant part:

```text
adjoint = antisymmetric-pair projector - determinant projector.
```

Call the first gluon's adjoint projector `A1` and the second one's `A2`. They
act on different groups of boxes and commute. Their product

```text
Q = A1 @ A2
```

keeps precisely the six-box states describing two adjoints. The rank of `Q`
must equal the number of fusion paths in the block. The implementation checks
this for every block.

## 8. Middle states and prefix projectors

A path contains a middle label:

```text
left -> middle -> right.
```

After boxes `1, 2, 3` have been added, the current Young diagram identifies
that middle label. A diagonal **prefix projector** selects tableaux with the
requested intermediate diagram.

Restricting `Q` with this prefix projector splits the complete path space into
independent middle-label sectors. The rank of each sector says how many paths
run through that middle label.

The tables currently require only these patterns:

| sector rank | path multiplicities | meaning |
|---:|---|---|
| 1 | `(0,0)` | a unique path |
| 2 | `(0,0), (1,0)` | the first edge has two copies |
| 2 | `(0,0), (0,1)` | the second edge has two copies |
| 4 | `(0,0), (0,1), (1,0), (1,1)` | both edges have two copies |

Here a pair `(a,b)` gives the copy number on the first and second edge. An
unexpected pattern is treated as a structural error rather than guessed.

## 9. What does “multiplicity” mean?

Sometimes two genuinely different states have the same outer SU(3) label.
Labels alone cannot distinguish them. This is called **multiplicity**.

Imagine two roads that start and end in the same towns. Saying only the town
names does not tell you which road was taken, so the roads receive numbers `0`
and `1`.

Choosing any rotated pair of roads would be mathematically valid, but the
stored table and its Fortran consumer require one fixed convention. Direct
Specht constructs that convention locally.

### 9.1 Generic repeated self-edge

A repeated edge occurs when adding a gluon can leave the SU(3) label unchanged
in two independent ways. For a generic label:

- multiplicity `0` is the ordinary color-action direction;
- multiplicity `1` is its oriented orthogonal complement.

The action direction is selected with a **Jucys--Murphy operator**. The name is
more intimidating than the operation: in the tableau basis it is just a
diagonal matrix whose entries are box contents. The code sandwiches this
diagonal matrix between determinant and adjoint projectors. The result singles
out one direction without consulting color tensors or another backend.

The second direction is obtained by projecting coordinate vectors in a fixed
order, subtracting the first direction, normalizing, and making the first
nonzero component positive. That is the repository's deterministic
lexicographic rule.

### 9.2 The special `8 x 8 -> 8` edge

When the old state and the added gluon are both adjoints, the two copies have a
familiar matrix interpretation:

- multiplicity `0` is the even, symmetric Jordan-product channel;
- multiplicity `1` is the odd, antisymmetric Lie-bracket channel.

“Even” means unchanged under exchanging the two adjoints; “odd” means the
state changes sign. The bracket's operand order supplies one additional minus
sign. This special handling is necessary for exact agreement with the
repository convention.

## 10. Why signs need an explicit convention

A one-dimensional direction is physically unchanged if its vector `v` is
replaced by `-v`. A projector cannot see the difference because

```text
v v^T = (-v) (-v)^T.
```

Therefore projectors determine path spaces and absolute coefficient sizes,
but not every displayed sign. This freedom is often called a **gauge** or
**phase convention**. Here all numbers are real, so the phase is just `+1` or
`-1`.

The project has one established convention used by all table builders and by
the Fortran consumer. Direct Specht fixes it using:

1. the Jucys--Murphy action direction;
2. the lexicographically oriented complement;
3. the Jordan/bracket convention for `8 x 8 -> 8`;
4. a deterministic Pieri/determinant sign relating lifted three-row diagrams
   to the repository's Young basis.

This last sign is bookkeeping caused by representing the same SU(3) label
with diagrams that differ by full three-box columns. It does not change the
physical subspace, but it does change the coordinates used to print a matrix.

## 11. Constructing the path basis

For one local edge, the code constructs a three-box frame. A unique edge has
one normalized vector; a repeated self-edge has the two convention-fixed
vectors described above.

For a complete path

```text
left -> middle -> right
```

the first three-box edge vector is multiplied by the second three-box edge
vector. Every six-box tableau splits naturally into its first-three and
last-three row histories, so this product gives one column of the full path
basis matrix `V`.

The columns are placed in the exact `Path` order used by the table format. The
implementation verifies

```text
V.T @ V = identity
V @ V.T = Q
```

up to floating-point tolerance. The first check says that the path vectors are
orthonormal. The second says that they fill the complete two-adjoint subspace:
nothing is missing and nothing extra was included.

## 12. The literal swap

The six box labels begin in this order:

```text
1 2 3 | 4 5 6
```

Exchanging the two gluons means

```text
4 5 6 | 1 2 3.
```

The implementation builds this permutation, called `tau`, from nine adjacent
box exchanges. It checks that

```text
tau.T = tau
tau @ tau = identity
tau @ Q = Q @ tau.
```

The first two statements say that the exchange is a symmetric involution: do
it twice and you return to the starting point. The last says that swapping the
two three-box groups preserves the two-gluon subspace.

Finally, the desired swap matrix is simply

```text
W = V.T @ tau @ V
```

Read this from right to left:

1. `V` converts path coefficients into six-box tableau coordinates;
2. `tau` exchanges the two three-box groups;
3. `V.T` converts the result back into path coordinates.

That one line is the heart of Direct Specht.

## 13. One block, step by step

Given `left` and `right`, `DirectSpechtSwapOracle.block()` performs these
operations:

1. Enumerate all legal `left -> middle -> right` paths.
2. Convert the outer SU(3) labels to compatible three-row partitions.
3. Enumerate only the standard skew tableaux of the six new boxes.
4. Build Young's adjacent-transposition matrices from tableau contents.
5. Construct the two local adjoint projectors and their product `Q`.
6. Construct and validate the literal block swap `tau`.
7. Group paths by middle label and validate each observed rank pattern.
8. Build convention-fixed three-box frames for both edges of every path.
9. Assemble the orthonormal, path-ordered matrix `V`.
10. Compute `W = V.T @ tau @ V`.
11. Remove numerical near-zero noise and validate symmetry and involution.
12. Return a `SwapBlock` containing labels, paths, and `W`.

The full `DirectSpechtSwapTableBuilder` first finds all left labels reachable
within the requested depth, then repeats this block-local procedure for every
reachable right label.

## 14. What the method deliberately does not build

At runtime, Direct Specht does not use:

- explicit SU(3) representation matrices;
- explicit color-index tensors;
- Clebsch--Gordan embeddings;
- highest-weight nullspaces in carrier space;
- lower-complexity Wigner matrices;
- results from `dimension-only` or another table builder.

It uses only fusion-label rules, skew tableaux, permutation matrices,
projectors, contents, and deterministic local sign conventions.

## 15. Built-in checks

Many errors could still produce plausible-looking numbers, so the method
checks structural identities before accepting a result.

For the six-box kernel it checks:

- `Q` is symmetric and idempotent;
- `rank(Q)` equals the number of paths;
- `tau` is symmetric and squares to the identity;
- `tau` commutes with `Q`;
- the restricted swap has the correct norm;
- its trace gives an integral number of `+1` and `-1` eigenvalues.

For each middle sector it checks:

- the prefix restriction is a projector;
- its rank equals the number of paths through that middle label;
- its multiplicity pattern is one of the supported structural patterns.

For the final path basis and matrix it checks:

- the columns of `V` are orthonormal;
- the columns close exactly onto `Q`;
- `W` is finite and symmetric;
- `W @ W` is the identity.

The normal table generator can additionally check three-gluon braid
relations. Those relations verify that different sequences of neighboring
swaps produce the same final permutation.

## 16. Using the method

Generate a table from the `Wigner6j` directory:

```sh
python3 generate_table.py \
  --method direct-specht \
  --max-prefix-gluons 6 \
  --output /tmp/su3_adjoint_swap_prefix_6.tbl
```

Cross-check it against every other backend at an affordable depth:

```sh
python3 generate_table.py \
  --method direct-specht \
  --cross-check \
  --max-prefix-gluons 2 \
  --braid-check-depth 2 \
  --output /tmp/su3_adjoint_swap_checked.tbl
```

Build from Python:

```python
from su3wigner import DirectSpechtSwapTableBuilder

blocks = DirectSpechtSwapTableBuilder().build(6)
```

Inspect one block:

```python
from su3wigner import DirectSpechtSwapOracle

oracle = DirectSpechtSwapOracle()
block = oracle.block((1, 1), (1, 1))

print(block.paths)
print(block.matrix)
```

Request gauge-independent diagnostics:

```python
diagnostics = oracle.diagnostics((1, 1), (1, 1))

print(diagnostics.skew_dimension)
print(diagnostics.product_rank)
print(diagnostics.prefix_ranks)
print(diagnostics.swap_signature)
```

## 17. Performance intuition

The main cost is controlled by the number of legal arrangements of the six
new boxes, not directly by the dimension of the full SU(3) representation.
Matrices inside one block are dense but small, and bounded caches retain
recent tableau and local-frame data.

For the depth-32 comparison performed during development:

- 545 prefix labels were included;
- 9,766 blocks were built;
- 32,807 paths were present;
- 175,859 matrix entries were compared;
- Direct Specht agreed with `dimension-only` to a worst absolute difference
  of about `3.1e-15`.

The standard serialized depth-32 table is about 6.3 MiB. These measurements
are examples, not hard complexity guarantees.

## 18. A small glossary

- **Adjoint / `8`:** The eight-component color representation carried by a
  gluon. Direct Specht represents it with a three-box `(2,1)` projector.

- **Braid relation:** The statement that two sequences of neighboring swaps
  which perform the same overall permutation give the same matrix.

- **Content:** `column - row` for a box in a tableau. Contents determine the
  seminormal exchange matrices and Jucys--Murphy eigenvalues.

- **Fusion path:** One allowed sequence of intermediate labels when gluons are
  added one at a time.

- **Gauge / phase convention:** A deterministic choice between equally valid
  basis vectors such as `v` and `-v`. It changes printed matrix signs but not
  the underlying subspace.

- **Jucys--Murphy operator:** A sum of exchanges involving one selected box.
  In the tableau basis used here it is diagonal, with box contents as
  eigenvalues.

- **Multiplicity:** The occurrence of more than one independent state with
  the same outer label. An extra copy number distinguishes them.

- **Partition / Young diagram:** A left-aligned arrangement of boxes, stored
  here as three nonincreasing row lengths.

- **Projector:** A matrix `P` with `P @ P = P`. It filters vectors down to a
  chosen subspace.

- **Specht representation:** A way to represent permutations as matrices
  acting on standard tableaux. The code uses only the small skew version
  associated with the six new boxes.

- **Standard skew tableau:** A legal ordering for adding only the boxes between
  a starting and ending Young diagram.

- **Swap block:** The square matrix that exchanges two adjacent gluons for
  fixed left and right outer labels.

## 19. Where to go next

After this guide, the most useful references are:

- [`README.md`](README.md) for commands, table format, and cutoff rules;
- [`DERIVATION.md`](DERIVATION.md) for the full mathematical derivation and
  convention audit;
- [`PERFORMANCE.md`](PERFORMANCE.md) for scaling measurements;
- [`su3wigner/direct_specht.py`](su3wigner/direct_specht.py) for the executable
  implementation;
- [`tests/test_direct_specht.py`](tests/test_direct_specht.py) for compact,
  runnable examples of the invariants and reference comparisons.
