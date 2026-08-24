# SU(3) Wigner-6j / adjacent-gluon swap table

This directory generates the local coefficient used in Eq. (8) of
[*All-gluon amplitudes with off-shell recursion in multiplet bases*](https://arxiv.org/abs/2507.22636)
([JHEP 10 (2025) 168](https://doi.org/10.1007/JHEP10(2025)168)). The stored
object is the **complete normalized adjacent-adjoint swap matrix**, not a bare
closed tetrahedral graph. Dimensions, three-vertex normalizations,
multiplicity labels, vertex reversals, phases, and the antisymmetric
three-gluon sign are already included.

The generated ASCII tables are consumed by
[`AmpliGluonMultiplet`](../AmpliGluonMultiplet/README.md), whose Fortran
reader is [`src/wigner_table.f90`](../AmpliGluonMultiplet/src/wigner_table.f90).
The mathematical derivation and convention audit are in
[DERIVATION.md](DERIVATION.md); measured scaling and backend tradeoffs are in
[PERFORMANCE.md](PERFORMANCE.md). A beginner-oriented explanation of the
permutation-based backend is in [DIRECT_SPECHT.md](DIRECT_SPECHT.md).

## Generate a table

The code is source-tree Python rather than an installed package. The commands
below run from this directory; invoking `/path/to/Wigner6j/generate_table.py`
also works from another directory, while importing `su3wigner` there requires
adding `Wigner6j` to `PYTHONPATH`.

```sh
cd Wigner6j
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 generate_table.py \
  --method dimension-only \
  --max-prefix-gluons 2 \
  --output /tmp/su3_adjoint_swap_prefix_2.tbl
```

`dimension-only` is the recommended production builder. The explicit output
path keeps this example from replacing a checked-in reference table.

The CLI defaults are:

| option | default | meaning |
|---|---|---|
| `--method` | `explicit-index` | audit-oriented coefficient backend |
| `--max-prefix-gluons` | `2` | largest minimum prefix depth included |
| `--braid-check-depth` | `1` | check prefix labels whose minimum depth is at most 1 |
| `--zero-tolerance` | `5e-14` | omit serialized entries with absolute value at most this value |
| `--cross-check` | off | do not build the other algorithms |
| `--output` | `data/su3_adjoint_swap_prefix_N.tbl` beside this script | destination table |

Use `--braid-check-depth -1` to skip braid checks; values below `-1` are
rejected. A braid depth smaller than `--max-prefix-gluons` is a deliberate
partial check, not validation of every generated prefix layer.

The generator writes a temporary sibling of the destination, reads it back,
validates the restored blocks, and atomically replaces the destination only
after every check succeeds. Any build, cross-check, braid, serialization, or
round-trip validation failure before that replacement preserves an existing
table; normal exception cleanup removes the temporary directory. Calling the
lower-level `write_table()` function directly does not provide that
transaction.

## The five constructions

All five methods return complete signed tables in the repository gauge:

| CLI method | construction | intended use |
|---|---|---|
| `explicit-index` | contracts the complete adjoint-index coupling trees and checks Eq. (8) on every explicit tensor index | strongest local carrier-space audit |
| `fundamental-split` | resolves each adjoint into traceless `3 x 3-bar`, cancels the common split isometry, and evaluates one highest-weight column | faster independent coefficient contraction |
| `recursive-reduction` | expands adjoint vertices into elementary fusion paths and composes terminal two-line recouplings | independent reduction hierarchy |
| `dimension-only` | uses analytic fusion labels, dimensions, Casimirs, and Young-content swaps without irreps, generators, or CG tensors | deep production tables |
| `direct-specht` | works in exact six-box skew-tableau spaces and applies the literal block permutation, with local Jucys--Murphy convention frames | independent symmetric-group backend |

`DirectSpechtSwapOracle` provides the same construction one block at a time,
plus gauge-independent projector diagnostics for invariant audits. See
[Direct Specht, from zero group theory to the swap table](DIRECT_SPECHT.md) for
a step-by-step introduction.

The three carrier-space methods use different coefficient algorithms, but a
CLI `--cross-check` intentionally shares one already validated
`AdjointCouplings` object among them. This saves the cost of constructing the
same irrep and CG input three times. `dimension-only` and Direct Specht remain
structurally independent of that carrier data.

Run all affordable checks at a small cutoff with:

```sh
python3 generate_table.py \
  --method dimension-only \
  --cross-check \
  --max-prefix-gluons 2 \
  --braid-check-depth 2 \
  --output /tmp/su3_adjoint_swap_cross_checked_2.tbl
```

`--cross-check` builds every other full-table method and rejects any block-key,
path, or coefficient mismatch above `5e-11`. It becomes carrier-space limited
at deep cutoffs; it is not intended for a large production run.

## Choosing the cutoff

`--max-prefix-gluons N` includes every left irrep reachable from the singlet
after at most `N` adjoint fusions. Each local block also includes every right
irrep reachable after the two adjoints being swapped, so its right label may
have minimum depth `N+2`.

For the current Fortran consumer, a table used with `G` total gluons must obey

```text
MAX_PREFIX_GLUONS >= floor(G/2).
```

A prefix-`p` table therefore supports the consumer through `2*p+1` total
gluons. See the consumer's [table-depth explanation](../AmpliGluonMultiplet/README.md#table-depth)
for the distinction between minimum fusion depth and position in one path.

A fully braid-checked depth-6 production build is:

```sh
python3 generate_table.py \
  --method dimension-only \
  --max-prefix-gluons 6 \
  --braid-check-depth 6 \
  --output /tmp/su3_adjoint_swap_prefix_6.tbl
```

Checked-in reference tables for depths 2 through 6 are in [`data/`](data/).
To intentionally regenerate one, point `--output` at that tracked path; the
atomic replacement rule still applies. Numerical agreement is expected, but
different NumPy/BLAS eigensolvers can change the last few binary64 bits, so a
regenerated file need not be byte-identical.

## What is validated

Every full-table builder rejects non-finite, nonsquare, nonsymmetric, or
nonorthogonal blocks and checks the adjacent swap is an involution. Additional
checks depend on the selected construction:

- `explicit-index` checks normalized coupling-tree bases and the relative
  Frobenius residual of Eq. (8) over every explicit tensor index;
- `fundamental-split` checks the normalized highest-weight equation, promoted
  to the target irrep by already validated intertwiners and Schur's lemma;
- `recursive-reduction` checks terminal swaps, highest-weight vertex
  projections, elementary pair permutations, and final path-space closure;
- `dimension-only` checks exact fusion dimension sums, analytic elementary
  frames, vertex normalization, pair permutations, and projected closure; and
- Direct-Specht diagnostics check projector/prefix ranks, symmetry,
  idempotence, literal-swap involution and commutation, restricted norm, and
  swap signature.

When enabled, braid validation checks `s1*s2*s1 = s2*s1*s2` and both local
involutions in every selected three-adjoint `(left,right)` sector; the CLI
rejects a worst residual above `5e-8`. Serialization preserves block keys and
ordered paths exactly. Values omitted by `--zero-tolerance` return as zero, so
the accepted round-trip matrix difference is
`max(zero_tolerance, 1e-13)`; restored blocks are revalidated with tolerance
`5e-10`.

The committed regression suite currently:

- executes the all-backend CLI cross-check through cutoff 3;
- compares the complete `dimension-only` result with every checked-in table
  from cutoff 2 through 6;
- checks all prefix sectors through depth 16 for braid and involution
  residuals; and
- compares all 354 signed Direct-Specht cutoff-6 blocks with the checked-in
  table and audits their projector diagnostics.

Manual deep validation runs have additionally reached depth 8 with
`fundamental-split` and depth 7 with `explicit-index` and
`recursive-reduction`. These are empirical reach results, not routine CI.
See [PERFORMANCE.md](PERFORMANCE.md#correctness-reach-versus-committed-regression-coverage)
for counts and residuals.

## Memory behavior

Carrier-space construction restricts lowering to exact weight sectors and
retains validated CG maps as packed weight blocks. Non-adjoint ambient
embeddings are lazy and their transforms are packed. The public
`Coupling.embedding`, `ElementaryCoupling.embedding`, and
`Irrep.ambient_embedding` dataclass fields remain ordinary NumPy arrays, but
reading them—including indirectly through `repr`, equality, `asdict`,
`replace`, pattern matching, or pickle—materializes and caches the dense form.

Backend working sets are reduced, not globally constant:

- `explicit-index` targets 32 MiB per complete-chain slab, but decoded operands
  and overlap/reconstruction temporaries are additional;
- `fundamental-split` contracts only relevant entries of one highest-weight
  column;
- `recursive-reduction` packs elementary CGs and bounds its elementary LRU,
  while completed reductions and blocks remain cached; and
- `dimension-only` bounds support caches and avoids dense path-swap chains,
  while its requested output blocks necessarily grow with the table.

The practical asymptotic limit of the three carrier-space methods is dense
irrep/generator construction and validation. `dimension-only` is the deep
production method. Exact measurements, cache tradeoffs, and reproduction
commands are in [PERFORMANCE.md](PERFORMANCE.md).

## Run tests and use the Python API

The declared runtime requirement is Python 3.10 or newer with the only
non-standard-library Python dependency `numpy>=1.26,<3`. This declared range
is not a bit-reproducible lock file or a promise that every version has a
dedicated CI job. The reported benchmarks use Python 3.12.3 and NumPy 1.26.4.

```sh
cd Wigner6j
python3 -m unittest discover -v
```

Programmatic construction and loading return ordered tuples of `SwapBlock`s:

```python
from su3wigner import DimensionOnlySwapTableBuilder, read_table

generated = DimensionOnlySwapTableBuilder().build(6)
loaded = read_table("data/su3_adjoint_swap_prefix_6.tbl")
```

`DirectSpechtSwapOracle().diagnostics(left, right)` returns the gauge-
independent diagnostic record for one valid block;
`DirectSpechtSwapOracle().block(left, right)` returns the signed `SwapBlock` in
the repository gauge for every valid key.

## Table lookup contract

A block is identified by its outer irreps

```text
left=(p_L,q_L), right=(p_R,q_R).
```

Its ordered paths are

```text
left --8[left_mult]--> middle --8[right_mult]--> right.
```

The path order defines both matrix axes. For fixed outer labels, stored value
rows follow

```text
block  out_path  in_path  coefficient
```

with the convention

```text
B(g1,g2)[in_path] = sum_out W[out_path,in_path] B(g2,g1)[out_path].
```

Indices and multiplicities are zero-based. Dynkin labels identify SU(3)
irreps: `8=(1,1)`, `10=(3,0)`, `10bar=(0,3)`, and `27=(2,2)`. A parity of zero
means “not applicable”; `+1/-1` is recorded when both inputs to a vertex are
adjoints. In particular, `8 x 8 -> 8` multiplicity 0 is the symmetric Jordan
product (`+1`) and multiplicity 1 is the antisymmetric Lie bracket (`-1`). For
other `R x 8 -> R` self channels, multiplicity 0 is the normalized generator
action; multiplicity 1, when present, is its canonically oriented complement.

### ASCII V1 layout

The canonical file order is:

```text
SU3_ADJOINT_SWAP_TABLE_V1
NC 3
ADJOINT_PQ 1 1
MAX_PREFIX_GLUONS N
NBLOCKS number_of_block_rows
NPATHS number_of_path_rows
NVALUES number_of_stored_value_rows
BEGIN_BLOCKS
block left_p left_q right_p right_q size
...
END_BLOCKS
BEGIN_PATHS
block path middle_p middle_q left_mult right_mult left_parity right_parity
...
END_PATHS
BEGIN_VALUES
block out_path in_path coefficient
...
END_VALUES
END_TABLE
```

The magic string must occupy physical line 1. After that line, blank lines and
lines whose first nonspace character is `#` are ignored. Block IDs and
per-block path indices are contiguous and zero-based. Each block has exactly
`size` path rows. Values are sparse coordinates: an omitted
`(out_path,in_path)` entry is exactly zero, and `NVALUES` counts only rows that
were written after applying `zero_tolerance`. Coefficients are formatted as
`%.17e`, i.e. 17 digits after the decimal point and 18 significant decimal
digits, which is more than sufficient to round-trip binary64.

Python `read_table()` validates the V1 syntax, required metadata, declared row
counts, index bounds, and finite coefficients, then returns dense
`SwapBlock`s. It does not return `MAX_PREFIX_GLUONS` and does not rediscover
fusion reachability, enforce semantic uniqueness of outer/path labels, or
rerun braid/coherence checks. Generated tables satisfy those stronger
conditions; use generator output or the regression validators when those
semantic guarantees are required. The Fortran reader loads the cutoff and
blocks and performs its own symmetry/orthogonality checks; its exact parser
contract is defined in
[`wigner_table.f90`](../AmpliGluonMultiplet/src/wigner_table.f90).
