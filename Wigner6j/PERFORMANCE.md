# Wigner-table scaling and backend tradeoffs

This document describes the performance of the five SU(3) adjacent-adjoint
swap constructions in this directory. All five are complete signed table
builders in the repository gauge. The mathematical constructions are derived in
[DERIVATION.md](DERIVATION.md).

`MAX_PREFIX_GLUONS`, called *depth* below, includes every left irrep whose
minimum number of adjoint fusions from the singlet is at most that value. The
right irrep of a local two-adjoint block can consequently have minimum depth
up to `depth + 2`. See the [README cutoff definition](README.md#choosing-the-cutoff)
and the [Fortran-facing discussion](DERIVATION.md#8-finite-table-cutoff-and-fortran-facing-format).

## Choosing a method

Use `dimension-only` for production tables and use carrier-space methods at
affordable cutoffs as independent coefficient audits.

| construction | complete signed table | principal retained data | main limit |
|---|---:|---|---|
| `explicit-index` | yes | dense irrep generators, packed CG sectors, completed blocks, 64 MiB decode LRU | full two-adjoint chain contractions and replay |
| `fundamental-split` | yes | dense irrep generators, packed CG sectors, completed blocks | carrier-space irrep/CG construction and validation |
| `recursive-reduction` | yes | dense irrep generators, packed adjoint and elementary CG sectors, completed reductions/blocks | carrier hierarchy plus elementary reductions |
| `dimension-only` | yes | label graph, bounded support caches, completed blocks | growth of the requested output table |
| `direct-specht` | yes | bounded skew-tableau, local-frame, projector, and recent-block caches | block-local eigensystems in the six-box skew space |

The three carrier-space coefficient algorithms differ, but they share the
same phase-fixed `AdjointCouplings` object during a CLI `--cross-check`. This
avoids rebuilding identical CG input and means the cross-check is not an
independent repetition of representation construction. `dimension-only` and
Direct Specht do not consume those carrier-space objects.

Carrier builders intentionally retain validated dense irrep generators,
packed adjoint decompositions, and requested results in unbounded instance-
method caches. Recursive reduction also retains terminal, vertex, path-group,
and final-block results; only its elementary-decomposition LRU is bounded.
These policies make one increasing-cutoff build incremental, but transient
builders in a long-lived process can retain their graphs. Run large,
independent carrier audits in fresh processes.

The `dimension-only` and Direct-Specht *support* caches are bounded, but their
public result caches still retain requested blocks or diagnostics. The output
itself grows approximately quadratically with depth and cannot be bounded
without streaming it to a consumer.

## Why the carrier methods stop first

For a balanced label `(p,q)` of scale `n`, the irrep dimension is `O(n^3)`.
Each of its nine dense generator matrices therefore contains `O(n^6)` values.
Construction also uses a harmonic basis whose ambient dimension is `O(n^4)`,
so its transient ambient-to-irrep array contains `O(n^7)` values. The cached
source set covers an approximately two-dimensional region of labels.

Packed CG storage removes zero-filled weight-forbidden entries, but it cannot
remove dense irrep construction, generator validation, or the dense temporary
work used to validate newly constructed CG maps. These shared operations are
the eventual limit for all three carrier backends. At smaller cutoffs the
explicit backend additionally spends substantial time contracting and
replaying its full two-adjoint chains.

By contrast, one local adjoint swap has a fixed fusion neighborhood.
`dimension-only` performs sparse actions on a few path columns, and the
Direct-Specht skew representation had dimension at most 90 in every block
tested through depth 64. Their observed work therefore follows the number of
requested labels and blocks rather than carrier dimension.

## Implemented reductions

Shared carrier-space changes are:

- descendant lowering actions and orthogonalization restricted to exact
  coordinate-weight sectors, while retaining the complete output columns;
- an all-incoming-equations orthogonal Procrustes solve for multiply connected
  target weight spaces, preventing deep normalization drift;
- exact weight-sector packing of validated adjoint and elementary CG maps;
- lazy non-adjoint ambient embeddings and packed harmonic-to-irrep transforms;
  and
- target-local construction, validation, and packing so dense copies of
  unrelated target irreps do not coexist.

The exported `Coupling.embedding`, `ElementaryCoupling.embedding`, and
`Irrep.ambient_embedding` fields still produce ordinary mutable NumPy arrays.
Reading those fields materializes and caches the dense form. Dataclass
operations that read public fields—such as `repr`, equality, `asdict`,
`replace`, positional pattern matching, and pickle—can materialize them too.
Production code uses private packed slab/column accessors and does not trigger
that promotion.

Backend-specific changes are:

- `explicit-index` targets 32 MiB for each all-path outer-index chain slab. It
  retains a block tensor when the complete block fits one slab and replays
  multi-slab blocks for the index equation. This target bounds the chain
  payload, not peak RSS: decoded operands, overlap arrays, reconstruction
  arrays, and residuals are separate. A builder-local dense-decode LRU retains
  at most 64 MiB by default, although two decoded operands can coexist during
  one contraction.
- `fundamental-split` analytically cancels the common split isometry and
  contracts only nonzero packed entries contributing to one normalized
  highest-weight column.
- `recursive-reduction` uses highest-weight Schur projections for terminal and
  vertex overlaps, packed elementary CG maps, source slabs capped at 8 MiB,
  canonical reverse swaps, shared multiplicity reductions, and a 256-entry
  elementary-decomposition LRU.
- `dimension-only` applies four elementary swaps directly to expansion columns
  without forming or caching full path-space swap chains. Its structural,
  frame, vertex, elementary-block, and fusion-path support caches are bounded.
- Direct Specht applies sparse seminormal actions to build small dense local
  projector/swap matrices, uses boolean prefix masks, avoids a dense
  Jucys--Murphy matrix, and bounds its edge, path, tableau, signed-block, and
  diagnostics caches.

## Output growth

The following counts come from the final `dimension-only` implementation.
`Dense entries` is `sum(block.size**2)` before sparse serialization;
`in-memory nonzeros` is `np.count_nonzero` after the builder's roundoff
cleaning and before applying a user-selected output tolerance.

| depth | reachable irreps | blocks | paths | dense entries | in-memory nonzeros |
|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 5 | 6 | 8 | 6 |
| 1 | 2 | 13 | 35 | 153 | 119 |
| 2 | 5 | 46 | 137 | 629 | 557 |
| 3 | 8 | 89 | 291 | 1,453 | 1,317 |
| 4 | 13 | 162 | 523 | 2,623 | 2,431 |
| 5 | 18 | 243 | 805 | 4,139 | 3,875 |
| 6 | 25 | 354 | 1,165 | 6,001 | 5,673 |
| 7 | 32 | 473 | 1,575 | 8,209 | 7,801 |
| 8 | 41 | 622 | 2,063 | 10,763 | 10,283 |
| 16 | 145 | 2,454 | 8,215 | 43,651 | 42,403 |
| 32 | 545 | 9,766 | 32,807 | 175,859 | 172,307 |
| 64 | 2,113 | 38,982 | 131,143 | 706,003 | 694,771 |

Doubling depth from 32 to 64 multiplies blocks by 3.99 and dense entries by
4.01, matching the expected approximately quadratic output growth.

## Benchmark protocol

Unless a row says otherwise, the measurements below use code commit
`8fdf218`, were taken on 2026-08-14, and use:

- Intel Core i7-8700K, six physical cores / twelve hardware threads;
- 31 GiB RAM and 31 GiB configured swap, Linux 7.0.0-28-generic x86-64;
- Python 3.12.3, NumPy 1.26.4, system OpenBLAS;
- one fresh process per point, with Python import and process startup included;
- one observation per point, not a median or statistical sample;
- `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and
  `BLIS_NUM_THREADS=1`; and
- GNU `/usr/bin/time`, whose `%M` value is peak RSS in KiB. Tables divide it by
  1024 and round to MiB.

The standard workload is `builder.build(depth)`. It includes construction and
all production validation performed by that builder, but excludes table
serialization, braid validation, `--cross-check`, and any prewarming. Filesystem
caches, allocator state, dynamic CPU frequency, and host load were not pinned,
so wall times are representative engineering measurements rather than portable
performance guarantees.

Direct Specht is not a full table builder. Its rows first build the
`dimension-only` reference to enumerate valid keys, then call
`oracle.diagnostics(left, right)` for every block. Its reported wall/RSS thus
includes that reference table and is deliberately reproducible through public
APIs; it is not an isolated kernel-only measurement.

## Current final-code measurements

All five paths at two common cutoffs:

| construction | depth 4 wall / peak RSS | depth 6 wall / peak RSS |
|---|---:|---:|
| `explicit-index` | 5.36 s / 293 MiB | 89.03 s / 619 MiB |
| `fundamental-split` | 3.18 s / 104 MiB | 30.86 s / 419 MiB |
| `recursive-reduction` | 6.35 s / 116 MiB | 54.60 s / 569 MiB |
| `dimension-only` | 0.16 s / 35 MiB | 0.22 s / 35 MiB |
| Direct-Specht diagnostic traversal | 0.33 s / 38 MiB | 0.60 s / 38 MiB |

The complete production curve measured in fresh `dimension-only` processes is:

| depth | wall | peak RSS |
|---:|---:|---:|
| 6 | 0.22 s | 35 MiB |
| 8 | 0.34 s | 36 MiB |
| 16 | 0.94 s | 40 MiB |
| 32 | 3.40 s | 49 MiB |
| 64 | 13.45 s | 73 MiB |

The public-API Direct-Specht traversal at depth 64 checked 38,982 blocks and
131,143 paths, with maximum skew dimension 90, in 56.33 s / 77 MiB. About
13.45 s of that wall time is the included `dimension-only` reference build.

### Explicit decode-cache tradeoff

The current depth-5 explicit builder was also measured with different private
dense-decode budgets. `Retained` is the cache payload at completion; it can be
below the configured cap when the working set is smaller.

| configured cap | wall | peak RSS | retained at completion |
|---:|---:|---:|---:|
| 0 MiB | 24.96 s | 336 MiB | 0 MiB |
| 64 MiB (default) | 23.39 s | 432 MiB | 63 MiB |
| 128 MiB | 21.92 s | 506 MiB | 128 MiB |
| 256 MiB | 22.41 s | 537 MiB | 223 MiB |

This single-run sweep shows a genuine speed/memory tradeoff rather than a hard
optimum. The 64 MiB default recovers some repeated decoding while avoiding the
roughly 74 MiB additional peak of the fastest observed 128 MiB point.

The sweep used the standard environment above and this loop; assigning the
private budget is benchmark instrumentation, not a supported public setting:

```sh
for cap_mib in 0 64 128 256; do
  env \
    PYTHONDONTWRITEBYTECODE=1 \
    OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    BLIS_NUM_THREADS=1 \
    WIGNER_CAP_MIB="$cap_mib" \
    /usr/bin/time -f "cap_mib=$cap_mib wall=%e peak_kib=%M" \
    python3 - <<'PY'
import os
from su3wigner import SwapTableBuilder

builder = SwapTableBuilder()
builder._dense_decode_budget = int(os.environ["WIGNER_CAP_MIB"]) * 1024**2
blocks = builder.build(5)
print(len(blocks), builder._dense_decode_bytes)
PY
done
```

## Historical comparisons

The three carrier baselines below were measured from commit `a638a1c`, before
the packed/streamed changes, with the same fresh-process depth-4 builder
boundary. The final column is the current measurement above.

| construction | baseline wall / peak RSS | final wall / peak RSS |
|---|---:|---:|
| `explicit-index` | 6.00 s / 593 MiB | 5.36 s / 293 MiB |
| `fundamental-split` | 4.24 s / 169 MiB | 3.18 s / 104 MiB |
| `recursive-reduction` | 7.76 s / 256 MiB | 6.35 s / 116 MiB |

The original CG-free `dimension-only` implementation on commit `ee9149f`
took 19.06 s / 330 MiB at depth 64. The final 13.45 s / 73 MiB result is 29%
faster and uses 78% less peak memory. Historical figures cannot remove normal
run-to-run variation, but the large memory reductions are well outside it.

## Correctness reach versus committed regression coverage

These are distinct:

- The committed tests cross-check all four complete builders and Direct Specht
  through cutoff 3, compare `dimension-only` with every checked-in table from
  cutoff 2 through 6, run Direct-Specht raw diagnostics over all 354 cutoff-6
  blocks, and check all prefix sectors through depth 16 for the braid identity
  and both adjacent-swap involutions. See
  [test_generate_table.py](tests/test_generate_table.py),
  [test_dimension_only.py](tests/test_dimension_only.py), and
  [test_direct_specht.py](tests/test_direct_specht.py).
- One-off late precommit optimization worktrees reached depth 8 with
  `fundamental-split` and depth 7 with both `explicit-index` and
  `recursive-reduction`. They used the final coefficient formulae but preceded
  some storage/cache integration in commit `8fdf218`, so they are not claims of
  final-commit runtime coverage. All 473 depth-7 block keys and 1,575 ordered
  paths matched; the largest reported coefficient differences were
  `5.662e-15` (explicit versus split) and `1.998e-15` (recursive versus split).
  At depth 8, `fundamental-split` agreed with `dimension-only` within
  `2.970e-15` over 622 blocks and 2,063 paths. These were manual correctness
  runs, not routine CI benchmarks.
- A one-off `dimension-only` depth-24 run covered 313 prefix labels and 10,359
  `(left,right)` three-adjoint sectors. Building took 1.834 s and braid plus
  both involution checks took 3.864 s; the worst combined residual was
  `3.154e-15`.

Numerical agreement is expected, not byte-identical output across BLAS builds.
The broader validation rationale and binary64 scope are in
[DERIVATION.md](DERIVATION.md#9-numerical-scope).

## Reproducing the measurements

Run from `Wigner6j`. Change `WIGNER_BUILDER` and `WIGNER_DEPTH` for any full
builder row:

```sh
env \
  PYTHONDONTWRITEBYTECODE=1 \
  OPENBLAS_NUM_THREADS=1 \
  OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  BLIS_NUM_THREADS=1 \
  WIGNER_BUILDER=DimensionOnlySwapTableBuilder \
  WIGNER_DEPTH=64 \
  /usr/bin/time -f 'wall=%e peak_kib=%M' \
  python3 - <<'PY'
import os
import su3wigner

builder_class = getattr(su3wigner, os.environ["WIGNER_BUILDER"])
blocks = builder_class().build(int(os.environ["WIGNER_DEPTH"]))
print(
    len(blocks),
    sum(block.size for block in blocks),
    sum(block.matrix.size for block in blocks),
)
PY
```

Full builder names are `SwapTableBuilder`,
`FundamentalSplitSwapTableBuilder`, `RecursiveReductionSwapTableBuilder`,
`DimensionOnlySwapTableBuilder`, and `DirectSpechtSwapTableBuilder`. Do not add
`--cross-check` to a timing run: it
builds multiple algorithms and intentionally shares carrier CG input.

The Direct-Specht full-builder traversal is:

```sh
env \
  PYTHONDONTWRITEBYTECODE=1 \
  OPENBLAS_NUM_THREADS=1 \
  OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  BLIS_NUM_THREADS=1 \
  /usr/bin/time -f 'wall=%e peak_kib=%M' \
  python3 - <<'PY'
from su3wigner import DirectSpechtSwapTableBuilder

depth = 64
builder = DirectSpechtSwapTableBuilder()
blocks = builder.build(depth)
print(
    len(blocks),
    sum(block.size for block in blocks),
    sum(block.matrix.size for block in blocks),
)
PY
```

For an end-to-end file-generation measurement, include serialization and the
desired braid depth explicitly:

```sh
env \
  PYTHONDONTWRITEBYTECODE=1 \
  OPENBLAS_NUM_THREADS=1 \
  OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 \
  BLIS_NUM_THREADS=1 \
  /usr/bin/time -f 'wall=%e peak_kib=%M' \
  python3 generate_table.py \
    --method dimension-only \
    --max-prefix-gluons 16 \
    --braid-check-depth 16 \
    --output /tmp/su3_adjoint_swap_prefix_16.tbl
```

For correctness rather than timing, run the complete regression suite or the
CLI cross-check described in [README.md](README.md#the-five-constructions).
