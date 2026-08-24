# Lightweight MadGraph5_aMC@NLO all-gluon standalones

This directory vendors the generated Fortran matrix-element kernels for

- `g g > g g` (`N=4`),
- `g g > g g g` (`N=5`),
- `g g > g g g g` (`N=6`).

They were generated with MadGraph5_aMC@NLO 3.6.0 and the Standard Model using
`QED=0`.  The original process cards, reported generator version, source Git
commit and Standard Model tree, per-process metadata, and SHA-256 hashes of the
generated `matrix.f` files are retained.
Only the generated matrix kernels, their required pure-gluon HELAS routines,
and benchmark inputs are stored; parameter cards, HTML, diagrams, binaries,
object files, and unused Standard Model code are omitted.  The result is about
0.8 MiB.  `N=7` and above are deliberately not included: the generated `N=7`
matrix source alone is about 18 MiB, and compiling it with the common `-O3`
benchmark flags proved too costly in both time and memory for a convenient
reference benchmark.

## Rebuild and benchmark

From the repository root, run

```sh
make -C Benchmark/MadGraph5
```

This needs only Python 3 and a current `gfortran`; it does not invoke MadGraph.
The result is written to `Benchmark/build/results-madgraph-fixed.md`.  Generated
objects and executables stay in the ignored `Benchmark/build/madgraph` tree.
For a quick compile-and-run check through `N=6`, use

```sh
make -C Benchmark/MadGraph5 smoke
```

Extra options can be passed through `OPTIONS`, for example

```sh
make -C Benchmark/MadGraph5 OPTIONS='--target-seconds 0.5 --batches 5'
```

The three checked-in event files are the exact deterministic, proxy-weighted
sample-1 rows from the three-point native N=4--11 benchmark.  They allow the
comparison to be reproduced without regenerating RAMBO points or re-running
the sampling stage.  Their native DDM values and hashes are retained in
`events/reference.json`, so the self-contained Make target also detects a
normalization or coupling regression.  To time a different sample set, pass
`--events-dir`; the
runner accepts either this directory's `N4.event` layout or the main harness's
`N4/sampled/gg_to_2g_sample_000001.event` layout.

The optional numerical cross-check against the native DDM result can be enabled
after building the main benchmark:

```sh
python3 Benchmark/run_madgraph_benchmark.py
```

## Single-helicity guarantee

The timing driver calls MadGraph's generated

```fortran
MATRIX(P, NHEL, IC)
```

function directly, passing one explicit physical `-1`/`+1` helicity vector.
It never calls the helicity-summing `SMATRIX` routine or the `SMATRIXHEL`
wrapper, which still scans the helicity-index array before evaluating one row.
For the checked-in kernels, the extracted `MATRIX` body contains no helicity
loop.  The runner also rejects the standard MadGraph `IHEL`/`NCOMB` loop form
if it occurs inside that body.

Incoming momenta retain positive energies and physical helicity labels are
passed unchanged; HELAS performs the crossing internally.  The direct kernel
returns the full colour sum with no incoming-colour or incoming-helicity
average and no identical-final-state symmetry factor.  This matches the native
benchmark convention.  The lightweight coupling adapter sets `g_s` from the
event file once, outside the timed region.

## Regenerating the vendored sources

Regeneration is normally unnecessary.  Given an MG5_aMC executable, the exact
layout can be recreated with

```sh
python3 Benchmark/vendor_madgraph_processes.py \
  --mg5 /absolute/path/to/bin/mg5_aMC \
  --min-gluons 4 --max-gluons 6
```

The script generates each standalone in a temporary directory, verifies that
there is exactly one subprocess, retains the provenance metadata, and checks
that shared HELAS sources are identical across multiplicities.
