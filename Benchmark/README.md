# Proxy-weighted single-helicity all-gluon benchmark

This directory contains the timing drivers and orchestration code for comparing

- `AmpliGluonMultipletOptimizedMHV`, using normalized SU(3) multiplets and the
  optimized analytic MHV path;
- `AmpliGluonMultipletDefaultBG`, using the same multiplet basis but the default
  Berends-Giele reduction for every nonzero helicity configuration;
- `AmpliGluonTraceOptimizedMHV` and `AmpliGluonTraceDefaultBG`, using
  `(N-1)!` fundamental trace orders, respectively the analytic MHV or default
  Berends-Giele kinematic path, and the symmetric-group FFT colour
  contraction;
- `AmpliGluonTraceOptimizedMHVDirectColour` and
  `AmpliGluonTraceDefaultBGDirectColour`, using the same trace kinematic paths
  with the normal direct colour-matrix multiplication;
- `AmpliGluonAdjointOptimizedMHV` and `AmpliGluonAdjointDefaultBG`, using the
  `(N-2)!` DDM adjoint basis, respectively the analytic MHV or default
  Berends-Giele kinematic path, and the symmetric-group FFT colour
  contraction;
- `AmpliGluonAdjointOptimizedMHVDirectColour` and
  `AmpliGluonAdjointDefaultBGDirectColour`, using the same adjoint kinematic
  paths with the normal direct colour-matrix multiplication.

The harness estimates the cost of one matrix-element evaluation in a realistic
helicity-assigned workload. Such events are expected to be distributed like
the underlying helicity-summed leading-colour result. Evaluating every
full-colour backend for all `2^N` helicities would reproduce that distribution
directly, but becomes prohibitively expensive at ten and eleven gluons.

Instead, the event generator still writes every physical helicity row and a
cheap proxy produces a weight for every row. It evaluates one canonical trace
order for each nonzero row and assigns analytic zeros directly. The squared
canonical partial amplitude is used only as a leading-colour sampling weight.
The benchmark then selects deterministic, globally proxy-weighted samples from
the complete set of phase-space-point/helicity pairs. By default it selects
three MHV/anti-MHV rows, one general-helicity row when that class is nonempty,
and one analytically vanishing row as a numerical-validation sentinel.
Optimized backends receive the full bounded pool. Default-BG backends receive
only its representative nonzero row and the sentinel, because additional
helicities would repeat the same nonzero control flow.

This preserves the important MHV-versus-general-helicity mixture without
turning the benchmark into a full-colour helicity sweep. The proxy is an
approximation to the desired event distribution, not a reported full-colour
matrix element or helicity sum.

For high multiplicities, implementations can diverge in feasibility. The
script supports backend selection (`--backend`) and per-multiplicity preflight
checks, so you can run only the useful overlap instead of forcing every
implementation at every multiplicity.

## Run

From the repository root:

```sh
make -C Benchmark
```

The report is written to `Benchmark/build/results.md` and also printed to
standard output. The default run covers four through eight total gluons, ten
phase-space points per multiplicity, and requests three MHV samples and one
non-MHV sample from each nonempty class. The exhaustive canonical-order proxy
remains inexpensive enough to determine the sampling weights and exact proxy
MHV fraction at each multiplicity. The default per-backend timeout is 600
seconds.

For a quick smoke test:

```sh
make -C Benchmark test
```

### MadGraph5_aMC@NLO fixed-helicity reference

The low-multiplicity MadGraph5_aMC@NLO comparison is kept separate from the
ten native backends because it calls a generated diagram-based kernel rather
than an `OptimizedMHV` or `DefaultBG` path.  Vendored pure-gluon standalones and
the exact proxy-weighted sample-1 events are provided for `N=4` through `N=6`.
From the repository root, rebuild and time them with

```sh
make -C Benchmark/MadGraph5
```

This writes `Benchmark/build/results-madgraph-fixed.md`; it uses the direct
`MATRIX(P,NHEL,IC)` entry point for one explicit helicity vector and does not
call the helicity-summing `SMATRIX` routine or the index-scanning `SMATRIXHEL`
wrapper.  The default Make target is self-contained and skips the optional
numerical comparison with the native DDM backend.  A short-timing smoke check is

```sh
make -C Benchmark/MadGraph5 smoke
```

See [`MadGraph5/README.md`](MadGraph5/README.md) for the precise timing and
normalization conventions, optional DDM cross-check, source provenance, and
regeneration command.  No MadGraph installation is needed unless the vendored
sources themselves are regenerated.

Run the Python entry point directly to select the sample size and timing
precision:

```sh
python3 Benchmark/run_benchmark.py \
  --points 12 \
  --sqrt-s 1000 \
  --seed 12345 \
  --mhv-samples 3 \
  --non-mhv-samples 1 \
  --target-seconds 0.5 \
  --batches 5 \
  --initialization-runs 3 \
  --backend AmpliGluonTraceOptimizedMHV \
  --backend AmpliGluonTraceOptimizedMHVDirectColour \
  --max-memory-gib 8 \
  --backend-timeout 600 \
  --output timings.md
```

Use `--help` for the complete option list. In particular, `--fflags` applies
the same flags to all three matrix-element implementations, the proxy, and the
RAMBO driver; the default deliberately does not enable `-ffast-math`.

For maximum single-core performance on the current machine, use the explicit
native target:

```sh
make -C Benchmark benchmark-native
```

This builds every selected backend, the proxy, and the event generator with
`-Ofast -march=native`. The resulting executables are host-specific and may
reassociate floating-point expressions, so this mode remains opt-in; the
benchmark still checks sampled numerical agreement and records the exact
flags.

Important runtime controls:

- `--mhv-samples` sets the number of deterministic proxy-weighted
  MHV/anti-MHV timing samples; `--non-mhv-samples` does the same for the
  general-helicity sector. Their defaults are three and one, respectively.
- `--points` controls the phase-space ensemble over which the proxy weights
  and samples are constructed. The seed makes both the generated ensemble and
  the resulting sample selection reproducible.
- `--backend NAME` can be repeated. The multiplet choices are
  `AmpliGluonMultipletOptimizedMHV` and
  `AmpliGluonMultipletDefaultBG`. Trace and Adjoint each provide optimized-MHV
  and default-BG variants with the FFT contraction under the existing names;
  append `DirectColour` to select the corresponding normal colour-matrix
  multiplication. Omit both multiplet variants when a higher multiplicity
  exceeds Wigner support or when you only want the trace/adjoint paths. The
  direct Trace cost grows especially quickly because its colour dimension is
  `(N-1)!`; the default range through eight gluons keeps that comparison
  practical.
- `--max-memory-gib` is enforced as a per-process address-space limit and is
  also used by the adjoint preflight estimate; `--max-target-seconds` bounds
  timing calibration work.
- `--backend-timeout` aborts slow invocations per backend without aborting the
  full run.
- `--skip-initialization-preflight` avoids a separate initialization-only
  process when a high-multiplicity point must be attempted exactly once.
- Results are written incrementally while sweeping multiplicities.

The current one-sided multiplet consumer only needs representation labels
reachable within half of each closed colour chain. Consequently, a table with
`MAX_PREFIX_GLUONS = p` supports total multiplicities through `2p+1`: the
required depth is `floor(total_gluons/2)`. For example, the checked-in prefix-6
table covers the group-theory data through thirteen total gluons. Runtime and
memory requirements of an individual backend can impose substantially lower
practical limits.

## Timing estimate

Let `w_p,h` be the squared canonical-order proxy amplitude at phase-space point
`p` and helicity `h`. Analytic-zero rows have weight zero. The exhaustive but
cheap proxy calculation determines

```text
f_MHV = sum_(p,h in MHV) w_p,h / sum_(p,h nonzero) w_p,h.
```

The sums are global over all generated phase-space points and physical
helicities. The deterministic samples are drawn from the same globally
weighted population, separately within the MHV and non-MHV classes. At four
and five gluons all nonzero configurations are MHV/anti-MHV, so there is no
general-helicity sample or `T_BG` term for an optimized backend.

For an `OptimizedMHV` backend, the harness measures the analytic MHV path and
the general Berends-Giele path separately. In each calibrated batch, if their
sample-mean timings are `T_MHV` and `T_BG`, the combined prediction is

```text
T_predicted = f_MHV * T_MHV + (1 - f_MHV) * T_BG.
```

The headline entry is the median of these batch-by-batch predictions; the
separate path tables report the median batch timing for each path.

For a `DefaultBG` backend every nonzero helicity follows the same control flow.
Its timing therefore needs only one representative, proxy-weighted nonzero
configuration; the MHV fraction does not alter its prediction. Phase-space
dependence is still represented by choosing the row from the global weighted
population. The analytic-zero sentinel participates in the first validation
pass but is excluded from the calibrated timing cells and production estimate.

## Report contents

The generated report contains:

1. the colour-space dimension reported by every selected backend;
2. median fresh-process initialization timings;
3. measured peak resident memory of each fresh-process run;
4. the exact canonical-proxy MHV fraction and the deterministic sample rows;
5. sampled MHV-path and BG-path timings and the combined single-helicity
   prediction;
6. numerical agreement on the sampled nonzero configurations and the absolute
   result for the analytic-zero sentinel;
7. compiler, flags, seeds, and SHA-256 hashes of the event files and Wigner
   table.

The numerical checks deliberately cover the timing samples rather than every
full-colour helicity row. The report consequently does not claim a full-colour
helicity sum or a measured complete-helicity sweep. Initialization is kept
separate, and backend execution order is rotated with multiplicity and
initialization sample to reduce ordering bias.

RamboOnDiet writes every physical `-1`/`+1` assignment, from all minus through
all plus in canonical binary order. Configurations that vanish analytically
after crossing the two incoming legs are identified from the helicity row and
checked with an absolute threshold; sampled nonzero matrix elements are
checked relatively. For Multiplet, Trace, and Adjoint, the `OptimizedMHV`
variant uses its analytic MHV/anti-MHV path, while the `DefaultBG` variant
routes every nonzero sector through the general Berends-Giele recursion. Exact
zero sectors return before recursion in either mode. The Trace and Adjoint
variants additionally exercise both the FFT and direct colour contractions on
the same selected amplitudes, and the numerical table checks their agreement.

Fortran's intrinsic random-number stream is reproducible for the same compiler
runtime but is not standardized bit-for-bit across compilers. The report
therefore records hashes of the actual event files in addition to their seeds.

## Layout

- `run_benchmark.py` builds, generates, samples, runs, validates, and renders
  the report;
- `src/benchmark_events.f90` provides the shared `AMPLIGLUON_EVENT_V1` I/O;
- `src/benchmark_helicity_proxy.f90` exhaustively evaluates the cheap
  canonical-order leading-colour proxy;
- the three `src/benchmark_ampligluon_*.f90` sources are the isolated
  full-colour timing drivers.

All generated executables, module files, phase-space points, sampled event
files, and reports live under the ignored `Benchmark/build` directory.
