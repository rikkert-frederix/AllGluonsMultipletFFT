# RamboOnDiet event generator

This directory contains the RAMBO-on-diet phase-space mapping and a small
driver that writes random massless centre-of-mass events in the shared
`AMPLIGLUON_EVENT_V1` format. The same physical momentum and helicity
convention is accepted by AmpliGluonMultiplet, AmpliGluonTrace, and
AmpliGluonAdjoint: the first two momenta are incoming, all remaining momenta
are outgoing, and every energy is positive.

## Build and test

From the repository root:

```sh
make -C RamboOnDiet
make -C RamboOnDiet test
make -C RamboOnDiet debug
```

The normal build creates `RamboOnDiet/build/generate_ampligluon_events`. The
test checks seeded reproducibility and compares every generated helicity row
with AmpliGluonMultiplet and AmpliGluonTrace.

## Generate events

```text
generate_ampligluon_events FINAL_GLUONS POINTS OUTPUT_PREFIX [SQRT_S [SEED]]
```

For example, this generates five independent `g g -> 4 g` points at 1 TeV:

```sh
mkdir -p random
RamboOnDiet/build/generate_ampligluon_events \
  4 5 random/gg_to_4g 1000 12345
```

The parent directory of the output prefix must already exist. The files are
named `random/gg_to_4g_000001.event` through
`random/gg_to_4g_000005.event`. Existing files with those names are replaced.
`SQRT_S` defaults to `1000` and `SEED` defaults to `1`. Supplying the same seed
to the same executable reproduces the same sequence.

For `N = FINAL_GLUONS + 2` total gluons, each file contains one phase-space
point and all `2^N` physical helicity configurations. Every leg independently
takes helicity `-1` or `+1`; rows use canonical binary order, starting with all
minus and ending with all plus. The first two entries are the physical incoming
helicities and the remaining entries are outgoing, exactly as for the
momenta. Consequently, a consumer must apply its usual incoming-leg crossing
convention when classifying analytic-zero, MHV, and other sectors.

Exhaustive enumeration currently supports at most 20 total gluons, and event
size doubles with every additional leg. The phase-space weight is printed by
the generator but is not written into the event because a matrix-element
comparison does not use it. The strong coupling in each file is `1.0`.

Run a generated file through AmpliGluonMultiplet with:

```sh
AmpliGluonMultiplet/build/ampligluon_multiplet \
  Wigner6j/data/su3_adjoint_swap_prefix_6.tbl \
  random/gg_to_4g_000001.event
```

Its momentum and helicity rows can be passed unchanged to
`ampligluon_trace_t%evaluate`, making the file a common input for
point-by-point AmpliGluonMultiplet--AmpliGluonTrace comparisons.

For convenience, build and run the comparison driver directly:

```sh
make -C RamboOnDiet checker
RamboOnDiet/build/compare_ampligluon_multiplet_trace \
  Wigner6j/data/su3_adjoint_swap_prefix_6.tbl \
  random/gg_to_4g_000001.event
```

The checker evaluates every helicity row with both implementations, prints the
two colour sums and their relative difference, and exits unsuccessfully if a
nonzero result differs by more than `1e-10` relative. Results below `1e-24`
are treated as numerical zero and compared with an absolute tolerance instead.
It disables the MHV/anti-MHV optimization in both implementations so that
every nonzero row is evaluated with their general Berends-Giele recursions.

The generator and comparison checker remain exhaustive. The top-level
`Benchmark` harness, however, no longer evaluates every row with every
full-colour backend. It evaluates one cheap canonical trace order for every
nonzero row, assigns analytic zeros directly, and uses the squared partial
amplitude as a leading-colour proxy weight.
Across all generated phase-space points it then selects deterministic,
proxy-weighted timing samples: by default three MHV/anti-MHV configurations
and one general-helicity configuration when that class is nonempty. One
analytically vanishing row is also carried as an unweighted
numerical-validation sentinel.

The harness exposes six selectable backend entries: optimized-MHV and
default-BG variants for AmpliGluonMultiplet, AmpliGluonTrace, and
AmpliGluonAdjoint. Each optimized variant times its analytic MHV/anti-MHV path
and its general Berends-Giele path separately. If `f_MHV` is the exact MHV
fraction of the exhaustively evaluated proxy distribution, its predicted
single-helicity time is

```text
T_predicted = f_MHV * T_MHV + (1 - f_MHV) * T_BG.
```

Each default-BG variant sends every nonzero configuration through the same
control flow, so one representative proxy-weighted nonzero row is sufficient
for its timing prediction. The sentinel is checked in the first validation
pass but excluded from calibrated timing cells and weighting. The proxy is used
for sampling and for `f_MHV`; it is not reported as a
full-colour helicity sum. Full-colour numerical agreement is checked only on
the bounded sample. The benchmark defaults to ten phase-space points, three
MHV samples, one non-MHV sample, and a 600-second timeout per backend
invocation. Use `--mhv-samples` and `--non-mhv-samples` to change the timing
sample counts without changing this generator's exhaustive event format.
