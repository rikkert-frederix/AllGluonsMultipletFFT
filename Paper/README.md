# Paper build and numerical inputs

Build the manuscript from this directory with

```sh
make
```

The main source is [`gluon_colour_sums.tex`](gluon_colour_sums.tex), and the
build produces [`gluon_colour_sums.pdf`](gluon_colour_sums.pdf).

The numerical source for the manuscript is [`results.md`](results.md).  The
reported CPU time is the time for a matrix-element call after initialization
and completion of the first call.  The paper uses the general
Berends--Giele calculation for all non-vanishing helicities; it does not use
the empty special-helicity timing table or the approximate mixed-helicity
time in the report.

The values printed in the manuscript are in
[`tables/evaluation_results.tex`](tables/evaluation_results.tex) and
[`tables/memory_results.tex`](tables/memory_results.tex), and the plot input
is [`data/colour_sum_times.dat`](data/colour_sum_times.dat).  The
three direct-colour times missing from `results.md` are marked by an asterisk.
They use

```text
t(N2) = t(N1) [P(N2) / P(N1)]^2,
P_trace(N) = (N-1)!,
P_adjoint(N) = (N-2)!.
```

This gives 2405.7 s and 240570 s for direct trace at ten and eleven gluons,
and 1692.9 s for direct adjoint at eleven gluons.  Memory values are never
extrapolated.

The bibliography uses InspireHEP texkeys for every reference indexed there.
The three mathematical references for fast symmetric-group transforms are
not indexed by InspireHEP and therefore use descriptive BibTeX keys.
