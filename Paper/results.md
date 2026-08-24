# Proxy-weighted sampled-helicity all-gluon matrix-element benchmark

RAMBO-on-diet phase-space points at `sqrt(s) = 1000`, 10 point(s) per multiplicity. A canonical leading-colour proxy weight is produced for all `2^N` physical helicities (with analytic zeros assigned directly), while the expensive full-colour backends use a bounded proxy-weighted sample pool. Optimized backends use every selected nonzero row; default-BG backends use its representative row and the common zero sentinel.

Selected backends: `Multiplet (default BG)`, `Trace (default BG, FFT colour)`, `Adjoint (default BG, FFT colour)`, `Trace (default BG, direct colour)`, `Adjoint (default BG, direct colour)`.

## Problem dimensions

| Total gluons | Multiplet (default BG) | Trace (default BG, FFT colour) | Adjoint (default BG, FFT colour) | Trace (default BG, direct colour) | Adjoint (default BG, direct colour) |
|---:|---:|---:|---:|---:|---:|
| 4 | 8 | 6 | 2 | 6 | 2 |
| 5 | 32 | 24 | 6 | 24 | 6 |
| 6 | 145 | 120 | 24 | 120 | 24 |
| 7 | 702 | 720 | 120 | 720 | 120 |
| 8 | 3598 | 5040 | 720 | 5040 | 720 |
| 9 | 19280 | 40320 | 5040 | 40320 | 5040 |
| 10 | 107160 | 362880 | 40320 | N/A | 40320 |
| 11 | 614000 | 3628800 | 362880 | N/A | N/A |

## Peak resident memory

Maximum resident set size of the full fresh-process run, including initialization and all requested evaluations.

| Total gluons | Multiplet (default BG) | Trace (default BG, FFT colour) | Adjoint (default BG, FFT colour) | Trace (default BG, direct colour) | Adjoint (default BG, direct colour) |
|---:|---:|---:|---:|---:|---:|
| 4 | 2.89 MiB | 2.9 MiB | 2.89 MiB | 2.86 MiB | 2.89 MiB |
| 5 | 3.02 MiB | 2.9 MiB | 2.84 MiB | 2.92 MiB | 2.81 MiB |
| 6 | 3.16 MiB | 2.98 MiB | 2.85 MiB | 2.86 MiB | 2.88 MiB |
| 7 | 3.98 MiB | 3.16 MiB | 3.09 MiB | 2.96 MiB | 2.98 MiB |
| 8 | 10.1 MiB | 4.43 MiB | 3.54 MiB | 4.24 MiB | 3.34 MiB |
| 9 | 55.4 MiB | 14.8 MiB | 7.86 MiB | 14 MiB | 7.66 MiB |
| 10 | 425 MiB | 111 MiB | 44.6 MiB | N/A | 43.9 MiB |
| 11 | 3.96 GiB | 1.09 GiB | 398 MiB | N/A | N/A |

## Initialization

Median of 1 fresh-process runs.

| Total gluons | Multiplet (default BG) | Trace (default BG, FFT colour) | Adjoint (default BG, FFT colour) | Trace (default BG, direct colour) | Adjoint (default BG, direct colour) | Fastest vs next |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 12.1 ms | 3 µs | 15 µs | 3 µs | 12 µs | tie |
| 5 | 12.1 ms | 4 µs | 26 µs | 4 µs | 17 µs | tie |
| 6 | 12.6 ms | 3 µs | 63 µs | 4 µs | 41 µs | Trace (default BG, FFT colour) 1.33× |
| 7 | 14.6 ms | 4 µs | 296 µs | 4 µs | 266 µs | tie |
| 8 | 27 ms | 3 µs | 2.4 ms | 4 µs | 2.17 ms | Trace (default BG, FFT colour) 1.33× |
| 9 | 108 ms | 4 µs | 23.4 ms | 4 µs | 20.2 ms | tie |
| 10 | 700 ms | 3 µs | 240 ms | N/A | 226 ms | Trace (default BG, FFT colour) 75418.33× |
| 11 | 5.28 s | 4 µs | 2.78 s | N/A | N/A | Trace (default BG, FFT colour) 695276.00× |

## First sampled pass

CPU time for one pass over the selected validation rows. This includes any lazy MHV or BG setup triggered by the sample.

| Total gluons | Multiplet (default BG) | Trace (default BG, FFT colour) | Adjoint (default BG, FFT colour) | Trace (default BG, direct colour) | Adjoint (default BG, direct colour) |
|---:|---:|---:|---:|---:|---:|
| 4 | 16 µs | 39 µs | 23 µs | 24 µs | 29 µs |
| 5 | 40 µs | 69 µs | 45 µs | 48 µs | 56 µs |
| 6 | 202 µs | 206 µs | 154 µs | 265 µs | 153 µs |
| 7 | 1.42 ms | 1.17 ms | 895 µs | 7.06 ms | 991 µs |
| 8 | 13.5 ms | 9.24 ms | 6.4 ms | 394 ms | 11.9 ms |
| 9 | 150 ms | 83.7 ms | 55.7 ms | 29.8 s | 349 ms |
| 10 | 1.94 s | 943 ms | 541 ms | N/A | 21.3 s |
| 11 | 24.7 s | 12 s | 6.47 s | N/A | N/A |

## Estimated production-workload evaluation after initialization

Median of 10 calibrated batches. In each batch, an optimized-backend estimate is `f_MHV t_MHV + (1-f_MHV) t_BG`; the table reports the median of those mixtures. Here `f_MHV` is the global canonical leading-colour proxy fraction. A default-BG backend uses its representative nonzero BG timing directly.

| Total gluons | Multiplet (default BG) | Trace (default BG, FFT colour) | Adjoint (default BG, FFT colour) | Trace (default BG, direct colour) | Adjoint (default BG, direct colour) | Fastest vs next |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 923 ns | 968 ns | 645 ns | 729 ns | 543 ns | Adjoint (default BG, direct colour) 1.19× |
| 5 | 6.38 µs | 4.41 µs | 2.62 µs | 6.23 µs | 2.44 µs | Adjoint (default BG, direct colour) 1.07× |
| 6 | 62.3 µs | 26.2 µs | 13.3 µs | 138 µs | 15.3 µs | Adjoint (default BG, FFT colour) 1.15× |
| 7 | 718 µs | 186 µs | 83.2 µs | 6.16 ms | 196 µs | Adjoint (default BG, FFT colour) 2.23× |
| 8 | 8.26 ms | 1.62 ms | 620 µs | 381 ms | 5.78 ms | Adjoint (default BG, FFT colour) 2.62× |
| 9 | 104 ms | 15.5 ms | 5.23 ms | 29.7 s | 303 ms | Adjoint (default BG, FFT colour) 2.95× |
| 10 | 1.49 s | 215 ms | 58.9 ms | N/A | 20.9 s | Adjoint (default BG, FFT colour) 3.65× |
| 11 | 19.8 s | 3.08 s | 814 ms | N/A | N/A | Adjoint (default BG, FFT colour) 3.79× |

## MHV/anti-MHV evaluation after initialization

Median of 10 batches over 3 deterministic proxy-weighted sample(s).

| Total gluons | Multiplet (default BG) | Trace (default BG, FFT colour) | Adjoint (default BG, FFT colour) | Trace (default BG, direct colour) | Adjoint (default BG, direct colour) | Fastest vs next |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | N/A | N/A | N/A | N/A | N/A | N/A |
| 5 | N/A | N/A | N/A | N/A | N/A | N/A |
| 6 | N/A | N/A | N/A | N/A | N/A | N/A |
| 7 | N/A | N/A | N/A | N/A | N/A | N/A |
| 8 | N/A | N/A | N/A | N/A | N/A | N/A |
| 9 | N/A | N/A | N/A | N/A | N/A | N/A |
| 10 | N/A | N/A | N/A | N/A | N/A | N/A |
| 11 | N/A | N/A | N/A | N/A | N/A | N/A |

## Berends-Giele evaluation after initialization

Median of 10 batches. Optimized backends use the selected general-helicity samples. Default-BG backends use one representative nonzero row because their operation count does not depend on the nonzero helicity; at N=4 and N=5 that representative is necessarily MHV.

| Total gluons | Multiplet (default BG) | Trace (default BG, FFT colour) | Adjoint (default BG, FFT colour) | Trace (default BG, direct colour) | Adjoint (default BG, direct colour) | Fastest vs next |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 923 ns | 968 ns | 645 ns | 729 ns | 543 ns | Adjoint (default BG, direct colour) 1.19× |
| 5 | 6.38 µs | 4.41 µs | 2.62 µs | 6.23 µs | 2.44 µs | Adjoint (default BG, direct colour) 1.07× |
| 6 | 62.3 µs | 26.2 µs | 13.3 µs | 138 µs | 15.3 µs | Adjoint (default BG, FFT colour) 1.15× |
| 7 | 718 µs | 186 µs | 83.2 µs | 6.16 ms | 196 µs | Adjoint (default BG, FFT colour) 2.23× |
| 8 | 8.26 ms | 1.62 ms | 620 µs | 381 ms | 5.78 ms | Adjoint (default BG, FFT colour) 2.62× |
| 9 | 104 ms | 15.5 ms | 5.23 ms | 29.7 s | 303 ms | Adjoint (default BG, FFT colour) 2.95× |
| 10 | 1.49 s | 215 ms | 58.9 ms | N/A | 20.9 s | Adjoint (default BG, FFT colour) 3.65× |
| 11 | 19.8 s | 3.08 s | 814 ms | N/A | N/A | Adjoint (default BG, FFT colour) 3.79× |

## Numerical agreement

| Total gluons | Compared backends | Max relative spread (nonzero) | Max absolute value (zero sector) |
|---:|:---|---:|---:|
| 4 | Multiplet (default BG), Trace (default BG, FFT colour), Adjoint (default BG, FFT colour), Trace (default BG, direct colour), Adjoint (default BG, direct colour) | 1.212e-15 | 0.000e+00 |
| 5 | Multiplet (default BG), Trace (default BG, FFT colour), Adjoint (default BG, FFT colour), Trace (default BG, direct colour), Adjoint (default BG, direct colour) | 3.382e-15 | 0.000e+00 |
| 6 | Multiplet (default BG), Trace (default BG, FFT colour), Adjoint (default BG, FFT colour), Trace (default BG, direct colour), Adjoint (default BG, direct colour) | 3.023e-15 | 0.000e+00 |
| 7 | Multiplet (default BG), Trace (default BG, FFT colour), Adjoint (default BG, FFT colour), Trace (default BG, direct colour), Adjoint (default BG, direct colour) | 3.558e-14 | 0.000e+00 |
| 8 | Multiplet (default BG), Trace (default BG, FFT colour), Adjoint (default BG, FFT colour), Trace (default BG, direct colour), Adjoint (default BG, direct colour) | 3.146e-13 | 0.000e+00 |
| 9 | Multiplet (default BG), Trace (default BG, FFT colour), Adjoint (default BG, FFT colour), Trace (default BG, direct colour), Adjoint (default BG, direct colour) | 5.146e-13 | 0.000e+00 |
| 10 | Multiplet (default BG), Trace (default BG, FFT colour), Adjoint (default BG, FFT colour), Adjoint (default BG, direct colour) | 8.054e-12 | 0.000e+00 |
| 11 | Multiplet (default BG), Trace (default BG, FFT colour), Adjoint (default BG, FFT colour) | 2.807e-14 | 0.000e+00 |

## Leading-colour helicity proxy

The fraction is computed from the exhaustive canonical-order proxy before sampling: `sum_p,h in MHV |A_p,h|^2 / sum_p,h |A_p,h|^2`. The global ratio represents points drawn according to the proxy helicity-summed result. Because this is one colour order rather than the full leading-colour sum, the combined timing is explicitly an approximation.

| Total gluons | Global MHV fraction | Point-to-point range |
|---:|---:|---:|
| 4 | 1.000000 | 1.000000 – 1.000000 |
| 5 | 1.000000 | 1.000000 – 1.000000 |
| 6 | 0.619140 | 0.596784 – 0.915393 |
| 7 | 0.583352 | 0.458650 – 0.754855 |
| 8 | 0.312354 | 0.310599 – 0.661647 |
| 9 | 0.298519 | 0.251972 – 0.634733 |
| 10 | 0.264340 | 0.172767 – 0.632302 |
| 11 | 0.168751 | 0.153397 – 0.377835 |

## Selected helicity samples

Original event/configuration indices, crossed-helicity path, and canonical leading-colour proxy weight in the selected pool. Sample 1 is drawn from the complete nonzero proxy distribution. Optimized backends use the full pool; default-BG backends use samples 1 and 2. Repeated draws are retained.

| Gluons | Sample | Source point | Source helicity | Path | Proxy weight |
|---:|---:|---:|---:|:---|---:|
| 4 | 1 | 6 | 16 | mhv | 1.460178244421e+03 |
| 4 | 2 | 4 | 15 | zero | 0.000000000000e+00 |
| 4 | 3 | 6 | 16 | mhv | 1.460178244421e+03 |
| 4 | 4 | 10 | 16 | mhv | 5.991301881266e+02 |
| 5 | 1 | 9 | 18 | mhv | 8.715940052174e-03 |
| 5 | 2 | 1 | 17 | zero | 0.000000000000e+00 |
| 5 | 3 | 9 | 1 | mhv | 1.287681490604e-02 |
| 5 | 4 | 9 | 18 | mhv | 8.715940052174e-03 |
| 6 | 1 | 9 | 41 | mhv | 2.854344792951e-08 |
| 6 | 2 | 2 | 8 | zero | 0.000000000000e+00 |
| 6 | 3 | 7 | 60 | bg | 1.479159823135e-06 |
| 6 | 4 | 7 | 1 | mhv | 2.552633869632e-06 |
| 6 | 5 | 7 | 1 | mhv | 2.552633869632e-06 |
| 7 | 1 | 5 | 3 | bg | 6.402394639729e-12 |
| 7 | 2 | 8 | 33 | zero | 0.000000000000e+00 |
| 7 | 3 | 2 | 1 | mhv | 6.946710728683e-12 |
| 7 | 4 | 5 | 1 | mhv | 3.623847345813e-11 |
| 7 | 5 | 1 | 1 | mhv | 8.684865727226e-11 |
| 8 | 1 | 10 | 7 | bg | 6.874584731924e-14 |
| 8 | 2 | 7 | 225 | zero | 0.000000000000e+00 |
| 8 | 3 | 10 | 256 | mhv | 1.233035081552e-12 |
| 8 | 4 | 10 | 256 | mhv | 1.233035081552e-12 |
| 8 | 5 | 10 | 256 | mhv | 1.233035081552e-12 |
| 9 | 1 | 5 | 448 | bg | 2.133125197535e-18 |
| 9 | 2 | 4 | 389 | zero | 0.000000000000e+00 |
| 9 | 3 | 6 | 1 | mhv | 6.317682487434e-19 |
| 9 | 4 | 5 | 512 | mhv | 6.579967788823e-18 |
| 9 | 5 | 5 | 1 | mhv | 6.579967788823e-18 |
| 10 | 1 | 10 | 1000 | bg | 2.893053260716e-22 |
| 10 | 2 | 3 | 512 | zero | 0.000000000000e+00 |
| 10 | 3 | 5 | 1024 | mhv | 1.108548336498e-20 |
| 10 | 4 | 5 | 1024 | mhv | 1.108548336498e-20 |
| 10 | 5 | 5 | 1024 | mhv | 1.108548336498e-20 |
| 11 | 1 | 7 | 41 | bg | 8.692339498354e-24 |
| 11 | 2 | 7 | 496 | zero | 0.000000000000e+00 |
| 11 | 3 | 4 | 1 | mhv | 2.774769123523e-22 |
| 11 | 4 | 4 | 1 | mhv | 2.774769123523e-22 |
| 11 | 5 | 4 | 2048 | mhv | 2.774769123523e-22 |

## Skipped backends

| Total gluons | Backend | Reason |
|---:|---|---|
| 10 | Trace (default BG, direct colour) | AmpliGluonTraceDefaultBGDirectColour benchmark failed with status 124:<br>/usr/bin/timeout --signal=TERM --kill-after=5s 1800s /usr/bin/time '--format=BENCHMARK_MAX_RSS_KIB %M' /usr/bin/prlimit --as=10737418240 -- /export/tmp/rikkert/git/MultipletRecursion2/Benchmark/build/benchmark_ampligluon_trace/benchmark_ampligluon_trace default-bg direct 0.25 10 /export/tmp/rikkert/git/MultipletRecursion2/Benchmark/build/events/N10/sampled/gg_to_8g_sample_000001.event /export/tmp/rikkert/git/MultipletRecursio… |
| 11 | Trace (default BG, direct colour) | AmpliGluonTraceDefaultBGDirectColour benchmark failed with status 124:<br>/usr/bin/timeout --signal=TERM --kill-after=5s 1800s /usr/bin/time '--format=BENCHMARK_MAX_RSS_KIB %M' /usr/bin/prlimit --as=10737418240 -- /export/tmp/rikkert/git/MultipletRecursion2/Benchmark/build/benchmark_ampligluon_trace/benchmark_ampligluon_trace default-bg direct 0.25 10 /export/tmp/rikkert/git/MultipletRecursion2/Benchmark/build/events/N11/sampled/gg_to_9g_sample_000001.event /export/tmp/rikkert/git/MultipletRecursio… |
| 11 | Adjoint (default BG, direct colour) | AmpliGluonAdjointDefaultBGDirectColour benchmark failed with status 124:<br>/usr/bin/timeout --signal=TERM --kill-after=5s 1800s /usr/bin/time '--format=BENCHMARK_MAX_RSS_KIB %M' /usr/bin/prlimit --as=10737418240 -- /export/tmp/rikkert/git/MultipletRecursion2/Benchmark/build/benchmark_ampligluon_adjoint/benchmark_ampligluon_adjoint default-bg direct 0.25 10 /export/tmp/rikkert/git/MultipletRecursion2/Benchmark/build/events/N11/sampled/gg_to_9g_sample_000001.event /export/tmp/rikkert/git/MultipletRe… |

## Reproducibility

| Total gluons | Generator seed | Combined event SHA-256 |
|---:|---:|:---|
| 4 | 1733 | `d777cc1433b1bc764684f98add195192d4adf57728248bc8ded69d4ed8d81ea4` |
| 5 | 1734 | `28d6d8421ae0a3b784257a6d9f7da2459eb2cfce2d260585a6050971218a7711` |
| 6 | 1735 | `74c5bc4edaf321abcdacc14c61d5b00cbf212d106222b6f0be95d9f4465366ce` |
| 7 | 1736 | `e95dabf8aae73f18ed9f7c58ebfa359faa50705ef5e0977797bd4dc25fbfe5a2` |
| 8 | 1737 | `dc5d086d288ddcebdee19cda75b61002dd62780e3a94ffa2d702c08f0065ee41` |
| 9 | 1738 | `b9c6bdfda81db611c5cac77c16093c49b2610a5ccc742c287364c3ce291e6984` |
| 10 | 1739 | `1eeed085ac9fdbad27d11c8980110e6ceb0ce3a83e21e601da2031865f7fbe38` |
| 11 | 1740 | `72c6173219fbc4c83c78a4413b5d8b8f260a51c7ed237714bbf038fb92940c18` |

Compiler: `GNU Fortran (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0`

Flags for selected backends: `-O3 -std=f2018 -Wall -Wextra`

Backend address-space cap: `10 GiB`; peak resident memory is measured with GNU time
Backend invocation timeout: `1.8e+03 s`
Separate initialization preflight: `enabled`
Minimum calibration target per sampled timing batch: `0.25 s`
Requested expensive samples per multiplicity: `3` MHV/anti-MHV and `1` general-helicity; one additional analytic-zero row is used only for validation.

The selected analytically vanishing row must be below 1.000e-24; selected nonzero rows are compared relatively. Zero sectors have zero proxy weight and are excluded from the calibrated timing cells and production estimate.
Wigner table: `/export/tmp/rikkert/git/MultipletRecursion2/Wigner6j/data/su3_adjoint_swap_prefix_6.tbl`
Wigner prefix depth: `6` (supports through `13` total gluons)
Wigner table SHA-256: `679875197c908b1ea0f2c4d2f44db1da0fb639b87108a391a6ecab47f9badbba`
