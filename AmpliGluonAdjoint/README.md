# AmpliGluonAdjoint

`AmpliGluonAdjoint` is the Del Duca--Dixon--Maltoni (DDM) colour-basis
counterpart of `AmpliGluonTrace`. It evaluates the same tree-level,
fixed-helicity

```text
g g -> n g,  n >= 2
```

matrix elements and returns the result summed over all initial and final SU(3)
colours. The physical momentum, helicity, coupling, and colour-averaging
conventions are identical to AmpliGluonTrace and AmpliGluonMultiplet.

## Adjoint basis

For `N` total gluons, the colour decomposition is

```text
M = sum_sigma C_sigma A(1, sigma, N),   sigma in S_{N-2},
```

where `C_sigma` is a chain of adjoint structure constants. In the shared
normalization, `F^{abc} = Tr([T^a,T^b]T^c)` and no additional phase or
normalization multiplies `A(1,sigma,N)`. Consequently only `(N-2)!` partial
amplitudes are evaluated, compared with the `(N-1)!` trace orders retained by
AmpliGluonTrace.

The non-orthogonal DDM Gram matrix is constructed exactly for SU(3). Each
structure-constant chain is a nested commutator of fundamental traces. The
commutators factor into finite differences of a U(3) trace kernel; their U(1)
parts cancel exactly. Since each Gram entry depends only on the relative DDM
permutation, the default implementation uses a symmetric-group Fourier
transform and stores the matrix in `(N-2)!` coefficients rather than
materializing its quadratic number of entries.

For comparison and benchmarking, initialization accepts
`use_colour_fft=.false.`. This selects ordinary multiplication in the original
DDM basis. It retains one exact real row of the Gram matrix and obtains every
other row from the relative permutation, so storage remains linear in the
basis size (apart from the lexicographic permutation labels). The matrix is real
symmetric: diagonal terms are accumulated once and each upper-triangle term is
doubled. This path performs quadratic contraction work but neither initializes
nor stores a symmetric-group Fourier plan.

The equivalent fundamental completeness relation is

```text
sum_a T^a_ij T^a_kl = delta_il delta_jk - delta_ij delta_kl / 3.
```

For example, the four-gluon Gram matrix is

```text
[ 288  144 ]
[ 144  288 ]
```

in the normalization shared by this repository. The basis dimensions for four
through eleven total gluons are
`2, 6, 24, 120, 720, 5040, 40320, 362880`.

The Berends--Giele current graph is likewise restricted to the subcurrents
needed by the fixed-first DDM orders; it does not construct the unused full
`(N-1)!` trace-order output set.  Through eleven total gluons these currents
are evaluated together in one shared graph.  From twelve gluons the library
falls back to a memory-bounded recursion over individual DDM orders.

## Build and test

From the repository root:

```sh
make -C AmpliGluonAdjoint
make -C AmpliGluonAdjoint test
make -C AmpliGluonAdjoint debug
```

The normal build creates
`AmpliGluonAdjoint/build/libampligluon_adjoint.a`. The implementation reuses
the common colour-ordered kinematic primitives from `AmpliGluonTrace/src`;
the Makefile compiles them into the adjoint library, so no prebuilt
AmpliGluonTrace library is required.

## Library interface

```fortran
use, intrinsic :: iso_fortran_env, only: real64
use ampligluon_adjoint, only: ampligluon_adjoint_t

type(ampligluon_adjoint_t) :: amplitude
real(real64) :: p(0:3, 6), matrix2
integer :: helicity(6)

call amplitude%initialize(final_gluons=4, use_colour_fft=.true.)
call amplitude%evaluate(p, helicity, matrix2, &
                        strong_coupling=1.0_real64, &
                        average_initial_colours=.false., &
                        use_analytic_mhv=.true.)
```

Input momenta are physical: rows one and two are incoming and the rest are
outgoing, all with positive energies. Helicity values are `-1` or `+1`.
The default result is colour-summed and not initial-colour averaged; averaging
divides it by `8^2`. Helicity is never summed or averaged.

`number_of_basis_amplitudes()` returns `(N-2)!`. `evaluate` can optionally
return the DDM partial amplitudes through `adjoint_amplitudes`; they follow
lexicographic permutations of legs `2,...,N-1`, with legs `1` and `N` fixed.

`use_colour_fft` is an optional initialization argument and defaults to
`.true.`. Set it to `.false.` to use the normal, untransformed Gram-matrix
multiplication described above. Reinitializing an amplitude with a different
setting releases the storage belonging to the previous contraction method.

For pure-gluon tree amplitudes, `evaluate` applies the exact helicity selection
rule before constructing currents. After crossing the two incoming gluons to
the all-outgoing convention, configurations with fewer than two helicities of
either sign return an exactly zero adjoint-basis vector and squared matrix
element. A zero strong coupling likewise returns immediately.

MHV and anti-MHV sectors use analytic Parke--Taylor partial amplitudes by
default. When `adjoint_amplitudes` is requested, the evaluator computes one
canonical DDM order recursively to determine the event-wide
external-wavefunction phase, then fills the remaining orders analytically.
Thus the returned complex vector retains the same phase convention as the
fully recursive result.

Pass `use_analytic_mhv=.false.` to `evaluate` to force Berends--Giele
recursion for every DDM order of the same event. The switch is per evaluation,
so analytic and recursive matrix elements or complex vectors can be compared
without reinitializing the amplitude object. NMHV and all other nonzero
sectors always use recursion; the exact-zero selection rule remains
independent of the MHV switch.

The shared Berends--Giele current graph is constructed lazily on the first
recursive evaluation. Initialization and workloads that remain on analytic
MHV or exact-zero paths therefore avoid its setup cost.

The top-level `Benchmark` directory generates common RamboOnDiet inputs and
compares this library directly with AmpliGluonMultiplet and AmpliGluonTrace:

```sh
make -C Benchmark
```
