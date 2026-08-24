# AmpliGluonTrace

`AmpliGluonTrace` is the fundamental-trace all-gluon, full-colour
matrix-element implementation used to cross-check `AmpliGluonMultiplet`. It
evaluates tree-level, fixed-helicity

```text
g g -> n g,  n >= 2
```

amplitudes and returns the squared matrix element summed over all initial and
final SU(3) colours. It is self-contained and requires only a Fortran 2018
compiler; no external linear-algebra library is required.

This is deliberately not a general event generator. The particle database,
Standard Model initialization, quarks, electroweak vertices, phase-space and
event generation, PDFs, cuts, scale choices, helicity sums, process parsing,
reweighting, and code-generation machinery from AmpliCol are absent.

## Method

Through eleven total gluons, AmpliGluonTrace builds all cyclic colour orders
together in one Berends--Giele current graph, with one external gluon fixed.
This is the pure-gluon `use_symmetry=.true.` construction from AmpliCol.  At
every current length an ordered word and its reverse share one stored current,

```text
J(reverse(w)) = (-1)^(length(w)-1) J(w).
```

Each ordered child pair stores its signed target list once.  Its three-gluon
and auxiliary-tensor contributions are combined before a single accumulation
into every symmetry-related target current.  The auxiliary tensor is the
non-propagating representation of the four-gluon vertex.  Terminal currents
fill the full `(N-1)!` amplitude vector in lexicographic order; in particular
reversed orders obey
`A(reverse(1,...,N-1),N) = (-1)^N A(1,...,N-1,N)`.

The fully materialized current graph becomes less memory-efficient at twelve
total gluons.  At that point the library automatically uses the earlier
memory-bounded per-order recursion with terminal reflection pairing.  This
preserves the previously validated high-multiplicity reach while the shared
graph is made more compact; it does not change the public ordering or result.
The kinematic graph or fallback recursion is initialized lazily on the first
helicity configuration that needs it. MHV-only or vanishing workloads therefore
avoid both its initialization time and storage.

MHV and anti-MHV configurations use the Parke--Taylor cyclic-denominator
ratios to fill the same full trace-order vector. One canonical Berends--Giele
amplitude fixes the library's phase and normalization convention. Reflection
then supplies the reversed half of the orders exactly; this is a trace identity,
not a Kleiss--Kuijf reduction. Configurations with fewer than two all-outgoing
helicities of either sign return zero before building wavefunctions or currents.

The trace Gram matrix is a convolution on `S_(N-1)`.  With `g` and `h` the
lexicographic colour orders, it has the form

```text
C(g,h) = k(g^-1 h).
```

For pure-gluon tree amplitudes the U(3) kernel gives exactly the same result
as the SU(3) trace sum because every U(1) component decouples. By default, the
first nonzero evaluation constructs the `O((N-1)!)` kernel and transforms it
with a symmetric-group Fourier transform; `initialize` only computes the
factorial order count. Evaluation transforms the ordered partial amplitudes and
contracts the resulting Young-irrep blocks. The colour storage is therefore
linear in the number of orders instead of a packed `O(((N-1)!)^2)` Gram
matrix. Each transformed kernel block is real symmetric, so the quadratic form
evaluates each off-diagonal column overlap only once. The contraction uses
internal loops, so BLAS is optional rather than a dependency.

For comparison, `initialize` accepts `use_colour_fft=.false.`. This retains the
untransformed first Gram-matrix row and the lexicographic permutation data, then
performs the ordinary `A^dagger G A` multiplication. It maps every matrix entry
from `G(g,h)=k(g^-1 h)` rather than materializing the quadratic Gram matrix, and
uses symmetry to evaluate only the diagonal and strict upper triangle. This
normal path also has linear colour storage, but quadratic contraction time; it
is intended as a transparent reference and benchmark alternative to the
symmetric-group Fourier method.

The number of colour orders is `(N-1)!` for `N` total gluons, so this reference
implementation intentionally retains AmpliCol's factorial high-multiplicity
scaling.  Sharing the ordered subcurrents avoids independently rebuilding a
complete recursion for every order.  The checked physical regressions cover
four through nine total gluons (`g g -> 2...7 g`), including reflection and
U(1) decoupling identities for the ordered amplitudes. Separate mathematical
regressions compare the Fourier quadratic form with direct dense convolution
for arbitrary complex vectors through `S_6`.

## Build and test

From the repository root:

```sh
make -C AmpliGluonTrace
make -C AmpliGluonTrace test
make -C AmpliGluonTrace debug
```

The normal build creates
`AmpliGluonTrace/build/libampligluon_trace.a`. The debug target enables bounds
checks and treats compiler warnings as errors.

## Library interface

Input momenta are physical: the first two rows are incoming and all remaining
rows are outgoing, with positive energies. Helicity values are `-1` or `+1`.

```fortran
use, intrinsic :: iso_fortran_env, only: real64
use ampligluon_trace, only: ampligluon_trace_t

type(ampligluon_trace_t) :: amplitude
real(real64) :: p(0:3, 5), matrix2
integer :: helicity(5)

call amplitude%initialize(final_gluons=3, use_colour_fft=.true.)
call amplitude%evaluate(p, helicity, matrix2, &
                        strong_coupling=1.0_real64, &
                        average_initial_colours=.false., &
                        use_mhv_optimization=.true.)
```

`evaluate` optionally returns every complex ordered partial amplitude through
the `ordered_amplitudes` argument. They are ordered by lexicographic
permutations of legs `1,...,N-1`, with leg `N` fixed last. The default result
is colour-summed but not initial-colour averaged; averaging divides it by
`8^2`. Helicity is never summed or averaged. Pure-gluon tree amplitudes with
fewer than two positive or two negative all-outgoing helicities vanish.
`evaluate` accounts for crossing of the two incoming gluons and returns these
forbidden sectors as exact zeros without constructing ordered currents or
performing the colour transform. Such zero-only workloads also avoid allocating
the factorial colour and internal amplitude storage. Requesting the optional
ordered-amplitude output necessarily allocates its `(N-1)!` zero entries.

`use_colour_fft` is an optional initialization switch and defaults to `.true.`.
Set it to `.false.` to use the normal colour-matrix multiplication described
above. Colour data remain lazy in either mode: exact-zero-only workloads do not
construct the Fourier plan, raw kernel, or permutation table.

When the crossed all-outgoing configuration is MHV or anti-MHV, the library
evaluates one canonical Berends--Giele order to preserve its complex phase and
coupling convention, then fills every other trace order from the corresponding
Parke--Taylor denominator ratio. All `(N-1)!` trace orders are retained; no
Kleiss--Kuijf reduction is used. More general helicity sectors continue to use
the shared-current graph.

`use_mhv_optimization` is an optional per-evaluation switch and defaults to
`.true.`. Set it to `.false.` to route MHV and anti-MHV configurations through
the general Berends--Giele recursion. The exact zero-helicity selection rule
remains active regardless of this switch. NMHV and all other nonzero sectors
already use the general recursion in either setting.

The tree amplitudes are homogeneous in the strong coupling. Currents and the
colour contraction are therefore evaluated at unit coupling, after which the
squared result is scaled by `g_s^(2N-4)`. The optional ordered-amplitude vector
is scaled by `g_s^(N-2)`. This avoids repeating the same coupling multiplication
at every attachment in the shared-current graph.

The lower-multiplicity numerical sentinels originate from the original
AmpliCol `imode=2` implementation.  The nine-gluon sentinel was independently
cross-checked with AmpliGluonMultiplet and AmpliGluonAdjoint.  A test-only copy
of the former per-order recursion also compares every ordered amplitude with
the shared-current evaluator through nine total gluons. The top-level
three-way timing harness compares random RamboOnDiet points with
AmpliGluonMultiplet, AmpliGluonTrace, and AmpliGluonAdjoint:

```sh
make -C Benchmark
```
