# First-principles derivation and convention audit

## 1. What is being tabulated

The paper's Eq. (8) moves two neighboring gluons past one another along a
multiplet chain.  A bare “Wigner 6j” is not by itself a safe software
interface: its numerical value changes when a three-vertex is rescaled, when a
multiplicity basis is rotated, or when a vertex orientation convention is
changed.  Eq. (8), on the other hand, needs one definite number after all of
those choices have been made.

The generator therefore stores the matrix `W` in the normalized index
identity

$$
 B^{S\mu\nu}_{r a b;t}
 =\sum_{U,\kappa,\lambda}
 W^{R,T}_{U\kappa\lambda,S\mu\nu}\,
 \widetilde B^{U\kappa\lambda}_{r a b;t}.                 \quad (1)
$$

Here

- $R=\alpha_i$ and $T=\alpha_k$ are the fixed outer representations;
- $S=\alpha_j$ is the old middle representation;
- $U=\gamma$ is the new middle representation;
- $a=g_1$ and $b=g_2$ are adjoint indices;
- $\mu,\nu,\kappa,\lambda$ distinguish repeated three-vertices; and
- the tilde means that the coupling tree sees $g_2$ first and $g_1$
  second.

Thus a table row is already the full $\mathcal W_\gamma$ of Eq. (8), with
all vertex labels restored.  It is not merely the tetrahedral numerator in
that equation.

## 2. SU(3) and representation bases

### 2.1 Lie-algebra convention

Colors are numbered $1,2,3$.  On the fundamental representation,

$$
 E_{ij}|k\rangle=\delta_{jk}|i\rangle,
 \qquad [E_{ij},E_{k\ell}]
 =\delta_{jk}E_{i\ell}-\delta_{i\ell}E_{kj}.              \quad (2)
$$

The scalar product is the usual positive Hermitian product, so
$E_{ij}^{\dagger}=E_{ji}$.  The anti-fundamental action is the dual action,

$$
 E_{ij}|\bar k\rangle=-\delta_{ik}|\bar j\rangle.          \quad (3)
$$

All matrices used by the code are real.  This is possible for these tensor
components even though Hermitian Gell-Mann generators would include imaginary
entries: the non-Hermitian step operators $E_{ij}$ form an equivalent real
basis of the complexified Lie algebra.

### 2.2 Constructing $(p,q)$ without representation tables

Introduce normalized bosonic occupation states in

$$
 \mathrm{Sym}^{p}(\mathbf 3)\otimes
 \mathrm{Sym}^{q}(\overline{\mathbf 3}).             \quad (4)
$$

If $a_i^\dagger$ creates an upper index and $b_i^\dagger$ a lower index,
then

$$
 E_{ij}=a_i^\dagger a_j-b_j^\dagger b_i.                  \quad (5)
$$

The contraction operator is

$$
 K=\sum_i a_i b_i.                                        \quad (6)
$$

The traceless tensors

$$
 V_{(p,q)}=\ker K                                          \quad (7)
$$

form the desired irrep.  The contraction is onto the corresponding
$(p-1,q-1)$ symmetric tensor space.  Subtracting dimensions gives

$$
 \begin{aligned}
 \dim V_{(p,q)}
 &=\binom{p+2}{2}\binom{q+2}{2}
   -\binom{p+1}{2}\binom{q+1}{2}\\
 &=\frac{(p+1)(q+1)(p+q+2)}{2}.                            \quad (8)
 \end{aligned}
$$

The code constructs the kernel separately in every Cartan-weight sector.  It
does not assume the dimension as a way of selecting states; instead, it checks
that the discovered kernel has the dimension in (8).

### 2.3 Phase and weight-multiplicity convention

The highest state of $(p,q)$ is

$$
 (a_1^\dagger)^p(b_3^\dagger)^q|0\rangle                 \quad (9)
$$

with positive coefficient.  All descendants are generated recursively.  For
each state already found, the code tries $E_{21}$ and then $E_{32}$, in
that order.  At a weight of multiplicity greater than one, components along
previous states of that weight are removed by two-pass modified
Gram--Schmidt.  The surviving vector is divided by its positive norm; its sign
is **not** changed afterward.  Relative phases are therefore inherited from
the lowering operators.

When a kernel itself has dimension greater than one, its basis is fixed as
follows.  Form the basis-independent orthogonal projector onto the kernel,
project ambient coordinate vectors in lexicographic order, and
Gram--Schmidt the nonzero results.  The first significant ambient component
of each initial kernel vector is positive.  This fixes weight and vertex
multiplicity gauges without relying on arbitrary singular vectors returned by
the linear-algebra library.

## 3. Three-vertices from first principles

Let $A=(1,1)$ denote the adjoint.  A coupling vertex is represented by an
isometric intertwiner

$$
 C^{R\to S,\mu}:V_S\longrightarrow V_R\otimes V_A,
 \qquad
 (C^{R\to S,\mu})^\dagger C^{R\to S,\nu}
 =\delta_{\mu\nu}I_{V_S}.                                 \quad (10)
$$

The code discovers the decomposition of $R\otimes A$.  In each product
weight sector it finds the simultaneous kernel of the raising operators

$$
 E_{12}^{R\otimes A},\qquad E_{23}^{R\otimes A}.           \quad (11)
$$

A kernel of Dynkin weight $(p,q)$ and dimension $m$ means that $(p,q)$
occurs with multiplicity $m$.  Starting from each such highest vector, the
abstract target irrep first records the weight order reached by the lowering
recursion of Section 2.3.  A product copy is then constructed one complete
weight space at a time.  If $X$ is the unknown isometry on a target weight,
all already-constructed parents provide simultaneous lowering equations

$$
 A=XB.
$$

Here $A$ contains the product-representation lowerings and $B$ the matching
blocks of the fixed target generators.  The code takes the orthogonal
Procrustes factor (the polar factor of $AB^T$), which is the isometric
least-squares solution.  Combining both simple-root lowerings avoids the
path-dependent normalization drift of replaying one descendant at a time,
while preserving the target gauge.  The program then explicitly checks

$$
 (E^R\otimes I+I\otimes E^A)C=CE^S                  \quad (12)
$$

for all six step operators and both Cartan generators.

For $8\otimes8$, factor exchange commutes with SU(3).  The highest-weight
kernel is diagonalized under that exchange.  To fix more than just the parity,
identify an adjoint state with its traceless $3\times3$ matrix $X$.  The
symmetric octet is oriented along the traceless Jordan product

$$
 X\circ Y=XY+YX-\frac{2}{3}\mathrm{Tr}(XY)I,
$$

and the antisymmetric octet is oriented along the Lie bracket $[X,Y]$.
After normalizing the adjoints of those product maps as in (10), the two octet
vertices are

$$
 \mu=0:8_s\quad(P=+1),\qquad
\mu=1:8_a\quad(P=-1).                                   \quad (13)
$$

There is also a representation-independent physical anchor for every
nontrivial self channel $R\otimes8\to R$:

$$
 J_R(v\otimes X)=\rho_R(X)v.                              \quad (14)
$$

The normalized adjoint $J_R^\dagger$ is multiplicity 0.  When a second copy
exists, multiplicity 1 is its orthogonal complement, with phase fixed by the
projector/lexicographic rule of Section 2.3.  Thus the Fortran consumer can
identify vertex 0 with the ordinary color-generator action; it need not
infer an unknown rotation in a repeated $R$ channel.  The adjoint case keeps
the explicit ordering (13), so its generator/Lie-bracket vertex is
multiplicity 1 and its projection is explicitly
$X_{\rm left}\otimes X_{\rm right}\mapsto [X_{\rm left},X_{\rm right}].$  (Writing the generic action (14) with the
right adjoint acting on the left adjoint would reverse that bracket and add a
minus sign; this is why the adjoint exception is stated explicitly.)  These
gauges are fully specified, but should not be confused with unnamed
conventions in another CG package.

## 4. Re-deriving Eq. (8)

### 4.1 Why the dimension and closed-three-vertex factors occur

First allow a vertex $C:V_S\to V_R\otimes V_A$ to have arbitrary norm

$$
 C^\dagger C=n_C I_{V_S}.                                 \quad (15)
$$

Closing its output and input gives the theta/three-vertex graph

$$
 \Theta_C=\mathrm{Tr}(C^\dagger C)=n_C d_S.         \quad (16)
$$

The orthogonal projector onto this copy of $S$ is therefore

$$
 P_S=\frac1{n_C}CC^\dagger
     =\frac{d_S}{\Theta_C}CC^\dagger.                     \quad (17)
$$

This is exactly the factor $d_\gamma/(\text{closed 3j})$ in the first line
of the paper's Eq. (8).  In the present isometric convention $n_C=1$, so
$\Theta_C=d_S$ and this factor is one.  It is nevertheless important to
derive it: omitting it while using non-isometric vertices would give a wrong
swap rule.

After inserting (17), the triangular loop is an intertwiner with the same
three external representations as the remaining vertex.  Schur's lemma says
that, within each multiplicity channel, it is a scalar times that vertex.
Contracting with the adjoint vertex and using (16) determines the scalar: its
numerator is the closed tetrahedral contraction (the 6j graph), and its
denominator is the appropriate closed three-vertex.  This gives the second
ratio in Eq. (8).  Thus both denominators in the paper are consequences of
ordinary orthogonal projection and Schur's lemma, not optional graphical
normalizations.

### 4.2 The normalized index equation

With the isometric vertices (10), define the original tree

$$
 B^{S\mu\nu}_{r a b;t}
 =\sum_s
 C^{R\to S,\mu}_{ra;s}
 C^{S\to T,\nu}_{sb;t},                                   \quad (18)
$$

and the tree after moving $g_1$ past $g_2$, written in the same explicit
index order $(r,a,b;t)$,

$$
 \widetilde B^{U\kappa\lambda}_{r a b;t}
 =\sum_u
 C^{R\to U,\kappa}_{rb;u}
 C^{U\to T,\lambda}_{ua;t}.                               \quad (19)
$$

Repeated use of (10) proves that both sets are orthonormal:

$$
 \frac1{d_T}\sum_{r,a,b,t}
 (B^x_{rab;t})^*B^y_{rab;t}=\delta_{xy},                  \quad (20)
$$

and likewise for the tilded trees.  Completeness of the CG decomposition says
that they span the same intertwiner space.  Taking the scalar product of (1)
with a tilded tree now gives the swap coefficient without any graphical
ambiguity:

$$
 \boxed{
 W^{R,T}_{U\kappa\lambda,S\mu\nu}
 =\frac1{d_T}\sum_{r,a,b,t,u,s}
 \left(
 C^{R\to U,\kappa}_{rb;u}
 C^{U\to T,\lambda}_{ua;t}
 \right)^*
 C^{R\to S,\mu}_{ra;s}
 C^{S\to T,\nu}_{sb;t}}
                                                               \quad (21)
$$

Closing the four vertices in (21) draws the tetrahedron in the middle of the
paper's Eq. (8).  Equation (21) is therefore that equation after its
completeness and vertex-correction normalization factors have been evaluated
in one consistent isometric convention.

### 4.3 The reversed central vertex and its minus sign

The star/adjoint on the first parenthesis in (21) reverses the orientation of
those trivalent vertices when the contraction is drawn as a tetrahedron.  In
addition, (18) uses the adjoint indices in the order $b,a$, not $a,b$.
These two facts are the index meaning of the minus/orientation marker at the
central vertex in the paper's 6j diagram.

No universal extra minus sign should be multiplied onto (21).  The result
depends on the vertex:

- for a symmetric $8\otimes8\to8_s$ vertex, exchanging its two adjoint
  arguments gives $+C$;
- for the antisymmetric triple-gluon vertex
  $8\otimes8\to8_a$, it gives $-C$; and
- with multiplicities, reversal can in general act as a matrix rather than a
  single sign.

The implementation performs the literal index exchange in (18), so all three
possibilities are included automatically.  The lowest nontrivial calculation
provides a sharp sign audit.  Set $R=1$.  The first coupling is uniquely
$1\otimes8\to8$, and (21) reduces to exchange parity in $8\otimes8$:

| final channel $T$ | $1$ | $8_s$ | $8_a$ | $10$ | $\overline{10}$ | $27$ |
|---|---:|---:|---:|---:|---:|---:|
| swap coefficient | +1 | +1 | -1 | -1 | -1 | +1 |

The dimensions also check the parity split:
$1+8+27=36=8\cdot9/2$ symmetric states and
$8+10+10=28=8\cdot7/2$ antisymmetric states.

### 4.4 Consequences used as checks

Because the two external particles are identical adjoints, both sides of (1)
use the same path list.  The physical exchange operator is real,
self-adjoint, and squares to one.  Hence every fixed-$(R,T)$ block must obey

$$
 W^T=W,\qquad W^TW=I,\qquad W^2=I.                         \quad (22)
$$

The `explicit-index` backend checks (1) on every explicit tensor index before
accepting a block; (22) is an additional check, not a replacement for Eq. (8).
The other three full-table builders establish the same signed identity through
the Schur or elementary-path reductions derived below.  Direct Specht
establishes it on its seven phase-bridged blocks and checks gauge-invariant
projector consequences elsewhere.

## 5. Independent construction by splitting adjoint lines

### 5.1 What was taken from the papers

The papers were used only for a structural suggestion.  In particular,
[*Wigner 6j symbols with gluon lines*](https://arxiv.org/abs/2312.16688)
suggests replacing an adjoint line by a fundamental and an antifundamental
line and removing the singlet component; [*Wigner 6j symbols for SU(N):
Symbols with at least two quark-lines*](https://arxiv.org/abs/2209.15013)
emphasizes completeness and orthogonality on the resulting fundamental-line
diagrams.  No coefficient, normalization, sign, or displayed formula from
either paper is input to the code.  The older explicit-projector construction
in [arXiv:1809.05002](https://arxiv.org/abs/1809.05002) provides a useful
contrast to this route.  The sources and convention differences are catalogued
in [LITERATURE_AUDIT.md](LITERATURE_AUDIT.md).

The implementation below re-derives the split directly in the concrete
unitary modules of Sections 2 and 3.  It then gives a second calculation of
every entry of (21).  Both calculations deliberately share the canonical
three-vertices (10): a numerical 6j table cannot agree entry by entry unless
the phase and multiplicity gauges of its vertices agree.

### 5.2 Deriving the traceless fundamental split

Let $M_3$ be the nine-dimensional space of $3\times 3$ matrices with
orthonormal coordinate matrices $e_{ij}$.  Its normalized invariant vector
is

$$
 s=I/\sqrt3.                                               \quad (23)
$$

The adjoint is the orthogonal complement of this singlet.  Its projector is
therefore obtained without a Fierz formula:

$$
 Q=I_{M_3}-|s\rangle\langle s|.                            \quad (24)
$$

Choose the same orthonormal adjoint basis as in Section 2 and write its
matrix components as $S_{aij}$.  Orthonormality and (24) immediately give

$$
 \sum_{ij}S_{aij}^*S_{bij}=\delta_{ab},                    \quad (25)
$$

$$
 \sum_a S_{aij}S_{ak\ell}^*
 =\delta_{ik}\delta_{j\ell}
  -\frac13\delta_{ij}\delta_{k\ell}.                      \quad (26)
$$

Thus (26), often drawn as the color completeness or Fierz relation, is just
the matrix of an orthogonal projector.  The program constructs $S$ from
the traceless-kernel basis and independently checks both (25) and (26).
Moreover, (2)--(3) act on a matrix by

$$
 E_{mn}\mathbin{\cdot}M=E_{mn}M-ME_{mn}.                  \quad (27)
$$

This preserves the trace, proving that $S:V_8\to V_3\otimes V_{\bar3}$ is
an SU(3) intertwiner as well as an isometry.

### 5.3 A second formula for the swap matrix

Fix a path $x=(S,\mu,\nu)$ and any normalized vector $h\in V_T$.  Splitting
both adjoints in its coupling tree gives

$$
 F^x_{r i j k\ell}
 =\sum_{s,a,b,t}
 C^{R\to S,\mu}_{ra;s}C^{S\to T,\nu}_{sb;t}
 h_t S_{aij}S_{bk\ell}.                                   \quad (28)
$$

For a path $y=(U,\kappa,\lambda)$ after the exchange, expressed in the same
external index order, define

$$
 \widetilde F^y_{r i j k\ell}
 =\sum_{u,a,b,t}
 C^{R\to U,\kappa}_{rb;u}C^{U\to T,\lambda}_{ua;t}
 h_t S_{aij}S_{bk\ell}.                                   \quad (29)
$$

The alternative coefficient is simply

$$
 W^{\rm split}_{yx}=\langle\widetilde F^y,F^x\rangle.     \quad (30)
$$

Applying (25) twice reduces (30) to the four-vertex contraction (21), except
that it is evaluated on $h$ rather than traced over all of $V_T$.  To see
why this is sufficient, regard each unsplit tree as an isometry
$B_x:V_T\to V_R\otimes V_8\otimes V_8$.  The endomorphism
$\widetilde B_y^\dagger B_x$ commutes with SU(3), so Schur's lemma makes it
$W_{yx}I_{V_T}$.  Its normalized trace is (21), while its expectation value
in any normalized $h$ is (30).  They are consequently identical.  The code
chooses the unique normalized highest-weight state so it never has to assume
this equality numerically.  Moreover, the same fixed isometry $S\otimes S$
acts on every path and commutes with exchanging the two adjoint pairs.
Consequently it cancels exactly from both the overlap and reconstruction
equations.  The implementation validates $S$ independently but keeps the
equivalent compact $8\times8$ highest-weight column rather than materializing
its $9\times9$ image.

This produces a genuinely separate coefficient path:

- the original backend constructs every final-state column, swaps two
  eight-valued adjoint axes, and divides their full trace by $d_T$, streaming
  bounded slabs rather than retaining a complete block tensor;
- the split backend derives the exchange through the traceless
  fundamental--antifundamental resolution, cancels its common isometry, and
  contracts the packed nonzero entries of one final-state column; and
- each backend reconstructs its own coupling-tree equation on all indices it
  retained before a block is accepted.

For a block with $k$ paths, the principal stored chain arrays scale as
a fixed slab budget for the full-index backend and $64k d_R$ numbers for the
split/highest-weight backend.  The latter removes the final-irrep dimension
from this contraction stage.  The slab budget bounds the explicit chain payload,
not decoded coupling operands or overlap and reconstruction temporaries, so it
is not a bound on total process RSS.  Shared phase-fixed three-vertices are
retained as exact weight-sector blocks rather than dense zero-filled matrices.

## 6. Recursive reduction to elementary coefficients

Section 5 still evaluates each requested two-adjoint coefficient directly:
it merely uses split indices instead of adjoint indices.  This section derives
a third construction which actually writes each coefficient as sums and
products of lower-complexity recouplings.

Only the *topological hierarchy* in the quark-line and gluon-line papers
audited in [LITERATURE_AUDIT.md](LITERATURE_AUDIT.md) motivated this route:
replace adjoint vertices by fundamental-line paths, reduce a graph with two
adjoint lines to graphs with one, and terminate at graphs with two fundamental
lines.  None of the paper equations, normalization factors, closed forms, or
sign prescriptions below were copied or evaluated.  The following derivation
starts again from isometries and orthogonal completeness.

### 6.1 Terminal two-line recouplings

Let $x$ be either $3$ or $\bar3$.  For every summand $S$ of
$R\otimes x$, construct an isometry

$$
 A^{R,x\to S}:V_S\longrightarrow V_R\otimes V_x.
$$

The code discovers its highest vector as the simultaneous raising-operator
kernel and generates all descendants with the convention of Section 2.3.
Tensoring by $3$ or $\bar3$ is multiplicity-free, and the independently
checked identities are

$$
 (A^{R,x\to S})^\dagger A^{R,x\to U}=\delta_{SU}I,
 \qquad
 \sum_S A^{R,x\to S}(A^{R,x\to S})^\dagger=I_{R\otimes x}. \quad (31)
$$

For two elementary lines $x,y\in\{3,\bar3\}$, fixed outer irreps $R,T$,
and an allowed middle irrep $S$, form the normalized coupling tree

$$
 B^S=(A^{R,x\to S}\otimes I_y)A^{S,y\to T}.              \quad (32)
$$

Define $\widetilde B^U$ in the same way with $x,y$ exchanged and its two
explicit elementary indices restored to the order $x,y$.  Schur's lemma
and (31) give the terminal recoupling matrix

$$
 K^{R,T}_{U,S}
 =\frac1{d_T}\mathrm{Tr}
   [ (\widetilde B^U)^\dagger B^S ],
 \qquad
 B^S=\sum_U K^{R,T}_{U,S}\widetilde B^U.                 \quad (33)
$$

These are the simplest coefficients used by the recursion: they contain two
elementary lines, no adjoint line, and no call to a more complicated swap
builder.  The implementation checks Equation (33), orthogonality of $K$, and
reconstruction on the unique normalized highest-weight state before a terminal
block is cached.  Because both trees were already validated as intertwiners,
Schur's lemma promotes that one-column check to the full target irrep.  As an
elementary regression independent of any paper formula, the permutation of the
last two factors in $3\otimes3\otimes3$ is checked to be $-1$ in its singlet,
$+1$ in its decuplet, and the canonical two-dimensional reflection in the
two octet paths.

### 6.2 Expanding a fixed adjoint vertex

Apply (31) to $3\otimes\bar3$.  It contains one singlet isometry $s$ and
one adjoint isometry $J$, which the elementary coupling construction finds
without a Fierz formula.  Direct orthogonal completeness gives

$$
 J^\dagger J=I_8,
 \qquad JJ^\dagger+ss^\dagger=I_{3\otimes\bar3}.          \quad (34)
$$

For a fixed convention-carrying adjoint vertex
$C^{R\to S,\mu}:V_S\to V_R\otimes V_8$, let
$D^{R,S}_\lambda$ be the two-step elementary tree

$$
 R\xrightarrow{3}\lambda\xrightarrow{\bar3}S.
$$

The $D_\lambda$ form an orthonormal basis of the copies of $S$ in
$R\otimes3\otimes\bar3$.  Hence (34), with no graphical normalization
assumption, implies the exact expansion

$$
 (I_R\otimes J)C^{R\to S,\mu}
   =\sum_\lambda E^{R,S}_{\lambda\mu}D^{R,S}_\lambda,     \quad (35)
$$

$$
 E^{R,S}_{\lambda\mu}
  =\frac1{d_S}\mathrm{Tr}
    [(D^{R,S}_\lambda)^\dagger
     (I_R\otimes J)C^{R\to S,\mu}].                     \quad (36)
$$

The program checks reconstruction (35) and unit norm for each expansion on the
unique normalized highest-weight state.  Prior intertwiner validation and
Schur's lemma extend those checks over the target irrep; the final expansion
matrix is checked again for orthonormal columns.  This is where the table's
unavoidable vertex gauge enters: the $C$'s are the phase- and
multiplicity-fixed vertices of Section 3.
No value of a two-adjoint coefficient is used.  In particular, the two
$8\otimes8\to8$ columns remain the Jordan/symmetric and
Lie-bracket/antisymmetric vertices rather than an arbitrary rotation.

### 6.3 Two adjoints as four elementary lines

For an adjoint path $z=(S,\mu,\nu)$, applying (35) at its two vertices gives
a four-line path

$$
 R\xrightarrow{3}\lambda\xrightarrow{\bar3}S
  \xrightarrow{3}\kappa\xrightarrow{\bar3}T
$$

with expansion matrix

$$
 \mathcal E_{(\lambda,S,\kappa),z}
   =E^{R,S}_{\lambda\mu}E^{S,T}_{\kappa\nu},
 \qquad \mathcal E^T\mathcal E=I.                        \quad (37)
$$

Now label the elementary factors as
$q_1,\bar q_1,q_2,\bar q_2$.  Exchanging the two adjoint pairs is the
literal permutation

$$
 (q_1,\bar q_1,q_2,\bar q_2)
 \longmapsto(q_2,\bar q_2,q_1,\bar q_1).
$$

If $\tau_i$ exchanges elementary positions $(i,i+1)$, its action on every
left-associated path is one of the terminal matrices (33), with all spectator
labels held fixed.  Moving the second pair to the front gives, with the
rightmost operation acting first,

$$
 P_{(12)(34)}=\tau_2\tau_3\tau_1\tau_2.                  \quad (38)
$$

The first two elementary moves reduce the exchange of one adjoint pair with
one elementary line; the last two perform the second such reduction.  Thus
the two-adjoint result is a sum of products of one-adjoint recouplings, each
of which is a product of the terminal $K$'s.  The recursion terminates after
these two levels: no step creates another adjoint coefficient.

Projecting the elementary path permutation back with (37) yields

$$
 W^{\rm recursive}=\mathcal E^T P_{(12)(34)}\mathcal E.  \quad (39)
$$

This formula follows solely by substituting the exact basis expansion (35)
on both sides of the pair permutation.  Since $P_{(12)(34)}\mathcal E$ is
again in the two-adjoint subspace, the code checks

$$
 P_{(12)(34)}\mathcal E
   =\mathcal E W^{\rm recursive}
$$

before accepting a block.  A separate audit constructs every four-line tree
as an explicit tensor and verifies that the four terminal moves in (38)
equal the literal permutation of its four elementary indices.  Finally, the
result is compared block by block with both direct constructions of Sections
4 and 5.  These checks separately test the terminal coefficients, the vertex
expansion, the recursion, and the final projection.

The distinction between the three backends is therefore sharp:

- `explicit-index` contracts complete two-adjoint trees;
- `fundamental-split` contracts one split highest-weight column directly;
- `recursive-reduction` never contracts a two-adjoint coefficient graph; it
  assembles it from the terminal two-elementary-line matrices (33).

### 6.4 Removing the CG tensors as well

The recursion above terminates, but it is not CG-free. Equations (31), (33),
and (36) are evaluated there with explicitly constructed carrier-space
embeddings. They can instead be replaced by closed representation-ring data.
The resulting `dimension-only` backend has no `Irrep`, generator,
highest-weight kernel, or `embedding` object.

The elementary fusion graph is supplied by the Pieri rules

$$
 (p,q)\otimes3=(p+1,q)\oplus(p-1,q+1)\oplus(p,q-1),
$$

with terms having negative labels omitted; conjugation gives tensoring with
$\bar3$. For two identical elementary lines, lift $(p,q)$ to the three-row
partition $(p+q,q,0)$. If the two added boxes have contents $c_1,c_2$, put
$r=c_2-c_1$. Young's real seminormal adjacent transposition is

$$
 K_{33}=
 \begin{pmatrix}
  1/r & \sqrt{1-1/r^2}\\
  \sqrt{1-1/r^2} & -1/r
 \end{pmatrix}.
$$

When only one path exists this is the scalar $+1$ or $-1$. Conjugating and
reordering paths gives $K_{\bar3\bar3}$. Thus both magnitude and sign are
fixed by box contents, without a CG phase recursion.

For a mixed return path

$$
 (p,q)\longrightarrow\lambda\longrightarrow(p,q),
$$

let $d_X$ and $C_X$ denote the Weyl dimension and quadratic Casimir, with
$C_3=4/3$. In the repository's Young gauge the normalized singlet and
representation-action directions have components

$$
 g_\lambda=\phi_\lambda
       \sqrt{\frac{d_\lambda}{3d_R}},
 \qquad
 a_\lambda\ \mathrel{\propto}\
 \phi_\lambda\sqrt{d_\lambda}
       \frac{C_\lambda-C_R-C_3}{2}.
$$

Here $\phi_\lambda=-1$ for the second-row branch of the relevant Pieri rule
and $+1$ otherwise. The vectors $g$ and $a$ are orthogonal; when a third path
exists, $b=g\mathbin\times a$ completes an oriented orthonormal frame. The
mixed adjacent swap is the change of frame between the $3\bar3$ and $\bar3 3$
orders, with eigenvalues $+1,-1,+1$ on the singlet, action, and complement
directions. A non-returning mixed space is one-dimensional and has coefficient
$+1$ in this gauge.

The same frames give the adjoint-vertex reductions without evaluating (36).
A unique non-self edge has coefficient $+1$; a generic repeated edge uses
$a$ for multiplicity zero and $g\times a$ for multiplicity one. For
$8\otimes8\to8$, the repository ordering is instead the even Jordan vertex
followed by the odd Lie-bracket vertex. These are respectively
$(-a)\times g$ and $-a$ in the split frame.

The two adjoints may now be resolved combinatorially as
$3,\bar3,3,\bar3$. The literal pair permutation is still (38), but every
factor is now an analytic dimension/content matrix. Projecting with the
analytic vertex frames gives the same expression (39), with no carrier-space
Clebsch--Gordan decomposition at any stage. The implementation also obtains
the complete adjoint fusion multigraph directly from $(p,q)\otimes(1,1)$,
including its two self edges when $p,q>0$, and verifies its dimension sum.
Fusion labels have not disappeared--they define the paths on which a 6j
matrix acts--but CG tensors and bases have.

The relationship to representation-only formulae in the literature is
audited in [LITERATURE_AUDIT.md](LITERATURE_AUDIT.md). In particular, the
magnitude theorem in Eq. (32) of arXiv:2209.15013 agrees with the seminormal
magnitude after normalization conventions are converted, but its Eq. (26d)
and Appendix B.3 sign argument are mutually inconsistent for a concrete
SU(3) example. Neither is used here.

### 6.5 A direct symmetric-group alternative

There is also a genuinely non-recursive CG-free construction. By
Schur--Weyl duality, realize $(p,q)$ at tensor degree $m$ by a partition

$$
 \lambda_m(p,q)=(p+q+c,q+c,c),\qquad
 c=\frac{m-p-2q}{3}.
$$

Columns of height three are determinant factors and therefore trivial on
SU(3). Represent each adjoint by a Hermitian primitive Young idempotent of
shape $(2,1)$. On one three-line block a convenient $3\bar3$-oriented copy is

$$
 e_A=\frac{1-(23)}2-e_{\det},\qquad
 e_{\det}=\frac1{6}\sum_{\pi\in S_3}\operatorname{sgn}(\pi)\pi.
$$

Both terms are Hermitian permutation-algebra projectors, and $e_A^2=e_A$.
For a fixed final partition $\Lambda$, work only in the Young-seminormal
skew-tableau space $\Lambda/\lambda_R$ of the six boxes added by the two
adjoints. Let $Q$ be the product of the left-irrep and two adjoint idempotents,
and let $Z_S$ select the prefix shape belonging to the middle irrep $S$.
Orthonormal, convention-fixed bases of the ranges of $QZ_S$ form the columns
of $V$.

If $\tau$ is the literal permutation that exchanges the two three-index
blocks, the answer is obtained in one step:

$$
 W=V^T\rho_\Lambda(\tau)V.
$$

There is no expansion into four elementary-line recouplings here. Repeated
self channels are split within the permutation algebra: the action ray is
selected by the Jucys--Murphy/color-exchange sandwich

$$
 (e_Re_A)J_{m_0+1}(e_Re_{\det}),\qquad
 J_{m_0+1}=\sum_{i=1}^{m_0}(i,m_0+1),
$$

and its orthogonal complement supplies the second copy. In a seminormal basis
$J_{m_0+1}$ is diagonal, with the content of the new box as its eigenvalue.
The $8\otimes8\to8$ rays are then oriented as the even Jordan and odd bracket
channels. Projectors determine absolute values but not all relative signs, so
matching this repository also requires the local vertex phase anchors above.

The executable `DirectSpechtSwapTableBuilder` implements this phase bridge for
every valid block, while `DirectSpechtSwapOracle` exposes the same construction
and its projector diagnostics block by block. The
exact six-box skew spaces have dimension at most 90 over every block in the
supplied prefix-six table, whose largest path rank is 10. This makes the
permutation-algebra route practical as a structurally independent oracle.

## 7. Braid coherence between prefix layers

Consider three neighboring adjoints and a fusion path

$$
 R_0\xrightarrow{8}R_1\xrightarrow{8}R_2
     \xrightarrow{8}R_3.                                  \quad (40)
$$

Let $s_1$ exchange the first two adjoints and $s_2$ the last two.  In the
path basis,

$$
 s_1\big|_{R_0,R_2}=W^{R_0,R_2},\qquad
 s_2\big|_{R_1,R_3}=W^{R_1,R_3}.                           \quad (41)
$$

Literal permutations of three tensor factors give

$$
 s_1^2=s_2^2=I,
 \qquad s_1s_2s_1=s_2s_1s_2.                              \quad (42)
$$

Expanding (42) over every allowed intermediate representation and every
vertex multiplicity relates products and sums of blocks with prefix $R_0$
to blocks with the higher prefix $R_1$.  This is an independent coherence
identity, distinct from the terminating fundamental-line reduction in
Section 6.  `su3wigner/validation.py` constructs these matrices from the
table paths and checks (42) for every final $R_3$.  Starting at $R_0=1$
tests the seed parities against all octet-prefix blocks; proceeding to
$R_0=8$, then to every irrep in $8\otimes8$, recursively checks the higher
layers.

There is an important limitation to using (42) as the *sole* generator.  If a
vertex occurs with multiplicity $m$, replacing its vertex basis by an
orthogonal matrix $O\in O(m)$ changes the numerical entries of nearby
$W$ blocks, while all equations (22) and (42) remain true.  Therefore the
simplest 6j values plus braid coherence do not uniquely determine a
convention-specific table.  A phase/multiplicity gauge must be supplied.

The tensor backends resolve that underdetermination from first principles:

1. construct each higher irrep as the traceless kernel (7);
2. discover each product highest-weight kernel (11);
3. fix its multiplicity gauge by the canonical projector prescription;
4. recursively lower to obtain all CG components; and
5. calculate the swap by its full-index, highest-weight/Schur, or elementary-
   path formula, then enforce (22) and the corresponding index or closure
   equation.

The dimension-only backend reaches the same fixed gauge with the oriented
dimension/Casimir frames of Section 6.4; the direct symmetric-group method
uses local morphism orientations in the relevant projector ranges.  This is
safer than numerically solving the nonlinear braid equations and silently
inheriting an arbitrary multiplicity rotation.  When enabled, the separately
implemented CLI braid validator applies the prefix-layer coherence identity to
any full-table backend.

## 8. Finite table cutoff and Fortran-facing format

Pure gluon chains start at $(0,0)$ and repeatedly tensor with $(1,1)$, so
only triality-zero irreps occur.  A finite amplitude calculation needs only a
finite fusion depth.  `--max-prefix-gluons N` includes every distinct left
irrep first reachable in at most $N$ steps.  For each such left irrep, the
table contains all fixed-right blocks after two additional adjoints.  Raising
the cutoff requires no new formula or hand-entered seed.

The ASCII table stores:

1. fixed outer labels $(R,T)$ for each block;
2. every path $(S,\mu,\nu)$, including optional adjoint-exchange parities;
3. sparse entries `W[out_path,in_path]` in the convention (1).

Coefficients use `%.17e`: 17 digits after the decimal point, hence 18
significant decimal digits.  This round-trips an IEEE binary64 value.  The
format is loaded by the existing
[`AmpliGluonMultiplet` Fortran consumer](../AmpliGluonMultiplet/README.md#table-depth)
and avoids coupling runtime evaluation to Python, NumPy, a CG library, or a
closed-graph normalization convention.

The dense shape of an adjoint coupling is $(8d_R)\times d_S$, but weight
conservation permits entries only when the product-row and target-column
weights agree.  After phase anchoring and dense validation, the generator
packs exactly those allowed rectangular sectors into one coefficient buffer.
Off-sector values must be bitwise zero; packing rejects rather than hides any
violation.  Slab and column accessors let all carrier-space backends consume
this form directly.  The exported `Coupling.embedding` interface is unchanged
and materializes a cached dense array on public field access.  Dataclass
operations that inspect that field can trigger the same materialization.

The internal harmonic-to-irrep transform has the same exact weight-sector
structure and is packed as well.  The exported `Irrep` constructor and
`ambient_embedding` ndarray dataclass field retain their dense interfaces;
factory objects reconstruct and cache that ambient array only if a caller
requests it, except for the small adjoint ambient embedding, which remains
eagerly available.

## 9. Numerical scope

The construction is algebraic but evaluated in binary64.  It does not try to
recognize values such as $\sqrt3/2$ symbolically: the Fortran consumer uses
binary floating-point values, and symbolic recognition could introduce a
second convention/error path.  Across the combined audit suite, results are
overconstrained by Lie-algebra, isometry, completeness, direct Eq. (8), parity,
involution, projector, cross-backend, and braid checks.  A standalone backend
performs the subset appropriate to its derivation, as summarized in
[README.md](README.md#what-is-validated).  Tested residuals for the supplied
reference cutoff are at ordinary floating-point roundoff; current scaling and
deeper empirical reach are recorded in [PERFORMANCE.md](PERFORMANCE.md).
