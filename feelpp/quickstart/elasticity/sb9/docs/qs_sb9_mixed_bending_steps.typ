// Step-by-step SB9 mixed-shell implementation note

#set page(
    paper: "a4",
    margin: (x: 18mm, y: 16mm),
    footer: context align(right)[
      #counter(page).display("1 / 1", both: true)
    ],
)

#set text(size: 10.5pt)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 16pt, weight: "bold")[Steps towards implementying SB9 into Feel++]

  #v(4pt)
  #text(size: 11pt)[Christophe Prud'homme (Cemosis)]
]

#v(8pt)

#outline()

#v(8pt)

= SB9 Mixed Shell Implementation Steps

This note summarizes an implementation path needed to build the mixed SB9
shell prototype in Feel++.
It is intentionally short and procedural: the goal is to explain the order of
the technical steps, why each one is needed, and which design decisions are
critical.

= Shell Geometry First

The first step is to add shell-specific geometric objects in `feelvf`:
`zeta()`, `shellThickness()`, `shellNormal()`, and `shellFrame()`.
SB9 is not a standard 3D solid model on stacked hexahedra: it is a shell model
carried by a 3D hexahedron with exactly one element through the thickness.
We therefore need a local shell thickness coordinate, a midsurface thickness,
and a consistent local shell frame before we can even define membrane, bending,
pinching, and transverse shear terms.

This step also leads naturally to a cached cell geometry layer:
center covariant and contravariant bases, local shell Jacobian, reduced SB9
Jacobians, and Hallquist data.
Without this cache, the later SB9 operators would recompute the same geometry
many times and become both fragile and expensive.

== DSEL proposition

At the DSEL level, shell geometry needs to appear as first-class geometric
terminals.
The quadrature-point reference coordinates and the cellwise midsurface data
need to be available directly inside variational expressions:

- `zeta()` for the reference thickness coordinate,
- `shellThickness()` and `shellArea0()` for cellwise shell measures,
- `shellNormal()` and `shellFrame()` for the local shell orientation,
- `shellCovariantBasis0()`, `shellContravariantBasis0()`,
  `shellMetric0()`, and `shellJacobian0()` for the reduced shell geometry.

This keeps the shell formulation readable.
Instead of rebuilding these objects by hand in each example, the formulation
can use them as native geometric ingredients of the language.

== Mathematical equivalent and corresponding C++ DSEL

#table(
  columns: (1fr, 1fr),
  stroke: 0.4pt,
  inset: 5pt,
  align: left + top,
  table.header(
    [*Mathematical equivalent*],
    [*Corresponding C++ DSEL*],
  ),
  [
    The proposed shell-geometric quantities are
    $
    zeta in [-1,1],
    $
    $
    h_K = frac( |K|, A_(0,K) ),
    $
    $
    n_K = frac( X_(K,xi)(0,0,0) times X_(K,eta)(0,0,0),
                || X_(K,xi)(0,0,0) times X_(K,eta)(0,0,0) || ),
    $
    $
    R_K = [ t_(1,K) quad t_(2,K) quad n_K ],
    $
    and, for a displacement field $u_h$, the shell-normal displacement is
    $
    u_(n,h)|_K = u_h dot n_K.
    $
  ],
  [
    ```cpp
    using namespace Feel;
    using namespace vf;

    auto zt = zeta();
    auto h = shellThickness();
    auto n = shellNormal();
    auto R = shellFrame();
    auto J0 = shellJacobian0();

    auto Sh = Pdh<0>( mesh );
    auto thicknessField = Sh->element();
    auto normalDisplacement = Sh->element();

    thicknessField.on(
        _range=elements( mesh ),
        _expr=shellThickness() );

    normalDisplacement.on(
        _range=elements( mesh ),
        _expr=inner( idv( uh ), shellNormal() ) );
    ```
  ],
)

= Symmetric Tensor Storage Is Necessary

The next step is to support symmetric tensor algebra in Voigt and Mandel
notations, together with basis vectors and component builders.
SB9 is naturally written as sums of membrane, bending, shear, and pinching
contributions inserted into a symmetric strain vector.
Trying to express that only with raw `3 times 3` tensors quickly becomes
verbose and error-prone.

This is why the implementation needs to introduce:
- tensor-basis support for symmetric storage,
- Voigt/Mandel vector builders,
- component insertion operators,
- higher-order constitutive tensors and contraction operators.

Mechanically, this matters because the constitutive law acts on specific strain
slots. A shell-normal correction, for example, must enter the normal slot and
not accidentally be mixed with a shear term.

== Mandel and Voigt notations

For a symmetric tensor
$
E = mat(
  E_11, E_12, E_13;
  E_12, E_22, E_23;
  E_13, E_23, E_33
),
$
the proposed Feel++ storage order is
$
[11, 12, 13, 22, 23, 33].
$

The two supported storage maps are:
- tensorial Voigt notation,
- Mandel notation.

The important point is that the proposed Feel++ Voigt notation is not the
engineering-strain convention.
Its off-diagonal entries stay equal to the tensor components themselves.
The difference between Mandel and Voigt therefore lies in the storage map and
in the inner product, not in a different isotropic stiffness matrix.

#table(
  columns: (1fr, 1fr),
  stroke: 0.4pt,
  inset: 5pt,
  align: left + top,
  table.header(
    [*Mathematical form*],
    [*Corresponding Feel++ DSEL*],
  ),

  [
    *Storage map*

    Tensorial Voigt:
    $
    V(E) = [ E_11, E_12, E_13, E_22, E_23, E_33 ]^T.
    $

    Mandel:
    $
    M(E) =
    [ E_11, sqrt(2) E_12, sqrt(2) E_13, E_22, sqrt(2) E_23, E_33 ]^T.
    $
  ],
  [
    ```cpp
    auto E_voigt = voigt( E );
    auto E_mandel = mandel( E );
    ```
  ],

  [
    *Inverse map*

    For
    $
    e_v = [ a, b, c, d, e, f ]^T,
    $
    $
    U_V(e_v) = mat(
      a, b, c;
      b, d, e;
      c, e, f
    ).
    $

    For
    $
    e_m = [ a, b, c, d, e, f ]^T,
    $
    $
    U_M(e_m) = mat(
      a, frac(b, sqrt(2)), frac(c, sqrt(2));
      frac(b, sqrt(2)), d, frac(e, sqrt(2));
      frac(c, sqrt(2)), frac(e, sqrt(2)), f
    ).
    $
  ],
  [
    ```cpp
    auto E_from_voigt = unvoigt( e_v );
    auto E_from_mandel = unmandel( e_m );
    ```
  ],

  [
    *Inner product*

    The Frobenius product satisfies
    $
    E : F = M(E)^T M(F).
    $

    In tensorial Voigt storage, the same quantity is
    $
    E : F =
      e_(v,11) f_(v,11)
      + 2 e_(v,12) f_(v,12)
      + 2 e_(v,13) f_(v,13) \
      + e_(v,22) f_(v,22)
      + 2 e_(v,23) f_(v,23)
      + e_(v,33) f_(v,33).
    $
  ],
  [
    ```cpp
    auto mf = inner( mandel( E ), mandel( F ) );
    auto vf = voigt_inner( voigt( E ), voigt( F ) );
    ```
  ],

  [
    *Constitutive action in Mandel notation*

    If
    $
    sigma = C epsilon,
    $
    then in Mandel storage
    $
    sigma_m = C e_m.
    $
  ],
  [
    ```cpp
    auto C = isotropic_stiffness<3>( lambda, mu );
    auto sigma_mandel = contract( C, mandel( E ) );
    ```
  ],

  [
    *Constitutive action in Voigt notation*

    If
    $
    sigma = C epsilon,
    $
    then in tensorial Voigt storage
    $
    sigma_v = C e_v.
    $
  ],
  [
    ```cpp
    auto C = isotropic_stiffness<3>( lambda, mu );
    auto sigma_voigt = voigt_contract( C, voigt( E ) );
    ```
  ],

  [
    *Constitutive bilinear form in Mandel notation*

    The constitutive bilinear form is
    $
    ( C e_m ) dot f_m.
    $
  ],
  [
    ```cpp
    auto C = isotropic_stiffness<3>( lambda, mu );
    auto a_mandel = contract( C, mandel( E ), mandel( F ) );
    ```
  ],

  [
    *Constitutive bilinear form in Voigt notation*

    The constitutive bilinear form is
    $
    ( C e_v )^T W f_v.
    $
  ],
  [
    ```cpp
    auto C = isotropic_stiffness<3>( lambda, mu );
    auto a_voigt = voigt_contract( C, voigt( E ), voigt( F ) );
    ```
  ],
)

== Operation dictionary

The previous table explains the storage maps.
The following table defines the main tensor operations used by the DSEL in
mathematical terms.

#table(
  columns: (1.1fr, 2.1fr, 1.8fr),
  stroke: 0.4pt,
  inset: 5pt,
  align: left + top,
  table.header(
    [*DSEL function*],
    [*Mathematical definition*],
    [*Typical DSEL usage*],
  ),

  [`inner(a,b)`],
  [
    Standard Euclidean inner product on vectors, or Frobenius product on
    tensors.
    In the Mandel setting,
    $
    M(E)^T M(F) = E : F.
    $
    This is why Mandel storage is convenient: the usual `inner` already
    reproduces the tensor contraction.
  ],
  [
    ```cpp
    auto s = inner( mandel( E ), mandel( F ) );
    auto t = inner( C*e_m, f_m );
    ```
  ],

  [`voigt_inner(a,b)`],
  [
    Weighted inner product on tensorial Voigt vectors, defined so that it still
    reproduces the Frobenius product:
    with
    $
    W = mat(
      1, 0, 0, 0, 0, 0;
      0, 2, 0, 0, 0, 0;
      0, 0, 2, 0, 0, 0;
      0, 0, 0, 1, 0, 0;
      0, 0, 0, 0, 2, 0;
      0, 0, 0, 0, 0, 1
    ),
    $
    $
    e_v^T W f_v =
      e_(v,11) f_(v,11)
      + 2 e_(v,12) f_(v,12)
      + 2 e_(v,13) f_(v,13) \
      + e_(v,22) f_(v,22)
      + 2 e_(v,23) f_(v,23)
      + e_(v,33) f_(v,33).
    $
  ],
  [
    ```cpp
    auto s = voigt_inner( voigt( E ), voigt( F ) );
    ```
  ],

  [`isotropic_stiffness<3>
                (lambda,mu)`],
  [
    Isotropic linear-elastic constitutive matrix associated with
    $
    sigma = lambda tr(epsilon) I + 2 mu epsilon.
    $
    In Mandel storage,
    $
    sigma_m = C_("iso") e_m,
    quad \
    e_m = [ epsilon_11, sqrt(2) epsilon_12, sqrt(2) epsilon_13,
            epsilon_22, sqrt(2) epsilon_23, epsilon_33 ]^T.
    $
    In the proposed Feel++ tensorial Voigt storage,
    $
    sigma_v = C_("iso") e_v,
    quad \
    e_v = [ epsilon_11, epsilon_12, epsilon_13,
            epsilon_22, epsilon_23, epsilon_33 ]^T.
    $
    With storage order
    $
    [11,12,13,22,23,33],
    $
    it is the matrix
    $
    C_("iso") = mat(
      lambda + 2 mu, 0, 0, lambda, 0, lambda;
      0, 2 mu, 0, 0, 0, 0;
      0, 0, 2 mu, 0, 0, 0;
      lambda, 0, 0, lambda + 2 mu, 0, lambda;
      0, 0, 0, 0, 2 mu, 0;
      lambda, 0, 0, lambda, 0, lambda + 2 mu
    ).
    $
    So in the proposed implementation the same matrix is used for Mandel and
    for Voigt.
    The difference is carried by the storage map and by the bilinear product:
    `inner` for Mandel, `voigt_inner` for tensorial Voigt.
  ],
  [
    ```cpp
    auto C = isotropic_stiffness<3>( lambda, mu );
    ```
  ],

  [`contract(C,eps)`],
  [
    Constitutive action in symmetric storage:
    $
    sigma_m = C e_m
    $
    in Mandel notation, or more generally the storage-compatible constitutive
    action selected by the overload.

    When used with three arguments,
    $
    a_m = (C e_m)^T f_m in RR,
    $
    which is a scalar and corresponds to the constitutive bilinear form.
  ],
  [
    ```cpp
    auto sigma = contract( C, e_m );
    auto a = contract( C, e_m, f_m );
    ```
  ],

  [`voigt_contract(C,eps)`],
  [
    Voigt-specific constitutive action.
    It is the shorthand for using the tensorial Voigt convention:
    $
    sigma_v = C e_v.
    $
    The associated bilinear form must then be interpreted with
    `voigt_inner`, or equivalently through the Voigt specialization of
    `contract`.
    In other words,
    $
    a_v = (C e_v)^T W f_v in RR,
    $
    so this row also produces a scalar.
  ],
  [
    ```cpp
    auto sigma_v = voigt_contract( C, voigt( E ) );
    auto a_v = voigt_contract( C, voigt( E ), voigt( F ) );
    ```
  ],
)

== DSEL proposition

The DSEL needs to provide symmetric-storage constructors and operators instead
of forcing the shell model to go back and forth between raw `3 times 3`
tensors and hand-written storage vectors.
The useful building blocks are:

- `mandel_vec<3>(...)` and `voigt_vec<3>(...)`,
- `mandel_basis<3,I,K>()` and `voigt_basis<3,I,K>()`,
- `mandel_component<3,I,K>(expr)` and `voigt_component<3,I,K>(expr)`,
- `isotropic_stiffness<3>(lambda, mu)` and `contract(C, eps)`,
- `contract(C, left, right)` for bilinear forms written directly in symmetric
  storage.

Mechanically, this is the right abstraction level because SB9 is written as a
sum of membrane, shear, pinching, and stabilization contributions that all
target precise symmetric-strain slots.

== Mathematical equivalent and corresponding C++ DSEL

#table(
  columns: (1fr, 1fr),
  stroke: 0.4pt,
  inset: 5pt,
  align: left + top,
  table.header(
    [*Mathematical equivalent*],
    [*Corresponding C++ DSEL*],
  ),
  [
    The symmetric shell strain is written in Mandel storage as
    $
    epsilon_("shell") = [ e_11, e_12, e_13, e_22, e_23, e_33 ]^T,
    $
    while the mixed SB9 scalar contributes only to the shell-normal slot:
    $
    epsilon_("alpha") = [ 0, 0, 0, 0, 0, gamma_alpha ]^T,
    quad
    gamma_alpha = -4 frac( zeta, h_K ) alpha_h.
    $

    With the isotropic stiffness tensor $C$, the bilinear contribution reads
    $
    a_("shell")( (u_h,alpha_h), v_h )
      = integral_(Omega_h) ( C ( epsilon_("shell") + epsilon_("alpha") ) ) : epsilon_("shell")(v_h).
    $

  ],
  [
    ```cpp
    using namespace Feel;
    using namespace vf;

    auto C = isotropic_stiffness<3>( lambda, mu );

    auto epsShell = mandel_vec<3>(
        e11,
        e12,
        e13,
        e22,
        e23,
        e33 );

    auto epsAlpha = mandel_component<3,2,2>(
        ( -4.0*zeta()/shellThickness() )*id( alpha ) );

    auto aShell = integrate(
        _range=elements( mesh ),
        _expr=contract( C, epsShell + epsAlpha, epsShell ) );
    ```
  ],
)

= Ordering Is a Structural Requirement

Two ordering conventions are critical.

First, the local displacement degrees of freedom of one Q1 hexahedron must be
read in the component-blocked order
$
[u_1,dots,u_8, v_1,dots,v_8, w_1,dots,w_8].
$
This is the order used when rebuilding the SB9 operators from the nodal data.

Second, the symmetric-storage order used by the shell operators is
$
[11, 12, 13, 22, 23, 33].
$
This point is not cosmetic.
A mismatch in the storage ordering produces mechanically wrong results while
the code still assembles and solves. In practice, a wrong ordering mixes
membrane, shear, and pinching contributions in the constitutive contraction.

== DSEL proposition

The DSEL needs to hide the ordering rules behind semantic SB9 operators.
The user should not have to remember local permutations every time a shell
strain is assembled.
In practice this means:

- the component-blocked local Q1 ordering stays encapsulated inside the SB9
  operator implementation,
- the symmetric-storage ordering stays encapsulated inside the Mandel and
  Voigt helpers,
- the user-facing formulation should work with semantic operators such as
  `sb9MembraneBending(u, zeta())`, `sb9Pinching(u)`, `sb9Shear(u, weight)`,
  and `sb9W9(alpha)`.

This is not only a convenience issue.
It is the safest way to make the formulation robust, because the ordering logic
is implemented once and reused everywhere.

== Mathematical equivalent and corresponding C++ DSEL

#table(
  columns: (1fr, 1fr),
  stroke: 0.4pt,
  inset: 5pt,
  align: left + top,
  table.header(
    [*Mathematical trial term*],
    [*Corresponding C++ DSEL term*],
  ),
  [
    Shared shell coordinate and shear weight:
    $
    zeta in [-1,1],
    quad
    w_s (zeta) = s_s ( 1 - zeta^2 ).
    $
  ],
  [
    ```cpp
    auto zt = zeta();
    auto shellShearWeight = cst( shellShearFactor ) * ( cst( 1.0 ) - zt*zt );
    ```
  ],

  [
    Membrane-bending trial term:
    $
    epsilon_("mb") (u_h; zeta) = B_("m0") (u_h) + zeta B_("b0") (u_h).
    $
  ],
  [
    ```cpp
    auto membraneTrial = sb9MembraneBending( u, zt );
    ```
  ],

  [
    Pinching trial term:
    $
    epsilon_("p") (u_h) = B_("pc") (u_h).
    $
  ],
  [
    ```cpp
    auto pinchingTrial = sb9Pinching( u );
    ```
  ],

  [
    Shear trial term:
    $
    epsilon_("s") (u_h; zeta) = w_s (zeta) B_("c0") (u_h).
    $
  ],
  [
    ```cpp
    auto shearTrial = sb9Shear( u, shellShearWeight );
    ```
  ],

  [
    Mixed scalar trial enrichment:
    $
    epsilon_("alpha") (alpha_h) = W_9 (alpha_h).
    $
  ],
  [
    ```cpp
    auto epsAlphaTrialMandel = sb9W9( alpha, cst( alphaScale ) );
    ```
  ],

  [
    Displacement shell strain in the fixed storage order
    $
    [11, 12, 13, 22, 23, 33]
    $
    is
    $
    epsilon_("shell") (u_h) \
      = [ epsilon_(( "mb", 11 )), epsilon_(( "mb", 12 )), epsilon_(( "s", 13 )), 
          epsilon_(( "mb", 22 )), epsilon_(( "s", 23 )), epsilon_(( "p", 33 )) ]^T.
    $
  ],
  [
    ```cpp
    auto epsShellTrial = vec(
        component<0,0>( membraneTrial ),
        component<1,0>( membraneTrial ),
        component<2,0>( shearTrial ),
        component<3,0>( membraneTrial ),
        component<4,0>( shearTrial ),
        component<5,0>( pinchingTrial ) );
    ```
  ],
)

= Build the Mixed Q1 x P0 Formulation

With the shell geometry and tensor tools in place, the next step is to build the mixed
finite element model:
the displacement lives in continuous vector Q1 and the extra SB9 scalar lives
in discontinuous P0, one scalar per cell.

The first version needs to be assembled monolithically with `blockform2` and
`blockform1`.
That choice is deliberate.
It keeps the formulation explicit and inspectable, which is essential while
checking signs, frame conventions, block couplings, and the effect of the extra
cellwise scalar.

This stage also exposes an important blockform issue:
strong Dirichlet elimination for a mixed product space acts on a row block, not
only on the `(0,0)` block.
That is why the row-view elimination path needs to be added for
`a.row(0_c) += on(...)`.

= Recover the Actual SB9 Mechanics

Once the mixed infrastructure is in place, the generic shell-gradient
prototype needs to be replaced by the actual SB9 kinematic operators.
This needs to introduce the dedicated operators for membrane-bending, pinching, and
assumed-shear terms, together with the Hallquist corrections and the
stabilization terms.

The reason is mechanical.
SB9 is not just "small strain plus one extra scalar".
Its good bending behavior comes from the specific split of the strain field,
the reduced shell geometry, and the stabilization of spurious modes.
Without these pieces, the model runs, but the bending response might not be(or is not) reliable.

== How the `B` operators act on the FE unknown

Hanna's question in section 1.4.1 of her note is the right one:
the SB9 matrices $B$ are built from the shape functions, the reduced Jacobian
evaluated at the element center, and the Hallquist vectors.
So what exactly should they multiply in the discrete formulation?

For one hexahedral shell element, let
$
  U_e = [ u_1, v_1, w_1, dots, u_8, v_8, w_8 ]^T
$
be the local vector of displacement coefficients, and let
$
  u_h (xi, eta, zeta) = sum_a N_a (xi, eta, zeta) u_a
$
be the Q1 displacement field on the element.
Then the SB9 strain operator acts as
$
  epsilon_h (xi, eta, zeta) = B_e (xi, eta, zeta) U_e
  = sum_a B_a (xi, eta, zeta) u_a.
$

At this stage, $U_e$ is only a *generic local coefficient vector* used to
write the element kinematics compactly.
When assembling the finite element matrix, one does not keep a generic $U_e$:
one applies the same operator column by column to the local trial basis
functions.
Equivalently, each column of $B_e$ is obtained by taking $U_e = e_j$, or by
evaluating the operator on one local basis function at a time.

The important point is that each block $B_a$ already contains the differentiated
shape-function contribution, together with the geometric transformation through
$bar(J)^(-T)$ and the Hallquist corrections.
Therefore one must *not* build a matrix-valued $B$ from those quantities and
then multiply it again by `idt(u)`.
In Feel++ notation, if `u` is the trial proxy associated with the displacement
space, `idt(u)` represents the trial-basis contribution carried by that proxy.
That is precisely where the Q1 interpolation enters the discrete form, namely
$
  u_h (xi, eta, zeta) = sum_a N_a (xi, eta, zeta) u_a
$
at the element level.
Therefore, if the matrix-valued $B$ has already been assembled from the
differentiated shape functions, the reduced Jacobian, and the Hallquist data,
doing `B * idt(u)` would reintroduce the FE basis a second time and would not
represent the intended SB9 kinematics.

There are therefore two equivalent viewpoints:

- *local matrix viewpoint*:
  build the element matrix $B_e$ from the shape functions, $bar(J)$, and the
  Hallquist vectors, then apply it to the local coefficient vector $U_e$;
- *DSEL/operator viewpoint*:
  define an operator that takes the trial or test proxy and internally
  contracts the SB9 coefficients with the corresponding FE basis
  contribution.

In Feel++, the `sb9*` operators should follow the second viewpoint.
That is why the right usage is
`sb9Bm0(u)`, `sb9Bb0(u)`, `sb9Pinching(u)`, `sb9Shear(u, weight)`, etc.,
where `u` is the trial or test FE proxy.
The operator hides the contraction with the local nodal unknowns and basis
functions attached to that proxy.
At the bilinear-form level, Feel++ then assembles the stiffness matrix from
expressions such as
$
  a_K (u_h, v_h) = integral_(K) (B (v_h))^T D B (u_h),
$
that is, by applying the operator to the actual trial and test basis functions.
So the matrix assembly is not "multiply an already built $B_e$ by `idt(u)`":
it is "evaluate the SB9 operator on the trial/test basis and integrate the
resulting bilinear form".

Said differently:

- if one writes a *manual local-matrix implementation*, $B$ multiplies the
  local coefficient vector $U_e$;
- if one writes a *Feel++ operator*, the operator itself must represent the
  action of the SB9 strain-displacement operator on the trial/test proxy, and
  one should not multiply the resulting expression by `idt(u)` again.

== Definitions of the `sb9*` keywords

The meaning of the `sb9*` keywords is aligned with the notation used in
`Implementation_SB9g25.pdf` and in the Matlab reference implementation,
especially `matrice_Ke.m`, `matrice_Bm0.m`, `matrice_Bb0.m`,
`matrices_Bp0.m`, `matrices_BcS.m`, and `matrices_Bs.m`.

#table(
  columns: (1.2fr, 1.8fr, 2.4fr),
  stroke: 0.4pt,
  inset: 5pt,
  align: left + top,
  table.header(
    [*Keyword*],
    [*Mathematical meaning*],
    [*Mechanical role in SB9*],
  ),

  [`sb9Bm0(u)`],
  [
    $
    B_("m0") (u_h)
    $
    as in
    $
    epsilon_("mb") (u_h; zeta) = B_("m0") (u_h) + zeta B_("b0") (u_h).
    $
  ],
  [
    Midsurface membrane part of the shell strain.
    It carries the in-plane terms at $zeta = 0$ and provides the base membrane
    contribution before the through-thickness bending correction is added.
  ],

  /*[`sb9Bb0(u)`],
  [
    $
    B_("b0") (u_h)
    $
    in the same split
    $
    epsilon_("mb") (u_h; zeta) = B_("m0") (u_h) + zeta B_("b0") (u_h).
    $
  ],
  [
    Bending part of the shell strain.
    It is multiplied by the thickness coordinate $zeta$ and therefore controls
    the linear variation of the membrane strain through the thickness.
  ],

  [`sb9Bpc(u)`],
  [
    $
    B_("pc") (u_h)
    $
    entering the pinching term
    $
    epsilon_("p") (u_h) = B_("pc") (u_h).
    $
  ],
  [
    Pinching or thickness-normal strain contribution.
    In the current mixed formulation it supplies the displacement-induced part
    of the $33$ strain component.
  ],

  [`sb9Bc0(u)`],
  [
    $
    B_("c0") (u_h)
    $
    used in
    $
    epsilon_("s") (u_h; zeta) = w_s (zeta) B_("c0") (u_h).
    $
  ],
  [
    Main assumed transverse-shear operator.
    It represents the physical shear part kept in the shell strain after
    weighting by the transverse shape function.
  ],

  [`sb9Bc1(u)`],
  [
    Auxiliary assumed-shear operator corresponding to the Matlab matrix
    $B_("c1")$, commented there as "a multiplier par eta".
  ],
  [
    First transverse-shear stabilization operator.
    It does not define the physical shell strain directly; it enters the
    additional shear-stabilization bilinear form.
  ],

  [`sb9Bc2(u)`],
  [
    Auxiliary assumed-shear operator corresponding to the Matlab matrix
    $B_("c2")$, commented there as "a multiplier par ksi".
  ],
  [
    Second transverse-shear stabilization operator.
    Like `sb9Bc1`, it acts in the stabilization part rather than in the main
    shell strain.
  ],

  [`sb9Bs1(u)`],
  [
    First stabilization mode
    $
    B_("s1") (u_h).
    $
  ],
  [
    Scalar stabilization mode associated with the first `Vgamma` family.
    In the Matlab reference it contributes to the extra stabilization matrix
    `Ks`.
  ],

  [`sb9Bs2(u)`],
  [
    Second stabilization mode
    $
    B_("s2") (u_h).
    $
  ],
  [
    Scalar stabilization mode associated with the second `Vgamma` family.
    It complements `sb9Bs1` for the shell-normal residual modes.
  ],

  [`sb9Bs3(u)`],
  [
    Third stabilization mode
    $
    B_("s3") (u_h).
    $
  ],
  [
    Two-component stabilization mode.
    In the current implementation it is built from the local in-plane
    directions and the third `Vgamma` family.
  ],

  [`sb9Bs4(u)`],
  [
    Fourth stabilization mode
    $
    B_("s4") (u_h).
    $
  ],
  [
    Three-component stabilization mode.
    It is the most general of the four `Bs*` operators and completes the
    stabilization space used for the additional `Ks` contribution.
  ],

  [`sb9MembraneBending(u,zeta)`],
  [
    Composite operator
    $
    epsilon_("mb") (u_h; zeta)
      = B_("m0") (u_h) + zeta B_("b0") (u_h).
    $
  ],
  [
    Direct DSEL access to the membrane-bending strain part in Mandel storage.
    It is the operator that is usually used in forms instead of manipulating
    `sb9Bm0` and `sb9Bb0` separately.
  ],

  [`sb9Pinching(u)`],
  [
    Composite operator
    $
    epsilon_("p") (u_h) = B_("pc") (u_h).
    $
  ],
  [
    Direct DSEL access to the pinching or $epsilon_33$ part generated by the
    displacement field.
  ],

  [`sb9Shear(u,weight)`],
  [
    Composite operator
    $
    epsilon_("s") (u_h; zeta) = w_s (zeta) B_("c0") (u_h).
    $
  ],
  [
    Direct DSEL access to the assumed transverse-shear part after weighting by
    the chosen thickness function.
  ],

  [`sb9ShellStrain(u,zeta,weight)`],
  [
    Final assembled shell strain
    $
    epsilon_("shell") (u_h)
      = \
      [ epsilon_(( "mb", 11 )), epsilon_(( "mb", 12 )), epsilon_(( "s", 13 )),
          epsilon_(( "mb", 22 )), epsilon_(( "s", 23 )), epsilon_(( "p", 33 )) ]^T.
    $
  ],
  [
    Ready-to-use displacement shell strain in the internal symmetric-storage
    order expected by the constitutive operators.
  ],*/

  [`sb9W9(alpha)`],
  [
    Mixed normal enrichment
    $
    W_9 (alpha_h) = -4 frac( zeta, h_K ) alpha_h
    $
    inserted in the $33$ strain slot.
  ],
  [
    Additional cellwise SB9 mode carried by the discontinuous scalar field.
    It plays the role of the extra normal mode in the mixed formulation.
    In the current Feel++ implementation this is the public mixed operator,
    while the Matlab reference also exposes a separate `Bpz` contribution.
  ],
)

= Then Add Static Condensation

After the monolithic mixed system is validated, static condensation needs to be added
for the `2 times 2` block system.
The monolithic path stays unchanged, and the solve strategy becomes a runtime
choice between monolithic and static condensation.

This order is important.
Condensing too early hides assembly and coupling bugs.
Doing it second lets us compare the condensed and monolithic responses on the
same test cases.

Here the right interpretation is mechanical as well as algebraic:
the P0 scalar is an internal cell mode, so it is natural to eliminate it
locally after its role in the element stiffness has been accounted for.

= Validate Only After the Framework Was Stable

Only after the shell geometry, tensor algebra, SB9 operators, mixed assembly,
and condensation path are in place does it make sense to extract the cases from
`validation.pdf` into JSON specifications.
The benchmarks, loads, constraints, and meshes are encoded in json specs, including marked points for nodal constraints and nodal loads.


= Recap

So the actual implementation path is:
1. shell geometry and cell cache,
2. symmetric tensor storage and contractions,
3. SB9-specific operator layer,
4. mixed $Q_1 times P_0$ monolithic assembly,
5. row-block Dirichlet treatment and benchmark plumbing,
6. static condensation for the $2 times 2$ system,
7. JSON validation cases extracted from `validation.pdf`.

= Benchmark Matrix

The benchmark layer needs to stay outside the C++ code and be driven by JSON
specifications.
This keeps the quickstart focused on the formulation and makes it possible to
change meshes, loads, markers, and references without editing the solver.

#table(
  columns: (1.7fr, 0.9fr, 2.3fr, 3.1fr),
  inset: 4pt,
  stroke: 0.4pt,
  align: left,
  table.header(
    [*Case*],
    [*Kind*],
    [*Setup*],
    [*What it tests or checks*],
  ),

  [`square-plate`],
  [benchmark],
  [Rectangular shell, all outer edges clamped, pressure on `ZPlus`.],
  [Global bending response of a thin plate and center transverse deflection under pressure.],

  [`bending-patch`],
  [benchmark],
  [Rectangular patch, all outer edges clamped, pressure on `ZPlus`.],
  [Patch-type bending behavior on a coarse elongated mesh and sensitivity of the mixed shell response.],

  [`circular-plate`],
  [benchmark],
  [Circular shell, outer boundary clamped, pressure on `ZPlus`.],
  [Axisymmetric bending behavior, center deflection, and robustness of the circular shell mesh and markers.],

  [`validation-traction-x`],
  [validation],
  [Clamp `XMoins`, apply a unit total force on `XPlus`.],
  [Exact axial response in the first material direction and consistency of shell-frame and constitutive slots.],

  [`validation-traction-y`],
  [validation],
  [Clamp `YMoins`, apply a unit total force on `YPlus`.],
  [Exact axial response in the second material direction and symmetry of the in-plane formulation.],

  [`validation-traction-z`],
  [validation],
  [Clamp `ZMoins`, apply a unit total force on `ZPlus`.],
  [Exact shell-normal response and correct treatment of the thickness direction.],

  [`validation-traction-x-partial`],
  [validation],
  [PDF point constraints on `P1,P3,P7,P5`, unit total force on `XPlus`.],
  [Exact axial response with pointwise blocking, marked points, and partial kinematic constraints.],

  [`validation-traction-y-partial`],
  [validation],
  [PDF point constraints on `P2,P6,P5,P1`, unit total force on `YPlus`.],
  [Exact axial response with pointwise blocking in the second in-plane direction.],

  [`validation-traction-z-partial`],
  [validation],
  [PDF point constraints on `P1,P2,P4,P7`, unit total force on `ZPlus`.],
  [Exact shell-normal response with pointwise blocking and one-layer shell kinematics.],

  [`validation-moment-x`],
  [validation],
  [Clamp `XMoins`, apply an equivalent unit moment on `XPlus` by tangential traction.],
  [Consistency of the moment loading path and bending/shear coupling about the first axis.],

  [`validation-moment-y`],
  [validation],
  [Clamp `YMoins`, apply an equivalent unit moment on `YPlus` by tangential traction.],
  [Consistency of the moment loading path and bending/shear coupling about the second axis.],

  [`validation-moment-z`],
  [validation],
  [Clamp `ZMoins`, apply an equivalent unit moment on `ZPlus` by tangential traction.],
  [Consistency of in-plane twisting actions and tangential traction moments.],

  [`validation-flexion-z-moment`],
  [validation],
  [Clamp `ZMoins`, apply a unit resultant moment on `ZPlus`.],
  [Flexural response generated by a top-face moment with a nonzero Poisson ratio.],

  [`validation-cantilever-z`],
  [validation],
  [Clamp `XMoins`, apply a total `-10` N force on `XPlus` in the `z` direction.],
  [Cantilever bending deflection against the `validation.pdf` reference on a `6 times 6 times 1` shell mesh.],

  [`validation-flexion-corner`],
  [validation],
  [PDF point constraints on `P1,P6,P8`, unit point load on `P3`.],
  [Corner bending response, marked-point loading, and pointwise essential constraints on one shell cell.],
)

= Public DSEL API

The table below is restricted to the general tensor-basis, Mandel/Voigt, and
contraction layer on which the SB9 implementation is built.
The SB9-specific public operators will be documented separately.

In the table below, we reuse the following objects:

- $d$ is the space dimension and `Dim` is its C++ counterpart.
- $A$ is a generic second-order tensor in dimension $d$.
- $S$ and $T$ are symmetric second-order tensors.
- $s_M$ and $s_V$ denote the Mandel and Voigt storage vectors associated with a symmetric tensor.
- $epsilon$ and $eta$ denote symmetric strain tensors, or their stored counterparts when the context is clear.
- $C$ denotes a constitutive matrix acting on symmetric-storage vectors.
- $lambda, mu$ are the Lamé coefficients.
- the code indices are zero-based, and the Feel++ 3D storage order is `(00,01,02,11,12,22)`.

#show raw: set text(size: 0.82em)

#table(
  columns: (2.45fr, 2.7fr, 3.25fr),
  inset: 4pt,
  stroke: 0.4pt,
  align: left,
  table.header(
    [*Name*],
    [*Description*],
    [*Example calls*],
  ),

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`delta<Dim,I,K>()`],
    [`symm_delta<Dim,I,K>()`],
    [`mandel_delta<Dim,I,K>()`],
  )],
  [Operation:
   `delta` gives $E_(i,k)$,
   `symm_delta` gives $E_(i,k)+E_(k,i)$,
   and, for $i != k$,
   `mandel_delta` gives $(E_(i,k)+E_(k,i))/sqrt(2)$.
   These are the basic tensor-basis expressions used to expose components and contractions explicitly.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`delta<3,0,2>()`],
    [`symm_delta<3,0,2>()`],
    [`mandel_delta<3,0,2>()`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`trans( delta<...>() )`],
    [`sym( delta<...>() )`],
    [`inner( A, delta<...>() )`],
  )],
  [Operation:
   `trans` maps $E_(i,k) -> E_(k,i)$,
   `sym` maps $E_(i,k) -> (E_(i,k)+E_(k,i))/2$,
   and `inner` extracts $A_(i,k)$ from $A$.
   This is the basic manipulation layer on top of the tensor basis.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`trans( delta<3,0,2>() )`],
    [`sym( delta<3,0,2>() )`],
    [`inner( A, mandel_delta<3,0,2>() )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`mandel( expr )`],
    [`voigt( expr )`],
  )],
  [Operation:
   `mandel` maps $S -> s_M$ with Mandel scaling on off-diagonal terms, and
   `voigt` maps $S -> s_V$ in classical Voigt storage.
   Mandel preserves the Euclidean tensor inner product directly, whereas Voigt uses the classical engineering storage.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`mandel( symm_grad( idv( u_h ) ) )`],
    [`voigt( symm_grad( idv( u_h ) ) )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`unmandel( expr )`],
    [`unvoigt( expr )`],
  )],
  [Operation:
   `unmandel` maps $s_M -> S$ and `unvoigt` maps $s_V -> S$.
   These inverse maps reconstruct a symmetric tensor from its Mandel or Voigt storage representation.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`unmandel( mandel( symm_grad( idv( u_h ) ) ) )`],
    [`unvoigt( voigt( symm_grad( idv( u_h ) ) ) )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`symm_storage_vec<Dim>( ... )`],
    [`mandel_vec<Dim>( ... )`],
    [`voigt_vec<Dim>( ... )`],
  )],
  [Operation:
   build $s = [s_0, dots, s_(m-1)]^T$ directly in Feel++ symmetric-storage order, either generically or with the Mandel and Voigt aliases.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`mandel_vec<3>( e00, e11, e22, e01, e12, e02 )`],
    [`voigt_vec<3>( e00, e11, e22, e01, e12, e02 )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`symm_storage_basis<Dim,I,K>()`],
    [`mandel_basis<Dim,I,K>()`],
    [`voigt_basis<Dim,I,K>()`],
  )],
  [Operation:
   build the storage-space basis vector associated with the $(i,k)$ slot.
   These are the symmetric-storage analogs of `delta`, `symm_delta`, and `mandel_delta`.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`mandel_basis<3,0,2>()`],
    [`voigt_basis<3,0,2>()`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`symm_storage_component<...>( value )`],
    [`mandel_component( ... )`],
    [`voigt_component( ... )`],
  )],
  [Operation:
   `mandel_component` maps $c -> c e_(i,k)^M$ and
   `voigt_component` maps $c -> c e_(i,k)^V$.
   A scalar field is injected into a single symmetric-storage slot while all other slots stay zero.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`mandel_component<3,2,2>( gamma33 )`],
    [`voigt_component<3,0,2>( gamma13 )`],
  )],

  [`scale_symm_storage( scale, expr )`],
  [Operation:
   $(c,s) -> c s$ inside symmetric-storage space.
   The scaling is done componentwise without leaving the stored representation.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`scale_symm_storage( cst( 0.5 ), eps )`],
    [`scale_symm_storage( shearWeight, eps )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`voigt_inner( left, right )`],
  )],
  [Operation:
   $(s_V,t_V) -> S:T$.
   This is the Voigt-storage inner product with the engineering weights required to match the tensor inner product.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`voigt_inner( voigt( S ), voigt( T ) )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`isotropic_stiffness<Dim>( lambda, mu )`],
    [`isotropic_stiffness<Dim,`],
    [`SymmetricTensorNotation::Voigt>( lambda, mu )`],
  )],
  [Operation:
   `isotropic_stiffness` maps $(lambda,mu)$ to the isotropic constitutive matrix $C$ in Mandel or Voigt storage.
   These are the standard isotropic constitutive matrices used by the contraction API.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`isotropic_stiffness<3>( lambda, mu )`],
    [`isotropic_stiffness<3,`],
    [`SymmetricTensorNotation::Voigt>( lambda, mu )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`contract( C, eps )`],
    [`contract( C, eps, eta )`],
  )],
  [Operation:
   $(C,epsilon) -> sigma = C:epsilon$,
   or
   $(C,epsilon,eta) -> epsilon:C:eta$.
   This is the generic constitutive contraction vocabulary.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`contract( C, eps )`],
    [`contract( C, eps, eta )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`ddot( C, eps_u, eps_v )`],
    [`ddot<SymmetricTensorNotation::Voigt>(`],
    [`  C, eps_u, eps_v )`],
  )],
  [Operation:
   $(C,epsilon,eta) -> epsilon:C:eta$,
   or
   $(C,epsilon) -> C:epsilon$.
   This is the most common DSEL spelling for the constitutive double contraction, on tensors or on stored vectors.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`ddot( C, symm_grad( u ), symm_grad( v ) )`],
    [`ddot<SymmetricTensorNotation::Voigt>(`],
    [`  Cvoigt, voigt( symm_grad( u ) ),`],
    [`  voigt( symm_grad( v ) ) )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`double_contract( C, eps )`],
    [`double_contract( C, eps, eta )`],
  )],
  [Operation:
   same map as `ddot`, namely
   $(C,epsilon) -> C:epsilon$
   and
   $(C,epsilon,eta) -> epsilon:C:eta$.
   This alias emphasizes the double-contraction viewpoint explicitly.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`double_contract( C, eps )`],
    [`double_contract( C, eps, eta )`],
  )],

  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`voigt_contract( Cvoigt, eps )`],
    [`voigt_contract( Cvoigt, eps, eta )`],
  )],
  [Operation:
   same constitutive maps as `contract`, but explicitly in Voigt notation.
   It makes the chosen storage convention visible in the form itself.],
  [#stack(
    dir: ttb,
    spacing: 0.18em,
    [`voigt_contract( Cvoigt, epsVoigt )`],
    [`voigt_contract( Cvoigt, epsVoigt, etaVoigt )`],
  )],
)
