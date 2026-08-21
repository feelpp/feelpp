// Mixed SB9 shell formulation used by qs_sb9_mixed_bending.cpp

#set page(
    paper: "a4",
    margin: (x: 18mm, y: 16mm),
)

#set text(size: 10.5pt)
#set heading(numbering: "1.")

= Mixed SB9 Shell Formulation on One-Layer Hexahedra

This note describes the mathematical formulation used in
`feelpp/quickstart/qs_sb9_mixed_bending.cpp`.
It gives the shell geometry, the mixed discrete unknowns, the SB9 kinematic
operators, the constitutive law, the stabilization terms, and the final mixed
weak form.

Only the formulation is described here. No implementation details, solver
choices, or numerical results are discussed.

= Geometric Setting

Let
$
Omega_h = union_(K in T_h) K
$
be a 3D hexahedral mesh of a shell-like body.
Each shell cell is a first-order hexahedron with exactly one element through the
thickness.

For each cell $K$, let the reference element be
$
hat(K) = [-1, 1]^3
$
with coordinates
$
(xi, eta, zeta).
$
The geometric map is
$
X_K : hat(K) -> K subset RR^3.
$

The reference coordinate $zeta in [-1,1]$ is the shell thickness coordinate.
The shell midsurface is the section $zeta = 0$.

This one-layer setting is essential.
The model is a shell model carried by a 3D hexahedral mesh, not a standard
3D solid discretization stacked through the thickness.
The dependence on $zeta$ already represents the through-thickness kinematics,
so refining the thickness direction by adding several element layers would
change the mechanical model instead of refining the same shell approximation.

= Cellwise Shell Geometry

For each cell $K$, the midsurface area is defined at the reference center by
$
A_(0,K) = 4 || X_(K,xi)(0,0,0) times X_(K,eta)(0,0,0) ||.
$

The shell thickness is the cellwise quantity
$
h_K = frac(V_K, A_(0,K)),
$
where $V_K = |K|$ is the cell volume.

The shell normal is
$
n_K =
  frac(
    X_(K,xi)(0,0,0) times X_(K,eta)(0,0,0),
    || X_(K,xi)(0,0,0) times X_(K,eta)(0,0,0) ||
  ),
$
with orientation chosen so that
$
n_K dot X_(K,zeta)(0,0,0) > 0.
$

The local orthonormal shell frame is
$
R_K = [ t_(1,K) quad t_(2,K) quad n_K ],
$
with
$
t_(1,K) = frac( X_(K,xi)(0,0,0), || X_(K,xi)(0,0,0) || ),
$
$
t_(2,K) = frac( n_K times t_(1,K), || n_K times t_(1,K) || ).
$

Define the center covariant basis
$
G_(0,K) = [ X_(K,xi)(0,0,0) quad X_(K,eta)(0,0,0) quad X_(K,zeta)(0,0,0) ]
$
and the center contravariant basis
$
G_(0,K)^* = G_(0,K)^(-T).
$
The corresponding metric tensor is
$
g_(0,K) = G_(0,K)^T G_(0,K).
$

The local shell Jacobian is
$
J_(0,K) = R_K^T G_(0,K).
$
We also use
$
hat(J)_(0,K) = J_(0,K)^(-T), quad G_(0,K)^(-1) = J_(0,K)^(-1).
$

These quantities separate the different shell mechanisms.
The frame $R_K$ identifies the two local tangential directions and the shell
normal.
The thickness $h_K$ measures the true local shell thickness of a distorted
hexahedron.
The Jacobian $J_(0,K)$ transfers reference derivatives into the local shell
frame, where membrane, bending, shear, and pinching have a direct mechanical
interpretation.

= Local Nodal Coordinates and Reduced SB9 Geometry

Let the eight cell vertices be $X_a in RR^3$, $a = 1,dots,8$, and let
$
X_c = frac(1,8) sum_(a=1)^8 X_a
$
be the cell center.
Their local shell-frame coordinates are
$
hat(X)_a = R_K^T ( X_a - X_c ) = [ hat(x)_a, hat(y)_a, hat(z)_a ]^T.
$

Define the center derivative matrix
$
B_xi = frac(1,8) mat(
    -1,  1,  1, -1, -1,  1,  1, -1;
    -1, -1,  1,  1, -1, -1,  1,  1;
    -1, -1, -1, -1,  1,  1,  1,  1
).
$
The cellwise coefficients $b_x^a$, $b_y^a$, $b_z^a$ are defined by
$
mat(
    b_x^1, b_x^2, b_x^3, b_x^4, b_x^5, b_x^6, b_x^7, b_x^8;
    b_y^1, b_y^2, b_y^3, b_y^4, b_y^5, b_y^6, b_y^7, b_y^8;
    b_z^1, b_z^2, b_z^3, b_z^4, b_z^5, b_z^6, b_z^7, b_z^8
) = hat(J)_(0,K) B_xi.
$

Next define the four Hallquist vectors
$
H = mat(
     1,  1,  1, -1;
     1, -1, -1,  1;
    -1, -1,  1, -1;
    -1,  1, -1,  1;
    -1, -1,  1,  1;
    -1,  1, -1, -1;
     1,  1,  1,  1;
     1, -1, -1, -1
).
$
For each $gamma = 1,dots,4$, define
$
h^(gamma) = sum_(a=1)^8 H_(a gamma) hat(X)_a
= [ h_1^(gamma), h_2^(gamma), h_3^(gamma) ]^T
$
and the corrected Hallquist scalars
$
v_(a gamma) =
frac(1,8) (
  H_(a gamma)
  - h_1^(gamma) b_x^a
  - h_2^(gamma) b_y^a
  - h_3^(gamma) b_z^a
),
quad a = 1,dots,8.
$

Using the local nodal in-plane coordinates
$
(hat(x)_a, hat(y)_a),
$
define the reduced SB9 Jacobians
$
J_a = mat(
    frac( hat(x)_2 - hat(x)_1 + hat(x)_6 - hat(x)_5, 4 ),
    frac( hat(y)_2 - hat(y)_1 + hat(y)_6 - hat(y)_5, 4 );
    frac( -hat(x)_1 - hat(x)_2 + hat(x)_3 + hat(x)_4 - hat(x)_5 - hat(x)_6 + hat(x)_7 + hat(x)_8, 8 ),
    frac( -hat(y)_1 - hat(y)_2 + hat(y)_3 + hat(y)_4 - hat(y)_5 - hat(y)_6 + hat(y)_7 + hat(y)_8, 8 )
),
$
$
J_b = mat(
    frac( -hat(x)_1 + hat(x)_2 + hat(x)_3 - hat(x)_4 - hat(x)_5 + hat(x)_6 + hat(x)_7 - hat(x)_8, 8 ),
    frac( -hat(y)_1 + hat(y)_2 + hat(y)_3 - hat(y)_4 - hat(y)_5 + hat(y)_6 + hat(y)_7 - hat(y)_8, 8 );
    frac( -hat(x)_2 + hat(x)_3 - hat(x)_6 + hat(x)_7, 4 ),
    frac( -hat(y)_2 + hat(y)_3 - hat(y)_6 + hat(y)_7, 4 )
),
$
$
J_c = mat(
    frac( hat(x)_3 - hat(x)_4 + hat(x)_7 - hat(x)_8, 4 ),
    frac( hat(y)_3 - hat(y)_4 + hat(y)_7 - hat(y)_8, 4 );
    frac( -hat(x)_1 - hat(x)_2 + hat(x)_3 + hat(x)_4 - hat(x)_5 - hat(x)_6 + hat(x)_7 + hat(x)_8, 8 ),
    frac( -hat(y)_1 - hat(y)_2 + hat(y)_3 + hat(y)_4 - hat(y)_5 - hat(y)_6 + hat(y)_7 + hat(y)_8, 8 )
),
$
$
J_d = mat(
    frac( -hat(x)_1 + hat(x)_2 + hat(x)_3 - hat(x)_4 - hat(x)_5 + hat(x)_6 + hat(x)_7 - hat(x)_8, 8 ),
    frac( -hat(y)_1 + hat(y)_2 + hat(y)_3 - hat(y)_4 - hat(y)_5 + hat(y)_6 + hat(y)_7 - hat(y)_8, 8 );
    frac( -hat(x)_1 + hat(x)_4 - hat(x)_5 + hat(x)_8, 4 ),
    frac( -hat(y)_1 + hat(y)_4 - hat(y)_5 + hat(y)_8, 4 )
).
$

Their reduced inverses are embedded as
$
bar(J)_s^(-1) = "diag"( J_s^(-1), frac(2, h_K) ),
quad s in { a, b, c, d }.
$

This reduced geometry is the core of the SB9 construction.
The vectors $b_x$, $b_y$, $b_z$ carry the center derivative information of the
Q1 hexahedron in the shell frame.
The corrected Hallquist coefficients $v_(a gamma)$ remove affine content from
the raw Hallquist patterns.
The reduced Jacobians $J_a$, $J_b$, $J_c$, $J_d$ are then used to build the
assumed-shear and stabilization terms that distinguish the shell formulation
from a plain small-strain 3D solid law.

= Discrete Spaces and Local Unknowns

The mixed discrete spaces are
$
U_h = [Q_1]^3,
quad
A_h = P_0^("disc"),
quad
X_h = U_h times A_h.
$

The unknown is
$
(u_h, alpha_h) in X_h,
$
with test pair
$
(v_h, beta_h) in X_h.
$

On each cell $K$, the displacement field is
$
u_h|_K(X) = sum_(a=1)^8 N_a(X) d_a,
quad d_a in RR^3,
$
and its local shell-frame nodal components are
$
q_a = R_K^T d_a = [ u_a, v_a, w_a ]^T.
$

The scalar field is cellwise constant:
$
alpha_h|_K = alpha_K,
quad beta_h|_K = beta_K.
$

This choice is mechanical, not merely algebraic.
The mixed scalar is an internal thickness-correction mode attached to the cell,
so it is represented by one scalar degree of freedom per shell cell rather than
by additional nodal displacement components.

= Ordering and Symmetric Storage

Two ordering conventions are fundamental to the mixed SB9 formulation.

First, the local displacement degrees of freedom are grouped componentwise in
the shell frame:
$
underline(q)_K =
[ u_1, dots, u_8, v_1, dots, v_8, w_1, dots, w_8 ]^T.
$
The SB9 operators act on this component-blocked ordering.
This matters because the local shell operators do not manipulate a node-major
vector of triples; they multiply separately the tangential and normal nodal
components.

Second, the symmetric strain is stored in the compressed order
$
(11, 12, 13, 22, 23, 33),
$
namely
$
underline(epsilon) =
[ epsilon_11,
  sqrt(2) epsilon_12,
  sqrt(2) epsilon_13,
  epsilon_22,
  sqrt(2) epsilon_23,
  epsilon_33 ]^T.
$

This ordering is mechanically critical.
The membrane terms occupy the $(11,12,22)$ slots, the transverse shear terms
occupy the $(13,23)$ slots, and the pinching and mixed $"W9"$ correction occupy
the $(33)$ slot.
If the shell strain pieces are inserted in a different compressed order, the
elastic operator couples the wrong mechanisms.
That ordering issue is not cosmetic; it changes the meaning of the element
stiffness.

= Linear Elasticity in Symmetric Storage

The isotropic linear elastic law is
$
underline(sigma) = C underline(epsilon),
$
with Lam'e coefficients
$
lambda = frac( E nu, (1 + nu)(1 - 2 nu) ),
quad
mu = frac( E, 2(1 + nu) ),
$
and stiffness matrix
$
C = mat(
    lambda + 2 mu, 0, 0, lambda, 0, lambda;
    0, 2 mu, 0, 0, 0, 0;
    0, 0, 2 mu, 0, 0, 0;
    lambda, 0, 0, lambda + 2 mu, 0, lambda;
    0, 0, 0, 0, 2 mu, 0;
    lambda, 0, 0, lambda, 0, lambda + 2 mu
).
$

The elastic energy density between two stored symmetric vectors is
$
underline(epsilon)^T C underline(eta).
$

The square-root factors on the shear slots ensure that this Euclidean product
matches the double contraction of the underlying symmetric tensors.

= SB9 Kinematic Operators

For each cell $K$, the shell strain is split into membrane-bending, pinching,
transverse shear, and one mixed cell-scalar correction.

== Membrane and Bending Part

Define the center membrane operator $B_m^0$ by
$
underline(epsilon)_(m,0)(u_h)
= sum_(a=1)^8 mat(
    b_x^a u_a;
    frac( b_y^a u_a + b_x^a v_a, sqrt(2) );
    0;
    b_y^a v_a;
    0;
    0
).
$

Define the curvature operator $B_b^0$ from
$
tilde(b)_x^a = v_(a,1) hat(J)_(0,K)_(1 2) + v_(a,2) hat(J)_(0,K)_(1 1),
$
$
tilde(b)_y^a = v_(a,1) hat(J)_(0,K)_(2 2) + v_(a,2) hat(J)_(0,K)_(2 1),
$
so that
$
underline(epsilon)_(b,0)(u_h)
= sum_(a=1)^8 mat(
    tilde(b)_x^a u_a;
    frac( tilde(b)_y^a u_a + tilde(b)_x^a v_a, sqrt(2) );
    0;
    tilde(b)_y^a v_a;
    0;
    0
).
$

The membrane-bending contribution at a thickness point $zeta$ is
$
underline(epsilon)_("mb")(u_h; zeta)
= underline(epsilon)_(m,0)(u_h)
+ zeta underline(epsilon)_(b,0)(u_h).
$

Mechanically, this is the standard shell split between midsurface stretching
and bending.
The first term is the membrane strain of the midsurface.
The factor $zeta$ in the second term produces the linear variation of
tangential strain through the thickness that characterizes bending.

== Pinching Part

The shell-normal pinching operator is
$
underline(epsilon)_("pc")(u_h)
= sum_(a=1)^8 mat(
    0;
    0;
    0;
    0;
    0;
    b_z^a w_a
).
$

This term is the displacement-driven part of the shell-normal strain.
It accounts for thickness stretching or pinching created by the normal
displacement components of the shell kinematics.

== Mixed Cellwise W9 Correction

The additional mixed scalar correction acts only on the shell-normal strain:
$
underline(epsilon)_("w9")(alpha_h; zeta)
= s_alpha mat(
    0;
    0;
    0;
    0;
    0;
    frac( -4 zeta alpha_K, h_K )
),
$
where $s_alpha > 0$ is the scalar scaling parameter.

The scalar $alpha_K$ is therefore an internal cellwise correction attached to
the shell-normal channel.
The factor $frac(-4 zeta, h_K)$ makes the associated mode odd across the thickness
and zero at the midsurface, which is exactly the type of thickness correction
that a one-layer displacement field cannot represent by itself.

== Transverse Shear Part

The transverse shear shape factor is
$
omega(zeta) = s_("sh") (1 - zeta^2),
$
where $s_("sh") > 0$ is the shell shear factor.

For $kappa in {0,1,2}$ define the two-row operator
$
B_c^(kappa)(u_h)
= mat(
    Gamma_(13)^(kappa)(u_h);
    Gamma_(23)^(kappa)(u_h)
),
$
with
$
Gamma_(13)^(kappa)(u_h)
= sum_(a=1)^8 (
    A_(11)^(kappa,a) u_a
  + A_(12)^(kappa,a) v_a
  + A_(13)^(kappa,a) w_a
),
$
$
Gamma_(23)^(kappa)(u_h)
= sum_(a=1)^8 (
    A_(21)^(kappa,a) u_a
  + A_(22)^(kappa,a) v_a
  + A_(23)^(kappa,a) w_a
).
$

The coefficients are
$
A_(11)^(kappa,a) = f_(11)^(b_z,kappa) b_z^a + f_(11)^(g_1,kappa) v_(a,1),
$
$
A_(12)^(kappa,a) = f_(12)^(b_z,kappa) b_z^a + f_(12)^(g_1,kappa) v_(a,1),
$
$
A_(13)^(kappa,a) = f_(13)^(b_x,kappa) b_x^a + f_(13)^(b_y,kappa) b_y^a + f_(13)^(g_3,kappa) v_(a,3),
$
$
A_(21)^(kappa,a) = f_(21)^(b_z,kappa) b_z^a + f_(21)^(g_1,kappa) v_(a,1) + f_(21)^(g_2,kappa) v_(a,2),
$
$
A_(22)^(kappa,a) = f_(22)^(b_z,kappa) b_z^a + f_(22)^(g_1,kappa) v_(a,1) + f_(22)^(g_2,kappa) v_(a,2),
$
$
A_(23)^(kappa,a) = f_(23)^(b_x,kappa) b_x^a + f_(23)^(b_y,kappa) b_y^a + f_(23)^(g_3,kappa) v_(a,3).
$

For $kappa = 0$, the implemented coefficients are
$
f_(11)^(b_z,0) = frac( hat(J)_(0,K)_(1 1) ( J_a(1,1) + J_c(1,1) ), 2 ),
$
$
f_(11)^(g_1,0) =
hat(J)_(0,K)_(1 1)
frac( bar(J)_c^(-1)(3,3) J_c(1,1) - bar(J)_a^(-1)(3,3) J_a(1,1), 2 ),
$
$
f_(12)^(b_z,0) = frac( hat(J)_(0,K)_(1 1) ( J_a(1,2) + J_c(1,2) ), 2 ),
$
$
f_(12)^(g_1,0) =
hat(J)_(0,K)_(1 1)
frac( bar(J)_c^(-1)(3,3) J_c(1,2) - bar(J)_a^(-1)(3,3) J_a(1,2), 2 ),
$
$
f_(13)^(b_x,0) = f_(11)^(b_z,0), quad
f_(13)^(b_y,0) = f_(12)^(b_z,0),
$
$
f_(13)^(g_3,0) =
hat(J)_(0,K)_(1 1)
frac(
(
  bar(J)_c^(-1)(1,1) J_c(1,1) + bar(J)_c^(-1)(2,1) J_c(1,2)
  - bar(J)_a^(-1)(1,1) J_a(1,1) - bar(J)_a^(-1)(2,1) J_a(1,2)
), 2 ),
$
$
f_(21)^(b_z,0) =
frac( hat(J)_(0,K)_(2 1) ( J_a(1,1) + J_c(1,1) ), 2 )
+ frac( hat(J)_(0,K)_(2 2) ( J_b(2,1) + J_d(2,1) ), 2 ),
$
$
f_(21)^(g_1,0) =
hat(J)_(0,K)_(2 1)
frac( bar(J)_c^(-1)(3,3) J_c(1,1) - bar(J)_a^(-1)(3,3) J_a(1,1), 2 ),
$
$
f_(21)^(g_2,0) =
hat(J)_(0,K)_(2 2)
frac( bar(J)_b^(-1)(3,3) J_b(2,1) - bar(J)_d^(-1)(3,3) J_d(2,1), 2 ),
$
$
f_(22)^(b_z,0) =
frac( hat(J)_(0,K)_(2 1) ( J_a(1,2) + J_c(1,2) ), 2 )
+ frac( hat(J)_(0,K)_(2 2) ( J_b(2,2) + J_d(2,2) ), 2 ),
$
$
f_(22)^(g_1,0) =
hat(J)_(0,K)_(2 1)
frac( bar(J)_c^(-1)(3,3) J_c(1,2) - bar(J)_a^(-1)(3,3) J_a(1,2), 2 ),
$
$
f_(22)^(g_2,0) =
hat(J)_(0,K)_(2 2)
frac( bar(J)_b^(-1)(3,3) J_b(2,2) - bar(J)_d^(-1)(3,3) J_d(2,2), 2 ),
$
$
f_(23)^(b_x,0) = f_(21)^(b_z,0), quad
f_(23)^(b_y,0) = f_(22)^(b_z,0),
$
$
f_(23)^(g_3,0) =
frac(
  hat(J)_(0,K)_(2 1)
  (
    bar(J)_c^(-1)(1,1) J_c(1,1) + bar(J)_c^(-1)(2,1) J_c(1,2)
    - bar(J)_a^(-1)(1,1) J_a(1,1) - bar(J)_a^(-1)(2,1) J_a(1,2)
  ),
  2
)
+ frac(
  hat(J)_(0,K)_(2 2)
  (
    bar(J)_b^(-1)(1,2) J_b(2,1) + bar(J)_b^(-1)(2,2) J_b(2,2)
    - bar(J)_d^(-1)(1,2) J_d(2,1) - bar(J)_d^(-1)(2,2) J_d(2,2)
  ),
  2
).
$

For $kappa = 1$, the coefficients are
$
f_(11)^(b_z,1) = frac( hat(J)_(0,K)_(1 1) ( J_c(1,1) - J_a(1,1) ), 2 ),
$
$
f_(11)^(g_1,1) =
hat(J)_(0,K)_(1 1)
frac( bar(J)_c^(-1)(3,3) J_c(1,1) + bar(J)_a^(-1)(3,3) J_a(1,1), 2 ),
$
$
f_(12)^(b_z,1) = frac( hat(J)_(0,K)_(1 1) ( J_c(1,2) - J_a(1,2) ), 2 ),
$
$
f_(12)^(g_1,1) =
hat(J)_(0,K)_(1 1)
frac( bar(J)_c^(-1)(3,3) J_c(1,2) + bar(J)_a^(-1)(3,3) J_a(1,2), 2 ),
$
$
f_(13)^(b_x,1) = f_(11)^(b_z,1), quad
f_(13)^(b_y,1) = f_(12)^(b_z,1),
$
$
f_(13)^(g_3,1) =
hat(J)_(0,K)_(1 1)
frac(
(
  bar(J)_c^(-1)(1,1) J_c(1,1) + bar(J)_c^(-1)(2,1) J_c(1,2)
  + bar(J)_a^(-1)(1,1) J_a(1,1) + bar(J)_a^(-1)(2,1) J_a(1,2)
), 2 ),
$
$
f_(21)^(b_z,1) = frac( hat(J)_(0,K)_(2 1) ( J_c(1,1) - J_a(1,1) ), 2 ),
$
$
f_(21)^(g_1,1) =
hat(J)_(0,K)_(2 1)
frac( bar(J)_c^(-1)(3,3) J_c(1,1) + bar(J)_a^(-1)(3,3) J_a(1,1), 2 ),
$
$
f_(22)^(b_z,1) = frac( hat(J)_(0,K)_(2 1) ( J_c(1,2) - J_a(1,2) ), 2 ),
$
$
f_(22)^(g_1,1) =
hat(J)_(0,K)_(2 1)
frac( bar(J)_c^(-1)(3,3) J_c(1,2) + bar(J)_a^(-1)(3,3) J_a(1,2), 2 ),
$
$
f_(23)^(b_x,1) = f_(21)^(b_z,1), quad
f_(23)^(b_y,1) = f_(22)^(b_z,1),
$
$
f_(23)^(g_3,1) =
hat(J)_(0,K)_(2 1)
frac(
(
  bar(J)_c^(-1)(1,1) J_c(1,1) + bar(J)_c^(-1)(2,1) J_c(1,2)
  + bar(J)_a^(-1)(1,1) J_a(1,1) + bar(J)_a^(-1)(2,1) J_a(1,2)
), 2 ).
$

For $kappa = 2$, the nonzero coefficients are
$
f_(21)^(b_z,2) = frac( hat(J)_(0,K)_(2 2) ( J_b(2,1) - J_d(2,1) ), 2 ),
$
$
f_(21)^(g_2,2) =
hat(J)_(0,K)_(2 2)
frac( bar(J)_b^(-1)(3,3) J_b(2,1) + bar(J)_d^(-1)(3,3) J_d(2,1), 2 ),
$
$
f_(22)^(b_z,2) = frac( hat(J)_(0,K)_(2 2) ( J_b(2,2) - J_d(2,2) ), 2 ),
$
$
f_(22)^(g_2,2) =
hat(J)_(0,K)_(2 2)
frac( bar(J)_b^(-1)(3,3) J_b(2,2) + bar(J)_d^(-1)(3,3) J_d(2,2), 2 ),
$
$
f_(23)^(b_x,2) = f_(21)^(b_z,2), quad
f_(23)^(b_y,2) = f_(22)^(b_z,2),
$
$
f_(23)^(g_3,2) =
hat(J)_(0,K)_(2 2)
frac(
(
  bar(J)_b^(-1)(1,2) J_b(2,1) + bar(J)_b^(-1)(2,2) J_b(2,2)
  + bar(J)_d^(-1)(1,2) J_d(2,1) + bar(J)_d^(-1)(2,2) J_d(2,2)
), 2 ).
$

The transverse shear contribution entering the shell strain is
$
underline(epsilon)_("sh")(u_h; zeta)
= mat(
    0;
    0;
    frac( omega(zeta) Gamma_(13)^(0)(u_h), sqrt(2) );
    0;
    frac( omega(zeta) Gamma_(23)^(0)(u_h), sqrt(2) );
    0
).
$

The factor $omega(zeta) = s_("sh")(1-zeta^2)$ is a parabolic shear profile.
It vanishes on the top and bottom faces and is largest at the midsurface.
This is consistent with the shell-mechanics requirement that transverse shear
stresses vanish on the outer faces.
The auxiliary operators $B_c^(1)$ and $B_c^(2)$ are used only for
stabilization; they do not define an additional physical strain component.

== Additional Stabilization Operators

The implemented SB9 stabilization operators are
$
B_(s,1)(u_h) = sum_(a=1)^8 v_(a,1) w_a,
$
$
B_(s,2)(u_h) = sum_(a=1)^8 v_(a,2) w_a,
$
$
B_(s,3)(u_h) = mat(
    sum_(a=1)^8 v_(a,3) u_a;
    sum_(a=1)^8 v_(a,3) v_a
),
$
$
B_(s,4)(u_h) = mat(
    sum_(a=1)^8 v_(a,4) u_a;
    sum_(a=1)^8 v_(a,4) v_a;
    sum_(a=1)^8 v_(a,4) w_a
).
$

These operators are stabilization probes, not additional constitutive strains.
They detect residual Hallquist-type modes that would otherwise remain too soft.
Their role is to remove spurious low-energy deformation patterns while leaving
the intended membrane, bending, pinching, and shear mechanisms in charge of the
physical shell response.

= Final Shell Strain Split

For the displacement field, the SB9 shell strain on a cell $K$ is
$
underline(epsilon)^u_K(u_h; zeta)
= underline(epsilon)_("mb")(u_h; zeta)
+ underline(epsilon)_("pc")(u_h)
+ underline(epsilon)_("sh")(u_h; zeta).
$

For the mixed cellwise scalar,
$
underline(epsilon)^alpha_K(alpha_h; zeta)
= underline(epsilon)_("w9")(alpha_h; zeta).
$

Hence the total mixed strain is
$
underline(epsilon)_K(u_h, alpha_h; zeta)
= underline(epsilon)^u_K(u_h; zeta)
+ underline(epsilon)^alpha_K(alpha_h; zeta).
$

= Mixed Bilinear Form

Find
$
(u_h, alpha_h) in X_h
$
such that for all
$
(v_h, beta_h) in X_h,
$
$
a( (u_h, alpha_h), (v_h, beta_h) ) = ell( (v_h, beta_h) ).
$

The elastic part is split into four blocks:
$
a_(u,u)(u_h, v_h)
= sum_(K in T_h) integral_K
  underline(epsilon)^u_K(u_h; zeta)^T
  C
  underline(epsilon)^u_K(v_h; zeta),
$
$
a_(u,alpha)(alpha_h, v_h)
= sum_(K in T_h) integral_K
  underline(epsilon)^alpha_K(alpha_h; zeta)^T
  C
  underline(epsilon)^u_K(v_h; zeta),
$
$
a_(alpha,u)(u_h, beta_h)
= sum_(K in T_h) integral_K
  underline(epsilon)^u_K(u_h; zeta)^T
  C
  underline(epsilon)^alpha_K(beta_h; zeta),
$
$
a_(alpha,alpha)(alpha_h, beta_h)
= sum_(K in T_h) integral_K
  underline(epsilon)^alpha_K(alpha_h; zeta)^T
  C
  underline(epsilon)^alpha_K(beta_h; zeta).
$

The Hallquist transverse-shear stabilization is
$
a_c(u_h, v_h)
= c_c sum_(K in T_h) integral_K (
    B_c^(1)(u_h) dot B_c^(1)(v_h)
  + B_c^(2)(u_h) dot B_c^(2)(v_h)
),
$
with
$
c_c = frac( s_c mu 5, 18 ),
$
where $s_c > 0$ is the Hallquist shear stabilization scale.

This contribution acts only on the displacement block.
Its purpose is to control the special shear modes extracted by $B_c^(1)$ and
$B_c^(2)$ without replacing the physical shear strain carried by
$underline(epsilon)_("sh")$.

The additional SB9 stabilization uses
$
d_N = lambda + 2 mu
$
and the entries of $J_(0,K)^(-1)$:
$
d_(s,1,K) = d_(s,2,K) = frac( s_b d_N, 3 ) ( J_(0,K)^(-1)(3,3) )^2,
$
$
d_(s,3,x,K) = frac( s_b d_N, 3 ) ( J_(0,K)^(-1)(1,1) )^2,
$
$
d_(s,3,y,K) = frac( s_b d_N, 3 ) (
    ( J_(0,K)^(-1)(2,1) )^2 + ( J_(0,K)^(-1)(2,2) )^2
),
$
$
d_(s,4,x,K) = frac( s_b d_N 10^(-4), 9 ) ( J_(0,K)^(-1)(1,1) )^2,
$
$
d_(s,4,y,K) = frac( s_b d_N 10^(-4), 9 ) (
    ( J_(0,K)^(-1)(2,1) )^2 + ( J_(0,K)^(-1)(2,2) )^2
),
$
$
d_(s,4,z,K) = frac( s_b d_N, 9 ) ( J_(0,K)^(-1)(3,3) )^2,
$
where $s_b > 0$ is the SB9 stabilization scale.

The corresponding bilinear form is
$
a_s(u_h, v_h)
= sum_(K in T_h) integral_K A_(s,K)(u_h, v_h),
$
with
$
A_(s,K)(u_h, v_h)
= d_(s,1,K) B_(s,1)(u_h) B_(s,1)(v_h)
+ d_(s,2,K) B_(s,2)(u_h) B_(s,2)(v_h)
+ d_(s,3,x,K) B_(s,3)(u_h)_1 B_(s,3)(v_h)_1
+ d_(s,3,y,K) B_(s,3)(u_h)_2 B_(s,3)(v_h)_2
$
$
+ d_(s,4,x,K) B_(s,4)(u_h)_1 B_(s,4)(v_h)_1
+ d_(s,4,y,K) B_(s,4)(u_h)_2 B_(s,4)(v_h)_2
+ d_(s,4,z,K) B_(s,4)(u_h)_3 B_(s,4)(v_h)_3.
$

The stabilization coefficients depend on the inverse local shell Jacobian, so
the stabilization scales automatically with the local anisotropy of the shell
cell.
The three groups have distinct roles:

- $B_(s,1)$ and $B_(s,2)$ act on scalar shell-normal residual modes.
- $B_(s,3)$ acts on tangential residual modes.
- $B_(s,4)$ acts on the remaining fully coupled residual mode.

The final mixed bilinear form is
$
a = a_(u,u) + a_(u,alpha) + a_(alpha,u) + a_(alpha,alpha) + a_c + a_s.
$

The structure of the mixed form mirrors the mechanics.
The $("u,u")$ block carries the shell response of the displacement field itself.
The off-diagonal blocks couple the internal scalar correction to the
displacement shell strain through the same elastic law.
The $("alpha,alpha")$ block gives the internal thickness mode its own energy.

= Loads

Only the displacement field is directly loaded.
The scalar field has no direct right-hand side.

The body-force contribution is
$
ell_f(v_h) = sum_(K in T_h) integral_K f dot v_h.
$

For each pressure boundary $Gamma_p$, with outward unit normal $N$,
the pressure load is
$
ell_p(v_h) = integral_(Gamma_p) ( -p N ) dot v_h.
$
Thus a positive pressure acts inward through the traction $-p N$.

For each prescribed surface traction $t$ on $Gamma_t$,
$
ell_t(v_h) = integral_(Gamma_t) t dot v_h.
$

For a total resultant force $F_G in RR^3$ prescribed on a boundary piece
$Gamma_G$, the formulation uses the uniform equivalent traction
$
t_G = frac( F_G, |Gamma_G| ),
$
so that
$
ell_G(v_h) = integral_(Gamma_G) t_G dot v_h.
$

The full right-hand side is
$
ell( (v_h, beta_h) ) = ell_f(v_h) + ell_p(v_h) + ell_t(v_h) + ell_G(v_h).
$

All external work is applied to the displacement field.
This is consistent with the interpretation of $alpha_h$ as an internal shell
correction variable rather than a directly loaded primary field.

= Essential Boundary Conditions

The displacement field is clamped on the prescribed clamp boundary
$Gamma_C$:
$
u_h = 0 quad "on" Gamma_C.
$

The scalar field $alpha_h$ carries no Dirichlet condition.

At the discrete level, the essential condition is imposed strongly by
elimination on the displacement row block of the mixed system.

This is the correct mixed treatment.
The clamp prescribes the physical displacement and therefore acts on the whole
displacement test row block.
No essential condition is imposed on $alpha_h$, because the cellwise scalar is
an internal shell mode rather than an independent boundary-controlled field.
