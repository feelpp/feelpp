/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page Spaces Function Spaces



\li \b Previous: \ref Mesh
\li \b Next: \ref Integrals

<hr>

\tableofcontents


\section QRFunctionSpace Function Space

We now turn to the next crucial mathematical ingredient: the function space,
whose definition depends on \f$\Omega_h\f$ --- or more precisely its partitioning
\f$\calT_h\f$ --- and the choice of basis function. Function spaces in \feel
follow the same definition, see listing~\ref fspace1, and \feel provides
support for continuous and discontinuous Galerkin methods and in particular
approximations in \f$L^2\f$, \f$H^1\f$-conforming and \f$H^1\f$-nonconforming, \f$H^2\f$,
\f$H(\mathrm{div})\f$ and \f$H(\mathrm{curl})\f$\footnote{At the time of writing, \f$H^2\f$,
\f$H(\mathrm{div})\f$ and \f$H(\mathrm{curl})\f$ approximations are in experimental
support.

\anchor fspace1
\code
 // space of continuous piecewise
 // \f$\P_3\f$ functions defined on a mesh
 // of order 2 triangles in 3D
 FunctionSpace<Mesh<Simplex<2,2,3>,
               bases<Lagrange<3> > > Xh;
\endcode


The \c FunctionSpace class

 -  constructs the table of degrees of freedom which maps local (elementwise) degrees of
  freedom to the global ones with respect to the geometrical entities,
 -  embeds the definition of the elements of the function space allowing for a
  tight coupling between the elements and their function spaces,
 -  stores an interpolation data structure (\eg region tree) for rapid
  localisation of point sets (determining in which element they reside).


We introduce the following spaces
\f[
  \begin{aligned}
    \mathbb{W}_h &= \{v_h \in L^2(\Omega_h): \ \forall K \in \mathcal{T}_h, v_h|_K
    \in \mathbb{P}_K\},\\
    \mathbb{V}_h &= \mathbb{W}_h \cap C^0(\Omega_h)= \{ v_h \in \mathbb{W}_h: \ \forall F \in
    \mathcal{F}^i_h\ \jump{v_h}_F = 0\}\\
    \mathbb{H}_h &= \mathbb{W}_h \cap C^1(\Omega_h)= \{ v_h \in \mathbb{W}_h: \ \forall F \in
    \mathcal{F}^i_h\ \jump{v_h}_F = \jump{\nabla v_h}_F = 0\}\\
    \CR_h &= \{ v_h \in L^2(\Omega_h):\ \forall K \in \calT_h, v_h|_K \in
    \P_1; \forall F \in \calF^i_h\ \int_F \jump{v_h} = 0 \}\\
    \RaTu_h &= \{ v_h \in L^2(\Omega_h):\ \forall K \in \calT_h, v_h|_K \in
    \Span{1,x,y,x^2-y^2}; \forall F \in \calF^i_h\ \int_F \jump{v_h} = 0 \}\\
    \RT_h&=\{\bm{v}_h \in [L^2(\Omega_h)]^d:\ \forall K \in \calT_h, v_h|_K \in
    \RT_k; \forall F \in \calF^i_h\ \jump{\bm{v}_h \cdot \normal}_F = 0 \}\\
    \N_h&=\{\bm{v}_h \in [L^2(\Omega_h)]^d:\ \forall K \in \calT_h, v_h|_K \in
    \N_k; \forall F \in \calF^i_h\ \jump{\bm{v}_h \times \normal}_F = 0 \}
  \end{aligned}
\f]
where \f$\RT_k\f$ and \f$\N_k\f$ are respectively the Raviart-Thomas and N&eacute;d&eacute;lec finite
elements of degree \f$k\f$.

The Legrendre and Dubiner basis yield implicitely discontinuous
approximations, the Legendre and Dubiner boundary adapted basis,
see~\cite MR1696933, were designed to handle continuous approximations
whereas the Lagrange basis can yield either discontinuous or continuous
(default behavior) approximations.  \f$\RT_h\f$ and \f$\N_h\f$ are implicitely spaces
of vectorial functions \f$\bm{f}\f$ \st \f$\bm{f}: \Omega_h \subset \R^d \mapsto
\R^d\f$. As to the other basis functions, \ie Lagrange, Legrendre, Dubiner,
etc., they are parametrized by their values namely \c Scalar ,
\c Vectorial or \c Matricial.  Note that
\c FunctionSpace handles also products of function spaces.  This is
very powerful to describe complex multiphysics problems when coupled with
operators, functionals and forms described in the next section. Extracting
subspaces or component spaces are part of the interface.

\code
// continuous piecewise P3
// approximations
FunctionSpace<Mesh<Simplex<2> >,
  bases<Lagrange<3,Scalar,
                 Continuous> > > P3ch;
// discontinuous piecewise P3
// approximations
FunctionSpace<Mesh<Simplex<2> >,
  bases<Lagrange<3,Scalar,
                 Discontinuous> > > P3dh;
// mixed (P2 vectorial, P1 scalar,
// P1 Scalar) approximation
FunctionSpace<Mesh<Simplex<2>>,
 bases<Lagrange<2,Vectorial>,
       Lagrange<1,Scalar>,
       Lagrange<1,Scalar>>> P2P1P1;
\endcode

The most important feature in \c FunctionSpace is that it embeds the
definition of element which allows for the strict definition of an \c
Element of a \c FunctionSpace and thus ensures the correctness of the
code.  An element has its representation as a vector --- also in the
case of product of multiple spaces. --- The vector representation is
parametrized by one of the linear algebra backends. Other supported
operations are interpolation and extraction of components --- be it a
product of function spaces element or a vectorial/matricial element:

\code
FunctionSpace<Mesh<Simplex<2> >,
  bases<Lagrange<3,Scalar, Continuous> > > P3ch;
// get an element from P3ch
auto u = P3ch.element();
FunctionSpace<Mesh<Simplex<2> >,
 bases<Lagrange<2,Vectorial>, Lagrange<1,Scalar>,
       Lagrange<1,Scalar> > > P2P1P1;
auto U = P2P1P1.element();
// Views: changing a view changes U and vice versa
// view on element associated to P2
auto u = U.element<0>();
// extract view of first component
auto ux = u.comp(X);
// view on element associated to 1st P1
auto p = U.element<1>();
// view on element associated to 2nd P1
auto q = U.element<2>();
\endcode

Finally Feel++ provides the Lagrange, \f$\Iclag, \Idlag\f$, Crouzeix-Raviart, \f$\Icr\f$,
Raviart-Thomas, \f$\Irt\f$ and N&eacute;d&eacute;lec, \f$\In\f$ global interpolation operators.
In abstract form, they read
\f[
  \calI : \X \ni v \mapsto \sum_{i=1}^{\opdim\X} \ell_i(v) \phi_i
\f]
where \f$\X\f$ is the infinite dimensional space, \f$(\ell_i)_{i=1,...,\opdim\X}\f$ are
the linear forms and \f$(\phi_i)_{i=1...\opdim\X}\f$ the basis function associated
with the various approximations.

\section QRFunctionSpaceFn Function Space helper functions

 - \c Pch<N>(mesh) generates \f$\Pch[N](\Omega_h)\f$
 - \c Pchv<N>(mesh) generates \f$[\Pch[N](\Omega_h)]^d\f$
 - \c THch<N>(mesh) generates \f$[\Pch[N](\Omega_h)]^d\times \Pch[N](\Omega_h)\f$

<a href="#" class="top">top</a>
<hr>



*/
}
