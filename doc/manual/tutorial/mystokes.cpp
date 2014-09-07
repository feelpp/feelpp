/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>

  Date 2013-02-19

  Copyright (C) 2013 Université de Strasbourg

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <feel/feel.hpp>
using namespace Feel;

/**
 \page StokesTutorial Stokes Tutorial
\author Feel++ Consortium
\date 2013-02-19

\tableofcontents
<br>
<hr>
<br>

\section Stokes_Theory Theory
Let solve the stokes equation considering a Poiseuille flow profile. <br>
We have the following system of equations,
<br><center>\f$
\left\{
\begin{aligned}
-\mu\Delta \bf u + \nabla p & =  \bf f & \text{on}\; \Omega \;, \\
                \nabla\cdot\bf u & =  0 & \text{on}\; \Omega \;, \\
                \bf u  & =  g & \text{on}\; \Gamma \;, \\
\end{aligned}
\right.
\f$</center><br>
where \f$u\in [H_g^1(\Omega)]^d\f$ denotes the flow speed, \f$p\in [L_0^2(\Omega)]\f$ the fluid pressure, \f$\mu\f$ the fluid viscosity.<br>
The last boundary condition expresses a null pressure fixed on the outlet.<br>
The Poiseuille profile on the boundary is,
<br><center>\f$
g(x,y)=\left(
\begin{aligned}
 y(1-y) \\
 0      \\
\end{aligned}
\right)
\f$</center><br>
The method used to obtain the strong formulation is closed to the one used for the laplacian (see section \ref Laplacian ).<br>
We multiply the first equation by a test function \f$v\in H^1(\Omega)\f$ and we integrate on the domain \f$\Omega\f$,
<br><center>\f$
\begin{aligned}
 \left(
\int_\Omega \mu \nabla \mathbf u : \nabla \mathbf v
-\int_{\partial\Omega} \frac{\partial \mathbf u}{\partial \mathbf n} \cdot \mathbf v
\right)
+\int_\Omega ( \nabla\cdot(p \mathbf v) - \mathbf p \nabla\cdot v )
=\int_\Omega \mathbf f \cdot \mathbf v \;.
\end{aligned}
\f$</center><br>
where \f$n\f$ denotes a normal vector on the boundary.<br>
The divergence theorem (or Gauss's theorem) gives,
<br><center>\f$
\begin{aligned}
\int_\Omega \nabla\cdot(p \mathbf v) = \int_{\partial\Omega} p \mathbf v\cdot \mathbf n \;.
\end{aligned}
\f$</center><br>
We have to add a consistency terms to the equation to guaranty the symmetry of the bilinear form.<br>
This term is provided by the second equation. We multiply this equation by a test function \f$q\in L_2(\Omega)\f$ and we integrate on the domain \f$\Omega\f$,
<br><center>\f$
\begin{aligned}
\int_{\Omega} \nabla\cdot\mathbf u q = 0 \;,
\end{aligned}
\f$</center><br>
Finally, we deduce from the equations and after rearranging the integrals the variationnal formulation,
<br><center>\f$
\begin{aligned}
\int_\Omega \mu \nabla \mathbf u :\nabla \mathbf v
+\int_\Omega \left( \nabla\cdot\mathbf u q - p \nabla\cdot\mathbf v \right)
+
    \int_{\partial\Omega} \left( p \mathbf n -
   \frac{\partial \mathbf u}{\partial \mathbf n} \right)
     \cdot \mathbf v
=\int_\Omega \mathbf f \cdot \mathbf v
\end{aligned}
\f$</center><br>
Let us assume now that \f$(\mathbf v,q) \in [H_0^1(\Omega)]^d \times L_0^2(\Omega)\f$, the variationnal formulation leads to:
Find \f$(\mathbf u,p)\in [H_g^1(\Omega)]^d\times L_0^2(\Omega) \f$ such that for all \f$(\mathbf v,q) \in [H_0^1(\Omega)]^d \times L_0^2(\Omega)\f$
<br><center>\f$
\begin{aligned}
\int_\Omega \mu \nabla \mathbf u :\nabla \mathbf v
+\int_\Omega \left( \nabla\cdot\mathbf u q - p \nabla\cdot\mathbf v \right)
=\int_\Omega \mathbf f \cdot \mathbf v
\end{aligned}
\f$</center><br>
Or equivalently:
<br><center>\f$
\begin{aligned}
  a((\mathbf u,p),(\mathbf v,q)) = l((\mathbf v,q))
\end{aligned}
\f$</center><br>
where \f$a\f$ is a bilinear form, continuous, coercive and where \f$l\f$ is a linear form.

\section Stokes_Implementation Implementation
Let's see the \feel code corresponding to this mathematical statement (source \c "doc/manual/tutorial/mystokes.cpp").<br>
We suppose for this example the viscosity \f$\mu=1\f$ and \f$\mathbf f = 0\f$.

The procedure to create the mesh is very simple.<br>
You have to provide to the command line (or via the cfg file) the gmsh.filename option.<br>
You can provide a \c .geo or a \c .msh file (created via gmsh).

As for the laplacian problem, the code is very closed to the mathematical formulation.<br>
We define the product of function spaces for the flow speed and the flow pressur using \c THch<order>()(TH stands for Taylor-Hoods) function which is \c Pch<N+1> \f$\times\f$ \c Pch<N> for respectively flow speed and pressure spaces.<br>
We take an element
\f$U=\left(
    \begin{array}{c}
        u \\
        p \\
    \end{array}
\right)
\f$
in this space. Then we define the integrals of the variationnal formulation for the left and the right hand side. Finally, we apply the Poiseuille profile on the boundary.<br>
We call the solver to resolve the problem (\ref Solver).
\snippet mystokes.cpp marker_main

*/
/// [marker_main]
int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="mystokes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    // create the mesh
    auto mesh = loadMesh(_mesh=new Mesh<Simplex< 2 > > );


    // function space
    auto Vh = THch<2>( mesh );

    // element U=(u,p) in Vh
    auto U = Vh->element();
    auto u = U.element<0>();
    auto p = U.element<1>();

    // left hand side
    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=trace(gradt(u)*trans(grad(u))) );

    a+= integrate(_range=elements(mesh),
                  _expr=-div(u)*idt(p)-divt(u)*id(p));

    auto syms = symbols<2>();
    auto u1 = parse( option(_name="functions.alpha").as<std::string>(), syms );
    auto u2 = parse( option(_name="functions.beta").as<std::string>(), syms );
    matrix u_exact = matrix(2,1);
    u_exact = u1,u2;
    auto p_exact = parse( option(_name="functions.gamma").as<std::string>(), syms );
	auto f = -laplacian( u_exact, syms ) + grad( p_exact, syms ).transpose();
    LOG(INFO) << "rhs : " << f;

    // right hand side
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=trans(expr<2,1,5>( f, syms ))*id(u));
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=expr<2,1,5>(u_exact,syms));

    // solve a(u,v)=l(v)
    a.solve(_rhs=l,_solution=U);
    //# endmarker_main #

    double mean_p = mean(_range=elements(mesh),_expr=idv(p))(0,0);
    double mean_p_exact = mean(_range=elements(mesh),_expr=expr(p_exact,syms))(0,0);
    double l2error_u = normL2( _range=elements(mesh), _expr=idv(u)-expr<2,1,5>( u_exact, syms ) );
    double l2error_p = normL2( _range=elements(mesh), _expr=idv(p)-mean_p-(expr( p_exact, syms )-mean_p_exact) );
    LOG(INFO) << "L2 error norm u: " << l2error_u;
    LOG(INFO) << "L2 error norm p: " << l2error_p;

    // save results
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->add( "p", p );
    e->save();
}
/// [marker_main]
