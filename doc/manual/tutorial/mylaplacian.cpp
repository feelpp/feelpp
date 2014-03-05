/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>

  Date 2013-02-18

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

/** \page LaplacianTutorial Laplacian with homogeneous Dirichlet conditions
\author Feel++ Consortium
\date 2013-02-11

\tableofcontents
<br>
<hr>
<br>

\section Laplacian_Theory Theory
This part explains how to solve the Laplacian equation for homogeneous dirichlet conditions,
<br><center>\f$
\left\{
\begin{aligned}
   -\Delta u & =  f & \text{on}\;\Omega \;, \\
            u & =  0 & \text{on}\;\partial\Omega \;,\\
\end{aligned}
\right.
\f$</center><br>
where \f$u\in\Omega\f$ is the unknown "trial" function and \f$\Omega\f$ the domain.

We multiply each part of the first equation by a "test" function \f$v\in H_0^1(\Omega)\f$ and we integrate the resulting equation on the domain \f$\Omega\f$,
<br><center>\f$
\begin{aligned}
 -\int_\Omega \Delta u v = \int_\Omega f v \;.
\end{aligned}
\f$</center><br>
We can integrate by parts this equation (Green Theorem) to obtain the variationnal formulation,
\f[
\begin{aligned}
\int_\Omega \nabla u \nabla v
-\underbrace{ \int_{\partial\Omega} \frac{\partial u}{\partial n} v }_{= 0}
=\int_\Omega f v \;
\end{aligned}
\f]
where \f$n\f$ denotes a unit outward normal vector to the boundary. We can rewrite the problem as find \f$u\in H_0^1(\Omega)\f$ such that for all \f$v\in H_0^1(\Omega)\f$,
<br><center>\f$
\begin{aligned}
a(u,v)=l(v) \;,
\end{aligned}
\f$</center><br>
where \f$a\f$ is a bilinear form, continuous, coercive and \f$l\f$ a linear form.

\section Laplacian_Implementation Implementation
Let's take a look at the \feel code (source \c "doc/manual/tutorial/mylaplacian.cpp").<br>
We consider for this example \f$f=1\f$ constant.
\snippet mylaplacian.cpp marker1

As you can see, the program looks very close to the mathematical formulation.<br>
We use the \c form2() function to define the bilinear form and \c form1() for the linear one (see \ref Forms ).<br>
The gradient for the trial functions is declared with the \c gradt() expression where as \c grad() is used for the test functions (see \ref Keywords).
Note that we need to transpose the second vector to perform the scalar product.

To introduce the homogeneous dirichlet conditions on the boundary, we use the function \c on(). Once the variationnal formulation and the boundary conditions are set, we call
the solver with \c solve().

*/

/// [marker1]
int main(int argc, char**argv )
{
    // initialize feel++
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="mylaplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    // create mesh
    auto mesh = unitSquare();

    // function space
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    // left hand side
    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );

    // right hand side
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v));

    // apply the boundary condition
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=constant(0.) );

    // solve the equation a(u,v) = l(v)
    a.solve(_rhs=l,_solution=u);

    // export results
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
}
/// [marker1]
