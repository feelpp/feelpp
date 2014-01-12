/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
	     Guillaume Dollé <guillaume.dolle@math.unistra.fr>

  Date 2013-02-25

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

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/matvec.hpp>
#include <feel/feelvf/on.hpp>


using namespace Feel;

/**
\page Advection Diffusion Advection Reaction Problem
\author Feel++ Consortium
\date 2013-02-25

\tableofcontents
<br>
<hr>

The diffusion advection reaction equation is a classical partial differential equation which can be found in many processes for example in chemistry or biology. This can be described by an equation containing a diffusion, an advection and a reaction term as follows,
<br><center>\f$
\left\{
\begin{aligned}
      -\epsilon\Delta  u + \bbeta \cdot \nabla  u + \mu u & =  f & \text{on}\; \Omega \;, \\
      u  & =  0 & \text{on}\; \partial\Omega \;, \\
\end{aligned}
\right
\f$</center><br>
We use here homogeneous Dirichlet boundary conditions.


\section Advection_Theory Theory
To establish the variationnal formulation, as always we mutiply the first equation by a test function \f$v\in H_0^1(\Omega)\f$ such that,
\f[
    H_0^1(\Omega) = \{ v\in H^1(\Omega),\; v=0 \; \text{on} \; \partial\Omega \} \;.
\f]

Then we integrate on the domain \f$\Omega\f$,
<br><center>\f$
\begin{aligned}
- \int_\Omega \epsilon \Delta u\ v
+ \int_\Omega \bbeta \cdot \nabla u\ v
+ \int_\Omega \mu\ u\ v
= \int_\Omega f\ v \;.
\end{aligned}
\f$</center><br>
We establish the variationnal formulation from the previous equation and using the Green formula, find \f$u \in \in H_0^1(\Omega)\f$
<br><center>\f$
\begin{aligned}
  \int_\Omega \epsilon \nabla u \cdot \nabla v
  - \underbrace{\int_{\partial\Omega} \epsilon (\nabla u \cdot \mathbf n)\ v}_{=0}
  + \int_\Omega (\beta \cdot \nabla u)\ v
  + \int_\Omega \mu\ u\ v
  = \int_\Omega f v \; \quad \forall v \in H_0^1(\Omega),
\end{aligned}
\f$</center><br>
where \f$\mathbf n\f$ is a unit outward normal vector. We can rewrite the problem, find \f$u \in \in H_0^1(\Omega)\f$
\f[
a(u,v) = l(v) \quad \forall v \in H_0^1(\Omega),
\f]
where \f$a\f$ is a bilinear form, continuous, coercive and \f$l\f$ is a linear form.

\section Advection_Implementation Implementation
We choose for our example \f$\mu = 1 \f$, \f$\epsilon = 1\f$, \f$f=1\f$, and \f$\bbeta=(1,1)^T\f$.
\snippet myadvection.cpp marker_main

Again the implementation is close to the mathematical formulation.<br>
Here again, we create the mesh for an unit square geometry.<br>
Then we define the function space \f$X_h\f$ we choose as order 1 Lagrange basis function using \c Pch<Order>().  Note that here, the function space is the same for "trial" and "test" functions.<br>
We declare the left and the right hand side integrals expressions for the equation (\ref Advection_Math ). <br>
Finally we add the Dirichlet boundary condition and we use the default solver to solve. We export the solution \f$u\f$ for post processing.

\section Advection_Results Results
There are various solutions of this problem displayed with Paraview for \f$\epsilon= 1, 0.01, 0.0001\f$.<br>
Notice how the solution gets unstable for\f$\epsilon = 0.0001\f$, this is classical and requires stabilisation methods to handle this issue.

<center>
<table border=0px>
<tr>
  <td>\image html sol-dar-1.png</td>
  <td>\image html sol-dar-2.png</td>
  <td>\image html sol-dar-3.png</td>
</tr>
<tr>
  <td><center>\f$\epsilon=1\f$</center></td>
  <td><center>\f$\epsilon=0.01\f$</center></td>
  <td><center>\f$\epsilon=0.0001\f$</center></td>
</tr>
</table>
</center>
 */

/// [marker_main]
int
main( int argc, char** argv )
{
    po::options_description opts ( "Advection diffusion reaction options ");
    opts.add_options()
        ( "epsilon", po::value<double>()->default_value( 1 ), "diffusion term coefficient" )
        ( "betax", po::value<double>()->default_value( 1 ), "convection term coefficient in x-direction" )
        ( "betay", po::value<double>()->default_value( 1 ), "convection term coefficient in y-direction" )
        ( "mu", po::value<double>()->default_value( 1 ), "reaction term coefficient" );
    // Initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=opts.add( feel_options() ),
                     _about=about(_name="myadvection",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org") );
    // create mesh
    auto mesh = unitSquare();

    // function space
    auto Xh = Pch<1>( mesh );
    auto u = Xh->element( "u" );
    auto v = Xh->element( "v" );

    // diffusion coeff.
    double epsilon = option(_name="epsilon").as<double>();
    // reaction coeff.
    double mu = option(_name="mu").as<double>();
    auto beta = vec( cst(option(_name="betax").as<double>()),
                     cst(option(_name="betay").as<double>()) );
    auto f = cst(1.);

    // left hand side
    auto a = form2( _test=Xh, _trial=Xh );
    a += integrate( _range=elements( mesh ),
                    _expr=( epsilon*gradt( u )*trans( grad( v ) )
                         + ( gradt( u )*beta )*id(v)
                         + mu*idt( u )*id( v ) ) );

    // right hand side
    auto l = form1( _test=Xh );
    l+= integrate( _range=elements( mesh ), _expr=f*id( v ) );

    // boundary condition
    a += on( _range=boundaryfaces( mesh ), _rhs=l, _element=u,
             _expr=cst(0.) );

    // solve the system
    a.solve( _rhs=l, _solution=u );

    // export results
    auto e = exporter( _mesh=mesh );
    e->add("u",u);
    e->save();
} // end main
/// [marker_main]






