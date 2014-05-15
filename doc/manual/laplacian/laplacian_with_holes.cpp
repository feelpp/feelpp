// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
/*
  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-01

  Copyright (C) 2013-2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <feel/feel.hpp>

namespace Feel
{
/**
   \page LaplacianWithHoles Laplacian in a perforated domain
   \author Feel++ Consortium

   \tableofcontents

   <br>
   <br>

   \section LaplacianWithHoles_Problem Problem statement

   We solve for the laplacian with homogeneous Dirichlet conditions in a domain with a hole
   \f[
   \left\{
   \begin{aligned}
   -\Delta u & =  f & \text{on}\;\Omega \;, \\
   u & =  0 & \text{on}\;\partial\Omega \;,\\
   \end{aligned}
   \right.
   \f]
   where \f$u\in\Omega\f$ is the unknown "trial" function and \f$\Omega\f$ the domain.

   The variational formulation reads, find \f$u \in H^1_0(\Omega)\f$ such that
   \f$\forall v \in H^1_0(\Omega)\f$

   \f[
   \int_\Omega \nabla u \cdot  \nabla v -\underbrace{ \int_{\partial\Omega} \frac{\partial u}{\partial n} v }_{= 0}\ =\ \int_\Omega f v \;
   \f]

   where \f$n\f$ denotes a unit outward normal vector to the boundary.  We can
   rewrite the problem as to find \f$u\in H_0^1(\Omega)\f$ such that for all \f$v\in
   H_0^1(\Omega)\f$,

   \f[
   a(u,v)=l(v) \;,
   \f]

   where \f$a\f$ is a bilinear form, continuous, coercive and \f$l\f$ a linear form.

   \section Laplacian_Implementation Implementation

   We defined \f$\Omega\f$ as the unit square with a circle inside of radius \f$0.25\f$
   \snippet laplacian_with_holes.cpp marker1

   We consider for this example \f$f=1\f$ constant.
   \snippet laplacian_with_holes.cpp marker2

   The complete example is here
   \snippet laplacian_with_holes.cpp marker3

   As you can see, the program looks very close to the mathematical
   formulation.<br> We use the \c form2() function to define the bilinear form and
   \c form1() for the linear one (see \ref Forms ).<br> The gradient for the trial
   functions is declared with the \c gradt() expression where as \c grad() is used
   for the test functions (see \ref Keywords).  Note that we need to transpose the
   second vector to perform the scalar product.

   To introduce the homogeneous dirichlet conditions on the boundary, we use the
   function \c on(). Once the variationnal formulation and the boundary conditions
   are set, we call the solver with \c solve().

   \section LaplacianWithHoles_Results Results

   The program is named `feelpp_doc_laplacian_with_holes`.

*/
}
int main(int argc, char**argv )
{
    /// [marker3]
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="laplacian-with_holes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    /// [marker1]
    typedef Mesh<Simplex<2> > mesh_type;

    GeoTool::Rectangle R1( 0.05,"R1",GeoTool::Node( 0,0 ),GeoTool::Node( 1,1 ) );
    GeoTool::Circle C1( 0.05,"C1",GeoTool::Node( 0.5,0.5 ),GeoTool::Node( 0.75,0.75 ) );

    auto R1mesh = R1.createMesh(_mesh=new mesh_type,_name="R1" );
    auto C1mesh = C1.createMesh(_mesh=new mesh_type,_name="C1" );
    auto R1mC1mesh = ( R1-C1 ).createMesh(_mesh=new mesh_type,_name="R1-C1" );
    /// [marker1]

    // auto mesh=R1mesh;
	// auto mesh=C1mesh;
    auto mesh=R1mC1mesh;

    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    /// [marker2]
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v));
    /// [marker2]

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=constant(0.) );
    a.solve(_rhs=l,_solution=u);

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
    /// [marker3]
    return 0;
}
