/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-05-01

  Copyright (C) 2013 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <feel/feel.hpp>
using namespace Feel;

/**
   \page StokesCurl
   \author Christophe Prud'homme
   \date 2013-05-01

   \tableofcontents
   \section StokesCurl_Theory Theory

   We are interested in solving
   \f$ -\nu \Delta \bvec{u} + \nabla p  = \bvec{f},\quad \nabla \cdot \bvec{u} = 0 \f$
   in curl-curl formulation in 2D.

   First we introduce the following notations:
   \f[
   u \times n = u_1 n_2 - u_2 n_!, \quad u = (u_1, u_2),\ n = (n_1, n_2)
   \f]
   The curl of a vectorial field \f$\phi\f$
   \f[
   \mathrm{curl} \phi  = (-\frac{\partial \phi_2}{\partial x_1}, -\frac{\partial \phi_1}{\partial x_2} )
   \f]
   The curl of a scalar field $\f$\psi\f$
   \f[
   \nabla \times \psi  = \frac{\partial \psi}{\partial x_1} -\frac{\partial \psi}{\partial x_2}
   \f]

   \section StokesCurl_Implementation Implementation
   \snippet stokes/stokes_curl.cpp marker
*/
/// [marker]
int main(int argc, char**argv )
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="stokes_curl",
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
                  _expr=curlxt(u)*curlx(u)  + divt(u)*div(u) );

    a+= integrate(_range=elements(mesh),
                  _expr=-div(u)*idt(p)-divt(u)*id(p));

    a+= integrate(_range=markedfaces(mesh,"inlet"),
                  _expr=( curlxt(u)*cross(id(u),N()) +
                          curlx(u)*cross(idt(u),N()) +
                          -trans(divt(u)*N())*id(u)  +
                          30*cross(id(u),N())*cross(idt(u),N())/hFace() ) );
    a+= integrate(_range=markedfaces(mesh,"outlet"),
                  _expr=( curlxt(u)*cross(id(u),N()) +
                          curlx(u)*cross(idt(u),N()) +
                          -trans(divt(u)*N())*id(u)  +
                          30*cross(id(u),N())*cross(idt(u),N())/hFace() ) );

    // right hand side
    auto l = form1( _test=Vh );

    l += integrate(_range=markedfaces(mesh,"inlet"),
                   _expr=2*trans(id(u))*N() );
    l += integrate(_range=markedfaces(mesh,"outlet"),
                   _expr=1*trans(id(u))*N() );

    // Dirichlet condition
    a+=on(_range=markedfaces(mesh,"wall"), _rhs=l, _element=u,
          _expr=vec(cst(0.),cst(0.)));

    // solve a(u,v)=l(v)
    a.solve(_rhs=l,_solution=U);

    // save results
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->add( "p", p );
    e->save();
}
/// [marker]
