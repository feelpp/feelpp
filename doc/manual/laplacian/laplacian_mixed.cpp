/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-02-10

  Copyright (C) 2014 Feel++ Consortium

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
/**
   \file laplacian_mixed.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-02-10
 */

#include <feel/feel.hpp>

namespace Feel
{
/*
   \page LaplacianMixed Laplacian in mixed formulation
   \author Feel++ Consortium

   We wish to solve the following Poisson problem, find \f$p\f$, such that:

   \f[ - \Delta p = f \mbox{ on } \Omega,\f]
   with
   \f[ \partial p / \partial n = (g1d,g2d) \cdot n \mbox{ on } \Gamma_{123}, and
   p = gd   \mbox{ on } \Gamma_{1} \f].

   We are interested in the mixed formulation, find $p\in L^2(\Omega)$ and $u\in
   H(div) $ such than

   \f[ u - \nabla p = 0, \quad - \nabla \cdot \mathbf{u}  =  f \f]
   with
   \f[ \mathbf{u} \cdot  \mathbf{n} = (g1d,g2d) \cdot \mathbf{n}  \mbox{ on } \Gamma_{123}, \quad
   p = gd  \mbox{ on } \Gamma_{1}\f\]

   The variational formulation reads

   \f[ \forall v\in H(div)$;  \mathbf{ v } \cdot \mathbf{ n } = 0 \mbox{ on } \Gamma_{4}\] :

   \f[ \int_\Omega  \mathbf{u} \cdot \mathbf{ v } + p div v -\int_{\Gamma_{123}} gd* \mathbf{v} \cdot \mathbf{n}  = 0 \f]

   \f[ \forall q\in L^2$:   $  -\int_\Omega q div u = \int_Omega f q \f]
   and \f[ u \cdot \mathbf{ n } = (g1n,g2n) \cdot \mathbf{ n }  \mbox{ on } \Gamma_4\f]

*/
}
int main(int argc, char**argv )
{
    /// [marker1]
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="laplacian_mixed",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3>> );
    auto Mh = DhPdh<0>( mesh );
    auto U = Mh->element();
    auto V = Mh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();

    auto l = form1( _test=Mh );
    l = integrate(_range=elements(mesh), _expr=id(q) );
    for( int i : std::vector<int>{ 1, 2, 3 } )
        l += integrate(_range=markedfaces(mesh,(boost::any)i), _expr=trans(id(v))*N());

    auto a = form2( _trial=Mh, _test=Mh );

    a = integrate(_range=elements(mesh), _expr=trans(idt(u))*id(v));
    a += integrate(_range=elements(mesh), _expr=id(q)*divt(u) - idt(p)*div(v) );
    a += integrate(_range=elements(mesh), _expr=1e-6*id(q)*idt(p));
    a += on( markedfaces(mesh, (boost::any)4 ), _element=u, _rhs=l, _expr=vec(cst(1.0),cst(1.),cst(0.) ) );

    a.solve(_rhs=l,_solution=U);

    auto Xh = Pchv<1>( mesh );
    auto w = Xh->element();
    auto b = form2( _test=Xh, _trial=Xh );
    b = integrate(_range=elements(mesh), _expr=trans(idt(w))*id(w));
    auto ll = form1( _test=Xh );
    ll = integrate( _range=elements(mesh), _expr=trans(idv(u))*id(w));
    b.solve( _solution=w, _rhs=ll, _rebuild=true );

    auto Yh = Pdh<0>( mesh );
    auto d = Yh->element();
    auto c = form2( _test=Yh, _trial=Yh );
    c = integrate(_range=elements(mesh), _expr=(idt(d))*id(d));
    auto crhs = form1( _test=Yh );
    crhs = integrate( _range=elements(mesh), _expr=divv(u)*id(d));
    c.solve( _solution=d, _rhs=crhs, _rebuild=true );

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->add( "ul2", w );
    e->add( "p", p );
    e->add( "laplacian(p)", d );
    e->save();
    /// [marker1]
}
