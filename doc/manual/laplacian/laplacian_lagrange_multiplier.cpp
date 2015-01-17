/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme
   Date     : Fri Mar 21 18:10:49 2014

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
#include <feel/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/product.hpp>



int main(int argc, char**argv )
{
    using namespace Feel;
    // initialize feel++
    Environment env( _argc=argc, _argv=argv,
                     _desc_lib=feel_options().add( backend_options( "u" ) ).add ( backend_options( "gradu" ) ),
                     _about=about(_name="laplacian_lagrange_multiplier",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh = new Mesh<Simplex<2>> );
    typedef FunctionSpace<Mesh<Simplex<2> >, bases<Lagrange<1, Scalar>, Lagrange<0, Scalar> > > space_type;
#if 0
    auto X1 = Pch<1>(mesh);
    auto X2 = Pch<0>(mesh);
    //auto Vh = space_type::NewFromList( X1, X2 );
    //auto Vh = product( X1, X2 );
    //BOOST_MPL_ASSERT_MSG( (boost::is_same<space_type, decltype(*Vh)>::value),
    //INVALID_SPACE,
    //(space_type, decltype(*Vh)));
    std::cout << "X : " << X1->nLocalDof() << "\n";
    std::cout << "Vh->X1 : " << Vh->functionSpace<0>()->nLocalDof() << "\n";
    //std::cout << "Vh->X2 : " << Vh->functionSpace<1>()->nLocalDof() << "\n";
#else
    auto Vh = space_type::New( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>() ;
    auto lambda = U.element<1>() ;
    auto v = V.element<0>() ;
    auto nu = V.element<1>() ;

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements(mesh),
                   _expr=gradt(u)*trans(grad(v)));
    a += integrate( _range=elements(mesh),
                    _expr=id( v )*idt( lambda ) + idt( u )*id( nu ) );

    auto f = expr(soption(_name="functions.f"));
    auto l = form1( _test=Vh );
    l = integrate( _range=elements(mesh),
                   _expr=f*id(v) );
    auto p = expr(soption(_name="functions.p"));
    auto q = expr(soption(_name="functions.q"));
    l+= integrate( _range=boundaryfaces(mesh),
                   _expr=trans(vec(p,q))*N()*id(v) );
    a.solve( _name="u", _rhs=l, _solution=U );

    //auto Xh = Vh->functionSpace<0>();
    auto Xh = Pchv<1>(mesh);
    auto gradu = Xh->element();
    auto b = form2( _trial=Xh, _test=Xh );
    b = integrate( _range=elements(mesh), _expr=trans(idt(gradu))*id(gradu));
    auto F = form1( _test=Xh );
    F = integrate( _range=elements(mesh), _expr=gradv( U.element<0>() )*id(gradu));
    b.solve( _name="gradu", _rhs=F, _solution=gradu );

    auto e = exporter( _mesh=mesh );
    auto z = expr(soption(_name="functions.z"));
    v.on( _range=elements(mesh), _expr=z );
    e->add( "z", v );
    e->add( "u", U.element<0>() );
    e->add( "grad_u", gradu );
    v.on( _range=elements(mesh), _expr=p );
    e->add( "p", v );
    v.on( _range=elements(mesh), _expr=q );
    e->add( "q", v );
    gradu.on( _range=elements(mesh), _expr=vec(p,q) );
    e->add( "(p,q)", gradu );
    e->save();
#endif
}
