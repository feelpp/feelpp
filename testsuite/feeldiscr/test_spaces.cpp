/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-07-15

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
/**
   \file test_spaces.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-07-15
 */
#include <sstream>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>

int main( int argc, char** argv)
{
    using namespace Feel;
    Environment env( _argc=argc,
                     _argv=argv );

    typedef Mesh<Simplex<2>> mesh_type;
    //typedef FunctionSpace<mesh_type,basis_type> fspace_type;
    
    auto mesh = loadMesh( _mesh = new mesh_type );
    auto P2v = Pchv<2>( mesh );
    auto P1 = Pch<1>( mesh );
    auto P0c =  Pch<0>( mesh );
    auto Xh = productPtr( P2v, P1, P0c );
    //auto Xh = fspace_type::New( mesh );
    auto u = P2v->element();
    auto p = P1->element();
    auto l = P0c->element();

    //auto t = U.element<3>();
    auto a = blockform2( _test=Xh, _strategy=solve::strategy::monolithic, _backend=backend() );
    auto b = blockform1( _test=Xh, _strategy=solve::strategy::monolithic, _backend=backend() );

    a(0_c,0_c)  = integrate( _range=elements(mesh), _expr=trans(idt(u))*id(u) );
    a(1_c,1_c) += integrate( _range=elements(mesh), _expr=trans(idt(p))*id(p) );
    a(2_c,2_c) += integrate( _range=elements(mesh), _expr=trans(idt(l))*id(l) );
    //a += integrate( elements(mesh), trans(idt(t))*id(t) );

    b(0_c)  = integrate( _range=elements(mesh), _expr=trans(vec(cst(1.),cst(1.)))*id(u) );
    b(1_c) += integrate( _range=elements(mesh), _expr=id(p) );
    b(2_c) += integrate( _range=elements(mesh), _expr=id(l) );
    //b += integrate( elements(mesh), id(t) );
    auto U = Xh->element();
    a.solve( _solution=U, _rhs=b );
#if 0
    a.matrixPtr()->printMatlab( "a.m" );
    b.vectorPtr()->printMatlab( "b.m" );
    U.printMatlab( "U.m" );
#endif
    U.vector()->setOnes();
    CHECK( math::abs(a(U,U)- 4)< 1e-10 ) << "invalid result a(1,1)=" << a(U,U) << " != " << 4;
}
