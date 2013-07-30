/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
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
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-07-15
 */
#include <sstream>
#include <boost/timer.hpp>
#include <feel/feel.hpp>

int main( int argc, char** argv)
{
    using namespace Feel;
    Environment env( _argc=argc,
                     _argv=argv );

    typedef Lagrange<2,Vectorial> b1_type;
    typedef Lagrange<1,Scalar> b2_type;
    typedef Lagrange<0,Scalar,Continuous> b3_type;
    typedef Lagrange<1,Scalar> b4_type;
    //typedef bases<b1_type,b2_type,b3_type,b4_type> basis_type;
    typedef bases<b1_type,b2_type,b3_type> basis_type;
    typedef Mesh<Simplex<2>> mesh_type;
    typedef FunctionSpace<mesh_type,basis_type> fspace_type;
    auto mesh = loadMesh( _mesh = new mesh_type );
    auto Xh = fspace_type::New( mesh );
    auto v = backend()->newVector( Xh );
    auto M = backend()->newMatrix( Xh, Xh );
    M->graph()->showMe( std::cout );

    auto U = Xh->element();
    auto u = U.element<0>();
    auto p = U.element<1>();
    auto l = U.element<2>();
    //auto t = U.element<3>();
    auto a = form2( Xh, Xh );
    a  = integrate( elements(mesh), trans(idt(u))*id(u) );
    a += integrate( elements(mesh), trans(idt(p))*id(p) );
    a += integrate( elements(mesh), trans(idt(l))*id(l) );
    //a += integrate( elements(mesh), trans(idt(t))*id(t) );
    auto b = form1( Xh );
    b  = integrate( elements(mesh), trans(vec(cst(1.),cst(1.)))*id(u) );
    b += integrate( elements(mesh), id(p) );
    b += integrate( elements(mesh), id(l) );
    //b += integrate( elements(mesh), id(t) );
    a.solve( _solution=U, _rhs=b );
#if 0
    a.matrixPtr()->printMatlab( "a.m" );
    b.vectorPtr()->printMatlab( "b.m" );
    U.printMatlab( "U.m" );
#endif
    U.setOnes();
    CHECK( math::abs(a(U,U)- 4)< 1e-10 ) << "invalid result a(1,1)=" << a(U,U) << " != " << 4;
}
