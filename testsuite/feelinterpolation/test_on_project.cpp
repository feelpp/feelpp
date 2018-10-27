/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Guillaume Dollé <gdolle@unistra.fr>
 Date: 26 Mar 2015

 Copyright (C) 2015 Feel++ Consortium

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

#define BOOST_TEST_MODULE test_on_dofs
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
//#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

template<int Dim, int PolyOrder=1>
void runTestAssign()
{
    // Check the on keyword for different topological entities.

    using mesh_type = Mesh< Simplex<Dim> >;
    auto mesh = loadMesh( _mesh=new mesh_type );
    auto Xh = Pch<PolyOrder>(mesh);
    auto p = Xh->element();
    auto l = Xh->element();
    auto s = Xh->element();
    auto v = Xh->element();

    size_type np = nelements( markedpoints( mesh,"P"),true );
    BOOST_CHECK( np > 0 );
    p.on( _range=markedpoints( mesh,"P"), _expr=cst(42.) );
    sync( p );
    double sump = p.sum();
    BOOST_CHECK_SMALL( sump - 42, 1e-12 );

    size_type ne = nelements( markededges( mesh,"L"),true );
    BOOST_CHECK( ne > 0 );
    l.on( _range=markededges( mesh,"L"), _expr=cst(42.) );
    sync( l );
    double maxl = l.max();
    BOOST_CHECK_SMALL( maxl - 42, 1e-12 );

    s.on( _range=markedfaces( mesh,"S"), _expr=cst(42.) );
    double ints = integrate(_range=markedfaces( mesh,"S"),_expr=idv(s) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( ints - 42, 1e-12 );

    v.on( _range=elements( mesh ), _expr=cst(42.) );
    double intv = integrate(_range=elements( mesh ),_expr=idv(v) ).evaluate()(0,0);
    BOOST_CHECK_SMALL( intv - 42*27, 1e-10 );

#if 0
    auto e = exporter( _mesh=mesh );
    e->add( "p", p);
    e->add( "l", l);
    e->add( "s", s);
    e->add( "v", v);
    e->save();
#endif
}


template<int Dim, int PolyOrder=1>
void runTestElimination()
{
    using mesh_type = Mesh< Simplex<Dim> >;
    auto mesh = loadMesh( _mesh=new mesh_type );
    auto Xh = Pch<PolyOrder>( mesh );

    auto u = Xh->element();
    auto mat = backend()->newMatrix(_test=Xh,_trial=Xh);
    auto rhs = backend()->newVector( Xh );
    auto l = form1( _test=Xh,_vector=rhs );
    l = integrate(_range=elements(mesh),
                  _expr=id(u));
    auto a = form2( _trial=Xh, _test=Xh,_matrix=mat );
    a = integrate(_range=elements(mesh),
                  _expr=inner(gradt(u),grad(u)) );
    a+=on(_range=markedfaces( mesh,"S"), _rhs=rhs, _element=u, _expr=cst(12.) );
    a+=on(_range=markededges( mesh,"L"), _rhs=rhs, _element=u, _expr=cst(22.) );
    a+=on(_range=markedpoints( mesh,"P"), _rhs=rhs, _element=u, _expr=cst(32.) );

    backend(_rebuild=true)->solve(_matrix=mat,_rhs=rhs,_solution=u );

#if 0
    // not yet work because of idv(u) in markededges and markedpoints
    // require to tell to the geomap that we manipulate edge or points
    // can be tested by remove if ( fusion::at_key<key_type>( geom )->faceId() != invalid_uint16_type_value  ) in feel/feevf/operator.hpp line 781
    auto err = Xh->element();
    err.on( _range=markedfaces( mesh,"S"), _expr=abs(idv(u)-cst(12.)) );
    err.on( _range=markededges( mesh,"L"), _expr=abs(idv(u)-cst(22.)) );
    err.on( _range=markedpoints( mesh,"P"), _expr=abs(idv(u)-cst(32.)) );
    BOOST_CHECK_SMALL( err.sum(), 1e-10 );
#endif
}

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_on_dofs )

BOOST_AUTO_TEST_CASE( assign_3d )
{
    runTestAssign<3>();
}

BOOST_AUTO_TEST_CASE( elimination_3d )
{
    runTestElimination<3>();
}

BOOST_AUTO_TEST_SUITE_END()

