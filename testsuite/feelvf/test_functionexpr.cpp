/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Thibaut Metivet <thibaut.metivet@inria.fr>
 Date: 2021-06-03

 Copyright (C) 2021 Feel++ Consortium

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

// give a name to the testsuite
#define BOOST_TEST_MODULE vf functionexpr expr testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/functionexpr.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( functionexpr )

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace Feel;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2,1>>);
    auto rangeMeshElements = elements( mesh );
    auto Xh = Pch<1>( mesh );
    auto Xhv = Pchv<1>( mesh );
    auto Xht = Pchm<1>( mesh );
    auto u = Xh->element( inner(P()) );

    // Identity function expr
    auto f1expr = functionExpr( []( double p ) { return p; }, idv(u) );
    auto f1elt = vf::project( _space=Xh, _range=rangeMeshElements, _expr=f1expr );
    auto f1err = normL2( elements( mesh ), idv(u) - f1expr );
    BOOST_CHECK_SMALL( f1err, 1e-10 );

    // Vectorial function expr (with Order=4)
    auto vecExpr = vec( cos( M_PI * Px() ), sin( M_PI * Py() ) );
    auto f2expr = functionExpr<4>( []( auto const& x ) { return 2.*x; }, vecExpr );
    auto f2elt = vf::project( _space=Xhv, _range=rangeMeshElements, _expr=f2expr );
    auto f2err = normL2( elements( mesh ), cst(2.)*vecExpr - f2expr );
    BOOST_CHECK_SMALL( f2err, 1e-10 );

    // Mix scalar-gradient
    auto f3expr = functionExpr( []( auto const& p, auto const& x ) { return p*p + x.norm(); }, idv(u), vecExpr );
    auto f3elt = vf::project( _space=Xh, _range=rangeMeshElements, _expr=f3expr );
    auto f3err = normL2( elements( mesh ), idv(u)*idv(u) + norm2(vecExpr) - f3expr );
    BOOST_CHECK_SMALL( f3err, 1e-10 );

    // Tensor
    auto f4expr = functionExpr( []( auto const& x ) { return x * x.transpose(); }, vecExpr );
    auto f4elt = vf::project( _space=Xht, _range=rangeMeshElements, _expr=f4expr );
    auto f4err = normL2( elements( mesh ), vecExpr*trans(vecExpr) - f4expr );
    BOOST_CHECK_SMALL( f4err, 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END()
