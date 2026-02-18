/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-02-18

  Copyright (C) 2026 Feel++ Consortium

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
/**
 * @file test_hdiv_runtime_order_parity.cpp
 * @brief Runtime-order constructor parity checks for H(div)/H(curl) wrappers.
 */
#define BOOST_TEST_MODULE test_hdiv_runtime_order_parity
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/neh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelvf/vf.hpp>

#include <string>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using mesh_type = Mesh<Simplex<2,1>>;
using mesh_ptrtype = std::shared_ptr<mesh_type>;

mesh_ptrtype
makeOneElementSimplexMesh( std::string const& name )
{
    auto desc = std::make_shared<GmshSimplexDomain>( 2, 1, Gmsh::GMSH_REFERENCE_DOMAIN );
    desc->setPrefix( name );
    desc->usePhysicalNames( true );
    desc->setCharacteristicLength( 2.0 );
    return createGMSHMesh( _mesh = new mesh_type, _desc = desc, _h = desc->h() );
}

template<typename SpacePtrType, typename MeshPtrType, typename ExprType>
void
checkOnParity( SpacePtrType const& staticSpace,
               SpacePtrType const& runtimeSpace,
               MeshPtrType const& mesh,
               ExprType const& expr,
               double tol )
{
    BOOST_REQUIRE_EQUAL( staticSpace->nDof(), runtimeSpace->nDof() );
    BOOST_CHECK_EQUAL( staticSpace->runtimeOrder(), runtimeSpace->runtimeOrder() );

    auto uStatic = staticSpace->element();
    auto uRuntime = runtimeSpace->element();
    uStatic.on( _range = elements( mesh ), _expr = expr );
    uRuntime.on( _range = elements( mesh ), _expr = expr );

    auto const parityError = normL2( _range = elements( mesh ),
                                     _expr = idv( uStatic ) - idv( uRuntime ),
                                     _quad = _Q<12>() );
    BOOST_CHECK_SMALL( parityError, tol );
}

} // namespace

BOOST_AUTO_TEST_SUITE( hdiv_runtime_order_parity_suite )

BOOST_AUTO_TEST_CASE( rt_runtime_order_constructor_parity )
{
    auto mesh = makeOneElementSimplexMesh( "rt-runtime-parity" );
    auto expr = vec( cst( 1.0 ) + Px(), cst( -0.5 ) + 2.0 * Py() );

    auto rt0Static = RTh<0>( mesh );
    auto rt0Runtime = RTh<0>( mesh, RuntimeOrder{ 0 } );
    checkOnParity( rt0Static, rt0Runtime, mesh, expr, 1e-12 );
}

BOOST_AUTO_TEST_CASE( bdm_runtime_order_constructor_parity )
{
    auto mesh = makeOneElementSimplexMesh( "bdm-runtime-parity" );
    auto expr = vec( cst( 0.75 ) + Px() + Px() * Py(),
                     cst( -1.25 ) + Py() + Px() * Px() );

    auto bdm0Static = BDMh<0>( mesh );
    auto bdm0Runtime = BDMh<0>( mesh, RuntimeOrder{ 0 } );
    checkOnParity( bdm0Static, bdm0Runtime, mesh, expr, 1e-12 );
}

BOOST_AUTO_TEST_CASE( neh_runtime_order_constructor_parity )
{
    auto mesh = makeOneElementSimplexMesh( "neh-runtime-parity" );
    auto expr = vec( -Py(), Px() );

    auto neh0Static = Neh<0>( mesh );
    auto neh0Runtime = Neh<0>( mesh, RuntimeOrder{ 0 } );
    checkOnParity( neh0Static, neh0Runtime, mesh, expr, 1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
