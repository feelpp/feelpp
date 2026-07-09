/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE test_sb9mass
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/sb9_quadrature.hpp>
#include <feel/feelvf/vf.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>

#include "test_sb9_common.hpp"

using namespace Feel;
using namespace Feel::vf;
using namespace Feel::Tests::SB9;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( sb9_mass_suite )

BOOST_AUTO_TEST_CASE( sb9_lumped_mass_quadrature_is_full_vertex_lobatto )
{
    auto q = sb9LumpedMassLobatto();
    std::set<std::tuple<int, int, int>> vertices;

    BOOST_REQUIRE_EQUAL( q.nPoints(), 8 );
    for ( uint16_type k = 0; k < 8; ++k )
    {
        BOOST_CHECK_SMALL( q.weight( k ) - 1.0, 1e-14 );

        int const x = q.points()( 0, k ) > 0 ? 1 : -1;
        int const y = q.points()( 1, k ) > 0 ? 1 : -1;
        int const z = q.points()( 2, k ) > 0 ? 1 : -1;
        BOOST_CHECK_SMALL( std::abs( q.points()( 0, k ) ) - 1.0, 1e-14 );
        BOOST_CHECK_SMALL( std::abs( q.points()( 1, k ) ) - 1.0, 1e-14 );
        BOOST_CHECK_SMALL( std::abs( q.points()( 2, k ) ) - 1.0, 1e-14 );
        vertices.emplace( x, y, z );
    }

    BOOST_CHECK_EQUAL( vertices.size(), 8 );
    BOOST_CHECK_SMALL( q.weightsSum() - 8.0, 1e-14 );
}

BOOST_AUTO_TEST_CASE( sb9_lumped_mass_quadrature_diagonalizes_q1_displacement_mass )
{
    auto mesh = createAxisAlignedUnitPatch( "sb9_lumped_mass_axis_aligned_unit_patch" );
    auto Uh = Pchv<1>( mesh );
    auto u = Uh->element( "u" );
    auto v = Uh->element( "v" );
    auto massQuad = sb9LumpedMassLobatto();

    auto m = form2( _trial=Uh, _test=Uh );
    m = integrate( _range=elements( mesh ),
                   _quad=massQuad,
                   _quad1=massQuad,
                   _expr=inner( idt( u ), id( v ) ) );
    m.close();

    auto M = m.matrixPtr();
    double maxOffDiagonal = 0.0;
    double diagonalSum = 0.0;
    double minDiagonal = std::numeric_limits<double>::max();

    for ( size_type i = 0; i < M->size1(); ++i )
    {
        double const diagonal = ( *M )( i, i );
        diagonalSum += diagonal;
        minDiagonal = std::min( minDiagonal, diagonal );

        for ( size_type j = 0; j < M->size2(); ++j )
        {
            if ( i != j )
                maxOffDiagonal = std::max( maxOffDiagonal, std::abs( ( *M )( i, j ) ) );
        }
    }

    BOOST_CHECK_SMALL( maxOffDiagonal, 1e-13 );
    BOOST_CHECK_GT( minDiagonal, 0.0 );
    BOOST_CHECK_CLOSE( diagonalSum, 3.0, 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END()
