/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE test_sb9mass
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/sb9_inertia.hpp>
#include <feel/feelvf/sb9_quadrature.hpp>
#include <feel/feelvf/vf.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>
#include <Eigen/Cholesky>

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

BOOST_AUTO_TEST_CASE( sb9_consistent_internal_mass_has_analytic_blocks )
{
    auto mesh = createAxisAlignedUnitPatch( "sb9_consistent_internal_mass" );
    auto Uh = Pchv<Dynamic>( mesh, RuntimeOrder{1} );
    auto Ah = Pdh<Dynamic>( mesh, RuntimeOrder{0} );
    auto ps = product( Uh, Ah );
    auto u = trial( Uh );
    auto v = test( Uh );
    auto alpha = trial( Ah );
    auto beta = test( Ah );
    auto mass = blockform2( ps, solve::strategy::monolithic, backend() );
    assembleSb9Inertia( mass, elements( mesh ), u, v, alpha, beta, 1.0, 1.0, true );
    mass.close();

    auto M = mass.matrixPtr();
    auto const n = Uh->nDof();
    BOOST_REQUIRE_EQUAL( n, 24 );
    BOOST_REQUIRE_EQUAL( M->size1(), n + 1 );
    BOOST_CHECK_CLOSE( ( *M )( n, n ), 8.0 / 15.0, 1e-9 );
    int coupledDofs = 0;
    Eigen::MatrixXd dense( n + 1, n + 1 );
    for ( size_type i = 0; i <= n; ++i )
    {
        for ( size_type j = 0; j <= n; ++j )
            dense( i, j ) = ( *M )( i, j );
        if ( i < n )
        {
            BOOST_CHECK_SMALL( ( *M )( i, n ) - ( *M )( n, i ), 1e-13 );
            if ( std::abs( ( *M )( i, n ) ) > 1e-12 )
            {
                ++coupledDofs;
                BOOST_CHECK_CLOSE( ( *M )( i, n ), 1.0 / 12.0, 1e-9 );
            }
        }
    }
    BOOST_CHECK_EQUAL( coupledDofs, 8 );
    BOOST_CHECK( dense.llt().info() == Eigen::Success );

    auto bubble = sb9DisplacementBubble();
    double const vertexMass = integrate( _range=elements( mesh ),
                                         _quad=sb9LumpedMassLobatto(),
                                         _quad1=sb9LumpedMassLobatto(),
                                         _expr=inner( bubble, bubble ) ).evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( vertexMass, 1e-14 );
}

BOOST_AUTO_TEST_CASE( sb9_consistent_internal_mass_is_invariant_under_rotation )
{
    auto mesh = createRotatedShellPatch( "sb9_rotated_internal_mass" );
    auto Ah = Pdh<Dynamic>( mesh, RuntimeOrder{0} );
    auto alpha = trial( Ah );
    auto beta = test( Ah );
    auto bubble = sb9DisplacementBubble();
    auto mass = form2( _trial=Ah, _test=Ah );
    mass = integrate( _range=elements( mesh ), _quad=_Q<6>(), _quad1=_Q<6>(),
                      _expr=inner( bubble * id( alpha ), bubble * id( beta ) ) );
    mass.close();
    // This patch has area 2 and thickness 0.2.
    BOOST_CHECK_CLOSE( ( *mass.matrixPtr() )( 0, 0 ), ( 8.0 / 15.0 ) * 0.4, 1e-9 );
}

BOOST_AUTO_TEST_SUITE_END()
