/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file test_de_rham.cpp
 * @brief Discrete de Rham and curved Piola-map contracts.
 */
#define BOOST_TEST_MODULE test_de_rham
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/neh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using affine_mesh_type = Mesh<Simplex<2, 1>>;

std::shared_ptr<affine_mesh_type>
makeAffineMesh()
{
    auto desc = std::make_shared<GmshSimplexDomain>( 2, 1, Gmsh::GMSH_REFERENCE_DOMAIN );
    desc->setPrefix( "de-rham-affine" );
    desc->setCharacteristicLength( 0.35 );
    return createGMSHMesh( _mesh = new affine_mesh_type, _desc = desc, _h = desc->h() );
}

template<typename SpacePtrType, typename MeshPtrType, typename ExprType>
auto
interpolate( SpacePtrType const& space, MeshPtrType const& mesh, ExprType const& expr )
{
    auto result = space->element();
    result.on( _range = elements( mesh ), _expr = expr );
    return result;
}

template<typename MeshPtrType, typename ElementType>
double
divergenceTheoremResidual( MeshPtrType const& mesh, ElementType const& field )
{
    auto const volumeFlux = integrate( _range = elements( mesh ),
                                       _expr = divv( field ),
                                       _quad = _Q<16>(),
                                       _geomap = GeomapStrategyType::GEOMAP_HO ).evaluate()( 0, 0 );
    auto const boundaryFlux = integrate( _range = boundaryfaces( mesh ),
                                         _expr = trans( N() ) * idv( field ),
                                         _quad = _Q<16>(),
                                         _geomap = GeomapStrategyType::GEOMAP_HO ).evaluate()( 0, 0 );
    return std::abs( volumeFlux - boundaryFlux );
}

template<int Order>
auto
potentialExpression()
{
    if constexpr ( Order == 0 )
        return cst( 1.0 ) + 2.0 * Px() - 3.0 * Py();
    else if constexpr ( Order == 1 )
        return cst( 1.0 ) + Px() + 2.0 * Py() + Px() * Px() + Px() * Py() - Py() * Py();
    else if constexpr ( Order == 2 )
        return cst( 1.0 ) + Px() - Py() + Px() * Px() * Px() + Px() * Px() * Py() +
               Px() * Py() * Py() - Py() * Py() * Py();
    else if constexpr ( Order == 3 )
        return cst( 1.0 ) + Px() - Py() + Px() * Px() * Px() * Px() +
               Px() * Px() * Px() * Py() + Px() * Px() * Py() * Py() +
               Px() * Py() * Py() * Py() - Py() * Py() * Py() * Py();
    else
        return cst( 1.0 ) + Px() - Py() + Px() * Px() * Px() * Px() * Px() +
               Px() * Px() * Px() * Px() * Py() + Px() * Px() * Px() * Py() * Py() +
               Px() * Px() * Py() * Py() * Py() + Px() * Py() * Py() * Py() * Py() -
               Py() * Py() * Py() * Py() * Py();
}

template<int Order>
auto
homogeneousExpression()
{
    if constexpr ( Order == 0 )
        return cst( 1.0 );
    else if constexpr ( Order == 1 )
        return Px() + Py();
    else if constexpr ( Order == 2 )
        return Px() * Px() + Px() * Py() + Py() * Py();
    else if constexpr ( Order == 3 )
        return Px() * Px() * Px() + Px() * Px() * Py() + Px() * Py() * Py() +
               Py() * Py() * Py();
    else
        return Px() * Px() * Px() * Px() + Px() * Px() * Px() * Py() +
               Px() * Px() * Py() * Py() + Px() * Py() * Py() * Py() +
               Py() * Py() * Py() * Py();
}

/** Discrete operators assembled while checking one affine sequence. */
struct affine_sequence_operators
{
    Eigen::MatrixXd gradient;
    Eigen::MatrixXd curl;
};

/**
 * Check one static- or runtime-order two-dimensional NED1 de Rham sequence.
 *
 * The H(curl) space is injected so that both order policies exercise exactly
 * the same mathematical checks and operator construction.
 */
template<int Order, typename HCurlSpacePtrType>
affine_sequence_operators
checkAffineSequence( std::shared_ptr<affine_mesh_type> const& mesh,
                     HCurlSpacePtrType const& hcurl,
                     char const* orderPolicy )
{
    auto h1 = Pch<Order + 1>( mesh );
    auto l2 = Pdh<Order>( mesh );
    affine_sequence_operators operators;

    BOOST_TEST_CONTEXT( "2D affine " << orderPolicy << " de Rham order " << Order )
    {
        // Euler--Poincare identity for
        // 0 -> R -> P_{k+1} -> NED1_k -> P_k^disc -> 0.
        BOOST_CHECK_EQUAL( h1->nDof() - hcurl->nDof() + l2->nDof(), 1 );

        auto potential = interpolate( h1, mesh, potentialExpression<Order>() );
        auto gradient = interpolate( hcurl, mesh, trans( gradv( potential ) ) );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = idv( gradient ) - trans( gradv( potential ) ),
                                   _quad = _Q<14>() ), 2e-10 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = curlxv( gradient ),
                                   _quad = _Q<14>() ), 2e-10 );

        // (-y s_k, x s_k) is the homogeneous enrichment of NED1_k and
        // curl(-y s_k, x s_k) = (k+2) s_k for homogeneous s_k.
        auto homogeneous = homogeneousExpression<Order>();
        auto field = interpolate( hcurl, mesh, vec( -Py() * homogeneous, Px() * homogeneous ) );
        auto expectedCurl = static_cast<double>( Order + 2 ) * homogeneous;
        auto discreteCurl = interpolate( l2, mesh, curlxv( field ) );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = curlxv( field ) - expectedCurl,
                                   _quad = _Q<14>() ), 2e-10 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = idv( discreteCurl ) - expectedCurl,
                                   _quad = _Q<14>() ), 2e-10 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = curlxv( field ) - idv( discreteCurl ),
                                   _quad = _Q<14>() ), 2e-10 );

        // The canonical edge transform must make the tangential trace
        // single-valued independently of adjacent-cell edge orientation.
        auto traceProbe = hcurl->element();
        auto tangentialJump = form1( _test = hcurl );
        tangentialJump = integrate(
            _range = internalfaces( mesh ),
            _expr = leftface( trans( id( traceProbe ) ) * vec( Tx(), Ty() ) ) +
                    rightface( -trans( id( traceProbe ) ) * vec( Tx(), Ty() ) ) );
        BOOST_CHECK_SMALL( tangentialJump( field ), 2e-10 );

        // Full coefficient-by-coefficient incidence assembly grows quadratically
        // with the global number of dofs. Orders zero and one exercise the
        // matrix complex; the field-level contracts above retain orders 0--4.
        if ( Environment::worldComm().localSize() == 1 && Order <= 1 )
        {
            operators.gradient = Eigen::MatrixXd::Zero( hcurl->nDof(), h1->nDof() );
            for ( size_type column = 0; column < h1->nDof(); ++column )
            {
                auto h1BasisFunction = h1->element();
                h1BasisFunction.zero();
                h1BasisFunction( column ) = 1.0;
                auto image = interpolate( hcurl, mesh, trans( gradv( h1BasisFunction ) ) );
                for ( size_type row = 0; row < hcurl->nDof(); ++row )
                    operators.gradient( row, column ) = image( row );
            }

            operators.curl = Eigen::MatrixXd::Zero( l2->nDof(), hcurl->nDof() );
            for ( size_type column = 0; column < hcurl->nDof(); ++column )
            {
                auto hcurlBasisFunction = hcurl->element();
                hcurlBasisFunction.zero();
                hcurlBasisFunction( column ) = 1.0;
                auto image = interpolate( l2, mesh, curlxv( hcurlBasisFunction ) );
                for ( size_type row = 0; row < l2->nDof(); ++row )
                    operators.curl( row, column ) = image( row );
            }

            Eigen::FullPivLU<Eigen::MatrixXd> gradientRank( operators.gradient );
            Eigen::FullPivLU<Eigen::MatrixXd> curlRank( operators.curl );
            gradientRank.setThreshold( 1e-9 );
            curlRank.setThreshold( 1e-9 );
            BOOST_CHECK_EQUAL( gradientRank.rank(), h1->nDof() - 1 );
            BOOST_CHECK_EQUAL( curlRank.rank(), l2->nDof() );
            BOOST_CHECK_SMALL( ( operators.curl * operators.gradient ).norm(), 2e-9 );
            BOOST_CHECK_EQUAL( hcurl->nDof() - gradientRank.rank(), curlRank.rank() );
        }
    }
    return operators;
}

/** Check static/runtime NED1 sequence and discrete-operator parity. */
template<int Order>
void
checkAffineStaticRuntimeSequence( std::shared_ptr<affine_mesh_type> const& mesh )
{
    auto hcurlStatic = Neh<Order>( mesh );
    auto hcurlRuntime = Neh<Dynamic>( mesh, RuntimeOrder{ Order } );

    BOOST_TEST_CONTEXT( "2D affine NED1 static/runtime contract order " << Order )
    {
        BOOST_CHECK_EQUAL( hcurlRuntime->basis()->runtimeOrder(), Order );
        BOOST_CHECK_EQUAL( hcurlRuntime->basis()->localDofCount(),
                           hcurlStatic->basis()->localDofCount() );
        BOOST_CHECK_EQUAL( hcurlRuntime->nDof(), hcurlStatic->nDof() );
    }

    auto const staticOperators = checkAffineSequence<Order>( mesh, hcurlStatic, "static" );
    auto const runtimeOperators = checkAffineSequence<Order>( mesh, hcurlRuntime, "runtime" );
    if ( Environment::worldComm().localSize() == 1 && Order <= 1 )
    {
        BOOST_TEST_CONTEXT( "2D affine NED1 operator parity order " << Order )
        {
            BOOST_CHECK_EQUAL( runtimeOperators.gradient.rows(), staticOperators.gradient.rows() );
            BOOST_CHECK_EQUAL( runtimeOperators.gradient.cols(), staticOperators.gradient.cols() );
            BOOST_CHECK_EQUAL( runtimeOperators.curl.rows(), staticOperators.curl.rows() );
            BOOST_CHECK_EQUAL( runtimeOperators.curl.cols(), staticOperators.curl.cols() );
            BOOST_CHECK_SMALL( ( runtimeOperators.gradient - staticOperators.gradient ).norm(), 2e-9 );
            BOOST_CHECK_SMALL( ( runtimeOperators.curl - staticOperators.curl ).norm(), 2e-9 );
        }
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE( de_rham_suite )

/** Verify the affine two-dimensional static sequence for orders zero to four. */
BOOST_AUTO_TEST_CASE( affine_exact_sequence_and_commuting_interpolation )
{
    auto mesh = makeAffineMesh();
    checkAffineStaticRuntimeSequence<0>( mesh );
    checkAffineStaticRuntimeSequence<1>( mesh );
    checkAffineStaticRuntimeSequence<2>( mesh );
    checkAffineStaticRuntimeSequence<3>( mesh );
    checkAffineSequence<4>( mesh, Neh<4>( mesh ), "static" );

    auto rt = RTh<0>( mesh );
    auto rtRuntime = RTh<Dynamic>( mesh, RuntimeOrder{ 0 } );
    auto bdm = BDMh<0>( mesh );
    auto bdmRuntime = BDMh<Dynamic>( mesh, RuntimeOrder{ 0 } );
    auto fluxExpr = vec( Px(), Py() );
    auto rtFlux = interpolate( rt, mesh, fluxExpr );
    auto rtFluxRuntime = interpolate( rtRuntime, mesh, fluxExpr );
    auto bdmFlux = interpolate( bdm, mesh, fluxExpr );
    auto bdmFluxRuntime = interpolate( bdmRuntime, mesh, fluxExpr );

    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = divv( rtFlux ) - cst( 2.0 ), _quad = _Q<10>() ), 1e-11 );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = divv( bdmFlux ) - cst( 2.0 ), _quad = _Q<10>() ), 1e-11 );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = idv( rtFluxRuntime ) - idv( rtFlux ), _quad = _Q<10>() ), 1e-11 );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ), _expr = idv( bdmFluxRuntime ) - idv( bdmFlux ), _quad = _Q<10>() ), 1e-11 );
}

/** Verify that true runtime NED1 construction has no order-three dispatch ceiling. */
BOOST_AUTO_TEST_CASE( affine_runtime_ned1_construction_above_validation_envelope )
{
    auto mesh = makeAffineMesh();
    auto hcurl4 = Neh<Dynamic>( mesh, RuntimeOrder{ 4 } );
    auto hcurl6 = Neh<Dynamic>( mesh, RuntimeOrder{ 6 } );

    BOOST_REQUIRE( hcurl4 );
    BOOST_REQUIRE( hcurl6 );
    BOOST_CHECK_EQUAL( hcurl4->basis()->runtimeOrder(), 4 );
    BOOST_CHECK_EQUAL( hcurl6->basis()->runtimeOrder(), 6 );
    BOOST_CHECK_EQUAL( hcurl4->basis()->localDofCount(), 35 );
    BOOST_CHECK_EQUAL( hcurl6->basis()->localDofCount(), 63 );
    BOOST_CHECK_GT( hcurl4->nDof(), 0 );
    BOOST_CHECK_GT( hcurl6->nDof(), hcurl4->nDof() );
}

/** Verify the enabled contravariant Piola identities on a genuinely curved mesh. */
BOOST_AUTO_TEST_CASE( curved_simplex_piola_contracts )
{
    // This contract targets the curved geometric/Piola mapping itself. Keep it
    // serial: distributed curved H(div) conformity belongs to a dedicated MPI
    // contract and depends on the partition/interface geometry representation.
    if ( Environment::worldComm().localSize() > 1 )
    {
        BOOST_TEST_MESSAGE( "curved Piola contract is serial-only" );
        return;
    }

    auto mesh = unitCircle<2>( 0.35 );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_GT( mesh->numGlobalElements(), 0 );

    auto fluxExpr = vec( cst( 1.0 ) + Px() * Px() - Py(),
                         cst( -0.5 ) + Px() + Py() * Py() );
    auto rt = RTh<0>( mesh );
    auto rtRuntime = RTh<Dynamic>( mesh, RuntimeOrder{ 0 } );
    auto bdm = BDMh<0>( mesh );
    auto bdmRuntime = BDMh<Dynamic>( mesh, RuntimeOrder{ 0 } );
    auto rtFlux = interpolate( rt, mesh, fluxExpr );
    auto rtFluxRuntime = interpolate( rtRuntime, mesh, fluxExpr );
    auto bdmFlux = interpolate( bdm, mesh, fluxExpr );
    auto bdmFluxRuntime = interpolate( bdmRuntime, mesh, fluxExpr );
    auto tangentExpr = vec( cst( 0.25 ) - Py() + Px() * Py(),
                            cst( -0.75 ) + Px() + Px() * Px() );
    auto ned1 = Ned1h<0>( mesh );
    auto ned1Runtime = Ned1h<Dynamic>( mesh, RuntimeOrder{ 0 } );
    auto ned1Field = interpolate( ned1, mesh, tangentExpr );
    auto ned1FieldRuntime = interpolate( ned1Runtime, mesh, tangentExpr );
    auto h1 = Pch<1>( mesh );
    auto potential = interpolate(
        h1, mesh, cst( 0.5 ) + Px() - 0.75 * Py() );
    auto discreteGradient = interpolate(
        ned1Runtime, mesh, trans( gradv( potential ) ) );

    BOOST_CHECK_SMALL( divergenceTheoremResidual( mesh, rtFlux ), 2e-8 );
    BOOST_CHECK_SMALL( divergenceTheoremResidual( mesh, rtFluxRuntime ), 2e-8 );
    BOOST_CHECK_SMALL( divergenceTheoremResidual( mesh, bdmFlux ), 2e-8 );
    BOOST_CHECK_SMALL( divergenceTheoremResidual( mesh, bdmFluxRuntime ), 2e-8 );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = idv( rtFluxRuntime ) - idv( rtFlux ),
                               _quad = _Q<16>(),
                               _geomap = GeomapStrategyType::GEOMAP_HO ), 2e-9 );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = idv( bdmFluxRuntime ) - idv( bdmFlux ),
                               _quad = _Q<16>(),
                               _geomap = GeomapStrategyType::GEOMAP_HO ), 2e-9 );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = idv( ned1FieldRuntime ) - idv( ned1Field ),
                               _quad = _Q<16>(),
                               _geomap = GeomapStrategyType::GEOMAP_HO ), 2e-9 );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = curlxv( ned1FieldRuntime ) -
                                       curlxv( ned1Field ),
                               _quad = _Q<16>(),
                               _geomap = GeomapStrategyType::GEOMAP_HO ), 2e-9 );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = curlxv( discreteGradient ),
                               _quad = _Q<16>(),
                               _geomap = GeomapStrategyType::GEOMAP_HO ), 2e-9 );

    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = gradv( ned1FieldRuntime ) -
                                       gradv( ned1Field ),
                               _quad = _Q<16>(),
                               _geomap = GeomapStrategyType::GEOMAP_HO ), 2e-8 );
    BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                               _expr = gradv( rtFluxRuntime ) -
                                       gradv( rtFlux ),
                               _quad = _Q<16>(),
                               _geomap = GeomapStrategyType::GEOMAP_HO ), 2e-8 );
}

BOOST_AUTO_TEST_SUITE_END()
