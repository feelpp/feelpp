/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_de_rham_3d.cpp
 * @brief Three-dimensional affine de Rham contracts for static/runtime order.
 */
#define BOOST_TEST_MODULE test_de_rham_3d
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/neh.hpp>
#include <feel/feeldiscr/ned2h.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using affine_mesh_type = Mesh<Simplex<3, 1>>;

/** Build one affine tetrahedron for local incidence-matrix validation. */
std::shared_ptr<affine_mesh_type>
makeOneElementAffineMesh()
{
    auto description = std::make_shared<GmshSimplexDomain>(
        3, 1, Gmsh::GMSH_REFERENCE_DOMAIN );
    description->setPrefix( "de-rham-3d-one-element" );
    description->setCharacteristicLength( 2.0 );
    return createGMSHMesh( _mesh = new affine_mesh_type,
                           _desc = description,
                           _h = description->h() );
}

template<typename SpacePtrType, typename ExprType>
auto
interpolate( SpacePtrType const& space,
             std::shared_ptr<affine_mesh_type> const& mesh,
             ExprType const& expr )
{
    auto result = space->element();
    result.on( _range = elements( mesh ), _expr = expr );
    return result;
}

/** Return a scalar polynomial of exact total degree @p Order. */
template<int Order>
auto
homogeneousExpression()
{
    if constexpr ( Order == 0 )
        return cst( 1.0 );
    else if constexpr ( Order == 1 )
        return Px() + 2.0 * Py() - Pz();
    else if constexpr ( Order == 2 )
        return Px() * Px() + Px() * Py() + Py() * Pz() + Pz() * Pz();
    else
        return Px() * Px() * Px() + Px() * Px() * Py() + Px() * Py() * Pz() +
               Py() * Pz() * Pz() + Pz() * Pz() * Pz();
}

/** Return a scalar potential of exact total degree @p Order + 1. */
template<int Order>
auto
potentialExpression()
{
    if constexpr ( Order == 0 )
        return cst( 1.0 ) + 2.0 * Px() - 3.0 * Py() + 4.0 * Pz();
    else if constexpr ( Order == 1 )
        return cst( 1.0 ) + Px() - Py() + Pz() + Px() * Px() + Px() * Py() -
               Py() * Pz() + Pz() * Pz();
    else if constexpr ( Order == 2 )
        return cst( 1.0 ) + Px() - Py() + Pz() + Px() * Px() * Px() +
               Px() * Px() * Py() + Px() * Py() * Pz() - Py() * Pz() * Pz() +
               Pz() * Pz() * Pz();
    else
        return cst( 1.0 ) + Px() - Py() + Pz() + Px() * Px() * Px() * Px() +
               Px() * Px() * Px() * Py() + Px() * Px() * Py() * Pz() +
               Px() * Py() * Pz() * Pz() - Py() * Pz() * Pz() * Pz() +
               Pz() * Pz() * Pz() * Pz();
}

/**
 * Return a field in the homogeneous enrichment of NED1_Order.
 *
 * For homogeneous s_Order, x cross (s_Order,0,0) equals
 * (0,z s_Order,-y s_Order) and is orthogonal to x.
 */
template<int Order>
auto
nedelecEnrichmentExpression()
{
    auto s = homogeneousExpression<Order>();
    return vec( cst( 0.0 ), Pz() * s, -Py() * s );
}

/** Return the analytic curl of nedelecEnrichmentExpression<Order>(). */
template<int Order>
auto
nedelecEnrichmentCurlExpression()
{
    if constexpr ( Order == 0 )
        return vec( cst( -2.0 ), cst( 0.0 ), cst( 0.0 ) );
    else if constexpr ( Order == 1 )
        return vec( -2.0 * Px() - 6.0 * Py() + 3.0 * Pz(),
                    Py(),
                    Pz() );
    else if constexpr ( Order == 2 )
        return vec( -2.0 * Px() * Px() - 3.0 * Px() * Py() -
                        4.0 * Py() * Pz() - 4.0 * Pz() * Pz(),
                    2.0 * Px() * Py() + Py() * Py(),
                    2.0 * Px() * Pz() + Py() * Pz() );
    else
        return vec( -2.0 * Px() * Px() * Px() - 3.0 * Px() * Px() * Py() -
                        4.0 * Px() * Py() * Pz() - 5.0 * Py() * Pz() * Pz() -
                        5.0 * Pz() * Pz() * Pz(),
                    3.0 * Px() * Px() * Py() + 2.0 * Px() * Py() * Py() +
                        Py() * Py() * Pz(),
                    3.0 * Px() * Px() * Pz() + 2.0 * Px() * Py() * Pz() +
                        Py() * Pz() * Pz() );
}

/** Return the homogeneous Raviart--Thomas enrichment x s_Order. */
template<int Order>
auto
raviartThomasEnrichmentExpression()
{
    auto s = homogeneousExpression<Order>();
    return vec( Px() * s, Py() * s, Pz() * s );
}

/** Discrete 3D incidence operators for one order policy. */
struct affine_sequence_operators
{
    Eigen::MatrixXd gradient;
    Eigen::MatrixXd curl;
    Eigen::MatrixXd divergence;
};

/**
 * Assemble and validate the three local incidence operators of one sequence.
 *
 * A one-tetrahedron mesh makes coefficient-level checks affordable through
 * order three. Multi-tetrahedron conformity is checked independently by the
 * field-level kernel below.
 */
template<int H1Order, int L2Order,
         typename HCurlSpacePtrType, typename HDivSpacePtrType>
affine_sequence_operators
assembleAffineSequenceOperators( std::shared_ptr<affine_mesh_type> const& mesh,
                                 HCurlSpacePtrType const& hcurl,
                                 HDivSpacePtrType const& hdiv )
{
    auto h1 = Pch<Dynamic>( mesh, RuntimeOrder{ H1Order } );
    auto l2 = Pdh<Dynamic>( mesh, RuntimeOrder{ L2Order } );
    affine_sequence_operators operators;

    operators.gradient = Eigen::MatrixXd::Zero( hcurl->nDof(), h1->nDof() );
    for ( size_type column = 0; column < h1->nDof(); ++column )
    {
        auto basisFunction = h1->element();
        basisFunction.zero();
        basisFunction( column ) = 1.0;
        auto image = interpolate( hcurl, mesh, trans( gradv( basisFunction ) ) );
        for ( size_type row = 0; row < hcurl->nDof(); ++row )
            operators.gradient( row, column ) = image( row );
    }

    operators.curl = Eigen::MatrixXd::Zero( hdiv->nDof(), hcurl->nDof() );
    for ( size_type column = 0; column < hcurl->nDof(); ++column )
    {
        auto basisFunction = hcurl->element();
        basisFunction.zero();
        basisFunction( column ) = 1.0;
        auto image = interpolate( hdiv, mesh, curlv( basisFunction ) );
        for ( size_type row = 0; row < hdiv->nDof(); ++row )
            operators.curl( row, column ) = image( row );
    }

    operators.divergence = Eigen::MatrixXd::Zero( l2->nDof(), hdiv->nDof() );
    for ( size_type column = 0; column < hdiv->nDof(); ++column )
    {
        auto basisFunction = hdiv->element();
        basisFunction.zero();
        basisFunction( column ) = 1.0;
        auto image = interpolate( l2, mesh, divv( basisFunction ) );
        for ( size_type row = 0; row < l2->nDof(); ++row )
            operators.divergence( row, column ) = image( row );
    }

    Eigen::FullPivLU<Eigen::MatrixXd> gradientRank( operators.gradient );
    Eigen::FullPivLU<Eigen::MatrixXd> curlRank( operators.curl );
    Eigen::FullPivLU<Eigen::MatrixXd> divergenceRank( operators.divergence );
    gradientRank.setThreshold( 1e-9 );
    curlRank.setThreshold( 1e-9 );
    divergenceRank.setThreshold( 1e-9 );

    BOOST_CHECK_EQUAL( gradientRank.rank(), h1->nDof() - 1 );
    BOOST_CHECK_EQUAL( curlRank.rank(), hcurl->nDof() - gradientRank.rank() );
    BOOST_CHECK_EQUAL( divergenceRank.rank(), l2->nDof() );
    BOOST_CHECK_EQUAL( hdiv->nDof() - curlRank.rank(), divergenceRank.rank() );
    BOOST_CHECK_SMALL( ( operators.curl * operators.gradient ).norm(), 2e-8 );
    BOOST_CHECK_SMALL( ( operators.divergence * operators.curl ).norm(), 2e-8 );
    return operators;
}

/**
 * Validate one affine tetrahedral first-kind de Rham sequence.
 *
 * The injected H(curl) and H(div) spaces let compile-time and runtime order
 * exercise the same reproduction, commuting, trace, and operator kernels.
 */
template<int Order, typename HCurlSpacePtrType, typename HDivSpacePtrType>
affine_sequence_operators
checkAffineSequence( std::shared_ptr<affine_mesh_type> const& mesh,
                     HCurlSpacePtrType const& hcurl,
                     HDivSpacePtrType const& hdiv,
                     char const* orderPolicy )
{
    auto h1 = Pch<Dynamic>( mesh, RuntimeOrder{ Order+1 } );
    auto l2 = Pdh<Dynamic>( mesh, RuntimeOrder{ Order } );
    affine_sequence_operators operators;

    BOOST_TEST_CONTEXT( "3D affine " << orderPolicy << " de Rham order " << Order )
    {
        // 0 -> R -> P_{k+1} -> NED1_k -> RT_k -> P_k^disc -> 0.
        auto const eulerCharacteristic =
            static_cast<int64_type>( h1->nDof() ) -
            static_cast<int64_type>( hcurl->nDof() ) +
            static_cast<int64_type>( hdiv->nDof() ) -
            static_cast<int64_type>( l2->nDof() );
        BOOST_CHECK_EQUAL( eulerCharacteristic, 1 );

        auto potential = interpolate( h1, mesh, potentialExpression<Order>() );
        auto gradient = interpolate( hcurl, mesh, trans( gradv( potential ) ) );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = idv( gradient ) - trans( gradv( potential ) ),
                                   _quad = _Q<14>() ), 5e-10 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = curlv( gradient ),
                                   _quad = _Q<14>() ), 5e-10 );

        auto hcurlField = interpolate( hcurl, mesh, nedelecEnrichmentExpression<Order>() );
        auto expectedCurl = nedelecEnrichmentCurlExpression<Order>();
        auto hdivCurl = interpolate( hdiv, mesh, curlv( hcurlField ) );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = curlv( hcurlField ) - expectedCurl,
                                   _quad = _Q<14>() ), 2e-9 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = idv( hdivCurl ) - expectedCurl,
                                   _quad = _Q<14>() ), 2e-9 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = curlv( hcurlField ) - idv( hdivCurl ),
                                   _quad = _Q<14>() ), 2e-9 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = divv( hdivCurl ),
                                   _quad = _Q<14>() ), 2e-9 );

        auto hdivField = interpolate( hdiv, mesh, raviartThomasEnrichmentExpression<Order>() );
        auto expectedDivergence = static_cast<double>( Order + 3 ) *
                                  homogeneousExpression<Order>();
        auto l2Divergence = interpolate( l2, mesh, divv( hdivField ) );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = divv( hdivField ) - expectedDivergence,
                                   _quad = _Q<14>() ), 2e-9 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = idv( l2Divergence ) - expectedDivergence,
                                   _quad = _Q<14>() ), 2e-9 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = divv( hdivField ) - idv( l2Divergence ),
                                   _quad = _Q<14>() ), 2e-9 );

        // Multi-tetrahedron normal conformity exercises the RT face
        // transformations with the mesh's actual local numbering. Nedelec
        // edge and dense-face transforms are already covered above by exact
        // reproduction of a globally smooth field and its curl.
        auto hdivProbe = hdiv->element();
        auto normalJump = form1( _test = hdiv );
        normalJump = integrate(
            _range = internalfaces( mesh ),
            _expr = leftface( trans( id( hdivProbe ) ) * N() ) +
                    rightface( -trans( id( hdivProbe ) ) * N() ) );
        BOOST_CHECK_SMALL( normalJump( hdivField ), 2e-9 );

        // Keep one assembled multi-cell low-order check in addition to the
        // complete one-cell order matrix exercised below.
        if ( Environment::worldComm().localSize() == 1 && Order == 0 )
            operators = assembleAffineSequenceOperators<Order+1,Order>( mesh, hcurl, hdiv );
    }
    return operators;
}

/** Validate static/runtime spaces and coefficient-level incidence parity. */
template<int Order>
void
checkAffineStaticRuntimeSequence( std::shared_ptr<affine_mesh_type> const& mesh )
{
    auto hcurlStatic = Neh<Order>( mesh );
    auto hcurlRuntime = Neh<Dynamic>( mesh, RuntimeOrder{ Order } );
    auto hdivStatic = RTh<Order>( mesh );
    auto hdivRuntime = RTh<Dynamic>( mesh, RuntimeOrder{ Order } );

    BOOST_TEST_CONTEXT( "3D affine static/runtime space contract order " << Order )
    {
        BOOST_CHECK_EQUAL( hcurlRuntime->basis()->runtimeOrder(), Order );
        BOOST_CHECK_EQUAL( hdivRuntime->basis()->runtimeOrder(), Order );
        BOOST_CHECK_EQUAL( hcurlRuntime->basis()->localDofCount(),
                           hcurlStatic->basis()->localDofCount() );
        BOOST_CHECK_EQUAL( hdivRuntime->basis()->localDofCount(),
                           hdivStatic->basis()->localDofCount() );
        BOOST_CHECK_EQUAL( hcurlRuntime->nDof(), hcurlStatic->nDof() );
        BOOST_CHECK_EQUAL( hdivRuntime->nDof(), hdivStatic->nDof() );
    }

    auto const staticOperators =
        checkAffineSequence<Order>( mesh, hcurlStatic, hdivStatic, "static" );
    auto const runtimeOperators =
        checkAffineSequence<Order>( mesh, hcurlRuntime, hdivRuntime, "runtime" );

    if ( Environment::worldComm().localSize() == 1 && Order == 0 )
    {
        BOOST_TEST_CONTEXT( "3D affine operator parity order " << Order )
        {
            BOOST_REQUIRE_EQUAL( runtimeOperators.gradient.rows(),
                                 staticOperators.gradient.rows() );
            BOOST_REQUIRE_EQUAL( runtimeOperators.gradient.cols(),
                                 staticOperators.gradient.cols() );
            BOOST_REQUIRE_EQUAL( runtimeOperators.curl.rows(), staticOperators.curl.rows() );
            BOOST_REQUIRE_EQUAL( runtimeOperators.curl.cols(), staticOperators.curl.cols() );
            BOOST_REQUIRE_EQUAL( runtimeOperators.divergence.rows(),
                                 staticOperators.divergence.rows() );
            BOOST_REQUIRE_EQUAL( runtimeOperators.divergence.cols(),
                                 staticOperators.divergence.cols() );
            BOOST_CHECK_SMALL( ( runtimeOperators.gradient - staticOperators.gradient ).norm(),
                               2e-8 );
            BOOST_CHECK_SMALL( ( runtimeOperators.curl - staticOperators.curl ).norm(), 2e-8 );
            BOOST_CHECK_SMALL( ( runtimeOperators.divergence - staticOperators.divergence ).norm(),
                               2e-8 );
        }
    }
}

/** Run the global sequence through true runtime-order spaces only. */
template<int Order>
void
checkAffineRuntimeSequence( std::shared_ptr<affine_mesh_type> const& mesh )
{
    auto hcurl = Neh<Dynamic>( mesh, RuntimeOrder{ Order } );
    auto hdiv = RTh<Dynamic>( mesh, RuntimeOrder{ Order } );
    BOOST_CHECK_EQUAL( hcurl->basis()->runtimeOrder(), Order );
    BOOST_CHECK_EQUAL( hdiv->basis()->runtimeOrder(), Order );
    checkAffineSequence<Order>( mesh, hcurl, hdiv, "runtime" );
}

/** Validate one-cell runtime incidence matrices through order three. */
template<int Order>
void
checkLocalRuntimeOperators( std::shared_ptr<affine_mesh_type> const& mesh )
{
    auto const operators = assembleAffineSequenceOperators<Order+1,Order>(
        mesh,
        Neh<Dynamic>( mesh, RuntimeOrder{ Order } ),
        RTh<Dynamic>( mesh, RuntimeOrder{ Order } ) );
    BOOST_CHECK_GT( operators.gradient.size(), 0 );
    BOOST_CHECK_GT( operators.curl.size(), 0 );
    BOOST_CHECK_GT( operators.divergence.size(), 0 );
}

/** Validate one runtime second-kind Nedelec sequence. */
template<int Order, typename HDivSpacePtrType>
void
checkAffineNedelec2RuntimeSequence(
    std::shared_ptr<affine_mesh_type> const& mesh,
    HDivSpacePtrType const& hdiv,
    uint16_type l2Order,
    bool assembleOperators )
{
    auto h1 = Pch<Dynamic>( mesh, RuntimeOrder{ Order+2 } );
    auto hcurl = Ned2h<Dynamic>( mesh, RuntimeOrder{ Order } );
    auto l2 = Pdh<Dynamic>( mesh, RuntimeOrder{ l2Order } );
    BOOST_TEST_CONTEXT( "3D affine runtime NED2 sequence order " << Order )
    {
        auto const eulerCharacteristic =
            static_cast<int64_type>( h1->nDof() )-
            static_cast<int64_type>( hcurl->nDof() )+
            static_cast<int64_type>( hdiv->nDof() )-
            static_cast<int64_type>( l2->nDof() );
        BOOST_CHECK_EQUAL( eulerCharacteristic, 1 );

        auto potential = interpolate( h1, mesh, potentialExpression<Order+1>() );
        auto gradient = interpolate( hcurl, mesh, trans( gradv( potential ) ) );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = idv( gradient )-trans( gradv( potential ) ),
                                   _quad = _Q<14>() ), 2e-9 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = curlv( gradient ),
                                   _quad = _Q<14>() ), 2e-9 );

        auto hcurlField = interpolate(
            hcurl, mesh, nedelecEnrichmentExpression<Order>() );
        auto expectedCurl = nedelecEnrichmentCurlExpression<Order>();
        auto hdivCurl = interpolate( hdiv, mesh, curlv( hcurlField ) );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = curlv( hcurlField )-expectedCurl,
                                   _quad = _Q<14>() ), 5e-9 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = idv( hdivCurl )-expectedCurl,
                                   _quad = _Q<14>() ), 5e-9 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = divv( hdivCurl ),
                                   _quad = _Q<14>() ), 5e-9 );

        auto hdivField = [&]()
        {
            if constexpr ( Order == 0 )
                return interpolate( hdiv, mesh, raviartThomasEnrichmentExpression<0>() );
            else
                return interpolate( hdiv, mesh,
                                    raviartThomasEnrichmentExpression<Order-1>() );
        }();
        auto expectedDivergence = [&]()
        {
            if constexpr ( Order == 0 )
                return 3.0*homogeneousExpression<0>();
            else
                return static_cast<double>( Order+2 )*
                       homogeneousExpression<Order-1>();
        }();
        auto l2Divergence = interpolate( l2, mesh, divv( hdivField ) );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = divv( hdivField )-expectedDivergence,
                                   _quad = _Q<14>() ), 5e-9 );
        BOOST_CHECK_SMALL( normL2( _range = elements( mesh ),
                                   _expr = idv( l2Divergence )-expectedDivergence,
                                   _quad = _Q<14>() ), 5e-9 );

        // The preceding exact multi-cell reproduction already exercises BDM
        // face sharing and orientation: an incorrect global face action
        // changes the reconstructed field in at least one adjacent cell.
        // Also assert the lower-level coefficient operators consumed by
        // FunctionSpace are nonsingular.  The generic internal-face VF
        // context is validated separately once it supports runtime BDM
        // cardinalities.
        for ( auto const& elementWrapper : elements( mesh ) )
        {
            auto const& element = unwrap_ref( elementWrapper );
            auto const transform =
                hdiv->dof()->localToGlobalCoefficientTransform( element.id() );
            BOOST_CHECK_EQUAL( transform.fullPivLu().rank(), transform.rows() );
        }

        if ( assembleOperators && Environment::worldComm().localSize() == 1 )
        {
            auto const operators = assembleAffineSequenceOperators<
                Order+2, ( Order == 0 ? 0 : Order-1 )>( mesh, hcurl, hdiv );
            BOOST_CHECK_GT( operators.gradient.size(), 0 );
        }
    }
}

/** Dispatch the compatible H(div) partner of runtime NED2_Order. */
template<int Order>
void
checkAffineNedelec2RuntimeSequence(
    std::shared_ptr<affine_mesh_type> const& mesh,
    bool assembleOperators = false )
{
    if constexpr ( Order == 0 )
        checkAffineNedelec2RuntimeSequence<Order>(
            mesh, RTh<Dynamic>( mesh, RuntimeOrder{ 0 } ), 0,
            assembleOperators );
    else
        checkAffineNedelec2RuntimeSequence<Order>(
            mesh, BDMh<Dynamic>( mesh, RuntimeOrder{ Order-1 } ), Order-1,
            assembleOperators );
}
} // namespace

BOOST_AUTO_TEST_SUITE( de_rham_3d_suite )

/** Validate affine tetrahedral static/runtime sequences for k=0,...,3. */
BOOST_AUTO_TEST_CASE( affine_exact_sequence_and_runtime_parity )
{
    auto mesh = unitCube( 0.8 );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_GT( mesh->numGlobalElements(), 1 );

    // One static assembled spot check anchors the separately verified local
    // static/runtime FE parity. The order sweep itself is runtime-only.
    checkAffineStaticRuntimeSequence<0>( mesh );
    checkAffineRuntimeSequence<1>( mesh );
    checkAffineRuntimeSequence<2>( mesh );
    checkAffineRuntimeSequence<3>( mesh );

    if ( Environment::worldComm().localSize() == 1 )
    {
        auto oneElementMesh = makeOneElementAffineMesh();
        BOOST_REQUIRE_EQUAL( oneElementMesh->numGlobalElements(), 1 );
        checkLocalRuntimeOperators<0>( oneElementMesh );
        checkLocalRuntimeOperators<1>( oneElementMesh );
        checkLocalRuntimeOperators<2>( oneElementMesh );
        checkLocalRuntimeOperators<3>( oneElementMesh );
    }
}

/** Validate runtime NED2/RT-BDM affine sequences for k=0,...,3. */
BOOST_AUTO_TEST_CASE( affine_nedelec2_exact_sequence_runtime_orders )
{
    auto mesh = unitCube( 0.8 );
    BOOST_REQUIRE_GT( mesh->numGlobalElements(), 1 );
    checkAffineNedelec2RuntimeSequence<0>( mesh );
    checkAffineNedelec2RuntimeSequence<1>( mesh );
    checkAffineNedelec2RuntimeSequence<2>( mesh );
    checkAffineNedelec2RuntimeSequence<3>( mesh );

    if ( Environment::worldComm().localSize() == 1 )
    {
        auto oneElementMesh = makeOneElementAffineMesh();
        checkAffineNedelec2RuntimeSequence<0>( oneElementMesh, true );
        checkAffineNedelec2RuntimeSequence<1>( oneElementMesh, true );
        checkAffineNedelec2RuntimeSequence<2>( oneElementMesh, true );
        checkAffineNedelec2RuntimeSequence<3>( oneElementMesh, true );
    }
}

/** Verify that true runtime tetrahedral construction has no order-three ceiling. */
BOOST_AUTO_TEST_CASE( affine_runtime_construction_above_validation_envelope )
{
    auto mesh = unitCube( 0.8 );
    auto hcurl = Neh<Dynamic>( mesh, RuntimeOrder{ 4 } );
    auto hcurl2 = Ned2h<Dynamic>( mesh, RuntimeOrder{ 4 } );
    auto hdiv = RTh<Dynamic>( mesh, RuntimeOrder{ 4 } );

    BOOST_REQUIRE( hcurl );
    BOOST_REQUIRE( hcurl2 );
    BOOST_REQUIRE( hdiv );
    BOOST_CHECK_EQUAL( hcurl->basis()->runtimeOrder(), 4 );
    BOOST_CHECK_EQUAL( hcurl2->basis()->runtimeOrder(), 4 );
    BOOST_CHECK_EQUAL( hdiv->basis()->runtimeOrder(), 4 );
    BOOST_CHECK_EQUAL( hcurl->basis()->localDofCount(), 140 );
    BOOST_CHECK_EQUAL( hcurl2->basis()->localDofCount(), 168 );
    BOOST_CHECK_EQUAL( hdiv->basis()->localDofCount(), 120 );
    BOOST_CHECK_GT( hcurl->nDof(), 0 );
    BOOST_CHECK_GT( hdiv->nDof(), 0 );
}

BOOST_AUTO_TEST_SUITE_END()
