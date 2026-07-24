/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_curved_piola_2d.cpp
 * @brief Curved-mesh Piola contracts for two-dimensional compatible spaces.
 */
#define BOOST_TEST_MODULE test_curved_piola_2d
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/bdmh.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/ned2h.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelvf/vf.hpp>

#include <string_view>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Interpolate @p expression in @p space on all mesh cells. */
template <typename SpacePtrType, typename MeshPtrType, typename ExprType>
auto interpolateField( SpacePtrType const& space,
                       MeshPtrType const& mesh,
                       ExprType const& expression )
{
    auto field = space->element();
    field.on( _range = elements( mesh ), _expr = expression );
    return field;
}

/** Return the absolute residual in the divergence theorem. */
template <typename MeshPtrType, typename ElementType>
double
divergenceTheoremResidual( MeshPtrType const& mesh, ElementType const& field )
{
    auto const volumeFlux =
        integrate( _range = elements( mesh ),
                   _expr = divv( field ),
                   _quad = _Q<18>(),
                   _geomap = GeomapStrategyType::GEOMAP_HO )
            .evaluate()( 0, 0 );
    auto const boundaryFlux =
        integrate( _range = boundaryfaces( mesh ),
                   _expr = trans( N() ) * idv( field ),
                   _quad = _Q<18>(),
                   _geomap = GeomapStrategyType::GEOMAP_HO )
            .evaluate()( 0, 0 );
    return std::abs( volumeFlux - boundaryFlux );
}

/** Check the normal trace shared by adjacent cells. */
template <typename SpacePtrType, typename ElementType, typename MeshPtrType>
void checkNormalTrace( SpacePtrType const& space,
                       ElementType const& field,
                       MeshPtrType const& mesh )
{
    auto probe = space->element();
    auto traceForm = form1( _test = space );
    traceForm = integrate(
        _range = internalfaces( mesh ),
        _expr = leftface( trans( id( probe ) ) * N() ) +
                rightface( -trans( id( probe ) ) * N() ),
        _quad = _Q<18>(),
        _geomap = GeomapStrategyType::GEOMAP_HO );
    BOOST_CHECK_SMALL( traceForm( field ), 5e-8 );
}

/** Check the scalar tangential trace shared by adjacent triangles. */
template <typename SpacePtrType, typename ElementType, typename MeshPtrType>
void checkTangentialTrace( SpacePtrType const& space,
                           ElementType const& field,
                           MeshPtrType const& mesh )
{
    auto probe = space->element();
    auto traceForm = form1( _test = space );
    traceForm = integrate(
        _range = internalfaces( mesh ),
        _expr = leftface( trans( id( probe ) ) * vec( Tx(), Ty() ) ) +
                rightface( -trans( id( probe ) ) * vec( Tx(), Ty() ) ),
        _quad = _Q<18>(),
        _geomap = GeomapStrategyType::GEOMAP_HO );
    BOOST_CHECK_SMALL( traceForm( field ), 5e-8 );
}

/** Assemble the H(div) graph inner product and require positive field energy. */
template <typename SpacePtrType, typename ElementType, typename MeshPtrType>
void checkHdivAssembly( SpacePtrType const& space,
                        ElementType const& field,
                        MeshPtrType const& mesh )
{
    auto trial = space->element();
    auto test = space->element();
    auto graph = form2( _trial = space, _test = space );
    graph = integrate( _range = elements( mesh ),
                       _expr = inner( idt( trial ), id( test ) ) +
                               divt( trial ) * div( test ),
                       _quad = _Q<18>(),
                       _geomap = GeomapStrategyType::GEOMAP_HO );
    BOOST_CHECK_GT( graph( field, field ), 0.0 );
}

/** Assemble the H(curl) graph inner product and require positive field energy. */
template <typename SpacePtrType, typename ElementType, typename MeshPtrType>
void checkHcurlAssembly( SpacePtrType const& space,
                         ElementType const& field,
                         MeshPtrType const& mesh )
{
    auto trial = space->element();
    auto test = space->element();
    auto graph = form2( _trial = space, _test = space );
    graph = integrate( _range = elements( mesh ),
                       _expr = inner( idt( trial ), id( test ) ) +
                               curlt( trial ) * curl( test ),
                       _quad = _Q<18>(),
                       _geomap = GeomapStrategyType::GEOMAP_HO );
    BOOST_CHECK_GT( graph( field, field ), 0.0 );
}

/** Validate one static/runtime H(div) family pair on a curved mesh. */
template <int Order, typename StaticFactory, typename RuntimeFactory,
          typename MeshPtrType>
void checkHdivFamily( StaticFactory&& makeStatic,
                      RuntimeFactory&& makeRuntime,
                      MeshPtrType const& mesh,
                      std::string_view family )
{
    auto staticSpace = makeStatic();
    auto runtimeSpace = makeRuntime();
    auto expression = vec( cst( 1.0 ) + Px() * Px() - Py() + Px() * Py(),
                           cst( -0.5 ) + Px() + Py() * Py() - Px() * Py() );
    auto staticField = interpolateField( staticSpace, mesh, expression );
    auto runtimeField = interpolateField( runtimeSpace, mesh, expression );

    BOOST_TEST_CONTEXT( family << " order " << Order )
    {
        BOOST_CHECK_EQUAL( runtimeSpace->basis()->runtimeOrder(), Order );
        BOOST_CHECK_EQUAL( runtimeSpace->basis()->localDofCount(),
                           staticSpace->basis()->localDofCount() );
        BOOST_CHECK_EQUAL( runtimeSpace->nDof(), staticSpace->nDof() );
        BOOST_CHECK_SMALL(
            normL2( _range = elements( mesh ),
                    _expr = idv( runtimeField ) - idv( staticField ),
                    _quad = _Q<18>(),
                    _geomap = GeomapStrategyType::GEOMAP_HO ),
            2e-8 );
        BOOST_CHECK_SMALL(
            normL2( _range = elements( mesh ),
                    _expr = divv( runtimeField ) - divv( staticField ),
                    _quad = _Q<18>(),
                    _geomap = GeomapStrategyType::GEOMAP_HO ),
            2e-8 );
        BOOST_CHECK_SMALL(
            normL2( _range = elements( mesh ),
                    _expr = gradv( runtimeField ) - gradv( staticField ),
                    _quad = _Q<18>(),
                    _geomap = GeomapStrategyType::GEOMAP_HO ),
            5e-8 );
        BOOST_CHECK_SMALL( divergenceTheoremResidual( mesh, staticField ), 5e-8 );
        BOOST_CHECK_SMALL( divergenceTheoremResidual( mesh, runtimeField ), 5e-8 );
        checkNormalTrace( runtimeSpace, runtimeField, mesh );
        checkHdivAssembly( runtimeSpace, runtimeField, mesh );
    }
}

/** Validate one static/runtime H(curl) family pair on a curved mesh. */
template <int Order, typename StaticFactory, typename RuntimeFactory,
          typename MeshPtrType>
void checkHcurlFamily( StaticFactory&& makeStatic,
                       RuntimeFactory&& makeRuntime,
                       MeshPtrType const& mesh,
                       std::string_view family )
{
    auto staticSpace = makeStatic();
    auto runtimeSpace = makeRuntime();
    auto expression = vec( cst( 0.25 ) - Py() + Px() * Py(),
                           cst( -0.75 ) + Px() + Px() * Px() );
    auto staticField = interpolateField( staticSpace, mesh, expression );
    auto runtimeField = interpolateField( runtimeSpace, mesh, expression );
    auto h1 = Pch<Dynamic>( mesh, RuntimeOrder{ Order + 1 } );
    auto potential = interpolateField(
        h1, mesh, cst( 0.5 ) + Px() - 0.75 * Py() + 0.2 * Px() * Py() );
    auto discreteGradient = interpolateField(
        runtimeSpace, mesh, trans( gradv( potential ) ) );

    BOOST_TEST_CONTEXT( family << " order " << Order )
    {
        BOOST_CHECK_EQUAL( runtimeSpace->basis()->runtimeOrder(), Order );
        BOOST_CHECK_EQUAL( runtimeSpace->basis()->localDofCount(),
                           staticSpace->basis()->localDofCount() );
        BOOST_CHECK_EQUAL( runtimeSpace->nDof(), staticSpace->nDof() );
        BOOST_CHECK_SMALL(
            normL2( _range = elements( mesh ),
                    _expr = idv( runtimeField ) - idv( staticField ),
                    _quad = _Q<18>(),
                    _geomap = GeomapStrategyType::GEOMAP_HO ),
            2e-8 );
        BOOST_CHECK_SMALL(
            normL2( _range = elements( mesh ),
                    _expr = curlxv( runtimeField ) - curlxv( staticField ),
                    _quad = _Q<18>(),
                    _geomap = GeomapStrategyType::GEOMAP_HO ),
            2e-8 );
        BOOST_CHECK_SMALL(
            normL2( _range = elements( mesh ),
                    _expr = gradv( runtimeField ) - gradv( staticField ),
                    _quad = _Q<18>(),
                    _geomap = GeomapStrategyType::GEOMAP_HO ),
            5e-8 );
        BOOST_CHECK_SMALL(
            normL2( _range = elements( mesh ),
                    _expr = curlxv( discreteGradient ),
                    _quad = _Q<18>(),
                    _geomap = GeomapStrategyType::GEOMAP_HO ),
            5e-8 );
        checkTangentialTrace( runtimeSpace, runtimeField, mesh );
        checkHcurlAssembly( runtimeSpace, runtimeField, mesh );
    }
}

/** Exercise all compatible families for one FE order on geometry order G. */
template <int G, int Order, typename MeshPtrType>
void checkOrder( MeshPtrType const& mesh )
{
    checkHdivFamily<Order>(
        [&]
        { return RTh<Order>( mesh ); },
        [&]
        { return RTh<Dynamic>( mesh, RuntimeOrder{ Order } ); },
        mesh, "RT" );
    checkHdivFamily<Order>(
        [&]
        { return BDMh<Order>( mesh ); },
        [&]
        { return BDMh<Dynamic>( mesh, RuntimeOrder{ Order } ); },
        mesh, "BDM" );
    checkHcurlFamily<Order>(
        [&]
        { return Ned1h<Order>( mesh ); },
        [&]
        { return Ned1h<Dynamic>( mesh, RuntimeOrder{ Order } ); },
        mesh, "NED1" );
    checkHcurlFamily<Order>(
        [&]
        { return Ned2h<Order>( mesh ); },
        [&]
        { return Ned2h<Dynamic>( mesh, RuntimeOrder{ Order } ); },
        mesh, "NED2" );
}

/** Exercise FE orders zero through three on one curved geometry order. */
template <int G>
void checkGeometryOrder()
{
    auto mesh = unitCircle<G>( 0.55 );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_GT( mesh->numGlobalElements(), 0 );
    checkOrder<G, 0>( mesh );
    checkOrder<G, 1>( mesh );
    checkOrder<G, 2>( mesh );
    checkOrder<G, 3>( mesh );
}
} // namespace

BOOST_AUTO_TEST_SUITE( curved_piola_2d )

/** Validate the complete compatible-space matrix on quadratic triangles. */
BOOST_AUTO_TEST_CASE( quadratic_geometry )
{
    checkGeometryOrder<2>();
}

/** Validate the complete compatible-space matrix on cubic triangles. */
BOOST_AUTO_TEST_CASE( cubic_geometry )
{
#if FEELPP_MESH_MAX_ORDER >= 3
    checkGeometryOrder<3>();
#else
    BOOST_TEST_MESSAGE(
        "cubic mesh validation requires FEELPP_MESH_MAX_ORDER >= 3" );
#endif
}

BOOST_AUTO_TEST_SUITE_END()
