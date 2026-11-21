/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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

#include "fe_parity_checks.hpp"

#include <string_view>
#include <type_traits>

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
    field.on( _range = elements( mesh ), _expr = expression,
              _geomap = GeomapStrategyType::GEOMAP_HO );
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
        BOOST_CHECK_EQUAL( runtimeSpace->basis()->order(), Order );
        BOOST_CHECK_EQUAL( runtimeSpace->basis()->localDofCount(),
                           staticSpace->basis()->localDofCount() );
        BOOST_CHECK_EQUAL( runtimeSpace->nDof(), staticSpace->nDof() );
        Feel::test_detail::checkFiniteElementTransformParity(
            staticSpace, runtimeSpace, mesh );
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
    auto potential = [&]
    {
        if constexpr ( Order == 0 )
            return interpolateField(
                h1, mesh, cst( 0.5 ) + Px() - 0.75 * Py() );
        else
            return interpolateField(
                h1, mesh,
                cst( 0.5 ) + Px() - 0.75 * Py() +
                    0.2 * Px() * Py() );
    }();
    auto discreteGradient = interpolateField(
        runtimeSpace, mesh, trans( gradv( potential ) ) );

    BOOST_TEST_CONTEXT( family << " order " << Order )
    {
        BOOST_CHECK_EQUAL( runtimeSpace->basis()->order(), Order );
        BOOST_CHECK_EQUAL( runtimeSpace->basis()->localDofCount(),
                           staticSpace->basis()->localDofCount() );
        BOOST_CHECK_EQUAL( runtimeSpace->nDof(), staticSpace->nDof() );
        Feel::test_detail::checkFiniteElementTransformParity(
            staticSpace, runtimeSpace, mesh );
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
template <int Order, typename MeshPtrType>
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
    checkOrder<0>( mesh );
    checkOrder<1>( mesh );
    checkOrder<2>( mesh );
    checkOrder<3>( mesh );
}

/**
 * @brief Exercise compatible spaces on a truly runtime-order curved mesh.
 * @param order positive geometric polynomial order
 */
void checkRuntimeGeometryOrder( RuntimeOrder order )
{
    auto mesh = unitCircle( order, 0.55 );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_EQUAL(
        mesh->order(), order.value );
    BOOST_REQUIRE_GT(
        mesh->numGlobalElements(), 0 );
    checkOrder<0>( mesh );
    checkOrder<1>( mesh );
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

/**
 * @test Validate FunctionSpace, DoFTable, interpolation, and Piola operations
 *       on curved geometry whose polynomial order is selected at runtime.
 */
BOOST_AUTO_TEST_CASE( runtime_quadratic_and_cubic_geometry )
{
    checkRuntimeGeometryOrder( RuntimeOrder{ 2 } );
    checkRuntimeGeometryOrder( RuntimeOrder{ 3 } );
}

/**
 * @test Require equivalent nodal/Bezier storage in scalar PDE assembly and
 *       compatible Piola transforms on one unchanged curved mesh.
 */
BOOST_AUTO_TEST_CASE( nodal_bezier_pde_and_piola_parity )
{
    auto mesh = unitCircle<2>( 0.55 );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_GT( mesh->numGlobalElements(), 0 );
    auto const nodalCoordinates = mesh->nodalGeometryCoordinates();

    auto space = Pch<2>( mesh );
    auto trial = space->element();
    auto test = space->element();
    auto nodalMatrix = form2( _trial = space, _test = space );
    auto nodalVector = form1( _test = space );
    nodalMatrix = integrate(
        _range = elements( mesh ),
        _expr = inner( gradt( trial ), grad( test ) ) +
                idt( trial ) * id( test ),
        _quad = _Q<18>(),
        _geomap = GeomapStrategyType::GEOMAP_HO );
    nodalVector = integrate(
        _range = elements( mesh ),
        _expr = ( cst( 1.0 ) + Px() - 0.5 * Py() ) * id( test ),
        _quad = _Q<18>(),
        _geomap = GeomapStrategyType::GEOMAP_HO );
    nodalMatrix.close();
    nodalVector.close();

    mesh->setGeometryRepresentation(
        geometry::GeometryRepresentationKind::BernsteinBezier );
    BOOST_CHECK(
        mesh->geometryRepresentation() ==
        geometry::GeometryRepresentationKind::BernsteinBezier );
    BOOST_CHECK_SMALL(
        ( mesh->nodalGeometryCoordinates() - nodalCoordinates ).norm(),
        2e-11 );

    auto const elementRange = elements( mesh );
    auto const elementIterator = elementRange.begin();
    if ( elementIterator != elementRange.end() )
    {
        auto const& element = boost::unwrap_ref( *elementIterator );
        using gm_type =
            std::remove_reference_t<decltype( *mesh->gm() )>;
        typename gm_type::Inverse inverse( mesh->gm(), element );
        typename gm_type::node_t_type referencePoint( 2 );
        referencePoint( 0 ) = -0.35;
        referencePoint( 1 ) = -0.25;
        auto const physicalPoint = mesh->gm()->transform(
            referencePoint, inverse.G() );
        auto const inverseResult = inverse.inverse( physicalPoint );
        BOOST_CHECK( inverseResult.converged() );
        BOOST_CHECK_SMALL(
            inverseResult.reference_point( 0 ) - referencePoint( 0 ),
            2e-10 );
        BOOST_CHECK_SMALL(
            inverseResult.reference_point( 1 ) - referencePoint( 1 ),
            2e-10 );
    }

    auto bezierMatrix = form2( _trial = space, _test = space );
    auto bezierVector = form1( _test = space );
    bezierMatrix = integrate(
        _range = elements( mesh ),
        _expr = inner( gradt( trial ), grad( test ) ) +
                idt( trial ) * id( test ),
        _quad = _Q<18>(),
        _geomap = GeomapStrategyType::GEOMAP_HO );
    bezierVector = integrate(
        _range = elements( mesh ),
        _expr = ( cst( 1.0 ) + Px() - 0.5 * Py() ) * id( test ),
        _quad = _Q<18>(),
        _geomap = GeomapStrategyType::GEOMAP_HO );
    bezierMatrix.close();
    bezierVector.close();

    auto const matrixScale = std::max(
        { 1., std::abs( nodalMatrix.matrixPtr()->l1Norm() ),
          std::abs( bezierMatrix.matrixPtr()->l1Norm() ) } );
    bezierMatrix.matrixPtr()->addMatrix(
        -1., nodalMatrix.matrixPtr() );
    bezierMatrix.matrixPtr()->close();
    BOOST_CHECK_SMALL(
        bezierMatrix.matrixPtr()->l1Norm() / matrixScale, 2e-11 );
    BOOST_CHECK_SMALL(
        bezierMatrix.matrixPtr()->linftyNorm() / matrixScale, 2e-11 );

    auto vectorDifference = bezierVector.vectorPtr()->clone();
    vectorDifference->add( 1., bezierVector.vectorPtr() );
    vectorDifference->add( -1., nodalVector.vectorPtr() );
    vectorDifference->close();
    auto const vectorScale = std::max(
        1., nodalVector.vectorPtr()->l2Norm() );
    BOOST_CHECK_SMALL(
        vectorDifference->l2Norm() / vectorScale, 2e-11 );

    checkOrder<0>( mesh );

    mesh->setGeometryRepresentation(
        geometry::GeometryRepresentationKind::NodalLagrange );
    BOOST_CHECK_SMALL(
        ( mesh->pointCoordinates() - nodalCoordinates ).norm(),
        2e-11 );
}

BOOST_AUTO_TEST_SUITE_END()
