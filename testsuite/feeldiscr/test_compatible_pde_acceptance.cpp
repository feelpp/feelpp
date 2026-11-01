/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file test_compatible_pde_acceptance.cpp
 * @brief Static/runtime Darcy, Maxwell, and Hodge solve acceptance checks.
 */
#define BOOST_TEST_MODULE test_compatible_pde_acceptance
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelpoly/fecapabilities.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>

#include "fe_parity_checks.hpp"

#include <boost/hana/integral_constant.hpp>

#include <utility>

using namespace Feel;
using namespace boost::hana::literals;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template<typename MeshPtrType>
void
checkDarcyParity( MeshPtrType const& mesh )
{
    fe::requireCapability<fe::Family::RaviartThomas>(
        *mesh, RuntimeOrder{ 1 }, fe::Operation::LinearSolve );

    auto staticFluxSpace = RTh<1>( mesh );
    auto runtimeFluxSpace = RTh<Dynamic>(
        mesh, RuntimeOrder{ 1 } );
    auto staticPressureSpace = Pdh<1>( mesh );
    auto runtimePressureSpace = Pdh<Dynamic>(
        mesh, RuntimeOrder{ 1 } );
    auto staticProduct = product(
        staticFluxSpace, staticPressureSpace );
    auto runtimeProduct = product(
        runtimeFluxSpace, runtimePressureSpace );

    auto fluxExpression = vec(
        cst( 0.5 ) + Px() - 0.2*Py(),
        cst( -0.25 ) + 0.3*Px() + Py() );
    auto pressureExpression =
        cst( 0.75 ) + 0.2*Px() - 0.4*Py();
    auto staticFluxTarget = staticFluxSpace->element();
    auto runtimeFluxTarget = runtimeFluxSpace->element();
    auto staticPressureTarget = staticPressureSpace->element();
    auto runtimePressureTarget = runtimePressureSpace->element();
    staticFluxTarget.on(
        _range = elements( mesh ), _expr = fluxExpression );
    runtimeFluxTarget.on(
        _range = elements( mesh ), _expr = fluxExpression );
    staticPressureTarget.on(
        _range = elements( mesh ), _expr = pressureExpression );
    runtimePressureTarget.on(
        _range = elements( mesh ), _expr = pressureExpression );

    auto assembleSystem = [&]( auto const& productSpace,
                               auto const& fluxSpace,
                               auto const& pressureSpace,
                               auto const& fluxTarget,
                               auto const& pressureTarget )
    {
        auto flux = fluxSpace->element();
        auto pressure = pressureSpace->element();
        auto fluxTest = fluxSpace->element();
        auto pressureTest = pressureSpace->element();
        auto op = blockform2(
            productSpace, solve::strategy::monolithic,
            backend(), Pattern::COUPLED );
        auto rhs = blockform1(
            productSpace, solve::strategy::monolithic, backend() );
        op( 0_c, 0_c ) += integrate(
            _range = elements( mesh ),
            _expr = inner( idt( flux ), id( fluxTest ) ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        op( 0_c, 1_c ) += integrate(
            _range = elements( mesh ),
            _expr = -idt( pressure )*div( fluxTest ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        op( 1_c, 0_c ) += integrate(
            _range = elements( mesh ),
            _expr = divt( flux )*id( pressureTest ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        op( 1_c, 1_c ) += integrate(
            _range = elements( mesh ),
            _expr = idt( pressure )*id( pressureTest ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        rhs( 0_c ) += integrate(
            _range = elements( mesh ),
            _expr = inner( idv( fluxTarget ), id( fluxTest ) )-
                    idv( pressureTarget )*div( fluxTest ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        rhs( 1_c ) += integrate(
            _range = elements( mesh ),
            _expr = divv( fluxTarget )*id( pressureTest )+
                    idv( pressureTarget )*id( pressureTest ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        return std::make_pair( std::move( op ), std::move( rhs ) );
    };
    auto [staticOperator, staticRhs] = assembleSystem(
        staticProduct, staticFluxSpace, staticPressureSpace,
        staticFluxTarget, staticPressureTarget );
    auto [runtimeOperator, runtimeRhs] = assembleSystem(
        runtimeProduct, runtimeFluxSpace, runtimePressureSpace,
        runtimeFluxTarget, runtimePressureTarget );

    staticOperator.close();
    runtimeOperator.close();
    staticRhs.close();
    runtimeRhs.close();
    Feel::test_detail::checkAssembledMatrixNormParity(
        staticOperator, runtimeOperator, 5e-10 );
    Feel::test_detail::checkAssembledVectorParity(
        staticRhs, runtimeRhs, 5e-10 );

    auto staticSolution = staticProduct.element();
    auto runtimeSolution = runtimeProduct.element();
    auto const staticSolve = staticOperator.solve(
        _rhs = staticRhs, _solution = staticSolution,
        _rebuild = true );
    auto const runtimeSolve = runtimeOperator.solve(
        _rhs = runtimeRhs, _solution = runtimeSolution,
        _rebuild = true );
    Feel::test_detail::checkSolveParity( staticSolve, runtimeSolve );

    BOOST_CHECK_SMALL(
        normL2( _range = elements( mesh ),
                _expr = idv( staticSolution( 0_c ) )-
                        idv( runtimeSolution( 0_c ) ),
                _quad = _Q<14>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-8 );
    BOOST_CHECK_SMALL(
        normL2( _range = elements( mesh ),
                _expr = divv( staticSolution( 0_c ) )-
                        divv( runtimeSolution( 0_c ) ),
                _quad = _Q<14>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-8 );
    BOOST_CHECK_SMALL(
        normL2( _range = boundaryfaces( mesh ),
                _expr = trans( idv( staticSolution( 0_c ) )-
                               idv( runtimeSolution( 0_c ) ) )*N(),
                _quad = _Q<14>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-8 );
    BOOST_CHECK_SMALL(
        normL2( _range = elements( mesh ),
                _expr = idv( staticSolution( 1_c ) )-
                        idv( runtimeSolution( 1_c ) ),
                _quad = _Q<14>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-8 );
}

template<typename MeshPtrType, typename ExpressionType>
void
checkHcurlParity( MeshPtrType const& mesh,
                  ExpressionType const& targetExpression,
                  bool expectClosedForm )
{
    fe::requireCapability<fe::Family::NedelecFirstKind>(
        *mesh, RuntimeOrder{ 1 }, fe::Operation::LinearSolve );
    auto staticSpace = Ned1h<1>( mesh );
    auto runtimeSpace = Ned1h<Dynamic>(
        mesh, RuntimeOrder{ 1 } );
    auto staticTarget = staticSpace->element();
    auto runtimeTarget = runtimeSpace->element();
    staticTarget.on(
        _range = elements( mesh ), _expr = targetExpression );
    runtimeTarget.on(
        _range = elements( mesh ), _expr = targetExpression );

    auto assembleSystem = [&]( auto const& space, auto const& target )
    {
        auto trial = space->element();
        auto test = space->element();
        auto op = form2( _trial = space, _test = space );
        auto rhs = form1( _test = space );
        op = integrate(
            _range = elements( mesh ),
            _expr = inner( idt( trial ), id( test ) )+
                    curlt( trial )*curl( test ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        rhs = integrate(
            _range = elements( mesh ),
            _expr = inner( idv( target ), id( test ) )+
                    curlxv( target )*curl( test ),
            _quad = _Q<14>(),
            _geomap = GeomapStrategyType::GEOMAP_HO );
        return std::make_pair( std::move( op ), std::move( rhs ) );
    };
    auto [staticOperator, staticRhs] =
        assembleSystem( staticSpace, staticTarget );
    auto [runtimeOperator, runtimeRhs] =
        assembleSystem( runtimeSpace, runtimeTarget );

    Feel::test_detail::checkAssembledMatrixParity(
        staticOperator, runtimeOperator, 5e-10 );
    Feel::test_detail::checkAssembledVectorParity(
        staticRhs, runtimeRhs, 5e-10 );
    auto staticSolution = staticSpace->element();
    auto runtimeSolution = runtimeSpace->element();
    auto const staticSolve = staticOperator.solve(
        _rhs = staticRhs, _solution = staticSolution,
        _rebuild = true );
    auto const runtimeSolve = runtimeOperator.solve(
        _rhs = runtimeRhs, _solution = runtimeSolution,
        _rebuild = true );
    Feel::test_detail::checkSolveParity( staticSolve, runtimeSolve );

    BOOST_CHECK_SMALL(
        normL2( _range = elements( mesh ),
                _expr = idv( staticSolution )-idv( runtimeSolution ),
                _quad = _Q<14>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-8 );
    BOOST_CHECK_SMALL(
        normL2( _range = elements( mesh ),
                _expr = curlxv( staticSolution )-
                        curlxv( runtimeSolution ),
                _quad = _Q<14>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-8 );
    BOOST_CHECK_SMALL(
        normL2( _range = boundaryfaces( mesh ),
                _expr = trans( idv( staticSolution )-
                               idv( runtimeSolution ) )*
                        vec( Tx(), Ty() ),
                _quad = _Q<14>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-8 );
    if ( expectClosedForm )
        BOOST_CHECK_SMALL(
            normL2( _range = elements( mesh ),
                    _expr = curlxv( runtimeSolution ),
                    _quad = _Q<14>(),
                    _geomap = GeomapStrategyType::GEOMAP_HO ),
            2e-8 );
}
} // namespace

BOOST_AUTO_TEST_SUITE( compatible_pde_acceptance )

BOOST_AUTO_TEST_CASE( published_capability_matrix )
{
    constexpr auto tensorRt = fe::capability( {
        fe::Family::RaviartThomas,
        geometry::GeometryCellKind::Hypercube,
        2, 1, 1, fe::Operation::LinearSolve } );
    static_assert( !tensorRt.available() );
    constexpr auto curvedNed = fe::capability( {
        fe::Family::NedelecFirstKind,
        geometry::GeometryCellKind::Simplex,
        2, 1, 3, fe::Operation::LinearSolve } );
    static_assert( curvedNed.validated() );
    constexpr auto curvedExport = fe::capability( {
        fe::Family::Lagrange,
        geometry::GeometryCellKind::Simplex,
        2, 2, 3, fe::Operation::Export } );
    static_assert( curvedExport.available() && !curvedExport.validated() );
    constexpr auto highOrderNed = fe::capability( {
        fe::Family::NedelecFirstKind,
        geometry::GeometryCellKind::Simplex,
        3, 4, 2, fe::Operation::LinearSolve } );
    static_assert( highOrderNed.available() && !highOrderNed.validated() );
    constexpr auto bdmSolve = fe::capability( {
        fe::Family::BrezziDouglasMarini,
        geometry::GeometryCellKind::Simplex,
        2, 1, 2, fe::Operation::LinearSolve } );
    static_assert( bdmSolve.available() && !bdmSolve.validated() );

    BOOST_CHECK_NO_THROW(
        fe::requireCapability( {
            fe::Family::NedelecFirstKind,
            geometry::GeometryCellKind::Simplex,
            3, 4, 2, fe::Operation::LinearSolve } ) );
    BOOST_CHECK_THROW(
        fe::requireCapability( {
            fe::Family::BrezziDouglasMarini,
            geometry::GeometryCellKind::Hypercube,
            2, 1, 1, fe::Operation::SpaceConstruction } ),
        std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( unrelated_mesh_transfer_capabilities )
{
    // Existing aggregate queries continue to describe identity transfer
    // within one family. Distributed availability is not a validation claim.
    constexpr fe::CapabilityQuery nodalQuery{
        fe::Family::Lagrange, geometry::GeometryCellKind::Simplex,
        2, 1, 1, fe::Operation::UnrelatedMeshTransfer,
        fe::Execution::Distributed };
    static_assert( fe::capability( nodalQuery ).available() );
    static_assert( !fe::capability( nodalQuery ).validated() );
    BOOST_CHECK_NO_THROW( fe::requireCapability( nodalQuery ) );

    for ( auto domainFamily : { fe::Family::Lagrange,
                               fe::Family::RaviartThomas,
                               fe::Family::BrezziDouglasMarini,
                               fe::Family::NedelecFirstKind,
                               fe::Family::NedelecSecondKind } )
        for ( auto imageFamily : { fe::Family::Lagrange,
                                  fe::Family::RaviartThomas,
                                  fe::Family::BrezziDouglasMarini,
                                  fe::Family::NedelecFirstKind,
                                  fe::Family::NedelecSecondKind } )
        {
            auto query = nodalQuery;
            query.family = domainFamily;
            query.transfer_image_family = imageFamily;
            bool const supported = domainFamily == fe::Family::Lagrange &&
                                   imageFamily == fe::Family::Lagrange;
            BOOST_TEST_CONTEXT( fe::familyName( domainFamily ) << " -> " <<
                                fe::familyName( imageFamily ) )
            {
                BOOST_CHECK_EQUAL( fe::capability( query ).available(), supported );
                if ( supported )
                    BOOST_CHECK_NO_THROW( fe::requireCapability( query ) );
                else
                    BOOST_CHECK_THROW( fe::requireCapability( query ),
                                       std::invalid_argument );

                query.transfer_operand = fe::TransferOperand::Differential;
                BOOST_CHECK( !fe::capability( query ).available() );
                BOOST_CHECK_THROW( fe::requireCapability( query ),
                                   std::invalid_argument );
            }
        }

    auto serialQuery = nodalQuery;
    serialQuery.execution = fe::Execution::Serial;
    BOOST_CHECK( fe::capability( serialQuery ).validated() );
    serialQuery.family = fe::Family::RaviartThomas;
    BOOST_CHECK( fe::capability( serialQuery ).validated() );
    serialQuery.family = fe::Family::NedelecFirstKind;
    BOOST_CHECK( fe::capability( serialQuery ).validated() );
    serialQuery.transfer_operand = fe::TransferOperand::Differential;
    BOOST_CHECK( fe::capability( serialQuery ).available() );
    BOOST_CHECK( !fe::capability( serialQuery ).validated() );
    serialQuery.transfer_operand = fe::TransferOperand::Identity;
    serialQuery.transfer_image_family = fe::Family::Lagrange;
    BOOST_CHECK( fe::capability( serialQuery ).available() );
    BOOST_CHECK( !fe::capability( serialQuery ).validated() );
}

BOOST_AUTO_TEST_CASE( mixed_darcy_curved_geometry )
{
    auto mesh = unitCircle( RuntimeOrder{ 2 }, 0.8 );
    BOOST_REQUIRE( mesh );
    checkDarcyParity( mesh );
}

BOOST_AUTO_TEST_CASE( maxwell_and_hodge_geometry_envelope )
{
    auto affineMesh = unitCircle( RuntimeOrder{ 1 }, 0.8 );
    auto cubicMesh = unitCircle( RuntimeOrder{ 3 }, 0.8 );
    BOOST_REQUIRE( affineMesh );
    BOOST_REQUIRE( cubicMesh );
    checkHcurlParity(
        cubicMesh,
        vec( cst( 0.25 )-Py()+Px()*Py(),
             cst( -0.75 )+Px()+Px()*Px() ),
        false );
    checkHcurlParity(
        affineMesh,
        vec( cst( 1.0 )+0.4*Py(),
             cst( -0.75 )+0.4*Px() ),
        true );
}

BOOST_AUTO_TEST_SUITE_END()
