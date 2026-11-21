/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file test_interpolation_de_rham.cpp
 * @brief Discrete de Rham contracts assembled with OperatorInterpolation.
 */
#define BOOST_TEST_MODULE test_interpolation_de_rham
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelpoly/fecapabilities.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template<typename MeshPtrType, typename SpacePtrType, typename ExprType>
auto
interpolate( MeshPtrType const& mesh, SpacePtrType const& space, ExprType const& expr )
{
    auto result = space->element();
    result.on( _range = elements( mesh ), _expr = expr );
    return result;
}

template<typename MeshPtrType, typename LeftElementType, typename RightElementType>
double
differenceNorm( MeshPtrType const& mesh,
                LeftElementType const& left,
                RightElementType const& right )
{
    return normL2( _range = elements( mesh ),
                   _expr = idv( left ) - idv( right ),
                   _quad = _Q<16>(),
                   _geomap = GeomapStrategyType::GEOMAP_HO );
}

/**
 * Check the two-dimensional discrete de Rham complex assembled through the
 * public interpolation-operator API.
 */
template<int Order, typename MeshPtrType, typename HCurlSpacePtrType>
void
checkOperatorDeRhamComplex( MeshPtrType const& mesh,
                            HCurlSpacePtrType const& hcurlSpace )
{
    auto h1Space = Pch<Dynamic>(
        mesh, RuntimeOrder{ static_cast<uint16_type>( Order+1 ) } );
    auto l2Space = Pdh<Dynamic>(
        mesh, RuntimeOrder{ static_cast<uint16_type>( Order ) } );
    auto potentialExpression = [&]()
    {
        if constexpr ( Order == 0 )
            return cst( 0.5 ) + Px() - 0.75 * Py();
        else
            return cst( 0.5 ) + Px() - 0.75 * Py() +
                   Px()*Px() + 0.5*Px()*Py() - Py()*Py();
    }();
    auto potential = interpolate( mesh, h1Space, potentialExpression );

    auto discreteGradient = hcurlSpace->element();
    auto discreteCurl = l2Space->element();
    auto gradientOperator = Grad( _domainSpace = h1Space,
                                  _imageSpace = hcurlSpace );
    auto curlOperator = Curl( _domainSpace = hcurlSpace,
                              _imageSpace = l2Space );
    gradientOperator.apply( potential, discreteGradient );
    curlOperator.apply( discreteGradient, discreteCurl );

    BOOST_CHECK_SMALL(
        normL2( _range = elements( mesh ),
                _expr = idv( discreteGradient ) - trans( gradv( potential ) ),
                _quad = _Q<16>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-9 );
    BOOST_CHECK_SMALL(
        normL2( _range = elements( mesh ),
                _expr = idv( discreteCurl ),
                _quad = _Q<16>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-9 );

    auto hcurlField = interpolate(
        mesh, hcurlSpace, vec( -Py(), Px() ) );
    auto operatorCurl = l2Space->element();
    auto directCurl = interpolate(
        mesh, l2Space, curlxv( hcurlField ) );
    curlOperator.apply( hcurlField, operatorCurl );
    BOOST_CHECK_SMALL(
        normL2( _range = elements( mesh ),
                _expr = idv( operatorCurl ) - idv( directCurl ),
                _quad = _Q<16>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
        2e-9 );
}

/** Check ID, curl, and div assembly for static and true-runtime Piola spaces. */
template<typename MeshPtrType>
void
checkPiolaInterpolationOperators( MeshPtrType const& mesh )
{
    auto scalarImage = Pdh<1>( mesh );
    auto vectorImage = Pdhv<1>( mesh );

    auto nedStaticSpace = Ned1h<0>( mesh );
    auto nedRuntimeSpace = Ned1h<Dynamic>( mesh, RuntimeOrder{ 0 } );
    auto rtStaticSpace = RTh<0>( mesh );
    auto rtRuntimeSpace = RTh<Dynamic>( mesh, RuntimeOrder{ 0 } );

    auto tangentExpression = vec( cst( 0.25 ) - Py() + Px() * Py(),
                                  cst( -0.75 ) + Px() + Px() * Px() );
    auto fluxExpression = vec( cst( 1.0 ) + Px() * Px() - Py(),
                               cst( -0.5 ) + Px() + Py() * Py() );

    auto nedStatic = interpolate( mesh, nedStaticSpace, tangentExpression );
    auto nedRuntime = interpolate( mesh, nedRuntimeSpace, tangentExpression );
    auto rtStatic = interpolate( mesh, rtStaticSpace, fluxExpression );
    auto rtRuntime = interpolate( mesh, rtRuntimeSpace, fluxExpression );

    auto checkIdentity = [&]( auto const& domainSpace, auto const& field )
    {
        auto assembled = vectorImage->element();
        auto expected = interpolate( mesh, vectorImage, idv( field ) );
        auto identity = opInterpolation( _domainSpace = domainSpace,
                                         _imageSpace = vectorImage );
        identity->apply( field, assembled );
        BOOST_CHECK_SMALL( differenceNorm( mesh, assembled, expected ), 2e-10 );
    };

    auto checkCurl = [&]( auto const& domainSpace, auto const& field )
    {
        auto assembled = scalarImage->element();
        auto expected = interpolate( mesh, scalarImage, curlxv( field ) );
        auto curlOperator = Curl( _domainSpace = domainSpace,
                                  _imageSpace = scalarImage );
        curlOperator.apply( field, assembled );
        BOOST_CHECK_SMALL( differenceNorm( mesh, assembled, expected ), 2e-10 );
    };

    auto checkDiv = [&]( auto const& domainSpace, auto const& field )
    {
        auto assembled = scalarImage->element();
        auto expected = interpolate( mesh, scalarImage, divv( field ) );
        auto divOperator = Div( _domainSpace = domainSpace,
                                _imageSpace = scalarImage );
        divOperator.apply( field, assembled );
        BOOST_CHECK_SMALL( differenceNorm( mesh, assembled, expected ), 2e-10 );
    };

    checkIdentity( nedStaticSpace, nedStatic );
    checkIdentity( nedRuntimeSpace, nedRuntime );
    checkIdentity( rtStaticSpace, rtStatic );
    checkIdentity( rtRuntimeSpace, rtRuntime );
    checkCurl( nedStaticSpace, nedStatic );
    checkCurl( nedRuntimeSpace, nedRuntime );
    checkDiv( rtStaticSpace, rtStatic );
    checkDiv( rtRuntimeSpace, rtRuntime );

    auto checkHCurlGradient = [&]( auto const& field )
    {
        BOOST_CHECK_SMALL(
            normL2(
                _range = elements( mesh ),
                _expr = ( gradv( field )( 1, 0 ) -
                          gradv( field )( 0, 1 ) ) - curlxv( field ),
                _quad = _Q<16>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
            3e-9 );
    };
    auto checkHDivGradient = [&]( auto const& field )
    {
        BOOST_CHECK_SMALL(
            normL2(
                _range = elements( mesh ),
                _expr = gradv( field )( 0, 0 ) +
                            gradv( field )( 1, 1 ) - divv( field ),
                _quad = _Q<16>(),
                _geomap = GeomapStrategyType::GEOMAP_HO ),
            3e-9 );
    };
    checkHCurlGradient( nedStatic );
    checkHCurlGradient( nedRuntime );
    checkHDivGradient( rtStatic );
    checkHDivGradient( rtRuntime );

    // Image-functional path: a runtime Lagrange domain is evaluated by the
    // static/runtime RT dual through the same H(div) interpolation kernels.
    auto vectorDomain = Pchv<Dynamic>( mesh, RuntimeOrder{ 2 } );
    auto domainField = interpolate(
        mesh, vectorDomain, vec( cst( 1.0 ), cst( -0.5 ) ) );
    auto expectedRtStatic = interpolate(
        mesh, rtStaticSpace, idv( domainField ) );
    auto expectedRtRuntime = interpolate(
        mesh, rtRuntimeSpace, idv( domainField ) );
    auto operatorRtStatic = rtStaticSpace->element();
    auto operatorRtRuntime = rtRuntimeSpace->element();
    auto toRtStatic = opInterpolation( _domainSpace = vectorDomain,
                                       _imageSpace = rtStaticSpace );
    auto toRtRuntime = opInterpolation( _domainSpace = vectorDomain,
                                        _imageSpace = rtRuntimeSpace );
    toRtStatic->apply( domainField, operatorRtStatic );
    toRtRuntime->apply( domainField, operatorRtRuntime );
    BOOST_CHECK_SMALL(
        differenceNorm( mesh, operatorRtStatic, expectedRtStatic ), 2e-10 );
    BOOST_CHECK_SMALL(
        differenceNorm( mesh, operatorRtRuntime, expectedRtRuntime ), 2e-10 );
}

/** Check localized FE-context evaluation between independent meshes. */
template<typename MeshPtrType>
void
checkUnrelatedPiolaInterpolationOperators( MeshPtrType const& domainMesh,
                                           MeshPtrType const& imageMesh,
                                           bool checkAnalyticFields )
{
    BOOST_REQUIRE( !imageMesh->isRelatedTo( domainMesh ) );
    auto const execution = Environment::worldComm().globalSize() > 1
                               ? fe::Execution::Distributed
                               : fe::Execution::Serial;
    auto const rtCapability = fe::capability( {
        fe::Family::RaviartThomas,
        geometry::GeometryCellKind::Simplex,
        2, 1, imageMesh->order(),
        fe::Operation::UnrelatedMeshTransfer,
        execution } );
    auto const nedCapability = fe::capability( {
        fe::Family::NedelecFirstKind,
        geometry::GeometryCellKind::Simplex,
        2, 1, imageMesh->order(),
        fe::Operation::UnrelatedMeshTransfer,
        execution } );
    if ( execution == fe::Execution::Serial )
    {
        BOOST_REQUIRE( rtCapability.validated() );
        BOOST_REQUIRE( nedCapability.validated() );
    }
    else
    {
        BOOST_REQUIRE( !rtCapability.available() );
        BOOST_REQUIRE( !nedCapability.available() );
    }

    auto scalarImage = Pdh<1>( imageMesh );
    auto vectorImage = Pdhv<1>( imageMesh );
    auto nedStaticSpace = Ned1h<0>( domainMesh );
    auto nedRuntimeSpace = Ned1h<Dynamic>( domainMesh, RuntimeOrder{ 0 } );
    auto rtStaticSpace = RTh<0>( domainMesh );
    auto rtRuntimeSpace = RTh<Dynamic>( domainMesh, RuntimeOrder{ 0 } );

    auto nedExpression = vec( -Py(), Px() );
    auto rtExpression = vec( Px(), Py() );
    auto nedStatic = interpolate( domainMesh, nedStaticSpace, nedExpression );
    auto nedRuntime = interpolate( domainMesh, nedRuntimeSpace, nedExpression );
    auto rtStatic = interpolate( domainMesh, rtStaticSpace, rtExpression );
    auto rtRuntime = interpolate( domainMesh, rtRuntimeSpace, rtExpression );

    auto applyIdentity = [&]( auto const& domainSpace, auto const& field )
    {
        auto assembled = vectorImage->element();
        auto identity = opInterpolation( _domainSpace = domainSpace,
                                         _imageSpace = vectorImage );
        identity->apply( field, assembled );
        return assembled;
    };
    auto applyCurl = [&]( auto const& domainSpace, auto const& field )
    {
        auto assembled = scalarImage->element();
        auto curlOperator = Curl( _domainSpace = domainSpace,
                                  _imageSpace = scalarImage );
        curlOperator.apply( field, assembled );
        return assembled;
    };
    auto applyDiv = [&]( auto const& domainSpace, auto const& field )
    {
        auto assembled = scalarImage->element();
        auto divOperator = Div( _domainSpace = domainSpace,
                                _imageSpace = scalarImage );
        divOperator.apply( field, assembled );
        return assembled;
    };

    if ( Environment::worldComm().localSize() > 1 )
    {
        // The legacy unrelated-mesh MPI matrix protocol cannot yet transport
        // mapped Piola values or differential operands correctly. Keep this
        // contract executable in MPI and require a safe, explicit failure.
        BOOST_CHECK_THROW( applyIdentity( nedStaticSpace, nedStatic ),
                           std::logic_error );
        BOOST_CHECK_THROW( applyIdentity( rtStaticSpace, rtStatic ),
                           std::logic_error );
        BOOST_CHECK_THROW( applyCurl( nedStaticSpace, nedStatic ),
                           std::logic_error );
        BOOST_CHECK_THROW( applyDiv( rtStaticSpace, rtStatic ),
                           std::logic_error );
        return;
    }

    auto nedIdentityStatic = applyIdentity( nedStaticSpace, nedStatic );
    auto nedIdentityRuntime = applyIdentity( nedRuntimeSpace, nedRuntime );
    auto rtIdentityStatic = applyIdentity( rtStaticSpace, rtStatic );
    auto rtIdentityRuntime = applyIdentity( rtRuntimeSpace, rtRuntime );
    auto nedCurlStatic = applyCurl( nedStaticSpace, nedStatic );
    auto nedCurlRuntime = applyCurl( nedRuntimeSpace, nedRuntime );
    auto rtDivStatic = applyDiv( rtStaticSpace, rtStatic );
    auto rtDivRuntime = applyDiv( rtRuntimeSpace, rtRuntime );

    BOOST_CHECK_SMALL( differenceNorm( imageMesh, nedIdentityStatic,
                                       nedIdentityRuntime ), 5e-9 );
    BOOST_CHECK_SMALL( differenceNorm( imageMesh, rtIdentityStatic,
                                       rtIdentityRuntime ), 5e-9 );
    BOOST_CHECK_SMALL( differenceNorm( imageMesh, nedCurlStatic,
                                       nedCurlRuntime ), 5e-9 );
    BOOST_CHECK_SMALL( differenceNorm( imageMesh, rtDivStatic,
                                       rtDivRuntime ), 5e-9 );

    if ( checkAnalyticFields )
    {
        auto expectedNed = interpolate( imageMesh, vectorImage, nedExpression );
        auto expectedRt = interpolate( imageMesh, vectorImage, rtExpression );
        auto expectedDerivative = interpolate( imageMesh, scalarImage, cst( 2.0 ) );
        BOOST_CHECK_SMALL( differenceNorm( imageMesh, nedIdentityStatic,
                                           expectedNed ), 5e-9 );
        BOOST_CHECK_SMALL( differenceNorm( imageMesh, rtIdentityStatic,
                                           expectedRt ), 5e-9 );
        BOOST_CHECK_SMALL( differenceNorm( imageMesh, nedCurlStatic,
                                           expectedDerivative ), 5e-9 );
        BOOST_CHECK_SMALL( differenceNorm( imageMesh, rtDivStatic,
                                           expectedDerivative ), 5e-9 );
    }

    if ( Environment::worldComm().localSize() == 1 )
    {
        // Apply complete unrelated-mesh RT/Nedelec duals, including their
        // interior moments, to one runtime-order domain field. Distributed
        // tests exercise the independently transported point-functional path.
        auto vectorDomain = Pchv<Dynamic>( domainMesh, RuntimeOrder{ 2 } );
        auto vectorExpression = vec( cst( 0.6 ) + Px() - 0.2*Py(),
                                     cst( -0.4 ) + 0.3*Px() + Py() );
        auto vectorField = interpolate(
            domainMesh, vectorDomain, vectorExpression );
        auto rtMomentStaticSpace = RTh<1>( imageMesh );
        auto rtMomentRuntimeSpace = RTh<Dynamic>(
            imageMesh, RuntimeOrder{ 1 } );
        auto nedMomentStaticSpace = Ned1h<1>( imageMesh );
        auto nedMomentRuntimeSpace = Ned1h<Dynamic>(
            imageMesh, RuntimeOrder{ 1 } );

        auto checkMomentImage = [&]( auto const& imageSpace )
        {
            auto assembled = imageSpace->element();
            auto expected = interpolate(
                imageMesh, imageSpace, vectorExpression );
            auto interpolation = opInterpolation(
                _domainSpace = vectorDomain, _imageSpace = imageSpace );
            interpolation->apply( vectorField, assembled );
            BOOST_CHECK_SMALL(
                differenceNorm( imageMesh, assembled, expected ), 2e-8 );
        };
        checkMomentImage( rtMomentStaticSpace );
        checkMomentImage( rtMomentRuntimeSpace );
        checkMomentImage( nedMomentStaticSpace );
        checkMomentImage( nedMomentRuntimeSpace );
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE( interpolation_de_rham_suite )

/** Verify affine reduction of Piola-aware interpolation operators. */
BOOST_AUTO_TEST_CASE( affine_static_runtime_operator_contract )
{
    auto mesh = unitCircle<1>( 0.45 );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_GT( mesh->numGlobalElements(), 0 );
    checkPiolaInterpolationOperators( mesh );
    checkOperatorDeRhamComplex<0>( mesh, Ned1h<0>( mesh ) );
    checkOperatorDeRhamComplex<0>(
        mesh, Ned1h<Dynamic>( mesh, RuntimeOrder{ 0 } ) );
    checkOperatorDeRhamComplex<1>( mesh, Ned1h<1>( mesh ) );
    checkOperatorDeRhamComplex<1>(
        mesh, Ned1h<Dynamic>( mesh, RuntimeOrder{ 1 } ) );
}

/** Verify pointwise Piola operators on a genuinely curved quadratic mesh. */
BOOST_AUTO_TEST_CASE( curved_static_runtime_operator_contract )
{
    if ( Environment::worldComm().localSize() > 1 )
    {
        BOOST_TEST_MESSAGE( "curved OperatorInterpolation contract is serial-only" );
        return;
    }

    auto mesh = unitCircle<2>( 0.45 );
    BOOST_REQUIRE( mesh );
    BOOST_REQUIRE_GT( mesh->numGlobalElements(), 0 );
    checkPiolaInterpolationOperators( mesh );
    checkOperatorDeRhamComplex<0>( mesh, Ned1h<0>( mesh ) );
    checkOperatorDeRhamComplex<0>(
        mesh, Ned1h<Dynamic>( mesh, RuntimeOrder{ 0 } ) );
}

/** Verify localized mapped operands between independent affine meshes. */
BOOST_AUTO_TEST_CASE( unrelated_affine_static_runtime_operator_contract )
{
    auto domainMesh = unitCircle<1>( 0.42 );
    auto imageMesh = unitCircle<1>( 0.31 );
    BOOST_REQUIRE( domainMesh );
    BOOST_REQUIRE( imageMesh );
    checkUnrelatedPiolaInterpolationOperators( domainMesh, imageMesh, true );
}

/** Verify localized mapped operands between independent curved meshes. */
BOOST_AUTO_TEST_CASE( unrelated_curved_static_runtime_operator_contract )
{
    auto domainMesh = unitCircle<2>( 0.42 );
    auto imageMesh = unitCircle<2>( 0.31 );
    BOOST_REQUIRE( domainMesh );
    BOOST_REQUIRE( imageMesh );
    checkUnrelatedPiolaInterpolationOperators( domainMesh, imageMesh, false );
}

BOOST_AUTO_TEST_SUITE_END()
