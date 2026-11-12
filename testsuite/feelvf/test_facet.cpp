/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE feelvf_facet
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{

template <typename MatrixPtr>
void
checkMatrixClose( MatrixPtr const& actual, MatrixPtr const& expected )
{
    auto delta = actual->clone();
    delta->addMatrix( -1.0, expected );
    delta->close();
    BOOST_CHECK_SMALL( delta->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( delta->linftyNorm(), 1e-9 );
}

} // namespace

BOOST_AUTO_TEST_SUITE( feelvf_facet_suite )

BOOST_AUTO_TEST_CASE( facet_range_helpers_match_existing_ranges )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();

    auto const interiorMeasure =
        integrate( _range = facet::interior( mesh ),
                   _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
    auto const explicitInteriorMeasure =
        integrate( _range = internalfaces( mesh ),
                   _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( interiorMeasure - explicitInteriorMeasure, 1e-12 );

    auto const boundaryMeasure =
        integrate( _range = facet::boundary( mesh ),
                   _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
    auto const explicitBoundaryMeasure =
        integrate( _range = boundaryfaces( mesh ),
                   _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( boundaryMeasure - explicitBoundaryMeasure, 1e-12 );

    auto const skeletonMeasure =
        integrate( _range = facet::skeleton( mesh ),
                   _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
    auto const explicitSkeletonMeasure =
        integrate( _range = faces( mesh ),
                   _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( skeletonMeasure - explicitSkeletonMeasure, 1e-12 );
}

BOOST_AUTO_TEST_CASE( facet_test_trial_pairing_matches_explicit_left_right_form )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pdh<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto helperForm = form2( _test = Xh, _trial = Xh );
    helperForm +=
        integrate( _range = facet::interior( mesh ),
                   _expr = facet::testTrialSum( v, u ) )
        + integrate( _range = facet::boundary( mesh ),
                     _expr = v * u );
    helperForm.close();

    auto explicitForm = form2( _test = Xh, _trial = Xh );
    explicitForm +=
        integrate( _range = internalfaces( mesh ),
                   _expr = leftface( v ) * leftfacet( u ) +
                           rightface( v ) * rightfacet( u ) )
        + integrate( _range = boundaryfaces( mesh ),
                     _expr = v * u );
    explicitForm.close();

    checkMatrixClose( helperForm.matrixPtr(), explicitForm.matrixPtr() );
}

BOOST_AUTO_TEST_CASE( facet_role_inferred_helpers_match_explicit_side_helpers )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pdh<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto testCompound = cst( 2.0 ) * v;
    auto trialCompound = cst( 3.0 ) * u;
    auto ambiguous = u * v;

    static_assert( requires( decltype( v ) const& expr ) { facet::left( expr ); } );
    static_assert( requires( decltype( u ) const& expr ) { facet::left( expr ); } );
    static_assert( requires( decltype( testCompound ) const& expr ) { facet::sum( expr ); } );
    static_assert( requires( decltype( trialCompound ) const& expr ) { facet::sum( expr ); } );
    static_assert( facet::detail::FacetExpressionWithAmbiguousRole<decltype( ambiguous )> );

    auto helperForm = form2( _test = Xh, _trial = Xh );
    helperForm +=
        integrate( _range = facet::interior( mesh ),
                   _expr = facet::left( v ) * facet::left( u ) +
                           facet::right( v ) * facet::right( u ) +
                           0.25 * facet::sum( v ) * facet::sum( u ) +
                           facet::mean( testCompound ) * facet::mean( trialCompound ) );
    helperForm.close();

    auto explicitForm = form2( _test = Xh, _trial = Xh );
    explicitForm +=
        integrate( _range = internalfaces( mesh ),
                   _expr = leftface( v ) * leftfacet( u ) +
                           rightface( v ) * rightfacet( u ) +
                           0.25 * ( leftface( v ) + rightface( v ) ) *
                               ( leftfacet( u ) + rightfacet( u ) ) +
                           ( 0.5 * ( leftface( testCompound ) + rightface( testCompound ) ) ) *
                               ( 0.5 * ( leftfacet( trialCompound ) + rightfacet( trialCompound ) ) ) );
    explicitForm.close();

    checkMatrixClose( helperForm.matrixPtr(), explicitForm.matrixPtr() );
}

BOOST_AUTO_TEST_CASE( facet_trace_normal_helpers_match_explicit_hdg_notation )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto facetMesh = createSubmesh( _mesh = mesh, _range = faces( mesh ), _update = 0 );
    auto Vh = Pdhv<1>( mesh );
    auto Mh = Pdh<1>( facetMesh );

    auto phat = trial( Mh );
    auto mu = test( Mh );
    auto u = trial( Vh );
    auto v = test( Vh );

    auto helperFluxTrace = form2( _test = Vh, _trial = Mh );
    helperFluxTrace +=
        integrate( _range = facet::interior( mesh ),
                   _expr = facet::trace( phat ) * facet::normalSum( v ) )
        + integrate( _range = facet::boundary( mesh ),
                     _expr = facet::trace( phat ) * facet::normal( v ) );
    helperFluxTrace.close();

    auto explicitFluxTrace = form2( _test = Vh, _trial = Mh );
    explicitFluxTrace +=
        integrate( _range = internalfaces( mesh ),
                   _expr = phat * ( leftface( normal( v ) ) +
                                    rightface( normal( v ) ) ) )
        + integrate( _range = boundaryfaces( mesh ),
                     _expr = phat * normal( v ) );
    explicitFluxTrace.close();

    checkMatrixClose( helperFluxTrace.matrixPtr(), explicitFluxTrace.matrixPtr() );

    auto helperTraceFlux = form2( _test = Mh, _trial = Vh );
    helperTraceFlux +=
        integrate( _range = facet::interior( mesh ),
                   _expr = facet::trace( mu ) * facet::normalSum( u ) )
        + integrate( _range = facet::boundary( mesh ),
                     _expr = facet::trace( mu ) * facet::normal( u ) );
    helperTraceFlux.close();

    auto explicitTraceFlux = form2( _test = Mh, _trial = Vh );
    explicitTraceFlux +=
        integrate( _range = internalfaces( mesh ),
                   _expr = mu * ( leftfacet( normal( u ) ) +
                                  rightfacet( normal( u ) ) ) )
        + integrate( _range = boundaryfaces( mesh ),
                     _expr = mu * normal( u ) );
    explicitTraceFlux.close();

    checkMatrixClose( helperTraceFlux.matrixPtr(), explicitTraceFlux.matrixPtr() );
}

BOOST_AUTO_TEST_SUITE_END()
