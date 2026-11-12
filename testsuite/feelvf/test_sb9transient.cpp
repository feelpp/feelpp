/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE test_sb9transient
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelts/newmark.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/sb9_inertia.hpp>
#include <feel/feelvf/sb9_strain.hpp>

#include <Eigen/Cholesky>

#include "test_sb9_common.hpp"

using namespace Feel;
using namespace Feel::vf;
using namespace Feel::Tests::SB9;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template <typename Matrix>
Eigen::MatrixXd denseMatrix( Matrix const& matrix )
{
    Eigen::MatrixXd result( matrix.size1(), matrix.size2() );
    for ( size_type i = 0; i < matrix.size1(); ++i )
        for ( size_type j = 0; j < matrix.size2(); ++j )
            result( i, j ) = matrix( i, j );
    return result;
}

template <typename Vector>
Eigen::VectorXd denseVector( Vector const& vector )
{
    Eigen::VectorXd result( vector.size() );
    for ( size_type i = 0; i < vector.size(); ++i )
        result( i ) = vector( i );
    return result;
}

// Exercise the production inertia assembly and Newmark histories against an
// independent dense predictor/corrector update. A small elastic foundation
// removes rigid/hourglass modes from this one-cell time-integration test.
void checkTransientCondensation( bool consistent, bool pointLoaded = false )
{
    auto mesh = createAxisAlignedUnitPatch( consistent ? "sb9_transient_consistent" : "sb9_transient_massless" );
    if ( pointLoaded )
    {
        // Unmarked vertices do not acquire incident-element information in
        // the default mesh update; physical point markers normally supply it.
        for ( auto const& cell : elements( mesh ) )
        {
            auto const& elt = unwrap_ref( cell );
            for ( uint16_type p = 0; p < elt.nPoints(); ++p )
                mesh->point( elt.point( p ).id() ).addElement( elt.id(), p );
        }
    }
    auto Uh = Pchv<Dynamic>( mesh, RuntimeOrder{1} );
    auto Ah = Pdh<Dynamic>( mesh, RuntimeOrder{0} );
    auto ps = product( Uh, Ah );
    auto u = trial( Uh );
    auto v = test( Uh );
    auto alpha = trial( Ah );
    auto beta = test( Ah );
    auto C = isotropic_stiffness<3>( 1.0, 1.0 );
    auto epsU = sb9ShellStrain( u, zeta(), cst( 1.0 ) );
    auto epsV = sb9ShellStrain( v, zeta(), cst( 1.0 ) );
    auto epsAlpha = sb9W9( alpha );
    auto epsBeta = sb9W9( beta );
    auto assembleStiffness = [&]( auto& a )
    {
        a( 0_c, 0_c ) += integrate( _range=elements( mesh ), _quad=_Q<6>(),
                                    _expr=ddot( C, epsU, epsV ) + inner( id( u ), id( v ) ) );
        a( 0_c, 1_c ) += integrate( _range=elements( mesh ), _quad=_Q<6>(),
                                    _expr=ddot( C, epsAlpha, epsV ) );
        a( 1_c, 0_c ) += integrate( _range=elements( mesh ), _quad=_Q<6>(),
                                    _expr=ddot( C, epsU, epsBeta ) );
        a( 1_c, 1_c ) += integrate( _range=elements( mesh ), _quad=_Q<6>(),
                                    _expr=ddot( C, epsAlpha, epsBeta ) );
    };
    auto assembleLoad = [&]( auto& l, double time )
    {
        if ( pointLoaded )
            l( 0_c ) += integrate( _range=points( mesh ),
                                    _expr=time * inner( vec( cst( 0.2 ), cst( -0.1 ), cst( 1.0 ) ), id( v ) ) );
        else
            l( 0_c ) += integrate( _range=elements( mesh ),
                                    _expr=time * inner( vec( cst( 0.2 ), cst( -0.1 ), Pz() ), id( v ) ) );
        l( 1_c ) += integrate( _range=elements( mesh ), _expr=time * id( beta ) );
    };

    auto stiffness = blockform2( ps, solve::strategy::monolithic, backend() );
    assembleStiffness( stiffness );
    stiffness.close();
    auto K = denseMatrix( *stiffness.matrixPtr() );
    auto mass = blockform2( ps, solve::strategy::monolithic, backend() );
    assembleSb9Inertia( mass, elements( mesh ), u, v, alpha, beta, 2.0, 1.0, consistent );
    mass.close();
    auto M = denseMatrix( *mass.matrixPtr() );
    auto load = blockform1( ps, solve::strategy::monolithic, backend() );
    assembleLoad( load, 1.0 );
    load.close();
    auto F = denseVector( *load.vectorPtr()->getVector() );
    if ( pointLoaded )
        BOOST_CHECK_CLOSE( F.head( Uh->nDof() ).sum(), 8.8, 1e-8 );

    constexpr double dt = 0.025;
    constexpr double finalTime = 0.4;
    constexpr double loadingTime = 0.2;
    constexpr double newmarkBeta = 0.25;
    constexpr double newmarkGamma = 0.5;
    double const a0 = 1.0 / ( newmarkBeta * dt * dt );
    auto run = [&]( bool condense )
    {
        auto strategy = condense ? solve::strategy::static_condensation : solve::strategy::monolithic;
        auto a = blockform2( ps, strategy, backend() );
        auto l = blockform1( ps, strategy, backend() );
        auto X = ps.element();
        auto uh = Uh->element();
        auto alphah = Ah->element();
        uh.zero();
        alphah.zero();
        auto uTime = newmark( _space=Uh, _name="sb9-u-test", _initial_time=0.0,
                             _final_time=finalTime, _time_step=dt, _save=false, _restart=false,
                             _steady=false, _beta=newmarkBeta, _gamma=newmarkGamma );
        auto alphaTime = newmark( _space=Ah, _name="sb9-alpha-test", _initial_time=0.0,
                                 _final_time=finalTime, _time_step=dt, _save=false, _restart=false,
                                 _steady=false, _beta=newmarkBeta, _gamma=newmarkGamma );
        uTime->start( uh );
        alphaTime->start( alphah );
        assembleStiffness( a );
        assembleSb9Inertia( a, elements( mesh ), u, v, alpha, beta, 2.0, a0, consistent );
        Eigen::VectorXd q = Eigen::VectorXd::Zero( K.rows() );
        Eigen::VectorXd velocity = q;
        Eigen::VectorXd acceleration = q;
        auto effectiveMass = ( M + newmarkBeta * dt * dt * K ).eval();
        auto factorization = effectiveMass.ldlt();
        double maxAlpha = 0.0;
        double freeEnergy = -1.0;
        int steps = 0;
        while ( !uTime->isFinished() )
        {
            double const time = uTime->time();
            double const loadFactor = time < loadingTime
                                          ? std::sin( 2.0 * std::acos( -1.0 ) * time / loadingTime )
                                          : 0.0;
            l.zero();
            assembleSb9InertiaHistory( l, elements( mesh ), v, beta,
                                       uTime->polyDeriv(), alphaTime->polyDeriv(), 2.0, consistent );
            assembleLoad( l, loadFactor );
            l.close();
            a.close();
            auto status = a.solve( _solution=X, _rhs=l, _condense=condense, _condenser=condenser_sb9() );
            BOOST_REQUIRE( status.isConverged() );

            Eigen::VectorXd predictor = q + dt * velocity + dt * dt * ( 0.5 - newmarkBeta ) * acceleration;
            Eigen::VectorXd nextAcceleration = factorization.solve( loadFactor * F - K * predictor );
            Eigen::VectorXd nextQ = predictor + newmarkBeta * dt * dt * nextAcceleration;
            velocity += dt * ( ( 1.0 - newmarkGamma ) * acceleration + newmarkGamma * nextAcceleration );
            q = nextQ;
            acceleration = nextAcceleration;
            auto values = backend()->newBlockVector( _block=X, _copy_values=true );
            BOOST_CHECK_SMALL( ( denseVector( *values ) - q ).norm(), 2e-9 );

            uh = X( 0_c );
            alphah = X( 1_c );
            uTime->next( uh );
            alphaTime->next( alphah );
            auto derivatives = ps.element();
            derivatives( 0_c ) = uTime->currentAcceleration();
            derivatives( 1_c ) = alphaTime->currentAcceleration();
            auto accelerations = backend()->newBlockVector( _block=derivatives, _copy_values=true );
            if ( consistent )
                BOOST_CHECK_SMALL( ( denseVector( *accelerations ) - acceleration ).norm(), 2e-7 );
            derivatives( 0_c ) = uTime->currentVelocity();
            derivatives( 1_c ) = alphaTime->currentVelocity();
            auto velocities = backend()->newBlockVector( _block=derivatives, _copy_values=true );
            auto actualQ = denseVector( *values );
            auto actualVelocity = denseVector( *velocities );
            BOOST_CHECK_SMALL( ( actualVelocity - velocity ).norm(), 2e-8 );
            double const energy = 0.5 * ( actualVelocity.dot( M * actualVelocity ) + actualQ.dot( K * actualQ ) );
            if ( time >= loadingTime )
            {
                if ( freeEnergy < 0.0 )
                    freeEnergy = energy;
                BOOST_CHECK_SMALL( energy - freeEnergy, 1e-9 );
            }
            maxAlpha = std::max( maxAlpha, std::abs( q( q.size() - 1 ) ) );
            ++steps;
        }
        BOOST_CHECK_GE( steps, 15 );
        BOOST_CHECK_GT( maxAlpha, 1e-6 );
        BOOST_CHECK_GT( freeEnergy, 1e-8 );
        return q;
    };
    auto monolithic = run( false );
    auto condensed = run( true );
    BOOST_CHECK_SMALL( ( monolithic - condensed ).norm(), 2e-9 );
}
} // namespace

BOOST_AUTO_TEST_CASE( sb9_consistent_transient_condensation_preserves_internal_history )
{
    checkTransientCondensation( true );
}

BOOST_AUTO_TEST_CASE( sb9_massless_transient_condensation_matches_dense_reference )
{
    checkTransientCondensation( false );
}

BOOST_AUTO_TEST_CASE( sb9_point_loaded_transient_condensation_matches_dense_reference )
{
    checkTransientCondensation( true, true );
}
