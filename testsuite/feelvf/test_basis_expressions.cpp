/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE feelvf_basis_expressions
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( feelvf_basis_expressions_suite )

BOOST_AUTO_TEST_CASE( test_and_trial_basis_match_supported_id_operators )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    using space_type = std::decay_t<decltype( *Xh )>;
    using fe_type = typename space_type::reference_element_type;

    auto u = trial( Xh );
    auto v = test( Xh );
    auto uOld = Xh->element( "u_old" );
    auto vOld = Xh->element( "v_old" );

    static_assert( decltype( u )::template has_trial_basis<fe_type> );
    static_assert( !decltype( u )::template has_test_basis<fe_type> );
    static_assert( decltype( v )::template has_test_basis<fe_type> );
    static_assert( !decltype( v )::template has_trial_basis<fe_type> );
    static_assert( is_trial_basis_v<decltype( u )> );
    static_assert( is_test_basis_v<decltype( v )> );

    BOOST_CHECK( u.space() == Xh );
    BOOST_CHECK( v.space() == Xh );
    BOOST_CHECK_EQUAL( nelements( elements( u.support() ), true ),
                       nelements( elements( support( Xh ) ), true ) );
    BOOST_CHECK_EQUAL( nelements( elements( v.support() ), true ),
                       nelements( elements( support( Xh ) ), true ) );

    auto backend = Feel::backend( _worldcomm = Xh->worldCommPtr() );
    auto newMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto oldMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto newVector = backend->newVector( _test = Xh );
    auto oldVector = backend->newVector( _test = Xh );

    form2( _test = Xh, _trial = Xh, _matrix = newMatrix ) +=
        integrate( _range = elements( support( Xh ) ),
                   _expr = inner( grad( u ), grad( v ) ) + u * v );
    form2( _test = Xh, _trial = Xh, _matrix = oldMatrix ) +=
        integrate( _range = elements( support( Xh ) ),
                   _expr = inner( gradt( uOld ), grad( vOld ) ) + idt( uOld ) * id( vOld ) );

    form1( _test = Xh, _vector = newVector ) +=
        integrate( _range = elements( support( Xh ) ),
                   _expr = cst( 2.0 ) * v );
    form1( _test = Xh, _vector = oldVector ) +=
        integrate( _range = elements( support( Xh ) ),
                   _expr = cst( 2.0 ) * id( vOld ) );

    newMatrix->close();
    oldMatrix->close();
    newVector->close();
    oldVector->close();

    auto matrixDelta = newMatrix->clone();
    matrixDelta->addMatrix( -1.0, oldMatrix );
    matrixDelta->close();
    BOOST_CHECK_SMALL( matrixDelta->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( matrixDelta->linftyNorm(), 1e-9 );

    auto vectorDelta = vf::cloneVectorWithValues( newVector );
    vectorDelta->add( -1.0, *oldVector );
    vectorDelta->close();
    BOOST_CHECK_SMALL( vectorDelta->linftyNorm(), 1e-10 );
}

BOOST_AUTO_TEST_CASE( basis_expression_lifetime_is_independent_of_local_variables )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uOld = Xh->element( "u_old" );
    auto vOld = Xh->element( "v_old" );
    auto backend = Feel::backend( _worldcomm = Xh->worldCommPtr() );

    auto delayedTerm = [Xh]()
    {
        auto u = trial( Xh );
        auto v = test( Xh );
        return integrate( _range = elements( support( Xh ) ),
                          _expr = inner( grad( u ), grad( v ) ) + u * v );
    }();

    auto newMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto oldMatrix = backend->newMatrix( _test = Xh, _trial = Xh );

    form2( _test = Xh, _trial = Xh, _matrix = newMatrix ) += delayedTerm;
    form2( _test = Xh, _trial = Xh, _matrix = oldMatrix ) +=
        integrate( _range = elements( support( Xh ) ),
                   _expr = inner( gradt( uOld ), grad( vOld ) ) + idt( uOld ) * id( vOld ) );

    newMatrix->close();
    oldMatrix->close();

    auto matrixDelta = newMatrix->clone();
    matrixDelta->addMatrix( -1.0, oldMatrix );
    matrixDelta->close();
    BOOST_CHECK_SMALL( matrixDelta->l1Norm(), 1e-9 );
}

BOOST_AUTO_TEST_CASE( field_value_adapter_matches_supported_idv_operators )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "uh" );
    uh.on( _range = elements( support( Xh ) ),
           _expr = cst( 1.0 ) + Px() + Py() );
    uh.close();

    auto rawIntegral = integrate( _range = elements( support( Xh ) ),
                                  _expr = uh ).evaluate()( 0, 0 );
    auto idvIntegral = integrate( _range = elements( support( Xh ) ),
                                  _expr = idv( uh ) ).evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( rawIntegral - idvIntegral, 1e-12 );

    auto rawNorm = normL2( _range = elements( support( Xh ) ),
                           _expr = uh );
    auto idvNorm = normL2( _range = elements( support( Xh ) ),
                           _expr = idv( uh ) );
    BOOST_CHECK_SMALL( rawNorm - idvNorm, 1e-12 );

    auto projected = project( _space = Xh,
                              _range = elements( support( Xh ) ),
                              _expr = uh );
    auto projectionError = normL2( _range = elements( support( Xh ) ),
                                   _expr = projected - uh );
    BOOST_CHECK_SMALL( projectionError, 1e-12 );

    auto absError = normL2( _range = elements( support( Xh ) ),
                            _expr = abs( uh ) - abs( idv( uh ) ) );
    auto powError = normL2( _range = elements( support( Xh ) ),
                            _expr = pow( uh, 2 ) - pow( idv( uh ), 2 ) );
    BOOST_CHECK_SMALL( absError, 1e-12 );
    BOOST_CHECK_SMALL( powError, 1e-12 );

    auto valueError = normL2( _range = elements( support( Xh ) ),
                              _expr = uh - idv( uh ) );
    auto gradError = normL2( _range = elements( support( Xh ) ),
                             _expr = grad( value( uh ) ) - gradv( uh ) );

    BOOST_CHECK_SMALL( valueError, 1e-12 );
    BOOST_CHECK_SMALL( gradError, 1e-12 );
}

BOOST_AUTO_TEST_CASE( two_valued_operators_accept_first_class_basis )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pdh<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto uOld = Xh->element( "u_old" );
    auto vOld = Xh->element( "v_old" );

    auto backend = Feel::backend( _worldcomm = Xh->worldCommPtr() );
    auto newMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto oldMatrix = backend->newMatrix( _test = Xh, _trial = Xh );

    form2( _test = Xh, _trial = Xh, _matrix = newMatrix ) +=
        integrate( _range = internalfaces( mesh ),
                   _expr = leftface( v ) * leftfacet( u ) +
                           rightface( v ) * rightfacet( u ) )
        + integrate( _range = boundaryfaces( mesh ),
                     _expr = v * u );

    form2( _test = Xh, _trial = Xh, _matrix = oldMatrix ) +=
        integrate( _range = internalfaces( mesh ),
                   _expr = leftface( id( vOld ) ) * leftfacet( idt( uOld ) ) +
                           rightface( id( vOld ) ) * rightfacet( idt( uOld ) ) )
        + integrate( _range = boundaryfaces( mesh ),
                     _expr = id( vOld ) * idt( uOld ) );

    newMatrix->close();
    oldMatrix->close();

    auto matrixDelta = newMatrix->clone();
    matrixDelta->addMatrix( -1.0, oldMatrix );
    matrixDelta->close();
    BOOST_CHECK_SMALL( matrixDelta->l1Norm(), 1e-9 );
    BOOST_CHECK_SMALL( matrixDelta->linftyNorm(), 1e-9 );
}

BOOST_AUTO_TEST_CASE( two_valued_operators_promote_function_space_fields )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pdh<1>( mesh );
    auto uh = Xh->element( "uh" );
    uh.on( _range = elements( support( Xh ) ),
           _expr = cst( 1.0 ) + Px() + 2.0 * Py(),
           _close = true );

    auto leftRaw = integrate( _range = internalfaces( mesh ),
                              _expr = leftface( uh ) ).evaluate()( 0, 0 );
    auto leftValue = integrate( _range = internalfaces( mesh ),
                                _expr = leftface( value( uh ) ) ).evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( leftRaw - leftValue, 1e-12 );

    auto rightRaw = integrate( _range = internalfaces( mesh ),
                               _expr = rightface( uh ) ).evaluate()( 0, 0 );
    auto rightValue = integrate( _range = internalfaces( mesh ),
                                 _expr = rightface( value( uh ) ) ).evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( rightRaw - rightValue, 1e-12 );

    auto averageRaw = integrate( _range = internalfaces( mesh ),
                                 _expr = average( uh ) ).evaluate()( 0, 0 );
    auto averageValue = integrate( _range = internalfaces( mesh ),
                                   _expr = average( value( uh ) ) ).evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( averageRaw - averageValue, 1e-12 );

    auto leftValueBranchRaw = integrate( _range = internalfaces( mesh ),
                                         _expr = leftfacev( uh ) ).evaluate()( 0, 0 );
    auto leftValueBranch = integrate( _range = internalfaces( mesh ),
                                      _expr = leftfacev( value( uh ) ) ).evaluate()( 0, 0 );
    BOOST_CHECK_SMALL( leftValueBranchRaw - leftValueBranch, 1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
