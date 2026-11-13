/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file
 * @brief FeelSystem regression test.
 */

#define BOOST_TEST_MODULE feelsystem_nonlinear
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelpde/modeldescriptors.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/dirichletconstraints.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

namespace
{

Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "FeelSystem nonlinear test options" );
    options.add( Feel::backend_options( "nonlinear-reaction" ) )
           .add( Feel::backend_options( "nonlinear-affine" ) )
           .add( Feel::backend_options( "nonlinear-metadata" ) )
           .add( Feel::backend_options( "picard-reaction" ) );
    return options;
}

Feel::nl::json
makeNonlinearModelJson()
{
    return Feel::nl::json::parse(
        R"json(
        {
            "Name": "nonlinear-metadata",
            "Materials": {
                "Mat1": {
                    "markers": "Mat1",
                    "physics": [ "thermal" ],
                    "k": "2"
                }
            },
            "BoundaryConditions": {
                "u": {
                    "Neumann": {
                        "flux": {
                            "markers": "Neumann",
                            "expr": "1"
                        }
                    },
                    "Radiation": {
                        "outer": {
                            "markers": "Neumann",
                            "expr1": "5.67e-8",
                            "expr2": "1"
                        }
                    }
                }
            }
        }
        )json" );
}

template <typename ElementType>
void
setConstantClosed( ElementType& element, double value )
{
    element.setConstant( value );
    element.close();
}

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_nonlinear" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_nonlinear_suite )

BOOST_AUTO_TEST_CASE( newton_solves_state_dependent_residual )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    setConstantClosed( uh, 0.5 );

    auto problem =
        nonlinearProblem( _space = Xh, _name = "nonlinear-reaction" )
            .reaction( integrate( _range = elements( mesh ),
                                  _expr = u * v ) )
            .source( integrate( _range = elements( mesh ),
                                _expr = v ) )
            .residual( [&]( auto const& state, auto& residual )
            {
                residual += integrate( _range = elements( mesh ),
                                       _expr = 0.1 * pow( idv( state ), 3 ) * v );
            } )
            .jacobian( [&]( auto const& state, auto& jacobian )
            {
                jacobian += integrate( _range = elements( mesh ),
                                      _expr = 0.3 * pow( idv( state ), 2 ) * u * v );
            } );

    static_assert( systems::NonlinearProblemLike<decltype( problem )> );
    BOOST_CHECK_EQUAL( problem.residualTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.jacobianTerms().size(), 1 );

    auto system = problem.assemble();
    BOOST_CHECK( system.assembled() );
    BOOST_CHECK_EQUAL( system.backendName(), problem.name() );
    BOOST_CHECK( system.jacobianPtr() );
    BOOST_CHECK( system.residualPtr() );

    auto result = system.solve( _solution = uh, _maxit = 25, _rtolerance = 1e-12 );
    BOOST_CHECK( result.converged );

    auto error = Xh->element( "error" );
    error.on( _range = elements( mesh ), _expr = abs( uh - cst( 0.921698994204679 ) ) );
    BOOST_CHECK_SMALL( error.max(), 1e-6 );
}

BOOST_AUTO_TEST_CASE( affine_helpers_contribute_to_nonlinear_storage )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto reactionTerm = integrate( _range = elements( mesh ),
                                   _expr = u * v );
    auto sourceTerm = integrate( _range = elements( mesh ),
                                 _expr = v );
    auto neumannFlux = cst( 2.0 );
    auto robinAlpha = cst( 3.0 );
    auto robinValue = cst( 4.0 );
    auto fourierCoeff = cst( 5.0 );
    auto exteriorValue = cst( 6.0 );

    auto problem =
        nonlinearProblem( _space = Xh, _name = "nonlinear-affine" )
            .reaction( reactionTerm )
            .source( sourceTerm )
            .neumann( _range = boundaryfaces( mesh ),
                      _test = v,
                      _expr = neumannFlux )
            .robin( _range = boundaryfaces( mesh ),
                    _trial = u,
                    _test = v,
                    _alpha = robinAlpha,
                    _expr = robinValue )
            .fourier( _range = boundaryfaces( mesh ),
                      _trial = u,
                      _test = v,
                      _alpha = fourierCoeff,
                      _expr = exteriorValue );

    auto system = problem.assemble();
    auto backend = system.backend();
    auto directMatrix = backend->newMatrix( _test = Xh, _trial = Xh );
    auto directRhs = backend->newVector( _test = Xh );
    directMatrix->zero();
    directRhs->zero();
    form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) += reactionTerm;
    form1( _test = Xh, _vector = directRhs ) += sourceTerm;
    form1( _test = Xh, _vector = directRhs ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = neumannFlux * v );
    form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = robinAlpha * u * v );
    form1( _test = Xh, _vector = directRhs ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = robinValue * v );
    form2( _test = Xh, _trial = Xh, _matrix = directMatrix ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = fourierCoeff * u * v );
    form1( _test = Xh, _vector = directRhs ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = fourierCoeff * exteriorValue * v );
    directMatrix->close();
    directRhs->close();

    auto matrixDelta = system.affineMatrixPtr()->clone();
    matrixDelta->addMatrix( -1.0, directMatrix );
    matrixDelta->close();
    BOOST_CHECK_SMALL( matrixDelta->linftyNorm(), 1e-9 );

    auto rhsDelta = vf::cloneVectorWithValues( system.affineRhsPtr() );
    rhsDelta->add( -1.0, *directRhs );
    rhsDelta->close();
    BOOST_CHECK_SMALL( rhsDelta->linftyNorm(), 1e-10 );
}

BOOST_AUTO_TEST_CASE( nonlinear_affine_terms_solve_with_newton )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    setConstantClosed( uh, 0.8 );

    auto problem =
        nonlinearProblem( _space = Xh, _name = "nonlinear-reaction" )
            .reaction( integrate( _range = elements( mesh ),
                                  _expr = u * v ) )
            .source( integrate( _range = elements( mesh ),
                                _expr = v ) )
            .residual( [&]( auto const& state, auto& residual )
            {
                residual += integrate( _range = elements( mesh ),
                                       _expr = 0.1 * pow( idv( state ), 3 ) * v );
            } )
            .jacobian( [&]( auto const& state, auto& jacobian )
            {
                jacobian += integrate( _range = elements( mesh ),
                                      _expr = 0.3 * pow( idv( state ), 2 ) * u * v );
            } );

    auto result = problem.solve( _solution = uh, _maxit = 25, _rtolerance = 1e-12 );
    BOOST_CHECK( result.converged );

    auto error = Xh->element( "error" );
    error.on( _range = elements( mesh ), _expr = abs( uh - cst( 0.921698994204679 ) ) );
    BOOST_CHECK_SMALL( error.max(), 1e-6 );
}

BOOST_AUTO_TEST_CASE( picard_problem_reassembles_each_iteration )
{
    using namespace Feel;
    using namespace Feel::vf;

    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    setConstantClosed( uh, 0.5 );

    std::size_t assemblyCount = 0;
    auto problem =
        picardProblem( _space = Xh, _name = "picard-reaction" )
            .iteration( [&]( auto const& uk, auto& a, auto& l )
            {
                ++assemblyCount;
                a += integrate( _range = elements( mesh ),
                                _expr = ( 1.0 + 0.1 * pow( idv( uk ), 2 ) ) * u * v );
                l += integrate( _range = elements( mesh ),
                                _expr = v );
            } );

    static_assert( systems::PicardProblemLike<decltype( problem )> );
    BOOST_CHECK_EQUAL( problem.iterationTerms().size(), 1 );

    auto result = problem.solve( _solution = uh,
                                 _strategy = systems::picard( _maxit = 100,
                                                              _rtolerance = 1e-10 ) );
    BOOST_CHECK( result.converged );
    BOOST_CHECK_GT( assemblyCount, 1 );

    auto error = Xh->element( "error" );
    error.on( _range = elements( mesh ), _expr = abs( uh - cst( 0.921698994204679 ) ) );
    BOOST_CHECK_SMALL( error.max(), 1e-5 );
}

BOOST_AUTO_TEST_CASE( nonlinear_problem_uses_material_and_boundary_iteration )
{
    using namespace Feel;
    using namespace Feel::vf;

    using mesh_type = Mesh<Simplex<2>>;
    auto mesh = createGMSHMesh( _mesh = new mesh_type,
                                _desc = domain( _name = "feelsystem-nonlinear-metadata",
                                                 _shape = "hypercube",
                                                 _dim = 2,
                                                 _h = 0.5,
                                                 _usenames = true,
                                                 _addmidpoint = false ) );

    ModelProperties props;
    props.enableBoundaryConditions2();
    props.setup( makeNonlinearModelJson() );

    auto domains = pde::materials( _model = props, _mesh = mesh );
    auto bcs = pde::boundaryConditions( _model = props, _name = "u", _mesh = mesh );

    auto Xh = Pch<1>( mesh );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        nonlinearProblem( _space = Xh, _name = "nonlinear-metadata" )
            .forEach( domains.hasProperty( "k" ), [&]( auto& p, auto const& mat )
            {
                auto k = mat.property( "k" ).template expr<1, 1>();
                p.residual( [&, range = mat.elements(), k]( auto const& state, auto& residual )
                {
                    residual += integrate( _range = range,
                                           _expr = k * pow( idv( state ), 3 ) * v );
                } );
                p.jacobian( [&, range = mat.elements(), k]( auto const& state, auto& jacobian )
                {
                    jacobian += integrate( _range = range,
                                          _expr = 3.0 * k * pow( idv( state ), 2 ) * u * v );
                } );
            } )
            .forEach( bcs.neumann(), [&]( auto& p, auto const& bc )
            {
                p.neumann( _range = bc.faces(),
                           _test = v,
                           _expr = bc.template expr<1, 1>() );
            } )
            .forEach( bcs.type( "Radiation" ), [&]( auto& p, auto const& bc )
            {
                auto epssigma = bc.template expr1<1, 1>();
                auto exterior = bc.template expr2<1, 1>();
                p.residual( [&, range = bc.faces(), epssigma, exterior]( auto const& state, auto& residual )
                {
                    residual += integrate( _range = range,
                                           _expr = epssigma *
                                                   ( pow( idv( state ), 4 ) - pow( exterior, 4 ) ) *
                                                   v );
                } );
                p.jacobian( [&, range = bc.faces(), epssigma]( auto const& state, auto& jacobian )
                {
                    jacobian += integrate( _range = range,
                                          _expr = 4.0 * epssigma * pow( idv( state ), 3 ) *
                                                  u * v );
                } );
            } );

    BOOST_CHECK_EQUAL( domains.hasProperty( "k" ).size(), 1 );
    BOOST_CHECK_EQUAL( bcs.neumann().size(), 1 );
    BOOST_CHECK_EQUAL( bcs.type( "Radiation" ).size(), 1 );
    BOOST_CHECK_EQUAL( problem.residualTerms().size(), 2 );
    BOOST_CHECK_EQUAL( problem.jacobianTerms().size(), 2 );
    BOOST_CHECK_EQUAL( problem.rhsTerms().size(), 1 );
}

BOOST_AUTO_TEST_SUITE_END()
