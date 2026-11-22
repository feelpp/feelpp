/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file
 * @brief FeelSystem quickstart example.
 */

#include "../discretization.hpp"
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

int
main( int argc, char** argv )
{
    using namespace Feel;
    using namespace Feel::vf;

    po::options_description options( "FeelSystem nonlinear quickstart options" );
    options.add( backend_options( "qs-sys-nonlinear" ) )
           .add( backend_options( "qs-sys-picard" ) );
    options.add_options()
        ( "order", po::value<int>()->default_value( 1 ), "approximation order" )
        ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" );

    Environment env( _argc = argc, _argv = argv,
                     _desc = options,
                     _about = about( _name = "qs_sys_nonlinear",
                                     _author = "University of Strasbourg",
                                     _email = "feelpp-devel@feelpp.org" ) );

    using mesh_type = Mesh<Simplex<FEELPP_DIM, Dynamic>>;
    auto mesh = createGMSHMesh( _mesh = new mesh_type( quickstart::geometryOrder() ),
                                _desc = domain( _name = "qs-sys-nonlinear-domain",
                                                 _shape = "hypercube",
                                                 _dim = FEELPP_DIM,
                                                 _h = 0.5,
                                                 _usenames = true,
                                                 _addmidpoint = false ) );

    quickstart::checkMeshOrder( mesh );
    auto Xh = Pch<Dynamic>( mesh, quickstart::checkedOrder( "order" ) );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    uh.setConstant( 0.8 );
    uh.close();

    auto newtonProblem =
        nonlinearProblem( _space = Xh, _name = "qs-sys-nonlinear" )
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

    auto newtonResult = newtonProblem.solve( _solution = uh,
                                             _maxit = 25,
                                             _rtolerance = 1e-12,
                                             _strategy = systems::newton() );

    auto up = Xh->element( "up" );
    auto upTrial = trial( Xh );
    up.setConstant( 0.5 );
    up.close();

    auto picardProblemExample =
        picardProblem( _space = Xh, _name = "qs-sys-picard" )
            .iteration( [&]( auto const& uk, auto& a, auto& l )
            {
                a += integrate( _range = elements( mesh ),
                                _expr = ( 1.0 + 0.1 * pow( idv( uk ), 2 ) ) * upTrial * v );
                l += integrate( _range = elements( mesh ),
                                _expr = v );
            } );

    auto picardResult = picardProblemExample.solve( _solution = up,
                                                   _strategy = systems::picard( _maxit = 100,
                                                                                _rtolerance = 1e-10 ) );
    double const newtonMax = uh.max();
    double const picardMax = up.max();

    if ( Environment::isMasterRank() )
    {
        Feel::cout << "Newton converged: " << std::boolalpha << newtonResult.converged
                   << ", iterations: " << newtonResult.iterations
                   << ", max(u): " << newtonMax << "\n";
        Feel::cout << "Picard converged: " << std::boolalpha << picardResult.converged
                   << ", iterations: " << picardResult.iterations
                   << ", max(u): " << picardMax << "\n";
    }

    return ( newtonResult.converged && picardResult.converged ) ? 0 : 1;
}
