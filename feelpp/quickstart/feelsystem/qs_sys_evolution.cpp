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
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

int
main( int argc, char** argv )
{
    using namespace Feel;
    using namespace Feel::vf;

    po::options_description options( "FeelSystem evolution quickstart options" );
    options.add( backend_options( "qs-sys-heat-step" ) )
           .add( bdf_options( "qs-sys-heat-time" ) )
           .add( ts_options( "qs-sys-heat-time" ) );
    options.add_options()
        ( "order", po::value<int>()->default_value( 1 ), "approximation order" )
        ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" );

    Environment env( _argc = argc, _argv = argv,
                     _desc = options,
                     _about = about( _name = "qs_sys_evolution",
                                     _author = "University of Strasbourg",
                                     _email = "feelpp-devel@feelpp.org" ) );

    auto mesh = loadMesh( _mesh = new Mesh<Simplex<FEELPP_DIM, Dynamic>>( quickstart::geometryOrder() ) );
    quickstart::checkMeshOrder( mesh );
    auto Xh = Pch<Dynamic>( mesh, quickstart::checkedOrder( "order" ) );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    uh.zero();
    auto source = expr( "1+t:t" );

    auto heatTime = bdf( _space = Xh,
                         _name = "qs-sys-heat-time",
                         _prefix = "qs-sys-heat-time",
                         _order = 1,
                         _initial_time = 0.0,
                         _final_time = 0.2,
                         _time_step = 0.1,
                         _save = false );

    auto heat =
        parabolicProblem( _space = Xh, _name = "qs-sys-heat-step" )
            .timeScheme( heatTime )
            .mass( integrate( _range = elements( mesh ),
                              _expr = u * v ) )
            .reaction( integrate( _range = elements( mesh ),
                                  _expr = u * v ) )
            .rhs( _range = elements( mesh ),
                  _expr = source * v );

    auto evolution = heat.evolve();
    bool converged = true;
    int steps = 0;
    systems::solve_result lastResult;

    for ( evolution.start( uh ); !evolution.isFinished(); evolution.advance( uh ) )
    {
        lastResult = evolution.step().solve( _solution = uh );
        converged = converged && lastResult.converged;
        ++steps;
    }
    double const solutionMax = uh.max();

    if ( Environment::isMasterRank() )
    {
        Feel::cout << "Problem: " << heat.name() << "\n";
        Feel::cout << "time backend: " << heatTime->prefix() << "\n";
        Feel::cout << "step backend: " << evolution.backendName() << "\n";
        Feel::cout << "steps: " << steps << "\n";
        Feel::cout << "converged: " << std::boolalpha << converged << "\n";
        Feel::cout << "last iterations: " << lastResult.iterations << "\n";
        Feel::cout << "max(u): " << solutionMax << "\n";
    }

    return converged ? 0 : 1;
}
