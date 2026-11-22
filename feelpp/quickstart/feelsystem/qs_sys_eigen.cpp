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

#include <cstddef>

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

    po::options_description options( "FeelSystem eigen quickstart options" );
    options.add( backend_options( "laplacian-eigen" ) )
           .add( solvereigen_options( "laplacian-eigen" ) );
    options.add_options()
        ( "order", po::value<int>()->default_value( 1 ), "approximation order" )
        ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" );

    Environment env( _argc = argc, _argv = argv,
                     _desc = options,
                     _about = about( _name = "qs_sys_eigen",
                                     _author = "University of Strasbourg",
                                     _email = "feelpp-devel@feelpp.org" ) );

    auto mesh = loadMesh( _mesh = new Mesh<Simplex<FEELPP_DIM, Dynamic>>( quickstart::geometryOrder() ) );
    quickstart::checkMeshOrder( mesh );
    auto Xh = Pch<Dynamic>( mesh, quickstart::checkedOrder( "order" ) );
    auto u = trial( Xh );
    auto v = test( Xh );

    auto problem =
        generalizedEigenvalueProblem( _space = Xh, _name = "laplacian-eigen" )
            .operatorA( integrate( _range = elements( mesh ),
                                   _expr = inner( grad( u ), grad( v ) ) + u * v ) )
            .operatorB( integrate( _range = elements( mesh ),
                                   _expr = u * v ) )
            .nev( 5 )
            .ncv( 10 )
            .spectrum( systems::smallestMagnitude{} )
            .transform( systems::shiftInvert{} )
            .tolerance( 1e-10 );

    auto system = problem.assemble();
    auto result = system.solve( _strategy = systems::eigen() );

    if ( Environment::isMasterRank() )
    {
        Feel::cout << "Problem: " << problem.name() << "\n";
        Feel::cout << "A terms: " << problem.operatorATerms().size() << "\n";
        Feel::cout << "B terms: " << problem.operatorBTerms().size() << "\n";
        Feel::cout << "requested eigenpairs: " << problem.nev() << "\n";
        Feel::cout << "assembled: " << std::boolalpha << system.assembled() << "\n";
        Feel::cout << "converged eigenpairs: " << result.converged << "\n";
        for ( std::size_t i = 0; i != system.eigenvalues().size(); ++i )
            Feel::cout << "lambda[" << i << "]: " << system.eigenvalues()[i] << "\n";
        Feel::cout << result.message << "\n";
    }

    return result.converged >= problem.nev() ? 0 : 1;
}
