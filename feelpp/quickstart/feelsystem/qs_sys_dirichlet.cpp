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
#include <feel/feelvf/dirichletconstraints.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

int
main( int argc, char** argv )
{
    using namespace Feel;
    using namespace Feel::vf;

    po::options_description options( "FeelSystem Dirichlet quickstart options" );
    options.add( backend_options( "dirichlet-reaction-diffusion" ) );
    options.add_options()
        ( "order", po::value<int>()->default_value( 1 ), "approximation order" )
        ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" );

    Environment env( _argc = argc, _argv = argv,
                     _desc = options,
                     _about = about( _name = "qs_sys_dirichlet",
                                     _author = "University of Strasbourg",
                                     _email = "feelpp-devel@feelpp.org" ) );

    auto mesh = loadMesh( _mesh = new Mesh<Simplex<FEELPP_DIM, Dynamic>>( quickstart::geometryOrder() ) );
    quickstart::checkMeshOrder( mesh );
    auto Xh = Pch<Dynamic>( mesh, quickstart::checkedOrder( "order" ) );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );
    auto g = cst( 1.0 );

    auto problem =
        ellipticProblem( _space = Xh, _name = "dirichlet-reaction-diffusion" )
            .lhs( integrate( _range = elements( mesh ),
                             _expr = inner( grad( u ), grad( v ) ) + u * v ) )
            .rhs( integrate( _range = elements( mesh ),
                             _expr = v ) )
            .dirichlet( _range = boundaryfaces( mesh ),
                        _element = uh,
                        _expr = g );

    auto system = problem.assemble();
    auto result = system.solve( _solution = uh );

    auto boundaryError = Xh->element( "boundary-error" );
    boundaryError.on( _range = boundaryfaces( mesh ), _expr = abs( uh - g ) );
    double const boundaryMaxError = boundaryError.max();

    if ( Environment::isMasterRank() )
    {
        Feel::cout << "Problem: " << problem.name() << "\n";
        Feel::cout << "lhs terms: " << problem.lhsTerms().size() << "\n";
        Feel::cout << "rhs terms: " << problem.rhsTerms().size() << "\n";
        Feel::cout << "constraint terms: " << problem.constraintTerms().size() << "\n";
        Feel::cout << "backend name: " << system.backendName() << "\n";
        Feel::cout << "constrained active system: " << std::boolalpha << system.constrained() << "\n";
        Feel::cout << "converged: " << std::boolalpha << result.converged << "\n";
        Feel::cout << "iterations: " << result.iterations << "\n";
        Feel::cout << "residual: " << result.residual << "\n";
        Feel::cout << "boundary max error: " << boundaryMaxError << "\n";
    }

    return result.converged ? 0 : 1;
}
