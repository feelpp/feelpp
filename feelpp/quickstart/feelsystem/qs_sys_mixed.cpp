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

#include <boost/hana/integral_constant.hpp>

#include "../discretization.hpp"
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

int
main( int argc, char** argv )
{
    using namespace Feel;
    using namespace Feel::vf;
    using namespace boost::hana::literals;

    po::options_description options( "FeelSystem mixed quickstart options" );
    options.add( backend_options( "mixed-stokes" ) );
    options.add_options()
        ( "order", po::value<int>()->default_value( 2 ), "velocity order k >= 2 (pressure uses k-1)" )
        ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" );

    Environment env( _argc = argc, _argv = argv,
                     _desc = options,
                     _about = about( _name = "qs_sys_mixed",
                                     _author = "University of Strasbourg",
                                     _email = "feelpp-devel@feelpp.org" ) );

    auto const velocityOrder = quickstart::checkedOrder( "order", 2 );
    auto const pressureOrder = RuntimeOrder::checked( velocityOrder.value - 1 );
    auto mesh = loadMesh( _mesh = new Mesh<Simplex<FEELPP_DIM, Dynamic>>( quickstart::geometryOrder() ) );
    quickstart::checkMeshOrder( mesh );
    auto Vh = Pchv<Dynamic>( mesh, velocityOrder );
    auto Qh = Pch<Dynamic>( mesh, pressureOrder );
    auto Xh = product( Vh, Qh );
    auto solution = Xh.element();
    auto boundaryVelocity = Vh->element( "u_dirichlet" );

    auto u = trial( Vh );
    auto v = test( Vh );
    auto p = trial( Qh );
    auto q = test( Qh );

    auto uExact = vec( Px() * Px(), -2.0 * Px() * Py() );
    auto fExact = vec( cst( -2.0 ), cst( 0.0 ) );
    double const pressurePenalty = 1e-6;

    auto problem =
        mixedProblem( _space = Xh, _name = "mixed-stokes" )
            .block( 0_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = inner( grad( u ), grad( v ) ) ) )
            .block( 0_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = -p * div( v ) ) )
            .block( 1_c, 0_c,
                    integrate( _range = elements( mesh ),
                               _expr = -div( u ) * q ) )
            .block( 1_c, 1_c,
                    integrate( _range = elements( mesh ),
                               _expr = pressurePenalty * p * q ) )
            .row( 0_c,
                  integrate( _range = elements( mesh ),
                             _expr = inner( fExact, v ) ) )
            .constraint( [mesh, boundaryVelocity, uExact]( auto& a, auto& l )
            {
                a.row( 0_c ) += on( _range = boundaryfaces( mesh ),
                                    _rhs = l( 0_c ),
                                    _element = boundaryVelocity,
                                    _expr = uExact,
                                    _type = "elimination" );
            } );

    auto system = problem.assemble();
    auto result = system.solve( _solution = solution );

    double const velocityError = normL2( _range = elements( mesh ),
                                         _expr = solution( 0_c ) - uExact );
    double const pressureNorm = normL2( _range = elements( mesh ),
                                        _expr = solution( 1_c ) );

    if ( Environment::isMasterRank() )
    {
        Feel::cout << "Problem: " << problem.name() << "\n";
        Feel::cout << "lhs terms: " << problem.lhsTerms().size() << "\n";
        Feel::cout << "rhs terms: " << problem.rhsTerms().size() << "\n";
        Feel::cout << "constraint terms: " << problem.constraintTerms().size() << "\n";
        Feel::cout << "backend name: " << system.backendName() << "\n";
        Feel::cout << "converged: " << std::boolalpha << result.converged << "\n";
        Feel::cout << "iterations: " << result.iterations << "\n";
        Feel::cout << "residual: " << result.residual << "\n";
        Feel::cout << "velocity L2 error: " << velocityError << "\n";
        Feel::cout << "pressure L2 norm: " << pressureNorm << "\n";
    }

    return result.converged && velocityError < 1e-8 && pressureNorm < 1e-8 ? 0 : 1;
}
