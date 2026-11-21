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
#include <feel/feelpde/modeldescriptors.hpp>
#include <feel/feelsystem/system.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/options.hpp>

namespace
{

Feel::nl::json
modelJson()
{
    return Feel::nl::json::parse(
        R"json(
        {
            "Name": "qs_sys_model",
            "Materials": {
                "Mat1": {
                    "markers": "Mat1",
                    "physics": "elliptic",
                    "k": "2",
                    "f": "1"
                }
            },
            "BoundaryConditions": {
                "u": {
                    "Dirichlet": {
                        "wall": {
                            "markers": "Dirichlet",
                            "expr": "0"
                        }
                    },
                    "Neumann": {
                        "flux": {
                            "markers": "Neumann",
                            "expr": "1"
                        }
                    },
                    "Robin": {
                        "exchange": {
                            "markers": "Neumann",
                            "expr1": "2",
                            "expr2": "3"
                        }
                    }
                }
            }
        }
        )json" );
}

}

int
main( int argc, char** argv )
{
    using namespace Feel;
    using namespace Feel::vf;

    po::options_description options( "FeelSystem model quickstart options" );
    options.add( backend_options( "model-iteration" ) );
    options.add_options()
        ( "order", po::value<int>()->default_value( 1 ), "approximation order" )
        ( "geometry.order", po::value<int>()->default_value( 1 ), "mesh geometric transformation order" );

    Environment env( _argc = argc, _argv = argv,
                     _desc = options,
                     _about = about( _name = "qs_sys_model",
                                     _author = "University of Strasbourg",
                                     _email = "feelpp-devel@feelpp.org" ) );

    using mesh_type = Mesh<Simplex<FEELPP_DIM, Dynamic>>;
    auto mesh = createGMSHMesh( _mesh = new mesh_type( quickstart::geometryOrder() ),
                                _desc = domain( _name = "qs-sys-model-domain",
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

    ModelProperties props;
    props.enableBoundaryConditions2();
    props.setup( modelJson() );

    auto bcs = pde::boundaryConditions( _model = props, _name = "u", _mesh = mesh );
    auto domains = pde::materials( _model = props, _mesh = mesh );

    auto problem =
        ellipticProblem( _space = Xh, _name = "model-iteration" )
            .forEach( domains.hasProperty( "k" ), [&]( auto& p, auto const& mat )
            {
                auto k = mat.property( "k" ).template expr<1, 1>();
                p.diffusion( integrate( _range = mat.elements(),
                                        _expr = inner( k * grad( u ), grad( v ) ) ) );
            } )
            .forEach( domains.hasProperty( "f" ), [&]( auto& p, auto const& mat )
            {
                auto f = mat.property( "f" ).template expr<1, 1>();
                p.source( integrate( _range = mat.elements(),
                                     _expr = f * v ) );
            } )
            .forEach( bcs.neumann(), [&]( auto& p, auto const& bc )
            {
                p.neumann( _range = bc.faces(),
                           _test = v,
                           _expr = bc.template expr<1, 1>() );
            } )
            .forEach( bcs.robin(), [&]( auto& p, auto const& bc )
            {
                p.robin( _range = bc.faces(),
                         _trial = u,
                         _test = v,
                         _alpha = bc.template expr1<1, 1>(),
                         _expr = bc.template expr2<1, 1>() );
            } )
            .forEach( bcs.dirichlet(), [&]( auto& p, auto const& bc )
            {
                p.dirichlet( _range = bc.faces(),
                             _element = uh,
                             _expr = bc.template expr<1, 1>() );
            } );

    auto system = problem.assemble();
    auto result = system.solve( _solution = uh );

    if ( Environment::isMasterRank() )
    {
        Feel::cout << "Problem: " << problem.name() << "\n";
        Feel::cout << "materials: " << domains.size() << "\n";
        Feel::cout << "boundary conditions: " << bcs.size() << "\n";
        Feel::cout << "stiffness terms: " << problem.stiffnessTerms().size() << "\n";
        Feel::cout << "lhs terms: " << problem.lhsTerms().size() << "\n";
        Feel::cout << "rhs terms: " << problem.rhsTerms().size() << "\n";
        Feel::cout << "constraints: " << problem.constraintTerms().size() << "\n";
        Feel::cout << "converged: " << std::boolalpha << result.converged << "\n";
    }

    return result.converged ? 0 : 1;
}
