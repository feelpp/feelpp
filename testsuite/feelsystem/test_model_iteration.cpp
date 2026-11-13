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

#define BOOST_TEST_MODULE feelsystem_model_iteration
#include <feel/feelcore/testsuite.hpp>

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

Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "FeelSystem model iteration test options" );
    options.add( Feel::backend_options( "metadata-iteration" ) );
    return options;
}

Feel::nl::json
makeModelJson()
{
    return Feel::nl::json::parse(
        R"json(
        {
            "Name": "metadata-iteration",
            "Materials": {
                "Mat1": {
                    "markers": "Mat1",
                    "physics": [ "elliptic", "thermal" ],
                    "k": "2",
                    "f": "1"
                }
            },
            "BoundaryConditions": {
                "temperature": {
                    "Dirichlet": {
                        "temperature-wall": {
                            "markers": "Dirichlet",
                            "expr": "7"
                        }
                    }
                },
                "u": {
                    "Load": {
                        "bulk-load": {
                            "markers": "Mat1",
                            "expr": "6"
                        }
                    },
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
                    },
                    "Fourier": {
                        "ambient": {
                            "markers": "Neumann",
                            "expr1": "4",
                            "expr2": "5"
                        }
                    }
                }
            }
        }
        )json" );
}

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault( "feelpp_test_system_model_iteration" ), makeOptions() )

BOOST_AUTO_TEST_SUITE( feelsystem_model_iteration_suite )

BOOST_AUTO_TEST_CASE( model_descriptors_filter_and_assemble )
{
    using namespace Feel;
    using namespace Feel::vf;

    using mesh_type = Mesh<Simplex<2>>;
    auto mesh = createGMSHMesh( _mesh = new mesh_type,
                                _desc = domain( _name = "feelsystem-model-iteration",
                                                 _shape = "hypercube",
                                                 _dim = 2,
                                                 _h = 0.5,
                                                 _usenames = true,
                                                 _addmidpoint = false ) );

    ModelProperties props;
    props.enableBoundaryConditions2();
    props.setup( makeModelJson() );

    auto bcs = pde::boundaryConditions( _model = props, _mesh = mesh );
    auto uBcs = bcs.field( "u" );
    auto materialDomains = pde::materials( _model = props, _mesh = mesh );

    BOOST_CHECK_EQUAL( bcs.size(), 6 );
    BOOST_CHECK_EQUAL( uBcs.size(), 5 );
    BOOST_CHECK_EQUAL( uBcs.dirichlet().size(), 1 );
    BOOST_CHECK_EQUAL( uBcs.neumann().size(), 1 );
    BOOST_CHECK_EQUAL( uBcs.robin().size(), 1 );
    BOOST_CHECK_EQUAL( uBcs.fourier().size(), 1 );
    BOOST_CHECK_EQUAL( uBcs.load().size(), 1 );

    BOOST_CHECK_EQUAL( materialDomains.size(), 1 );
    BOOST_CHECK_EQUAL( materialDomains.hasProperty( "k" ).size(), 1 );
    BOOST_CHECK_EQUAL( materialDomains.hasProperty( "f" ).size(), 1 );
    BOOST_CHECK_EQUAL( materialDomains.withPhysic( "thermal" ).size(), 1 );
    BOOST_CHECK_EQUAL( materialDomains.withPhysics( { "mechanics", "elliptic" } ).size(), 1 );

    BOOST_CHECK_EQUAL( nelements( uBcs.dirichlet()[0].faces(), true ),
                       nelements( markedfaces( mesh, "Dirichlet" ), true ) );
    BOOST_CHECK_EQUAL( nelements( uBcs.neumann()[0].faces(), true ),
                       nelements( markedfaces( mesh, "Neumann" ), true ) );
    BOOST_CHECK_EQUAL( nelements( materialDomains[0].elements(), true ),
                       nelements( markedelements( mesh, "Mat1" ), true ) );

    auto Xh = Pch<1>( mesh );
    auto uh = Xh->element( "u" );
    auto u = trial( Xh );
    auto v = test( Xh );

    std::size_t countedBoundaryConditions = 0;
    auto problem =
        ellipticProblem( _space = Xh, _name = "metadata-iteration" )
            .forEach( uBcs, [&countedBoundaryConditions]( auto const& )
            {
                ++countedBoundaryConditions;
            } )
            .forEach( materialDomains.hasProperty( "k" ), [&]( auto& p, auto const& material )
            {
                auto k = material.property( "k" ).template expr<1, 1>();
                p.diffusion( integrate( _range = material.elements(),
                                        _expr = inner( k * grad( u ), grad( v ) ) ) );
            } )
            .forEach( materialDomains.hasProperty( "f" ), [&]( auto& p, auto const& material )
            {
                auto f = material.property( "f" ).template expr<1, 1>();
                p.source( integrate( _range = material.elements(),
                                     _expr = f * v ) );
            } )
            .forEach( uBcs.load(), [&]( auto& p, auto const& bc )
            {
                p.source( integrate( _range = bc.elements(),
                                     _expr = bc.template expr<1, 1>() * v ) );
            } )
            .forEach( uBcs.neumann(), [&]( auto& p, auto const& bc )
            {
                p.neumann( _range = bc.faces(),
                           _test = v,
                           _expr = bc.template expr<1, 1>() );
            } )
            .forEach( uBcs.robin(), [&]( auto& p, auto const& bc )
            {
                p.robin( _range = bc.faces(),
                         _trial = u,
                         _test = v,
                         _alpha = bc.template expr1<1, 1>(),
                         _expr = bc.template expr2<1, 1>() );
            } )
            .forEach( uBcs.fourier(), [&]( auto& p, auto const& bc )
            {
                p.fourier( _range = bc.faces(),
                           _trial = u,
                           _test = v,
                           _alpha = bc.template expr1<1, 1>(),
                           _expr = bc.template expr2<1, 1>() );
            } )
            .forEach( uBcs.dirichlet(), [&]( auto& p, auto const& bc )
            {
                p.dirichlet( _range = bc.faces(),
                             _element = uh,
                             _expr = bc.template expr<1, 1>() );
            } );

    BOOST_CHECK_EQUAL( countedBoundaryConditions, uBcs.size() );
    BOOST_CHECK_EQUAL( problem.lhsTerms().size(), 2 );
    BOOST_CHECK_EQUAL( problem.stiffnessTerms().size(), 1 );
    BOOST_CHECK_EQUAL( problem.rhsTerms().size(), 5 );
    BOOST_CHECK_EQUAL( problem.constraintTerms().size(), 1 );

    auto system = problem.assemble();
    BOOST_CHECK( system.assembled() );
    BOOST_CHECK( system.constrained() );
    BOOST_CHECK_EQUAL( system.backendName(), problem.name() );

    auto result = problem.solve( _solution = uh );
    BOOST_CHECK( result.converged );
}

BOOST_AUTO_TEST_SUITE_END()
