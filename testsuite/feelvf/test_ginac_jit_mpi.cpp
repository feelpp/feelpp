/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#define BOOST_TEST_MODULE ginac jit mpi testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelvf/vf.hpp>

#include <string>
#include <vector>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( ginac_jit_mpi )

BOOST_AUTO_TEST_CASE( persistent_expression_is_built_on_ensemble_worldcomm )
{
    auto& world = Feel::Environment::worldComm();
    int nsplit = ( world.globalSize() > 1 && world.globalSize() % 2 == 0 ) ? 2 : 1;
    auto [color, w, wglob] = Feel::Environment::worldCommPtr()->split( nsplit );

    std::vector<GiNaC::symbol> symbols{ GiNaC::symbol( "x" ), GiNaC::symbol( "y" ), GiNaC::symbol( "z" ) };
    std::string expression = "x+" + std::to_string( color + 1 );

    auto e = Feel::vf::expr<2>( expression, symbols, "", *w );
    e.setParameterValues( { { "x", 3.0 }, { "y", 0.0 }, { "z", 0.0 } } );

    BOOST_CHECK_CLOSE( e.evaluate()( 0, 0 ), 4.0 + color, 1e-12 );
}

BOOST_AUTO_TEST_SUITE_END()
