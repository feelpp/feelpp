#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_pow
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

using namespace Feel;
using namespace Feel::vf;

inline
AboutData
makeAbout()
{
    AboutData about( "test_pow" ,
                     "test_pow" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( pow_suite )

    BOOST_AUTO_TEST_CASE( test_0 )
{
    auto mesh = unitSphere();

    auto exprB = expr(soption("functions.g"));
    auto exprE = expr(soption("functions.h"));
    auto expo = doption("parameters.alpha");

    auto Vh = Pch<1>(mesh);
    auto u = Vh->element(exprB);
    auto v = Vh->element(exprE);

    u = project(_range=elements(mesh), _space=Vh, _expr=exprB);
    v = project(_range=elements(mesh), _space=Vh, _expr=exprE);

    auto err = normL2( _range=elements(mesh), _expr=pow(idv(u), expo) - idv(v));
    BOOST_CHECK_SMALL( err, 1e-1 );
}

BOOST_AUTO_TEST_SUITE_END()
