#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_inner_curl
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( inner_curl_suite )

BOOST_AUTO_TEST_CASE( test_3d )
{
    auto mesh = unitSphere();

    auto Vh = Pchv<2>(mesh);
    auto f = expr<3,1>("{2*x^2+3*y^2+4*z^2,3*x^2+4*y^2+2*z^2,4*x^2+2*y^2+3*z^2}:x:y:z");
    auto g = expr<3,1>("{4*y-4*z,8*z-8*x,6*x-6*y}:x:y:z");
    auto u = Vh->element(f);
    auto v = Vh->element(g);

    // check inner keyword
    auto a = form2( _test=Vh, _trial=Vh);
    a = integrate( _range=elements( mesh ), _expr=inner(curlt(u),curl(v)));

    auto err = normL2(_range=elements(mesh), _expr=curlv(u)-idv(v));
    BOOST_CHECK_SMALL( err, 1e-8 );
}

BOOST_AUTO_TEST_CASE( test_2d )
{
    auto mesh = unitCircle();

    auto Vh = Pchv<2>(mesh);
    auto Wh = Pch<1>(mesh);
    auto f = expr<2,1>("{2*x^2+3*y^2,3*x^2+2*y^2}:x:y");
    auto g = expr("6*x-6*y:x:y");
    auto u = Vh->element(f);
    auto v = Wh->element(g);

    // check inner keyword
    auto a = form2( _test=Vh, _trial=Vh);
    a = integrate( _range=elements( mesh ), _expr=inner(curlt(u),curl(u)) );

    auto nu = normL2(_range=elements(mesh), _expr=curlv(u));
    auto nv = normL2(_range=elements(mesh), _expr=idv(v));
    BOOST_CHECK_CLOSE( nu, nv, 1e-4 );

    auto err = normL2(_range=elements(mesh), _expr=curlv(u)-idv(v));
    BOOST_CHECK_SMALL( err, 1e-8 );
}

BOOST_AUTO_TEST_SUITE_END()
