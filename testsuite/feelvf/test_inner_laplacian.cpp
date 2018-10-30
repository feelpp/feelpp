#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_inner_laplacian
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( inner_laplacian_suite )

BOOST_AUTO_TEST_CASE( test_0 )
{
    auto mesh = unitSphere();

    auto Vh = Pchv<1>(mesh);
    auto u = Vh->element();
    auto v = Vh->element();

    auto a = form2( _test=Vh, _trial=Vh);
    a = integrate( _range=elements( mesh ), _expr=inner(laplaciant(u),laplacian(v)));

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
