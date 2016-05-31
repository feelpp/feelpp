#define BOOST_TEST_MODULE test_cross_disc
#include <testsuite/testsuite.hpp>


#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;


FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( cross_disc )


BOOST_AUTO_TEST_CASE( test_1 )
{
    auto mesh = unitSphere();
    auto Xh = Pdhv<1>(mesh);
    auto u = Xh->element();
    auto f = form1(_test=Xh);
    f = integrate(internalfaces(mesh), cross(leftface(id(u)),leftface(N()) ) );
    // f = integrate(internalfaces(mesh), leftface(trans(id(u))*N()));
}

BOOST_AUTO_TEST_SUITE_END()
