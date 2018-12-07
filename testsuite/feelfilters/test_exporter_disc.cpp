#define BOOST_TEST_MODULE test_exporter_disc
#include <feel/feelcore/testsuite.hpp>


#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pdh.hpp>

/** use Feel namespace */
using namespace Feel;


FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( exporter_disc )


BOOST_AUTO_TEST_CASE( test_1 )
{
    auto mesh = unitSquare();
    auto Xh = Pdhv<1>(mesh);
    auto Vh = Pdh<1>(mesh);
    auto u = Xh->element();
    auto v = Vh->element();
    auto e = exporter(_mesh=mesh);
    e->step(0)->add( "u", u );
    e->step(0)->add( "v", v );
    e->save();
    e->step(1)->add( "u", u );
    e->step(1)->add( "v", v );
    e->save();
}

BOOST_AUTO_TEST_SUITE_END()
