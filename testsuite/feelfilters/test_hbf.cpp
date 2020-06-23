#define BOOST_TEST_MODULE hbf

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/hbf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_hbf )

typedef boost::mpl::list<int8_t, int16_t,int32_t,int64_t, uint8_t, uint16_t,uint32_t,uint64_t,float,double> types;
using namespace Feel;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_hbf, T, types )
{
    holo3_image<T> x (2,2);
    x=holo3_image<T>::Ones(2,2);
    writeHBF("test.hbf", x);
    holo3_image<T> y = readHBF<T>("test.hbf");
    BOOST_CHECK_EQUAL(x,y);
}
BOOST_AUTO_TEST_SUITE_END()