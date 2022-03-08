#define BOOST_TEST_MODULE hbf

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/hbf.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_hbf )

typedef boost::mpl::list<int8_t, int16_t,int32_t,int64_t, uint8_t, uint16_t,uint32_t,uint64_t,float,double> types;
using namespace Feel;

BOOST_AUTO_TEST_CASE_TEMPLATE( test_hbf, T, types )
{
    for(int n = 2; n <= 6; n +=2)
    {
        holo3_image<T> x (n,n+1);
        x=holo3_image<T>::Random(n,n+1);

        std::string oname = (boost::format("test_%1%.hbf")%n).str();
        writeHBF(oname, x);
        Environment::worldComm().barrier();
        holo3_image<T> y = readHBF<T>(oname);
        Environment::worldComm().barrier();
        BOOST_CHECK_EQUAL(x,y);
    }
}
BOOST_AUTO_TEST_SUITE_END()
