#define BOOST_TEST_MODULE heat testsuite
#include "test_cavity_common.hpp"

using namespace Feel;

FEELPP_ENVIRONMENT_WITH_OPTIONS(
    Test::makeCavityAbout("test_cavity_triangular", "Triangular cavity test"),
    Test::makeCavityOptions()
);

BOOST_AUTO_TEST_SUITE( heatsuite )

BOOST_AUTO_TEST_CASE( test_triangular_cavity_radiation )
{
    Test::runCavityRadiationTest<FEELPP_DIM, FEELPP_ORDER>();
}

BOOST_AUTO_TEST_SUITE_END()