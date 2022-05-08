#define BOOST_TEST_MODULE crbmp testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelmor/crbmodelproperties.hpp>

using namespace Feel;


FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( crbmp_suite )

BOOST_AUTO_TEST_CASE( test_crbparameters )
{
    CRBModelProperties model_props;
    model_props.setupFromFilenameOption();

    auto params = model_props.parameters();
    BOOST_CHECK_EQUAL(params.size(), 2);
    auto um = params["Um"];
    BOOST_CHECK_CLOSE(um.min(), 1e-4, 1e-10);
    BOOST_CHECK_CLOSE(um.max(), 10, 1e-10);
    BOOST_CHECK_EQUAL(um.description(), "Um desc");
    BOOST_CHECK_EQUAL(um.sampling(), "log-random");
    BOOST_CHECK_EQUAL(um.samplingSize(), 5);
    auto H = params["H"];
    BOOST_CHECK_CLOSE(H.min(), -2, 1e-10);
    BOOST_CHECK_CLOSE(H.max(), 3, 1e-10);
    BOOST_CHECK_CLOSE(H.value(), 0.3, 1e-10);
    
}

BOOST_AUTO_TEST_CASE( test_crboutputs )
{
    CRBModelProperties model_props;
    model_props.setupFromFilenameOption();

    auto outputs = model_props.outputs();
    BOOST_CHECK_EQUAL(outputs.size(), 4);

    for( auto const& [name,output] : outputs )
    {
        if( output.type() == "mean" )
        {
            BOOST_CHECK_EQUAL(output.name(), "myoutput");
            BOOST_CHECK_EQUAL(output.dim(), 3);
            std::set<std::string> m {"marker1","marker2"};
            auto markers = output.markers();
            BOOST_CHECK_EQUAL_COLLECTIONS(markers.begin(), markers.end(), m.begin(), m.end());
        }
        else if( output.type() == "integrate" )
        {
            BOOST_CHECK_EQUAL(output.name(), "flux");
            BOOST_CHECK_EQUAL(output.dim(), 2);
            auto markers = output.markers();
            std::set<std::string> m {"face1"};
            BOOST_CHECK_EQUAL_COLLECTIONS(markers.begin(), markers.end(), m.begin(), m.end());
            BOOST_CHECK(output.hasExpression());
        }
        else if( output.type() == "sensor" )
        {
            BOOST_CHECK_EQUAL(output.name(), "sensor");
            BOOST_CHECK_CLOSE(output.radius(), 0.1, 1e-10);
            auto coord = output.coord();
            std::vector<double> c {0.2,0.1};
            BOOST_CHECK_EQUAL_COLLECTIONS(coord.begin(),coord.end(),c.begin(),c.end());
        }
        else if( output.type() == "point" )
        {
            BOOST_CHECK_EQUAL(output.name(), "point_A");
            BOOST_CHECK(output.hasExpression());
            auto coord = output.coord();
            std::vector<double> c {0.2,0.1};
            BOOST_CHECK_EQUAL_COLLECTIONS(coord.begin(),coord.end(),c.begin(),c.end());
        }
    }

}

BOOST_AUTO_TEST_SUITE_END()
