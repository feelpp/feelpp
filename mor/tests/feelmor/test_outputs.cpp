#define BOOST_TEST_MODULE crbmp testsuite

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelmor/crbmodelproperties.hpp>
#include <feel/feelmodels/modelproperties.hpp>

using namespace Feel;


FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( crbmp_suite )


BOOST_AUTO_TEST_CASE( test_valid_crboutputs )
{
    CRBModelProperties crb_model_props;
    crb_model_props.setupFromFilenameOption();

    auto crb_outputs = crb_model_props.outputs();
    BOOST_CHECK_EQUAL(crb_outputs.size(), 4);

    ModelProperties model_properties;
    model_properties.setupFromFilenameOption();

    auto outputs = model_properties.outputs();
    BOOST_CHECK_EQUAL(outputs.size(), 4);

    for( auto const& [name,crb_output] : crb_outputs )
    {
        Feel::cout << name << std::endl;

        auto output = outputs.find(name);
        std::cout << typeid(output).name() << std::endl;

        if( crb_output.type() == "mean" )
        {
            BOOST_CHECK_EQUAL(crb_output.name(), "myoutput");
            BOOST_CHECK_EQUAL(crb_output.dim(), 3);
            std::set<std::string> m {"marker1","marker2"};
            auto markers = crb_output.markers();
            BOOST_CHECK_EQUAL_COLLECTIONS(markers.begin(), markers.end(), m.begin(), m.end());
        }
        else if( crb_output.type() == "integrate" )
        {
            BOOST_CHECK_EQUAL(crb_output.name(), "flux");
            BOOST_CHECK_EQUAL(crb_output.dim(), 2);
            auto markers = crb_output.markers();
            std::set<std::string> m {"face1"};
            BOOST_CHECK_EQUAL_COLLECTIONS(markers.begin(), markers.end(), m.begin(), m.end());
            BOOST_CHECK(crb_output.hasExpression());
        }
        else if( crb_output.type() == "sensor" )
        {
            BOOST_CHECK_EQUAL(crb_output.name(), "sensor");
            BOOST_CHECK_CLOSE(crb_output.radius(), 0.1, 1e-10);
            auto coord = crb_output.coord();
            std::vector<double> c {0.2,0.1};
            BOOST_CHECK_EQUAL_COLLECTIONS(coord.begin(),coord.end(),c.begin(),c.end());
        }
        else if( crb_output.type() == "point" )
        {
            BOOST_CHECK_EQUAL(crb_output.name(), "point_A");
            BOOST_CHECK(crb_output.hasExpression());
            auto coord = crb_output.coord();
            std::vector<double> c {0.2,0.1};
            BOOST_CHECK_EQUAL_COLLECTIONS(coord.begin(),coord.end(),c.begin(),c.end());
        }
    }

}
BOOST_AUTO_TEST_SUITE_END()
