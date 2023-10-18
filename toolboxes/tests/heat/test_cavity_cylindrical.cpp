#define BOOST_TEST_MODULE heat testsuite
#include <feel/feelcore/testsuite.hpp>
#include "cavity_radiation_jacobian.hpp"

using namespace Feel;
inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "rht options" );
    options.add_options()
        ( "specs", Feel::po::value<std::string>(),"json spec file for rht" )
        ( "steady", Feel::po::value<bool>()->default_value( 1 ),"if 1: steady else unsteady" )
        ("deactivate-exporters",Feel::po::value<bool>()->default_value( false ),"deactivate exporters");
    options.add( Feel::backend_options("heatEq"));

    return options.add( Feel::feel_options() );

} 

inline
AboutData
makeAbout()
{
    AboutData about( "test_cavity_cylindrical" ,
                     "test_cavity_cylindrical" ,
                     "0.1",
                     "Cylindrical cavity test",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2023 Feel++ Consortium" );

    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );


BOOST_AUTO_TEST_SUITE( heatsuite )

BOOST_AUTO_TEST_CASE( test_cylindrical_cavity_radiation )
{
    using namespace Feel;

    // Read the json file associated to the heat transfer problem
    auto jsonfile = removeComments( readFromFile( Environment::expand( soption( "specs" ) ) ) );
    std::istringstream istr( jsonfile );
    json specs = json::parse( istr );    

    // Instantiate the class for the solution of the heat transfer problem
    RHT<FEELPP_DIM, FEELPP_ORDER> rht(specs);

    // Solve the heat transfer problem
    rht.executeNonLinear();

    // Checks are inside this function using CHECK()
    rht.checkResults();

}

BOOST_AUTO_TEST_SUITE_END()