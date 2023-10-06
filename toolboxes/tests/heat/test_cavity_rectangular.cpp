#include "cavity_radiation_jacobian.hpp"

int main( int argc, char** argv )
{
    using namespace Feel;
    try
    {
        Environment env( _argc = argc, _argv = argv,
                         _desc = makeOptions(),
                         _about = about( _name = fmt::format( "rht-{}dp{}", FEELPP_DIM, FEELPP_ORDER ),
                                         _author = "Feel++ Consortium",
                                         _email = "feelpp@cemosis.fr" ) );

        // Read the json file associated to the heat transfer problem
        auto jsonfile = removeComments( readFromFile( Environment::expand( soption( "specs" ) ) ) );
        std::istringstream istr( jsonfile );
        json specs = json::parse( istr );    

        // Instantiate the class for the solution of the heat transfer problem
        RHT<FEELPP_DIM, FEELPP_ORDER> rht(specs);

        // Solve the heat transfer problem
        rht.executeNonLinear();
        rht.checkResults();
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
} // end main()