#include "qs_active_elasticity.hpp"
#include "qs_active_elasticity_contact.hpp"


namespace Feel
{
extern template class ActiveContact<2, 1>;

void runModel( const nl::json& specs )
{
    int dimension = specs["/Models/HyperElasticity/dimension"_json_pointer];
    int order = specs["/Models/HyperElasticity/order"_json_pointer];
    int rigidmotion = specs["/Models/HyperElasticity/rigidMotion"_json_pointer];
    
    ActiveContact<2,1> model(specs);
    model.run();
}

} // namespace Feel


int main( int argc, char** argv )
{
    using namespace Feel;

    try
    {
        Environment env( _argc = argc, _argv = argv,
                         _desc = makeOptions(),
                         _about = about( _name = fmt::format( "active_hyper_elasticity_contact" ),
                                         _author = "Feel++ Consortium",
                                         _email = "feelpp@cemosis.fr" ) );

        auto jsonfile = removeComments( readFromFile( Environment::expand( soption( "specs" ) ) ) );
        std::istringstream istr( jsonfile );
        json specs = json::parse( istr );

        if ( specs["/Models/HyperElasticity"_json_pointer].empty() )
        {
            throw std::runtime_error( "No HyperElasticity model specified in the input file" );
        }

        runModel( specs );
    }
    catch ( ... )
    {
        handleExceptions();
    }

    return 0;
}