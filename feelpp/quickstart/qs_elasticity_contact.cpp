#include "qs_elasticity_contact.hpp"
#include "qs_elasticity_contact_dynamic.hpp"
#include "qs_elasticity_contact_static.hpp"
#include "qs_elasticity_rigid.hpp"

namespace Feel
{
extern template class ContactStatic<2, 1>;
extern template class ContactStatic<3, 1>;
extern template class ContactDynamic<2, 1,1>;
extern template class ContactDynamic<2, 2,2>;
extern template class ContactDynamic<3, 1,1>;
extern template class ContactDynamic<3, 2,2>;
extern template class ElasticRigid<2, 1>;
extern template class ElasticRigid<3, 1>;


void runModel( const nl::json& specs )
{
    int dimension = specs["/Models/LinearElasticity/dimension"_json_pointer];
    int order = specs["/Models/LinearElasticity/order"_json_pointer];
    int orderGeo = specs["/Models/LinearElasticity/orderGeo"_json_pointer];
    int rigidmotion = specs["/Models/LinearElasticity/rigidMotion"_json_pointer];
    bool steady = specs["/TimeStepping/LinearElasticity/steady"_json_pointer];
    
    if (rigidmotion == 0)
    {
        if (dimension == 2)
        {
            if (steady)
            {
                ContactStatic<2, 1> model( specs );
                model.run();
            }
            else 
            {
                if ( orderGeo == 1 )
                {
                    ContactDynamic<2, 1,1> model( specs );
                    model.run();
                }
                else if (orderGeo == 2) 
                {
                    ContactDynamic<2, 2, 2> model(specs);
                    model.run();
                }
            }
        }
        else if (dimension == 3)
        {
            if (steady)
            {
                ContactStatic<3, 1> model( specs );
                model.run();
            }
            else 
            {
                if ( orderGeo == 1 )
                {
                    ContactDynamic<3, 1,1> model( specs );
                    model.run();
                }
                else if (orderGeo == 2) 
                {
                    ContactDynamic<3, 2, 2> model(specs);
                    model.run();
                }
            }
        }
    }
    else if (rigidmotion == 1)
    {
        if (dimension == 2)
        {
            if ( orderGeo == 1 )
            {
                ElasticRigid<2, 1> model( specs );
                model.run();
            }
        }
        else if (dimension == 3)
        {
            if ( orderGeo == 1 )
            {
                ElasticRigid<2, 1> model( specs );
                model.run();
            }
        }
    }
    else
    {
        throw std::runtime_error( fmt::format( "Invalid dimension {} specified in the input file", dimension ) );
    }
}

} // namespace Feel


int main( int argc, char** argv )
{
    using namespace Feel;

    try
    {
        Environment env( _argc = argc, _argv = argv,
                         _desc = makeOptions(),
                         _about = about( _name = fmt::format( "elasticity_contact" ),
                                         _author = "Feel++ Consortium",
                                         _email = "feelpp@cemosis.fr" ) );

        auto jsonfile = removeComments( readFromFile( Environment::expand( soption( "specs" ) ) ) );
        std::istringstream istr( jsonfile );
        json specs = json::parse( istr );

        if ( specs["/Models/LinearElasticity"_json_pointer].empty() )
        {
            throw std::runtime_error( "No LinearElasticity model specified in the input file" );
        }

        runModel( specs );
    }
    catch ( ... )
    {
        handleExceptions();
    }

    return 0;
}