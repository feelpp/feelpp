#include "qs_elasticity_contact_static.hpp"
#include "qs_elasticity_contact_dynamic.hpp"
#include "qs_contact_lagrange.hpp"



int main(int argc, char** argv)
{
    using namespace Feel;
    int status;
    
    try
    {
        Environment env(_argc = argc, _argv = argv,
                        _desc = makeOptions(),
                        _about = about(_name = fmt::format("elasticity_contact"),
                                       _author = "Feel++ Consortium",
                                       _email = "feelpp@cemosis.fr"));
        
        auto jsonfile = removeComments(readFromFile(Environment::expand(soption("specs"))));
        std::istringstream istr(jsonfile);
        json specs = json::parse(istr);

        if ( specs["/Models/LinearElasticity"_json_pointer].empty() )
            throw std::runtime_error("No LinearElasticity model specified in the input file");

        int dimension = specs["/Models/LinearElasticity/dimension"_json_pointer];
        int order = specs["/Models/LinearElasticity/order"_json_pointer];
        int orderGeo = specs["/Models/LinearElasticity/orderGeo"_json_pointer];
        bool sta = specs["/TimeStepping/LinearElasticity/steady"_json_pointer];
        std::string method = specs["/Collision/LinearElasticity/method"_json_pointer];
        
        if ( dimension == 2 && order == 1)
        {
            if (sta)
            {
                ContactStatic<2, 1> ContactStatic(specs);
                ContactStatic.run();
            }
            else 
            {
                if (orderGeo == 1)
                {
                    ContactDynamic<2,1,1> ContactDynamic(specs);
                    ContactDynamic.run();
                }
                else if (orderGeo == 2)
                {
                    ContactDynamic<2,1,2> ContactDynamic(specs);
                    ContactDynamic.run();
                }
            }
        }
        else if ( dimension == 2 && order == 2)
        {
            if (sta)
            {
                ContactStatic<2, 2> ContactStatic(specs);
                ContactStatic.run();
            }
            else 
            {
                if (orderGeo == 1)
                {
                    if (method.compare("lagrange") == 0)
                    {
                        ContactDynamicLagrange<2,2,1> ContactDynamicLagrange(specs);
                        ContactDynamicLagrange.run();             
                    }
                    else 
                    {
                        ContactDynamic<2,2,1> ContactDynamic(specs);
                        ContactDynamic.run();
                    }
                }
                else if (orderGeo == 2)
                {
                    ContactDynamic<2,2,2> ContactDynamic(specs);
                    ContactDynamic.run();
                }
            }
        }
        else if ( dimension == 3 && order == 1)
        {
            if (sta)
            {
                ContactStatic<3, 1> ContactStatic(specs);
                ContactStatic.run();
            }
            else 
            {
                ContactDynamic<3,1,1> ContactDynamic(specs);
                ContactDynamic.run();
            }
        }
        else if ( dimension == 3 && order == 2)
        {
            if (sta)
            {
                ContactStatic<3, 2> ContactStatic(specs);
                ContactStatic.run();
            }
            else 
            {
                ContactDynamic<3,2,1> ContactDynamic(specs);
                ContactDynamic.run();
            }
        }
        else
            throw std::runtime_error(fmt::format("Invalid dimension {} specified in the input file", dimension));
    }
    catch (...)
    {
        handleExceptions();
    }
    return 0;
}
