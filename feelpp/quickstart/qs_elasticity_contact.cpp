#include "qs_elasticity_contact.hpp"

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

        if ( dimension == 2 && order == 1)
        {
            ElasticContact<2, 1> ElasticContact(specs);
            ElasticContact.run();
        }
        else if ( dimension == 2 && order == 2)
        {
            ElasticContact<2, 2> ElasticContact(specs);
            ElasticContact.run();
        }
        else if ( dimension == 3 && order == 1)
        {
            ElasticContact<3, 1> ElasticContact(specs);
            ElasticContact.run();
        }
        else if ( dimension == 3 && order == 2)
        {
            ElasticContact<3, 2> ElasticContact(specs);
            ElasticContact.run();
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
