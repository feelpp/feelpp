#pragma once

#include <feel/feelcore/testsuite.hpp>
#include "cavity_radiation_jacobian.hpp"

namespace Feel::Test {

/// @brief Create command-line options for radiative heat transfer cavity tests
inline po::options_description makeCavityOptions() 
{
    po::options_description options("rht options");
    options.add_options()
        ("specs", po::value<std::string>(), "json spec file for rht")
        ("steady", po::value<bool>()->default_value(1), "if 1: steady else unsteady")
        ("deactivate-exporters", po::value<bool>()->default_value(false), "deactivate exporters");
    options.add(backend_options("heatEq"));
    return options.add(feel_options());
}

/// @brief Create AboutData for cavity radiation tests
/// @param test_name Name of the test executable
/// @param description Brief description of the test case
inline AboutData makeCavityAbout(const std::string& test_name, const std::string& description) 
{
    return AboutData(test_name, test_name, "0.1", description, 
                     AboutData::License_GPL, 
                     "Copyright (c) 2023 Feel++ Consortium");
}

/// @brief Execute radiative heat transfer cavity test
/// @tparam Dim Spatial dimension
/// @tparam Order Polynomial order
template<int Dim, int Order>
void runCavityRadiationTest() 
{
    // Read the json file associated to the heat transfer problem
    auto jsonfile = removeComments(readFromFile(Environment::expand(soption("specs"))));
    std::istringstream istr(jsonfile);
    json specs = json::parse(istr);
    
    // Instantiate the class for the solution of the heat transfer problem
    RHT<Dim, Order> rht(specs);
    
    // Solve the heat transfer problem
    rht.executeNonLinear();
    
    // Checks are inside this function using CHECK()
    rht.checkResults();
}

} // namespace Feel::Test
