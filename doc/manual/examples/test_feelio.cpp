/* Test MasterStream behavior in different scenarios */

#include <feel/feelcore/feel.hpp>
#include <iostream>

using namespace Feel;

int main(int argc, char** argv)
{
    std::cout << "=== Testing Feel::cout behavior ===" << std::endl;
    
    // Test 1: Before Environment initialization (MPI not initialized)
    std::cout << "Test 1: Before MPI init" << std::endl;
    Feel::cout << "  Feel::cout before Environment (should print)" << std::endl;
    
    // Test 2: With Environment (MPI initialized)
    {
        std::cout << "\nTest 2: With Environment (MPI active)" << std::endl;
        po::options_description opts("Test options");
        Environment env( _argc=argc, _argv=argv, _desc=opts, _about=about(_name="test_feelio") );
        
        Feel::cout << "  Feel::cout with Environment (should print only on master)" << std::endl;
        Feel::cerr << "  Feel::cerr with Environment (should print only on master)" << std::endl;
        
        std::cout << "Test 2 complete, Environment will be destroyed..." << std::endl;
    }
    
    // Test 3: After Environment destruction (MPI finalized)
    std::cout << "\nTest 3: After Environment destruction (MPI finalized)" << std::endl;
    Feel::cout << "  Feel::cout after Environment (should print to avoid crash)" << std::endl;
    
    std::cout << "\n=== All tests completed successfully ===" << std::endl;
    return 0;
}
