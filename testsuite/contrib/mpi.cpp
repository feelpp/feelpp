// Taken from Boost Tutorial
#include <boost/mpi.hpp>
#include <iostream>
#include <string>
#include <boost/serialization/string.hpp>
namespace mpi = boost::mpi;

int main()
{
    mpi::environment env;
    mpi::communicator world;

    std::string value;
    if (world.rank() == 0) {
        value = "Hello, World!";
    }

    broadcast(world, value, 0);

    std::cout << "Process #" << world.rank() << " says " << value
              << std::endl;
    return 0;
}
