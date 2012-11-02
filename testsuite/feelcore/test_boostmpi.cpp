#include <cmath>

#include <boost/mpi.hpp>
#include <iostream>
#include <boost/serialization/string.hpp> // Needed to send/receive strings!
namespace mpi = boost::mpi;

#define BOOST_TEST_MODULE example
#include <boost/test/included/unit_test.hpp>

struct Feelpp {
    Feelpp()
        {
            BOOST_TEST_MESSAGE( "setup Feel++" );
            Feel::Environment env( boost::unit_test::framework::master_test_suite().argc,
                                   boost::unit_test::framework::master_test_suite().argv);
        }
    ~Feelpp()
        {
            BOOST_TEST_MESSAGE( "teardown Feel++" );
        }

};

//____________________________________________________________________________//

BOOST_FIXTURE_TEST_SUITE( s, Feelpp )

BOOST_AUTO_TEST_CASE( test_case1 )
{
    mpi::communicator world;

    if ( world.rank() == 0 )
    {
        world.send( 1, 0, std::string( "Hello" ) );
        std::string msg;
        world.recv( 1, 1, msg );
        std::cout << msg << "!" << std::endl;
    }

    else
    {
        std::string msg;
        world.recv( 0, 0, msg );
        std::cout << msg << ", ";
        std::cout.flush();
        world.send( 0, 1, std::string( "world" ) );
    }

    world.barrier();
    std::string value;

    if ( world.rank() == 0 )
    {
        value = "Hello, World!";
    }

    broadcast( world, value, 0 );

    std::cout << "Process #" << world.rank() << " says " << value << std::endl;

    world.barrier();
    double dvalue;

    if ( world.rank() == 1 )
    {
        dvalue = M_PI;
    }

    broadcast( world, dvalue, 1 );

    std::cout << "Process #" << world.rank() << " shows " << dvalue << std::endl;

    world.barrier();
    double res = 0;
    all_reduce( world, dvalue, res, std::plus<double>() );

    std::cout << "Process #" << world.rank() << " shows " << res << std::endl;

}


//____________________________________________________________________________//


BOOST_AUTO_TEST_SUITE_END()
