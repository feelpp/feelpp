#include <cmath>

#include <boost/mpi.hpp>
#include <iostream>
#include <boost/serialization/string.hpp> // Needed to send/receive strings!

#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE boostmpi
#include <testsuite/testsuite.hpp>

//____________________________________________________________________________//

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( boostmpi )

BOOST_AUTO_TEST_CASE( test_case1 )
{
	using namespace Feel;
	BOOST_TEST_MESSAGE( "test_case1" );
	mpi::communicator world;

	if ( world.rank() == 0 )
	{
		world.send( 1, 0, std::string( "Hello" ) );
		std::string msg;
		world.recv( 1, 1, msg );
		BOOST_TEST_MESSAGE( msg << "!" );
	}

	else
	{
		std::string msg;
		world.recv( 0, 0, msg );
		BOOST_TEST_MESSAGE( msg << ", " );
		world.send( 0, 1, std::string( "world" ) );
	}

	world.barrier();
	std::string value;

	if ( world.rank() == 0 )
	{
		value = "Hello, World!";
	}

	broadcast( world, value, 0 );

	BOOST_TEST_MESSAGE( "Process #" << world.rank() << " says " << value  );

	world.barrier();
	double dvalue;

	if ( world.rank() == 1 )
	{
		dvalue = M_PI;
	}

	broadcast( world, dvalue, 1 );

	BOOST_TEST_MESSAGE( "Process #" << world.rank() << " shows " << dvalue  );

	world.barrier();
	double res = 0;
	all_reduce( world, dvalue, res, std::plus<double>() );

	BOOST_TEST_MESSAGE( "Process #" << world.rank() << " shows " << res  );

}


//____________________________________________________________________________//


BOOST_AUTO_TEST_SUITE_END()
