/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_archive_derivation.cpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// should pass compilation and execution

#include <fstream>

#include <cstdio> // remove
#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::remove;
}
#endif

#include "test_tools.hpp"
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE(BOOST_ARCHIVE_TEST)

#include "D.hpp"

int test_main( int /* argc */, char* /* argv */[] )
{
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    const A a;
    A a1;

    {   
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os);
        oa << boost::serialization::make_nvp("a", a);
    }
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is);
        ia >> boost::serialization::make_nvp("a", a1);
    }
    BOOST_CHECK(a == a1);
    std::remove(testfile);
    return EXIT_SUCCESS;
}

// EOF
