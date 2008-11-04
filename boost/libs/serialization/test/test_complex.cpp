/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_complex.cpp

// (C) Copyright 2005 Matthias Troyer . 
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

#include <boost/serialization/complex.hpp>

int test_main( int /* argc */, char* /* argv */[] )
{
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    // test array of objects
    std::complex<float> a(std::rand(),std::rand());
    std::complex<double> b(std::rand(),std::rand());
    {   
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os);
        oa << boost::serialization::make_nvp("afloatcomplex", a);
        oa << boost::serialization::make_nvp("adoublecomplex", b);
    }
    std::complex<float> a1;
    std::complex<double> b1;
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is);
        ia >> boost::serialization::make_nvp("afloatcomplex", a1);
        ia >> boost::serialization::make_nvp("adoublecomplex", b1);
    }
    bool equal = (std::abs(a-a1) <= 2.*std::numeric_limits<float>::round_error()  
          && std::abs(b-b1) <= 2.*std::numeric_limits<double>::round_error() );
                  
    BOOST_CHECK(equal);
    std::remove(testfile);
    return EXIT_SUCCESS;
}

// EOF
