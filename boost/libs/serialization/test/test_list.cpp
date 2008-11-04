/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_list.cpp

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

#include <boost/archive/archive_exception.hpp>
#include "test_tools.hpp"
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE(BOOST_ARCHIVE_TEST)

#include <boost/serialization/list.hpp>
#ifdef BOOST_HAS_SLIST
#include <boost/serialization/slist.hpp>
#endif

#include "A.hpp"

int test_main( int /* argc */, char* /* argv */[] )
{
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    std::list<A> alist;
    alist.push_back(A());
    alist.push_back(A());
    {   
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os);
        oa << boost::serialization::make_nvp("alist",alist);
    }

    std::list<A> alist1;
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is);
        ia >> boost::serialization::make_nvp("alist",alist1);
    }
    BOOST_CHECK(alist == alist1);
    
    #ifdef BOOST_HAS_SLIST
    BOOST_STD_EXTENSION_NAMESPACE::slist<A> aslist;
    aslist.push_front(A());
    aslist.push_front(A());
    {   
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os);
        oa << boost::serialization::make_nvp("aslist", aslist);
    }
    BOOST_STD_EXTENSION_NAMESPACE::slist<A> aslist1;{
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is);
        ia >> boost::serialization::make_nvp("aslist", aslist1);
   }
    BOOST_CHECK(aslist == aslist1);
    
    #endif
    std::remove(testfile);
    return EXIT_SUCCESS;
}

// EOF

