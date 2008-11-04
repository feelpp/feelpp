#ifndef BOOST_SERIALIZATION_TEST_TOOLS_HPP
#define BOOST_SERIALIZATION_TEST_TOOLS_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_tools.hpp
//
// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/config.hpp>
#include <cstdio> // remove, tmpnam

// win32 has a brain-dead tmpnam implementation.
// which leaves temp files in root directory 
// regardless of environmental settings
#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)

#include <cstdlib>
#include <cstring>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::remove;
    using ::strcpy;
    using ::strcat;
    using ::tmpnam;
}
#endif // defined(BOOST_NO_STDC_NAMESPACE)

#include <direct.h>
#include <boost/archive/tmpdir.hpp>

#if defined(__COMO__)
    #define chdir _chdir
#endif

#if defined(NDEBUG) && defined(__BORLANDC__)
    #define STRCPY strcpy
#else
    #define STRCPY std::strcpy
#endif

namespace boost {
namespace archive {
    char * tmpnam(char * buffer){
        char old_dir[256];
        _getcwd(old_dir, sizeof(old_dir) - 1);

        char * temp_dir = boost::archive::tmpdir();
        chdir(temp_dir);

        char temp_name[256];
        std::tmpnam(temp_name);

        chdir(old_dir);
        static char ibuffer [512];

        if(NULL == buffer)
            buffer = ibuffer;

        STRCPY(buffer, temp_dir);
        std::strcat(buffer, temp_name);
        return buffer;
    }
} // archive
} // boost

#else // defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#if defined(__hpux)
// (C) Copyright 2006 Boris Gubenko.
// HP-UX has a restriction that for multi-thread applications, (i.e.
// the ones compiled -mt) if argument to tmpnam is a NULL pointer, then,
// citing the tmpnam(3S) manpage, "the operation is not performed and a
// NULL pointer is returned". tempnam does not have this restriction, so,
// let's use tempnam instead.
 
#define tmpnam(X) tempnam(NULL,X)
 
namespace boost {
namespace archive {
    using ::tempnam;
} // archive
} // boost

#else // defined(__hpux)

namespace boost {
namespace archive {
    using std::tmpnam;
} // archive
} // boost

#endif // defined(__hpux)
#endif // defined(_WIN32) || defined(__WIN32__) || defined(WIN32)

/////////////////////////////////////////////
// invoke header for a custom archive test.
#if ! defined(BOOST_ARCHIVE_TEST)
#define BOOST_ARCHIVE_TEST text_archive.hpp
#endif

//#include <boost/test/test_tools.hpp>
#include <boost/detail/lightweight_test.hpp>

#define BOOST_CHECK( P ) \
    BOOST_TEST( (P) )
#define BOOST_REQUIRE( P )  \
    BOOST_TEST( (P) )
#define BOOST_CHECK_MESSAGE( P, M )  \
    ((P)? (void)0 : ::boost::detail::error_impl( (M) , __FILE__, __LINE__, BOOST_CURRENT_FUNCTION))
#define BOOST_REQUIRE_MESSAGE( P, M ) \
    BOOST_CHECK_MESSAGE( (P), (M) )
#define BOOST_CHECK_EQUAL( A , B ) \
    BOOST_TEST( (A) == (B) )

namespace boost { namespace detail {
inline void msg_impl(char const * msg, char const * file, int line, char const * function)
{
    std::cerr << file << "(" << line << "): " << msg << " in function '" << function << "'" << std::endl;
}
} } // boost::detail

#define BOOST_WARN_MESSAGE( P, M )  \
    ((P)? (void)0 : ::boost::detail::msg_impl( (M) , __FILE__, __LINE__, BOOST_CURRENT_FUNCTION))
#define BOOST_MESSAGE( M ) \
    BOOST_WARN_MESSAGE( true , (M) )

#define BOOST_CHECKPOINT( M ) \
    BOOST_WARN_MESSAGE( true , (M) )

#define BOOST_TEST_DONT_PRINT_LOG_VALUE( T ) 

#define BOOST_FAIL( M ) BOOST_REQUIRE_MESSAGE( false, (M) )
#define EXIT_SUCCESS 0

int test_main(int argc, char * argv[]);

int
main(int argc, char * argv[]){
    test_main(argc, argv);
    return boost::report_errors();
}

// the following is to ensure that when one of the libraries changes
// BJAM rebuilds and relinks the test.
/*
#include "text_archive.hpp"
#include "text_warchive.hpp"
#include "binary_archive.hpp"
#include "xml_archive.hpp"
#include "xml_warchive.hpp"
*/

#endif // BOOST_SERIALIZATION_TEST_TOOLS_HPP
