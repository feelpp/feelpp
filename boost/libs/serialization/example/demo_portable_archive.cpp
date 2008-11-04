/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
//
// demo_portable_archive.cpp
//
// (C) Copyright 2002-4 Robert Ramey - http://www.rrsd.com .
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// should pass compilation and execution
#include <sstream>

#define BOOST_ARCHIVE_SOURCE
#include "portable_binary_oarchive.hpp"
#include "portable_binary_iarchive.hpp"

#include <cstdlib>
#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ using ::rand; }
#endif

// the following is required to be sure the "EXPORT" works if it is used
#define CUSTOM_ARCHIVE_TYPES portable_binary_oarchive,portable_binary_iarchive

class A
{
    friend class boost::serialization::access;
    int i;
    unsigned int ui;
    long l;
    unsigned long ul;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* version */){
        ar & i & ui & l & ul ;
    }
public:
    bool operator==(const A & rhs) const {
        return
            i == rhs.i && ui == rhs.ui && l == rhs.l && ul == rhs.ul
        ;
    }
    A() :
        i(std::rand()),
        ui(std::rand()),
        l(std::rand()),
        ul(std::rand())
    {}
};

int main( int /* argc */, char* /* argv */[] )
{
    const A a;
    A a1;

    std::stringstream ss;
    {   
        portable_binary_oarchive pboa(ss);
        pboa << a;
    }
    {
        portable_binary_iarchive pbia(ss);
        pbia >> a1;
    }
    return !(a == a1);
}


