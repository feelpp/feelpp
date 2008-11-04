/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_void_cast.cpp: test implementation of run-time casting of void pointers

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
// <gennadiy.rozental@tfn.com>

#include "test_tools.hpp"
#include <boost/serialization/extended_type_info_typeid.hpp>
#include <boost/serialization/void_cast.hpp>

class Base1
{
    char a;
};

class Base2
{
    int b;
};

class Derived : public Base1, public Base2
{
    long c;
};

class MostDerived : public Derived
{
    char d[32];
};

int
test_main( int /* argc */, char* /* argv */[] )
{
    MostDerived d;
    MostDerived* pd =& d;
    Derived* pc = static_cast<Derived*>(pd);

    Base2* pb = static_cast<Base2*>(pd);
    Base1* pa = static_cast<Base1*>(pd);

    void* vpd = static_cast<void*>(pd);
    void* vpc = static_cast<void*>(pc);
    void* vpb = static_cast<void*>(pb);
    void* vpa = static_cast<void*>(pa);

    // simple casts only requiring table lookup
    BOOST_CHECK(vpc == boost::serialization::void_downcast(
        * boost::serialization::extended_type_info_typeid<Derived>::get_instance(),
        * boost::serialization::extended_type_info_typeid<Base1>::get_instance(), 
        vpa
    ));
    BOOST_CHECK(vpa == boost::serialization::void_upcast(
        * boost::serialization::extended_type_info_typeid<Derived>::get_instance(),
        * boost::serialization::extended_type_info_typeid<Base1>::get_instance(), 
        vpc  
    ));
    BOOST_CHECK(vpc == boost::serialization::void_downcast(  
        * boost::serialization::extended_type_info_typeid<Derived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base2>::get_instance(), 
        vpb
    ));
    BOOST_CHECK(vpb == boost::serialization::void_upcast(
        * boost::serialization::extended_type_info_typeid<Derived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base2>::get_instance(), 
        vpc
    ));
    BOOST_CHECK(vpd == boost::serialization::void_downcast(  
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Derived>::get_instance(), 
        vpc
    ));
    BOOST_CHECK(vpc == boost::serialization::void_upcast(
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Derived>::get_instance(), 
        vpd
    ));

    // note relationship between MostDerived and Base1 is automatically derived
    BOOST_CHECK(vpd == boost::serialization::void_downcast(  
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base1>::get_instance(), 
        vpa
    ));
    BOOST_CHECK(vpa == boost::serialization::void_upcast(
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base1>::get_instance(), 
        vpd
    ));

    // note relationship between MostDerived and Base2 is automatically derived
    BOOST_CHECK(vpd == boost::serialization::void_downcast(  
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base2>::get_instance(), 
        vpb
    ));
    BOOST_CHECK(vpb == boost::serialization::void_upcast(
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base2>::get_instance(), 
        vpd
    ));

    // need to double check to validate speed up optimization of derivations
    BOOST_CHECK(vpd == boost::serialization::void_downcast(  
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base1>::get_instance(), 
        vpa
    ));
    BOOST_CHECK(vpa == boost::serialization::void_upcast(
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base1>::get_instance(), 
        vpd
    ));
    BOOST_CHECK(vpd == boost::serialization::void_downcast(
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base2>::get_instance(), 
        vpb
    ));
    BOOST_CHECK(vpb == boost::serialization::void_upcast(
        * boost::serialization::extended_type_info_typeid<MostDerived>::get_instance(), 
        * boost::serialization::extended_type_info_typeid<Base2>::get_instance(), 
        vpd
    ));

    // check things that should fail
    BOOST_CHECK(NULL == boost::serialization::void_downcast(
        * boost::serialization::extended_type_info_typeid<Base2>::get_instance(),
        * boost::serialization::extended_type_info_typeid<Base1>::get_instance(), 
        vpa
    ));

    // note that a fundamental feature is that derived/base pairs are created
    // at compiler time so that all are registered before the main program starts
    // so leave the registration here at the end to verify this. Note bogus arguments
    // to workaround msvc 6 bug
    boost::serialization::void_cast_register<Derived, Base1>(
        static_cast<Derived *>(NULL),
        static_cast<Base1 *>(NULL)
    );
    boost::serialization::void_cast_register<Derived, Base2>(
        static_cast<Derived *>(NULL),
        static_cast<Base2 *>(NULL)
    );
    boost::serialization::void_cast_register<MostDerived, Derived>(
        static_cast<MostDerived *>(NULL),
        static_cast<Derived *>(NULL)
    );
    return EXIT_SUCCESS;
}

// EOF
