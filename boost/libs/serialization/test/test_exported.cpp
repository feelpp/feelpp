/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_exported.cpp

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

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#include <boost/archive/archive_exception.hpp>
#include "test_tools.hpp"
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE(BOOST_ARCHIVE_TEST)

class polymorphic_base
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & /* ar */, const unsigned int /* file_version */){
    }
public:
    virtual ~polymorphic_base(){};
};

BOOST_IS_ABSTRACT(polymorphic_base)

class polymorphic_derived1 : public polymorphic_base
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int /* file_version */){
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(polymorphic_base);
    }
public:
    virtual ~polymorphic_derived1(){}
};

BOOST_CLASS_EXPORT(polymorphic_derived1)

// MWerks users can do this to make their code work
BOOST_SERIALIZATION_MWERKS_BASE_AND_DERIVED(polymorphic_base, polymorphic_derived1)

class polymorphic_derived2 : public polymorphic_base
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int /* file_version */){
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(polymorphic_base);
    }
public:
    virtual ~polymorphic_derived2(){}
};

BOOST_CLASS_EXPORT(polymorphic_derived2)

// MWerks users can do this to make their code work
BOOST_SERIALIZATION_MWERKS_BASE_AND_DERIVED(polymorphic_base, polymorphic_derived2)

// save exported polymorphic class
void save_exported(const char *testfile)
{
    test_ostream os(testfile, TEST_STREAM_FLAGS);
    test_oarchive oa(os);

    polymorphic_base *rb1 = new polymorphic_derived1;
    polymorphic_base *rb2 = new polymorphic_derived2;

    // export will permit correct serialization
    // through a pointer to a base class
    oa << BOOST_SERIALIZATION_NVP(rb1);
    oa << BOOST_SERIALIZATION_NVP(rb2);

    delete rb1;
    delete rb2;
}

// save exported polymorphic class
void load_exported(const char *testfile)
{
    test_istream is(testfile, TEST_STREAM_FLAGS);
    test_iarchive ia(is);

    polymorphic_base *rb1 = NULL;
    polymorphic_base *rb2 = NULL;

    // export will permit correct serialization
    // through a pointer to a base class
    ia >> BOOST_SERIALIZATION_NVP(rb1);
    BOOST_CHECK_MESSAGE(
        boost::serialization::type_info_implementation<polymorphic_derived1>
            ::type::get_instance()
        == boost::serialization::type_info_implementation<polymorphic_base>
            ::type::get_derived_extended_type_info(*rb1),
        "restored pointer b1 not of correct type"
    );

    ia >> BOOST_SERIALIZATION_NVP(rb2);
    BOOST_CHECK_MESSAGE(
        boost::serialization::type_info_implementation<polymorphic_derived2>
            ::type::get_instance()
        == boost::serialization::type_info_implementation<polymorphic_base>
            ::type::get_derived_extended_type_info(*rb2),
        "restored pointer b2 not of correct type"
    );

    delete rb1;
    delete rb2;
}

int
test_main( int /* argc */, char* /* argv */[] )
{
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    save_exported(testfile);
    load_exported(testfile);
    std::remove(testfile);
    return EXIT_SUCCESS;
}

// EOF
