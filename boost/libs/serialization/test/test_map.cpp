/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_map.cpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// should pass compilation and execution

#include <algorithm>
#include <vector>
#include <fstream>
#include <cstddef> // size_t

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

#include <cstdio>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::rand; 
    using ::size_t;
}
#endif

#include "test_tools.hpp"
#include <boost/preprocessor/stringize.hpp>
#include BOOST_PP_STRINGIZE(BOOST_ARCHIVE_TEST)

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/map.hpp>

#include "A.hpp"

///////////////////////////////////////////////////////
// a key value initialized with a random value for use
// in testing STL map serialization
struct random_key {
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(
        Archive & ar, 
        const unsigned int /* file_version */
    ){
        ar & boost::serialization::make_nvp("random_key", m_i);
    }
    int m_i;
    random_key() : m_i(std::rand()){};
    bool operator<(const random_key &rhs) const {
        return m_i < rhs.m_i;
    }
    bool operator==(const random_key &rhs) const {
        return m_i == rhs.m_i;
    }
    operator std::size_t () const {    // required by hash_map
        return m_i;
    }
};  

BOOST_BROKEN_COMPILER_TYPE_TRAITS_SPECIALIZATION(random_key)

void
test_map(){
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    BOOST_MESSAGE("map");
    // test map of objects
    std::map<random_key, A> amap;
    amap.insert(std::make_pair(random_key(), A()));
    amap.insert(std::make_pair(random_key(), A()));
    {   
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os);
        oa << boost::serialization::make_nvp("amap", amap);
    }
    std::map<random_key, A> amap1;
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is);
        ia >> boost::serialization::make_nvp("amap", amap1);
    }
    BOOST_CHECK(amap == amap1);
    std::remove(testfile);
}

void
test_map_2(){
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    BOOST_MESSAGE("map_2");
    std::pair<int, int> a(11, 22);
    std::map<int, int> b;
    b[0] = 0;
    b[-1] = -1;
    b[1] = 1;
    {
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        std::pair<int, int> * const pa = &a;
        std::map<int, int> * const pb = &b;
        test_oarchive oa(os);
        oa << BOOST_SERIALIZATION_NVP(pb);
        oa << BOOST_SERIALIZATION_NVP(pa);
    }
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        std::pair<int, int> *pa = 0;
        std::map<int, int> *pb = 0;
        test_iarchive ia(is);
        ia >> BOOST_SERIALIZATION_NVP(pb);
        ia >> BOOST_SERIALIZATION_NVP(pa);
        delete pa;
        delete pb;
    }
    std::remove(testfile);
}

void
test_multimap(){
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    BOOST_MESSAGE("multimap");
    std::multimap<random_key, A> amultimap;
    amultimap.insert(std::make_pair(random_key(), A()));
    amultimap.insert(std::make_pair(random_key(), A()));
    {   
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os);
        oa << boost::serialization::make_nvp("amultimap", amultimap);
    }
    std::multimap<random_key, A> amultimap1;
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is);
        ia >> boost::serialization::make_nvp("amultimap", amultimap1);
    }
    BOOST_CHECK(amultimap == amultimap1);
    std::remove(testfile);
}

#ifdef BOOST_HAS_HASH
#include <boost/serialization/hash_map.hpp>

namespace BOOST_STD_EXTENSION_NAMESPACE {
    template<>
    struct hash<random_key>{
        std::size_t operator()(const random_key& r) const {
            return (std::size_t)r;
        }
    };
} // namespace BOOST_STD_EXTENSION_NAMESPACE 

void
test_hash_map(){
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    BOOST_CHECKPOINT("hash_map");
    // test hash_map of objects
    BOOST_STD_EXTENSION_NAMESPACE::hash_map<random_key, A> ahash_map;
    ahash_map.insert(std::make_pair(random_key(), A()));
    ahash_map.insert(std::make_pair(random_key(), A()));
    {   
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os);
        oa << boost::serialization::make_nvp("ahashmap",ahash_map);
    }
    BOOST_STD_EXTENSION_NAMESPACE::hash_map<random_key, A> ahash_map1;
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is);
        ia >> boost::serialization::make_nvp("ahashmap",ahash_map1);
    }

    // at least one library - MSL notes: it doesn't make much sense
    // to implement the == operator for hash collections - but goes ahead
    // does it anyway even though it doesn't seem to work.  So sort into
    // vectors and then compare. Currently, it seems that STLPort versions
    // greater than 5.0 don't support the == operator on hashed sets
    #if ! defined(__SGI_STL_PORT) || (__SGI_STL_PORT < 0x500)
    BOOST_CHECK(ahash_map == ahash_map1);
    #else
        std::vector< std::pair<random_key, A> > tvec, tvec1;
        std::copy(ahash_map.begin(), ahash_map.end(), std::back_inserter(tvec));
        std::sort(tvec.begin(), tvec.end());
        std::copy(ahash_map1.begin(), ahash_map1.end(), std::back_inserter(tvec1));
        std::sort(tvec1.begin(), tvec1.end());
        BOOST_CHECK(tvec == tvec1);
    #endif
    std::remove(testfile);
}

void
test_hash_multimap(){
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(NULL != testfile);

    BOOST_CHECKPOINT("hash_multimap");
    BOOST_STD_EXTENSION_NAMESPACE::hash_multimap<random_key, A> ahash_multimap;
    ahash_multimap.insert(std::make_pair(random_key(), A()));
    ahash_multimap.insert(std::make_pair(random_key(), A()));
    {   
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os);
        oa << boost::serialization::make_nvp("ahash_multimap", ahash_multimap);
    }
    BOOST_STD_EXTENSION_NAMESPACE::hash_multimap<random_key, A> ahash_multimap1;
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is);
        ia >> boost::serialization::make_nvp("ahash_multimap", ahash_multimap1);
    }
    #if ! defined(__SGI_STL_PORT) || (__SGI_STL_PORT < 0x500)
    BOOST_CHECK(ahash_multimap == ahash_multimap1);
    #else
        std::vector< std::pair<random_key, A> > tvec, tvec1;
        tvec.clear();
        tvec1.clear();
        std::copy(ahash_multimap.begin(), ahash_multimap.end(), std::back_inserter(tvec));
        std::sort(tvec.begin(), tvec.end());
        std::copy(ahash_multimap1.begin(), ahash_multimap1.end(), std::back_inserter(tvec1));
        std::sort(tvec1.begin(), tvec1.end());
        BOOST_CHECK(tvec == tvec1);
    #endif
    std::remove(testfile);
}
#endif

int test_main( int /* argc */, char* /* argv */[] )
{
    test_map();
    test_map_2();
    test_multimap();
    #ifdef BOOST_HAS_HASH
    test_hash_map();
    test_hash_multimap();
    #endif
    return EXIT_SUCCESS;
}
