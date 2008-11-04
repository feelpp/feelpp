//  (C) Copyright Eric Niebler 2004.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 Revision history:
   13 December 2004 : Initial version.
*/

#include <list>
#include <vector>
#include <boost/test/minimal.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "../../../boost/foreach.hpp"

///////////////////////////////////////////////////////////////////////////////
// int_iterator
//
typedef boost::counting_iterator<int> int_iterator;

///////////////////////////////////////////////////////////////////////////////
// define come containers
//
char my_ntcs_buffer[] = "\1\2\3\4\5";
char *my_ntcs  = my_ntcs_buffer;
int my_array[] = { 1,2,3,4,5 };
std::list<int> my_list(int_iterator(1),int_iterator(6));
std::pair<int_iterator,int_iterator> my_pair(int_iterator(1),int_iterator(6));

int const (&my_const_array)[5] = my_array;
char const *my_const_ntcs  = my_ntcs;
std::list<int> const &my_const_list = my_list;
std::pair<int_iterator,int_iterator> const &my_const_pair = my_pair;

///////////////////////////////////////////////////////////////////////////////
// define a user-defined collection type and teach BOOST_FOREACH how to enumerate it
//
namespace mine
{
    struct dummy {};
}

namespace boost
{
    template<>
    struct range_iterator<mine::dummy>
    {
        typedef char * type;
    };
    template<>
    struct range_const_iterator<mine::dummy>
    {
        typedef char const * type;
    };
    char * begin(mine::dummy&) {return 0;}
    char const * begin(mine::dummy const&) {return 0;}
    char * end(mine::dummy&) {return 0;}
    char const * end(mine::dummy const&) {return 0;}
}

///////////////////////////////////////////////////////////////////////////////
// to_vector_for
//
template<typename Range>
std::vector<int> to_vector_for( Range & rng )
{
    std::vector<int> vect;
    typedef BOOST_DEDUCED_TYPENAME boost::range_result_iterator<Range>::type iterator;
    for(iterator begin = boost::begin(rng), end = boost::end(rng);
        begin != end; ++begin)
    {
        vect.push_back(*begin);
    }
    return vect;
}

///////////////////////////////////////////////////////////////////////////////
// to_vector_foreach_byval
//
template<typename Range>
std::vector<int> to_vector_foreach_byval( Range & rng )
{
    std::vector<int> vect;
    typedef BOOST_DEDUCED_TYPENAME boost::range_result_iterator<Range>::type iterator;
    typedef BOOST_DEDUCED_TYPENAME boost::iterator_value<iterator>::type value;
    BOOST_FOREACH( value i, rng )
    {
        vect.push_back(i);
    }
    return vect;
}

///////////////////////////////////////////////////////////////////////////////
// to_vector_foreach_byref
//
template<typename Range>
std::vector<int> to_vector_foreach_byref( Range & rng )
{
    std::vector<int> vect;
    typedef BOOST_DEDUCED_TYPENAME boost::range_result_iterator<Range>::type iterator;
    typedef BOOST_DEDUCED_TYPENAME boost::iterator_reference<iterator>::type reference;
    BOOST_FOREACH( reference i, rng )
    {
        vect.push_back(i);
    }
    return vect;
}

///////////////////////////////////////////////////////////////////////////////
// mutate_foreach_byref
//
template<typename Range>
void mutate_foreach_byref( Range & rng )
{
    typedef BOOST_DEDUCED_TYPENAME boost::range_result_iterator<Range>::type iterator;
    typedef BOOST_DEDUCED_TYPENAME boost::iterator_reference<iterator>::type reference;
    BOOST_FOREACH( reference i, rng )
    {
        ++i;
    }
}


///////////////////////////////////////////////////////////////////////////////
// test_main
//   
int test_main( int, char*[] )
{
    // non-const containers by value
    BOOST_CHECK(to_vector_foreach_byval(my_array) == to_vector_for(my_array));
    BOOST_CHECK(to_vector_foreach_byval(my_ntcs)  == to_vector_for(my_ntcs));
    BOOST_CHECK(to_vector_foreach_byval(my_list)  == to_vector_for(my_list));
    BOOST_CHECK(to_vector_foreach_byval(my_pair)  == to_vector_for(my_pair));

    // const containers by value
    BOOST_CHECK(to_vector_foreach_byval(my_const_array) == to_vector_for(my_const_array));
    BOOST_CHECK(to_vector_foreach_byval(my_const_ntcs)  == to_vector_for(my_const_ntcs));
    BOOST_CHECK(to_vector_foreach_byval(my_const_list)  == to_vector_for(my_const_list));
    BOOST_CHECK(to_vector_foreach_byval(my_const_pair)  == to_vector_for(my_const_pair));

    // non-const containers by reference
    BOOST_CHECK(to_vector_foreach_byref(my_array) == to_vector_for(my_array));
    BOOST_CHECK(to_vector_foreach_byref(my_ntcs)  == to_vector_for(my_ntcs));
    BOOST_CHECK(to_vector_foreach_byref(my_list)  == to_vector_for(my_list));
    BOOST_CHECK(to_vector_foreach_byref(my_pair)  == to_vector_for(my_pair));

    // const containers by reference
    BOOST_CHECK(to_vector_foreach_byref(my_const_array) == to_vector_for(my_const_array));
    BOOST_CHECK(to_vector_foreach_byref(my_const_ntcs)  == to_vector_for(my_const_ntcs));
    BOOST_CHECK(to_vector_foreach_byref(my_const_list)  == to_vector_for(my_const_list));
    BOOST_CHECK(to_vector_foreach_byref(my_const_pair)  == to_vector_for(my_const_pair));

    // mutate the mutable collections
    mutate_foreach_byref(my_array);
    mutate_foreach_byref(my_ntcs);
    mutate_foreach_byref(my_list);

    // compare the mutated collections to the actual results
    std::pair<int_iterator,int_iterator> results(int_iterator(2),int_iterator(7));
    BOOST_CHECK(to_vector_foreach_byval(my_array) == to_vector_for(results));
    BOOST_CHECK(to_vector_foreach_byval(my_ntcs)  == to_vector_for(results));
    BOOST_CHECK(to_vector_foreach_byval(my_list)  == to_vector_for(results));

    // loop over a user-defined type (just make sure this compiles)
    mine::dummy d;
    BOOST_FOREACH( char c, d )
    {
        ((void)c); // no-op
    }

    return 0;
}
