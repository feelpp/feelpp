/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman
    Copyright (c) 2005-2006 Dan Marsden

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/

#include <boost/fusion/adapted/variant.hpp>

#include <boost/fusion/support/is_sequence.hpp>
#include <boost/fusion/support/is_view.hpp>
#include <boost/fusion/support/category_of.hpp>

#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/sequence/intrinsic/begin.hpp>
#include <boost/fusion/sequence/intrinsic/end.hpp>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/iterator/prior.hpp>
#include <boost/fusion/iterator/deref.hpp>
#include <boost/fusion/iterator/distance.hpp>
#include <boost/fusion/iterator/value_of.hpp>

#include <boost/detail/lightweight_test.hpp>
#include <boost/variant.hpp>

#include <boost/mpl/assert.hpp>

#include <boost/type_traits/is_same.hpp>

#include <string>

int main()
{
    namespace fusion = boost::fusion;
    typedef boost::variant<double, std::string> var_type;
    var_type var = "hello";

    BOOST_MPL_ASSERT((fusion::traits::is_sequence<var_type>));
    BOOST_MPL_ASSERT_NOT((fusion::traits::is_view<var_type>));
    BOOST_MPL_ASSERT((boost::is_same<
                      fusion::traits::category_of<var_type>::type, 
                      fusion::forward_traversal_tag>));

    BOOST_TEST(fusion::size(var) == 2);
    BOOST_TEST(fusion::distance(fusion::begin(var), fusion::end(var)) == 2);
    BOOST_TEST(*fusion::next(fusion::begin(var)) == "hello");
    BOOST_TEST(fusion::next(fusion::next(fusion::begin(var))) == fusion::end(var));
    BOOST_MPL_ASSERT((boost::is_same<
                      fusion::result_of::value_of<fusion::result_of::begin<var_type>::type>::type, 
                      double>));
    BOOST_MPL_ASSERT((boost::is_same<
                      fusion::result_of::value_of<fusion::result_of::next<fusion::result_of::begin<var_type>::type>::type>::type, 
                      std::string>));
    return boost::report_errors();
}
