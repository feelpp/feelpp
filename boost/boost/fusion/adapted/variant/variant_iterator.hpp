/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman
    Copyright (c) 2005-2006 Dan Marsden

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_VARIANT_ITERATOR_12112006_1617)
#define BOOST_FUSION_VARIANT_ITERATOR_12112006_1617

#include <boost/fusion/iterator/iterator_facade.hpp>

#include <boost/mpl/next.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/distance.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <boost/variant/get.hpp>
#include <boost/detail/workaround.hpp>

namespace boost { namespace fusion
{
    struct forward_traversal_tag;

    template<typename Variant, typename MPLIterator>
    struct variant_iterator
        : iterator_facade<variant_iterator<Variant, MPLIterator>, forward_traversal_tag>
    {
        typedef Variant variant_type;
        typedef MPLIterator iterator;

        variant_iterator(Variant& var)
            : var_(var) {}
        Variant& var_;

        template<typename Iterator>
        struct next
        {
            typedef variant_iterator<
                typename Iterator::variant_type,
                typename mpl::next<typename Iterator::iterator>::type> type;

            static type
            call(Iterator const& i)
            {
                return type(i.var_);
            }
        };

        template<typename I1, typename I2>
        struct distance
            : mpl::distance<
            typename I1::iterator,
            typename I2::iterator>
        {
            typedef typename mpl::distance<
                typename I1::iterator,
                typename I2::iterator>::type type;

            static type call(I1 const& i1, I2 const& i2)
            {
                return type();
            }
        };

        template<typename Iterator>
        struct value_of
            : mpl::deref<typename Iterator::iterator>
        {};

        template <typename Iterator>
        struct deref
        {
            typedef typename
                mpl::eval_if<
                    is_const<typename Iterator::variant_type>
                  , add_const<typename mpl::deref<typename Iterator::iterator>::type>
                  , mpl::deref<typename Iterator::iterator>
                >::type
            value_type;

            typedef typename
                add_reference<value_type>::type
            type;

#if BOOST_WORKAROUND(BOOST_MSVC, < 1400)
// for some unknown reason (compiler bug) VC7.1 gets confused with
// variant and optional get functions.
            static type
            call(Iterator const & it)
            {
                boost::detail::variant::get_visitor<type> v;
                typedef typename mpl::deref<typename Iterator::iterator>::type type;
                if (type* result = it.var_.apply_visitor(v))
                    return *result;
                it.var_ = type(); // prime the variant
                return *it.var_.apply_visitor(v); // no-throw!
            }
#else
            static type
            call(Iterator const & it)
            {
                typedef typename mpl::deref<typename Iterator::iterator>::type type;
                if (type* result = boost::get<type>(&it.var_))
                    return *result;
                it.var_ = type(); // prime the variant
                return *boost::get<type>(&it.var_); // no-throw!
            }
#endif
        };
    };
}}

#endif
