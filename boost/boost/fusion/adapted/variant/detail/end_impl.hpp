/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman
    Copyright (c) 2005-2006 Dan Marsden

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_VARIANT_END_IMPL_12112006_2137)
#define BOOST_FUSION_VARIANT_END_IMPL_12112006_2137

#include <boost/mpl/end.hpp>

namespace boost { namespace fusion {

    struct variant_tag;

    template<typename Seq, typename Iterator>
    struct variant_iterator;

    namespace extension
    {
        template<typename T>
        struct end_impl;

        template <>
        struct end_impl<variant_tag>
        {
            template <typename Seq>
            struct apply 
            {
                typedef variant_iterator<
                    Seq, 
                    typename mpl::end<typename Seq::types>::type> type;
    
                static type
                call(Seq& v)
                {
                    return type(v);
                }
            };
        };
    }
}}

#endif
