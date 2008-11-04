/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman
    Copyright (c) 2005-2006 Dan Marsden

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_VARIANT_SIZE_IMPL_12112006_2115)
#define BOOST_FUSION_VARIANT_SIZE_IMPL_12112006_2115

#include <boost/mpl/size.hpp>

namespace boost { namespace fusion {

    struct variant_tag;

    namespace extension
    {
        template<typename T>
        struct size_impl;

        template<>
        struct size_impl<variant_tag>
        {
            template<typename Sequence>
            struct apply : mpl::size<typename Sequence::types> 
            {};
        };
    }
}}

#endif
