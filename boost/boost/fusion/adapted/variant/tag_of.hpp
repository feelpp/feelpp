/*=============================================================================
    Copyright (c) 2001-2006 Joel de Guzman
    Copyright (c) 2005-2006 Dan Marsden

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_FUSION_VARIANT_TAG_OF_12112006_1704)
#define BOOST_FUSION_VARIANT_TAG_OF_12112006_1704

#include <boost/fusion/support/tag_of_fwd.hpp>
#include <boost/variant/variant_fwd.hpp>

namespace boost { namespace fusion 
{
    struct variant_tag;

    namespace traits
    {
        template<BOOST_VARIANT_ENUM_PARAMS(typename T)>
#if defined(BOOST_PARTIAL_SPECIALIZATION_EXPLICT_ARGS)
        struct tag_of<boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>, void >
#else
        struct tag_of<boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)> >
#endif
        {
            typedef variant_tag type;
        };
    }
}}

#endif
