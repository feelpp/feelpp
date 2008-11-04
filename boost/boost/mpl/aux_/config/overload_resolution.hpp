
#ifndef BOOST_MPL_AUX_CONFIG_OVERLOAD_RESOLUTION_HPP_INCLUDED
#define BOOST_MPL_AUX_CONFIG_OVERLOAD_RESOLUTION_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2002-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Source$
// $Date: 2007-10-02 22:32:05 +0200 (mar, 02 oct 2007) $
// $Revision: 39663 $

#include <boost/mpl/aux_/config/workaround.hpp>

#if    !defined(BOOST_MPL_CFG_BROKEN_OVERLOAD_RESOLUTION) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE) \
    && (   BOOST_WORKAROUND(__BORLANDC__, < 0x590) \
        || BOOST_WORKAROUND(__MWERKS__, < 0x3001) \
        )

#   define BOOST_MPL_CFG_BROKEN_OVERLOAD_RESOLUTION

#endif

#endif // BOOST_MPL_AUX_CONFIG_OVERLOAD_RESOLUTION_HPP_INCLUDED
