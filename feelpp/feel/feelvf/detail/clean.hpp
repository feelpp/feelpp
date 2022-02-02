/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2016 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file clean.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined( FEELPP_DETAIL_CLEAN_HPP )
#define FEELPP_DETAIL_CLEAN_HPP 1

#if 0
#include <feel/feelcore/traits.hpp>

namespace Feel
{
namespace vf
{

/// \cond DETAIL
namespace detail
{
#if 0
template<typename TheArgs, typename Tag>
struct clean_type
{
    typedef typename boost::remove_pointer<
        typename boost::remove_const<
            typename boost::remove_reference<
                typename parameter::binding<TheArgs, Tag>::type
                >::type
            >::type
        >::type type;
};
template<typename TheArgs, typename Tag, typename Default>
struct clean2_type
{
    typedef typename boost::remove_pointer<
    typename boost::remove_const<
    typename boost::remove_reference<
    typename parameter::binding<TheArgs, Tag, Default>::type
    >::type
    >::type
    >::type type;
};
#else
template<typename TheArgs, typename Tag>
using clean_type = std::decay_t<typename parameter::binding<TheArgs, Tag>::type>;

template<typename TheArgs, typename Tag, typename Default>
using clean2_type = std::decay_t<typename parameter::binding<TheArgs, Tag, Default>::type>;


#endif
} } }

#endif

#endif /* FEELPP_DETAIL_CLEAN_HPP */
