/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013 Feel++ Consortium

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
   \file element_product.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_ELEMENT_PRODUCT_HPP)
#define FEELPP_ELEMENT_PRODUCT_HPP 1

#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

/**
   Computes the element wise product of two vectors and eventually in parallel
   \param v1 vector (eventually distributed)
   \param v2 vector (eventually distributed)

   \return the element wise product of \p v1 and \p v2
*/
template <typename ElementType>
ElementType
prod( ElementType const& v1, ElementType const& v2,  typename boost::enable_if<boost::is_base_of<FunctionSpaceBase::ElementBase,ElementType> >::type* dummy = 0  )
{
    FEELPP_ASSERT( v1.functionSpace() == v2.functionSpace() ).error( "incompatible function spaces" );

    typedef typename type_traits<typename ElementType::value_type>::real_type real_type;

    ElementType _t( v1.functionSpace() );
    size_type s = v1.localSize();
    size_type start = v1.firstLocalIndex();

    for ( size_type i = 0; i < s; ++i )
        _t.operator()( start+i ) = v1.operator()( start + i )* v2.operator()( start + i );

    return _t;
}

}

#endif /* FEELPP_ELEMENT_PRODUCT_HPP */
