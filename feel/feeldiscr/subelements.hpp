/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Tue Feb 25 06:38:19 2014

   Copyright (C) 2014 Feel++ Consortium

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
#ifndef FEELPP_SUBELEMENTS_HPP
#define FEELPP_SUBELEMENTS_HPP 1

#include <feel/feeldiscr/detail/createelementvector.hpp>

namespace Feel {

/**
   extract subelements from an element of a product of spaces

   @param e an element of a product of spaces
   @param n a vector of string to name the sub elements
   @return a vector of subelements
 */
template<typename EltType>
typename fusion::result_of::accumulate<typename EltType::functionspace_type::functionspace_vector_type,
                                       fusion::vector<>,
                                       Feel::detail::CreateElementVector<EltType> >::type
subelements( EltType const& e, std::vector<std::string> const& n )
{
    return fusion::accumulate( e.functionSpaces(), fusion::vector<>(), Feel::detail::CreateElementVector<EltType>( e, n ) );
}

} // Feel


#endif /* FEELPP_SUBELEMENTS_HPP */
