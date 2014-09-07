/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-11-28

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file lineordering.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-11-28
 */
#ifndef __LineOrdering_H
#define __LineOrdering_H 1

namespace Feel
{
/// \cond DETAIL
namespace details
{
struct vertex
{
    static uint16_type f2p( uint16_type /*f*/, uint16_type /*p*/ )
    {
        return invalid_uint16_type_value;
    }
    static uint16_type f2e( uint16_type /*f*/, uint16_type /*e*/ )
    {
        throw std::logic_error( "invalid call to line::f2e" );
        return invalid_uint16_type_value;
    }
    static uint16_type e2p( uint16_type /*e*/, uint16_type /*p*/ )
    {
        return invalid_uint16_type_value;
    }

    std::vector<uint16_type> entity( uint16_type /*topo_dim*/, uint16_type /*id*/ ) const
    {
        std::vector<uint16_type> __entity( 1 );
        __entity[0] = 0;
        return __entity;
    }
};

/**
 * \class line
 */
template<uint16_type Order>
struct line
{
    //static uint16_type f2p( uint16_type /*f*/, uint16_type /*p*/ ) { throw std::logic_error( "invalid call to line::f2p" ); return uint16_type(-1); }
    static uint16_type f2p( uint16_type f, uint16_type /*p*/ )
    {
        return f;
    }
    static uint16_type f2e( uint16_type /*f*/, uint16_type /*e*/ )
    {
        throw std::logic_error( "invalid call to line::f2e" );
        return uint16_type( -1 );
    }

    static uint16_type e2p( uint16_type /*e*/, uint16_type p )
    {
        return __e2p[p];
    }

    static const uint16_type __e2p[11];

    std::vector<uint16_type> entity( uint16_type /*topo_dim*/, uint16_type /*id*/ ) const
    {
        std::vector<uint16_type> __entity( 2 );
        __entity[0] = 0;
        __entity[1] = 1;
        return __entity;
    }
};
template<uint16_type Order> const uint16_type  line<Order>::__e2p[11] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
}
/// \endcond
}
#endif /* __LineOrdering_H */
