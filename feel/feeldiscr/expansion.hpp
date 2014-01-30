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
   \file expansion.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_EXPANSION_HPP)
#define FEELPP_EXPANSION_HPP 1

namespace Feel {

/**
   Given a set of coefficient \p c and a set of finite element function \p b, build
   \f[
   w = \sum_{i=1}^{M} c_i b_i
   \f]
   The last argument \p M allows to build only a subset of the expansion
 */
template<typename ElementType, typename CoeffType>
ElementType
expansion( std::vector<ElementType> const& b, CoeffType const& c, int M = -1 )
{
    auto res = b[0].functionSpace()->element();
    res.zero();
    if ( ( M == -1 ) || M > c.size() ) M = c.size() ;
    CHECK( (c.size() <= M) )
        << "Invalid coefficient or basis function elements "
        << "M=" << M << " coeff: " << b.size() << " elements: " << c.size();
    for( int i = 0; i < M; ++i )
    {
        res.add( c[i], b[i] );
    }

    return res;
}

}

#endif /* FEELPP_EXPANSION_HPP */
