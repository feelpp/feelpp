/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-06-11

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
/**
   \file functors.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-06-11
 */
#if !defined(FEELPP_CORE_FUNCTORS_HPP)
#define FEELPP_CORE_FUNCTORS_HPP 1

#include <vector>

namespace Feel
{
namespace Functor
{
/**
 * functors to add two vectors \p x and \p y.
 *
 * The structure must provide \c .size() and \c operator[] to access entries
 *
 * Requirements: x and y must be of the same size
 */
template<typename T>
struct AddStdVectors
{
    std::vector<T>
    operator()( std::vector<T> const& x, std::vector<T> const& y )
        {
            DCHECK( x.size() == y.size() ) << "incompatible vector size()";
            std::vector<T> z(x.size());
            const int s = x.size();
            for(int i = 0; i < s; ++i ) z[i] = x[i] + y[i];
            return z;
        }
};
} // Functor
} // Feel


#endif /* FEELPP_FUNCTORS_HPP */
