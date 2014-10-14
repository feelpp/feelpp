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
template<typename T, typename C>
static inline void AvgMinMax(const T* const in, T* const inout, const int* const len, C* type)
{
    for(int i = 0; i < *len; i += 3) {
        inout[0 + 3 * i] += in[0 + 3 * i];
        inout[1 + 3 * i] = std::min(in[1 + 3 * i], inout[1 + 3 * i]);
        inout[2 + 3 * i] = std::max(in[2 + 3 * i], inout[2 + 3 * i]);
    }
}
} // Functor
} // Feel


#endif /* FEELPP_FUNCTORS_HPP */
