/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 14 Aug 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
//!
//! @file
//! Defines `Feel::range`
//!
#ifndef FEELPP_RANGE_HPP
#define FEELPP_RANGE_HPP 1

#include <vector>

#include <boost/range/irange.hpp>

namespace Feel {

using boost::irange;

/**
 @ingroup Core
 Provides range of integers

 \c range(stop) 
 \code range(start, stop[, step]) \endcode
 
 This is a versatile function to create lists containing arithmetic
 progressions. It is most often used in for loops. The arguments must be plain
 integers. If the step argument is omitted, it defaults to 1. If the start
 argument is omitted, it defaults to 0. The full form returns a list of plain
 integers [start, start + step, start + 2 * step, ...]. If step is positive, the
 last element is the largest start + i * step less than stop; if step is
 negative, the last element is the smallest start + i * step greater than
 stop. step must not be zero.

 @code
 range(10)     : 0,1,2,3,4,5,6,7,8,9
 range(1,10)   : 1,2,3,4,5,6,7,8,9
 range(1,10,2) : 1,3,5,7,9
 @endcode
 @note Boost.irange should be preferred to range
 @deprecated
*/
template<typename T>
std::vector<T>
range( T && beg, T && end, T && step )
{
    std::vector<T> v( static_cast<size_type>((end-beg)/step) );

    T n{beg};
    std::generate(v.begin(), v.end(), [&n,&step]{ n+=step; return n-step; }); 

    return v;
}
//!
//! @ingroup Core
//! return a range of integers between \c beg and \c end with a step of 1
//! @deprecated
//!
template<typename T>
std::vector<T>
range( T && beg, T && end )
{
    return range( std::forward<T>(beg), std::forward<T>(end), T(1) );
}
//!
//! @ingroup Core
//! return a range of integers between 0 and \c end with a step of 1
//! @deprecated
//!
template<typename T>
std::vector<T>
range( T && end )
{
    return range( static_cast<T>(0), std::forward<T>(end) );
}


}

#endif
