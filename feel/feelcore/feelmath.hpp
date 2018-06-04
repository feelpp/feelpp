//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! Date: 27 Nov 2016
//!
//! Copyright (C) 2016 Feel++ Consortium
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!

#ifndef FEELPP_FEELMATH_HPP
#define FEELPP_FEELMATH_HPP 1

#include <cstdint>

namespace Feel
{

//!
//! compute integer power
//! \code
//! uint64_t i = ipow( 2, 3 );
//! \endcode
//!
constexpr int64_t ipow(int64_t base, int exp, int64_t result = 1)
{
    return exp < 1 ? result : ipow(base*base, exp/2, (exp % 2) ? result*base : result);
}


} // Feel
#endif
