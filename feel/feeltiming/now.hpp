/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-03-20

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file now.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-03-20
 */
#if !defined(FEELPP_TIMING_NOW_HPP)
#define FEELPP_TIMING_NOW_HPP 1

#include <chrono>

#include <sys/time.h>
#include <boost/cstdint.hpp>
#include <feel/feelcore/feel.hpp>


namespace Feel
{
using time_point = std::chrono::high_resolution_clock::time_point;
using seconds_t = double;
using microseconds_t = double;
using time_quantum_t = uint64_type;
using cycles_t = Feel::uint64_type;

namespace details
{
microseconds_t to_microseconds( time_quantum_t const t );
    
time_point now();

} // details
} // Feel


#endif /* FEELPP_TIMING_NOW_HPP */
