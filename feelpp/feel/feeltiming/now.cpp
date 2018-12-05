//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
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
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 29 Apr 2018
//! @copyright 2018 Feel++ Consortium
//!
#include <feel/feeltiming/now.hpp>

namespace Feel {

namespace details {

// TODO: Invesitgate boost/std::chrono. (18.10.2012.) (Domagoj Saric)
static double const timer_ticks_per_microsecond( CLOCKS_PER_SEC / 1000000.0 );

static microseconds_t const
inverse_timer_ticks_per_ms( static_cast<microseconds_t>(1)
                            / timer_ticks_per_microsecond
                            );

microseconds_t to_microseconds( time_quantum_t const t )
{
    return static_cast<microseconds_t>(t) * inverse_timer_ticks_per_ms;
}

//time_quantum_t now()
time_point now()
{
#if 0
    timeval tp;
    BOOST_VERIFY( ::gettimeofday( &tp, NULL ) == 0 );
    return static_cast<time_quantum_t>( tp.tv_sec ) * 1000000 + tp.tv_usec;
#endif
    return std::chrono::high_resolution_clock::now();
}


} // details
} // Feel
