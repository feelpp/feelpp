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
   \file tic.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-03-20
 */
#if !defined(FEELPP_TIMING_TIC_HPP)
#define FEELPP_TIMING_TIC_HPP 1

#include <feel/feelcore/environment.hpp>
#include <feel/feeltiming/timer.hpp>
#include <feel/feeltiming/now.hpp>

namespace Feel
{
namespace details
{
struct SecondBasedTimer
{
    static void print( std::string const& msg, const double& val )
    {
        if ( Environment::isMasterRank() )
        {
            if ( !msg.empty() )
                std::cout << "[" << msg << "] Time : " << val << "s\n";
            else
                std::cout << "Time : " << val << "s\n";
        }
    }
    static inline double  time()
    {
        return Feel::details::now();
    }
};

counter<double,SecondBasedTimer> const sec_timer = {};
}  // details
} // Feel

namespace Feel
{
namespace time
{
inline void tic()
{
    Feel::details::sec_timer.tic();
}

inline double  toc( std::string const& msg = "", bool display = true )
{
    return Feel::details::sec_timer.toc( msg, display );
}
} // time
} // Feel

namespace Feel
{
// Convenience namespace injection from time:: into Feel::
using time::tic;
using time::toc;
}

#endif /* FEELPP_TIMING_TIC_HPP */
