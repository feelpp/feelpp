/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-01-26

  Copyright (C) 2011-2012 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef __PATTERN_H
#define __PATTERN_H 1

#include <vector>
#include <feel/feelcore/feeltypes.hpp>



namespace Feel
{

//namespace vf
//{
enum  Pattern
{
    DEFAULT   = 1 << 0,
    EXTENDED  = 1 << 1,
    COUPLED   = 1 << 2,
    PATTERN_SYMMETRIC = 1 << 3,
    ZERO      = 1 << 4,
    HAS_NO_BLOCK_PATTERN = 1 << 5,
    HDG = 1 << 6,
    NONE= 1<<7
};

//} // namespace vf

namespace pattern
{
//!
//! @returns a pattern vector of the size of @param p using the ZERO pattern
//!
inline std::vector<uint32_type> toZero( std::vector<uint32_type> const& p )
{
    return std::vector( p.size(), uint32_type(Pattern::ZERO) );
}
}
} // namespace Feel


#endif // __PATTERN_H
