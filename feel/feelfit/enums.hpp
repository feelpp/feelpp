/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Vincent Huber <vincent.huber@cemosis.Fr>
Date: 01-04-2016

Copyright (C) 2007-2008 Universite Joseph Fourier (Grenoble I)

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

#ifndef FEELPP_INTERPOLATOR_ENUMS_HPP
#define FEELPP_INTERPOLATOR_ENUMS_HPP 1

namespace Feel
{

enum InterpolationType
{
  P0=0,
  P1,
  Spline,
  Akima
};

extern std::map<std::string, InterpolationType> InterpolationTypeMap;

/*
 * left : f(x_{i-1} < x < x_i     ) = x_i
 * right: f(x_i     < x < x_{i+1} ) = x_i
 * center: f((x_{i-1} + x_{i})/2 < x < (x_{i+1} + x_{i})/2) = x_i
 */
enum InterpolationType_P0
{
  left=0,
  right,
  center
};

/*
 * zero     : P1(x ! \in{x_0, x_n} ) = 0
 * constant : P1(x = x0) = x0 || P1(x = xn) = xn
 * extrapol : linear extrapolation
 */
enum ExtrapolationType_P1
{
  zero=0,
  constant,
  extrapol
};

/*
 * natural : f''(x0) = f''(xn) = 0
 * clamped : f' (x0) = f' (xn) = 0
 */
enum ExtrapolationType_spline
{
  natural=0,
  clamped
};

} // namespace Feel

#endif //FEELPP_INTERPOLATOR_ENUMS_HPP
