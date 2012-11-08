/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-16

  Copyright (C) 2009 Université de Grenoble 1
  Copyright (C) 2005,2006 EPFL

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
/**
   \file constants.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-16
 */
#ifndef __Constants_H
#define __Constants_H 1

#include <cmath>
#include <feel/feelcore/traits.hpp>

namespace Feel
{
/**
 * \namespace Feel::math
 *
 * The \c math namespace unifies calls to functions and constants with
 * respect to numerical types.
 */
namespace math
{
/**
 * \class Constant
 *\ingroup Core
 *\brief Structure that holds a constant with different numerical representation
 *
 * The \p Constant class/struct provides a place to put * the actual
 * definition of the value of each constant.  * It also turns an
 * implicit conversion operation * (which does not need any
 * parentheses) into an explicit call to the * correct function that
 * actually knows (and returns) the right value.
 */
template <typename Tag, typename Rep = double>
struct Constant
{
    Constant() {} // Returns the value of the constant.
    operator Rep() const; // Fully specialized for each Rep/Tag pair.
};


struct pi_tag {};
namespace float_constants
{
Constant<pi_tag, float> const pi;
}
namespace double_constants
{
Constant<pi_tag, double> const pi;
}
namespace long_double_constants
{
Constant<pi_tag, long double> const pi;
}
#if defined(FEELPP_HAS_QD_REAL)
namespace dd_real_constants
{
Constant<pi_tag, dd_real> const pi;
}
namespace qd_real_constants
{
Constant<pi_tag, qd_real> const pi;
}
#endif /*FEELPP_HAS_QD_REAL*/
template<> inline Constant<pi_tag, long double>::operator long double() const
{
    return 3.141592653589793238462643383279502884197L;
}
template<> inline Constant<pi_tag, double>::operator double() const
{
    return 3.141592653589793238462643383279502884197;
}
template<> inline Constant<pi_tag, float>::operator float() const
{
    return 3.141592653589793238462643383279502884197F;
}
#if defined(FEELPP_HAS_QD_REAL)
template<> inline Constant<pi_tag, dd_real>::operator dd_real() const
{
    return dd_real::_pi;
}
template<> inline Constant<pi_tag, qd_real>::operator qd_real() const
{
    return qd_real::_pi;
}
#endif /*FEELPP_HAS_QD_REAL*/
} // namespace math
} // namespace Feel
#endif /* __Constants_H */
