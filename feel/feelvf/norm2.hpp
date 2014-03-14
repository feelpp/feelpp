/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-01-13

  Copyright (C) 2014 Feel++ Consortium

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
#ifndef FEELPP_FEELVF_NORM2_HPP
#define FEELPP_FEELVF_NORM2_HPP 1

#include <feel/feelvf/stdmathfunctors.hpp>
#include <feel/feelvf/inner.hpp>

namespace Feel
{
namespace vf
{
/**
 * Norm-2 of the expression \p v
 *
 * \note it handles scalar, vectorial or matricial cases
 *
 * \return the norm 2 of the expression \f$\sqrt{\operatorname{tr}(v^T * v)}\f$
 */
template<typename ExprT>
inline
auto
norm2( ExprT v ) -> decltype( inner( v, v, mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() ) )
{
    return inner( v, v, mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
}


} // vf


} // Feel
#endif /* FEELPP_FEELVF_NORM2_HPP */
