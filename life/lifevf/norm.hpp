/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-04-17

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file norm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-04-17
 */
#ifndef __norm_H
#define __norm_H 1

namespace Life
{
namespace vf
{
/**
 * inner product of the left and right expressions
 *
 * \note it handles scalar, vectorial or matricial cases
 *
 * \return the inner product \f$\operatorname{tr}(l * r)\f$
 */
template<typename ExprL, typename ExprR>
inline
//Expr<Trace<Expr<vf_mul<Expr<Trans<ExprL> >, ExprR > > > >
__typeof__( trace( trans(ExprL()) * ExprR() ) )
inner_prod( ExprL l, ExprR r )
{
    return trace( trans(l) * r );
}

/**
 * Norm-2 of the expression \p v
 *
 * \note it handles scalar, vectorial or matricial cases
 *
 * \return the norm 2 of the expression \f$\sqrt{\operatorname{tr}(v^T * v)}\f$
 */
template<typename ExprT>
inline
//Expr<Sqrt<Expr<Trace<Expr<vf_mul<Expr<Trans<ExprT> >, ExprT > > > > > >
__typeof__( sqrt(  inner_prod( ExprT(), ExprT() ) ) )
norm2( ExprT v )
{
    return sqrt(  inner_prod( v, v ) );
}

} // vf


} // life
#endif /* __norm_H */
