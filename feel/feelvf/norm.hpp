/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-04-17

  Copyright (C) 2010 Universite Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-17
 */
#ifndef __norm_H
#define __norm_H 1

namespace Feel
{
namespace vf
{
/**
 * inner product of the left and right expressions
 *
 * \note it handles vectorial expressions
 *
 * \return the inner product \f$(l * r^T)\f$
 */
template<typename ExprL, typename ExprR>
inline
auto
outer_prod( ExprL l, ExprR r ) -> decltype( l * trans( r ) )
{
    return l * trans( r );
}

/**
 * double dot product of the left and right expressions
 *
 * \note it handles matricial(rank 2 tensors) expressions
 *
 * \return the double dot product \f$(l : r)\f$
 */
template<typename ExprL, typename ExprR>
inline
auto
ddot( ExprL l, ExprR r ) -> decltype( trace( trans( l ) * r ) )
{
    return trace( trans( l ) * r );
}


} // vf


} // feel
#endif /* __norm_H */
