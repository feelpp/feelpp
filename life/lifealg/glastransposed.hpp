/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-08-17

  Copyright (C) 2005,2006 EPFL

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
   \file glastransposed.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-08-17
 */
#ifndef __LIFE_GLAS_TRANSPOSED_HPP
#define __LIFE_GLAS_TRANSPOSED_HPP 1
namespace Life
{
namespace glas
{
template <typename Matrix>
inline
typename transposed_return<L *>::return_type
transposed( Matrix &m )
{
    return typename transposed_return<L *>::return_type(linalg_cast(l));
}

} // namespace glas
} // namespace Life

#endif /*  __LIFE_GLAS_TRANSPOSED_HPP */
