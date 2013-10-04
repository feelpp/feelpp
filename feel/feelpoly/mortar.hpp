/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2012-06-26

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
   \file mortar.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-06-26
 */
#ifndef __Mortar_H
#define __Mortar_H 1

namespace Feel
{
namespace detail
{
template<typename P>
class Mortar
{
public:
};
} // detail
/**
 * \class Mortar
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class Mortar
{
public:

    template<typename P>
    class apply
    {
        typedef typename detail::Mortar<P> type;
    };

};

template<typename P>
struct is_mortar : public boost::is_base_of<Mortar,P>
{};

} // Feel
#endif /* __Mortar_H */
