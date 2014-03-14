/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-24

  Copyright (C) 2013 Université de Strasbourg

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
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-10-24
 */
#ifndef FEELPP_Mortar_H
#define FEELPP_Mortar_H 1

namespace Feel
{
namespace detail{
class mortar_base {};
}
class Mortar : public Feel::detail::mortar_base
{
public:
    static const bool is_mortar = true;

};
class NoMortar : public Feel::detail::mortar_base
{
public:
    static const bool is_mortar = false;

};

}
#endif /* FEELPP_Mortar_H */
