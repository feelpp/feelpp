/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2021-04-13

  Copyright (C) 2021-present Feel++ Consortium

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
#pragma once
#include <type_traits>

namespace Feel
{
    enum class condenser_type
    {
        poisson = 1,
        stokes = 2
    };
    struct condenser_base
    {    
    };
    template <typename T, T V>
    struct condenser : public std::integral_constant<T, V>, public condenser_base
    {
    };
    using condenser_poisson = condenser<condenser_type, condenser_type::poisson>;
    using condenser_stokes = condenser<condenser_type, condenser_type::stokes>;

    template <typename T>
    using is_condenser = std::is_base_of<condenser_base, T>;
    template <typename T>
    inline constexpr bool is_condenser_v = std::is_base_of_v<condenser_base, T>;
    template <typename T>
    using is_condenser_poisson = std::is_same<condenser_poisson, T>;
    template <typename T>
    inline constexpr bool is_condenser_poisson_v = std::is_same_v<condenser_poisson, T>;
    template<typename T>
    using is_condenser_stokes = std::is_same<condenser_stokes,T>;
    template <typename T>
    inline constexpr bool is_condenser_stokes_v = std::is_same_v<condenser_stokes,T>;
} // namespace Feel

