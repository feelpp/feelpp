
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
      Date: 2022-02-23

 Copyright (C) 2022 Feel++ Consortium


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
#if !defined(FEELPP_FEELCORE_HANA_HPP)
#define FEELPP_FEELCORE_HANA_HPP

#include <boost/hana.hpp>

namespace std
{
template <std::size_t n, typename... Types>
struct tuple_element<n, boost::hana::tuple<Types...>>
{
    using type = typename decltype( +boost::hana::tuple_t<Types...>[boost::hana::size_c<n>] )::type;
};

template <typename... Types>
struct tuple_size<boost::hana::tuple<Types...>> : public integral_constant<std::size_t, sizeof...( Types )>
{
};
} // namespace std

namespace boost
{
namespace hana
{
template <std::size_t n, typename... Types>
constexpr decltype( auto ) get( hana::tuple<Types...>& t )
{
    return t[hana::size_c<n>];
}

template <std::size_t n, typename... Types>
constexpr decltype( auto ) get( const hana::tuple<Types...>& t )
{
    return t[hana::size_c<n>];
}

template <std::size_t n, typename... Types>
constexpr decltype( auto ) get( hana::tuple<Types...>&& t )
{
    return static_cast<hana::tuple<Types...>&&>( t )[hana::size_c<n>];
}

template <std::size_t n, typename... Types>
constexpr decltype( auto ) get( const hana::tuple<Types...>&& t )
{
    return static_cast<const hana::tuple<Types...>&&>( t )[hana::size_c<n>];
}
} // namespace hana
} // namespace boost
#endif