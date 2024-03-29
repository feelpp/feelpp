/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannnes@feelpp.org>
 Date: 12 dec. 2016

 Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_UNWRAPPTR_HPP
#define FEELPP_UNWRAPPTR_HPP 1

#include <feel/feelcore/traits.hpp>


namespace Feel {

template<typename T>
using unwrap_ptr_t = remove_shared_ptr_type<std::remove_pointer_t<T>>;

/**
 * @return the c++ objet pointed by c++ pointer or shared_ptr or same object otherwise
 */
template<typename C>
unwrap_ptr_t<C> const&
unwrap_ptr( C const& c )
{
    if constexpr ( is_ptr_or_shared_ptr<C>() )
        return *c;
    else
        return c;
}

template<typename C>
unwrap_ptr_t<C> &
unwrap_ptr( C & c )
{
    if constexpr ( is_ptr_or_shared_ptr<C>() )
        return *c;
    else
        return c;
}

template<typename C>
unwrap_ptr_t<C> &
unwrap_ptr( C * c )
{
    return *c;
}
#if 0
template<typename C>
unwrap_ptr_t<C>&&
unwrap_ptr( C && c )
{
    if constexpr ( is_ptr_or_shared_ptr<std::decay_t<C>>() )
        return *std::forward<C>(c);
    else
        return std::forward<C>(c);
}
#endif


} // Feel
#endif
