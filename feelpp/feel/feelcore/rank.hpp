/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 30 juil. 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_RANK_HPP
#define FEELPP_RANK_HPP 1

#include <feel/feelcore/traits.hpp>


namespace Feel {


/**
 * @return the local MPI rank of the data structure @p c 
 */
template<typename C>
auto rank ( C const& c ) 
{
    if constexpr ( is_ptr_or_shared_ptr<C>() )
        return c->worldComm().localRank();
    else
        return c.worldComm().localRank();
}

template<typename C>
FEELPP_DEPRECATED 
auto meshrank ( C const& c )
{
    if constexpr ( is_ptr_or_shared_ptr<C>() )
        return c->worldComm().globalRank();
    else
        return c.worldComm().globalRank();
}

/**
 * @return the global MPI rank of the data structure @p c 
 */
template<typename C>
auto globalRank ( C const& c )
{
    if constexpr ( is_ptr_or_shared_ptr<C>() )
        return c->worldComm().globalRank();
    else
        return c.worldComm().globalRank();
}

} // Feel
#endif
