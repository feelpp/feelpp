/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 25 Apr 2017

 Copyright (C) 2017 Feel++ Consortium

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
#ifndef FEELPP_INTERSECT_HPP
#define FEELPP_INTERSECT_HPP 1

#include <utility>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelmesh/detail/intersect.hpp>



namespace Feel {

template<typename IteratorType, typename ...Args>
using intersect_t = std::remove_const<decay_type<IteratorType>>;

/**
 * this function takes a set of iterators over geometrical entities of the same
 * type (e.g. faces, volumes, edges, points) and generate a data structure that
 * will hold the corresponding reference of mesh elements. This data structure
 * can then be used to iterate over meshes like other mesh filters.
 *
 * \code
 * // the following code below intersection of internal elements with marked elements omega
 * auto subset = intersect( internalelements(mesh), markedelements(mesh,"omega") );
 * // it can then be used eg to create a submesh
 * auto submesh = createSubmesh( mesh, subset );
 * \endcode
 */
template<typename RangeType, typename... Args>
FEELPP_EXPORT auto
intersect( RangeType&& it, Args&&... args )
{
    return Feel::detail::intersect_impl( std::forward<RangeType>(it), std::forward<Args>(args)... );
}





}
#endif
