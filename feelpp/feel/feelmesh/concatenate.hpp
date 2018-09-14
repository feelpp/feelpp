/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 21 Apr 2016

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
#ifndef FEELPP_CONCATENATE_HPP
#define FEELPP_CONCATENATE_HPP 1

#include <utility>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelmesh/detail/concatenate.hpp>



namespace Feel {

template<typename IteratorType, typename ...Args>
using concatenate_t = ext_entities_from_iterator_t<IteratorType>;


/**
 * this function takes a set of iterators over geometrical entities of the same
 * type (e.g. faces, volumes, edges, points) and generate a data structure that
 * will hold the corresponding reference of mesh elements. This data structure
 * can then be used to iterate over meshes like other mesh filters.
 *
 * \code
 * // the following code below merges internal faces with inlet faces
 * auto subset = concatenate( internalfaces(mesh), markedfaces(mesh,"inlet") );
 * // it can then be used eg to create a submesh
 * auto submesh = createSubmesh( mesh, subset );
 * \endcode
 */
template<typename IteratorType, typename... Args>
FEELPP_EXPORT Feel::concatenate_t<IteratorType>
concatenate( IteratorType&& it, Args&&... args )
{
    return Feel::detail::concatenate_impl( std::forward<IteratorType>(it), std::forward<Args>(args)... );
}





}
#endif
