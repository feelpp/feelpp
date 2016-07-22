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
#ifndef FEELPP_FEELMESH_ITERATOR_HPP
#define FEELPP_FEELMESH_ITERATOR_HPP 1

namespace Feel {

enum class EntityProcessType {LOCAL_ONLY,GHOST_ONLY,ALL,IGNORE_ENTITY_ON_INTERPROCESS_FACE};
using entity_process_t = EntityProcessType;

template<size_t S, class ITERATOR>
ITERATOR begin( boost::tuple<mpl::size_t<S>,ITERATOR,ITERATOR> &range )
{
    return range.template get<1>();
}

template<size_t S, class ITERATOR>
ITERATOR end( boost::tuple<mpl::size_t<S>,ITERATOR,ITERATOR> &range )
{
    return range.template get<2>();
}
template<size_t S, class ITERATOR>
ITERATOR begin( boost::tuple<mpl::size_t<S>,ITERATOR,ITERATOR> const& range )
{
    return range.template get<1>();
}

template<size_t S, class ITERATOR>
ITERATOR end( boost::tuple<mpl::size_t<S>,ITERATOR,ITERATOR> const& range )
{
    return range.template get<2>();
}
template<size_t S, class ITERATOR, class CONTAINER>
ITERATOR begin( boost::tuple<mpl::size_t<S>,ITERATOR,ITERATOR,CONTAINER> &range )
{
    return range.template get<1>();
}

template<size_t S, class ITERATOR, class CONTAINER>
ITERATOR end( boost::tuple<mpl::size_t<S>,ITERATOR,ITERATOR,CONTAINER> &range )
{
    return range.template get<2>();
}

template<size_t S, class ITERATOR, class CONTAINER>
ITERATOR begin( boost::tuple<mpl::size_t<S>,ITERATOR,ITERATOR,CONTAINER> const& range )
{
    return range.template get<1>();
}

template<size_t S, class ITERATOR, class CONTAINER>
ITERATOR end( boost::tuple<mpl::size_t<S>,ITERATOR,ITERATOR,CONTAINER> const& range )
{
    return range.template get<2>();
}

/**
 * entity identifier  for mesh iterators
 */
enum ElementsType
{
    MESH_ELEMENTS = 0,           /**< elements */
    MESH_FACES = 1,              /**< faces */
    MESH_INTERNAL_FACES = 2,     /**< internal faces */
    MESH_EDGES = 3,              /**< edges */
    MESH_INTERNAL_EDGES = 4,     /**< internal edges */
    MESH_POINTS = 5              /**< points */
};


}
#endif
