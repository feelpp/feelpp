/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-27

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2010 Universite de Grenoble 1 (Joseph Fourier)

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
/**
   \file filters.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-27
 */
#ifndef __FEELPP_FILTERS_HPP
#define __FEELPP_FILTERS_HPP 1

#include <utility>

#include <feel/feelcore/environment.hpp>
#include <feel/feelmesh/traits.hpp>

namespace Feel
{
enum class EntityProcessType {LOCAL_ONLY,GHOST_ONLY,ALL,IGNORE_ENTITY_ON_INTERPROCESS_FACE};

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

enum ElementsType
{
    MESH_ELEMENTS = 0,           /**< elements */
    MESH_FACES = 1,              /**< faces */
    MESH_INTERNAL_FACES = 2,     /**< internal faces */
    MESH_EDGES = 3,              /**< edges */
    MESH_INTERNAL_EDGES = 4,     /**< internal edges */
    MESH_POINTS = 5              /**< points */
};

namespace meta
{

template<typename MeshType>
struct elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
            typename MeshTraits<MeshType>::element_const_iterator,
            typename MeshTraits<MeshType>::element_const_iterator> type;
};

template<typename MeshType>
struct markedelements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
            typename MeshTraits<MeshType>::marker_element_const_iterator,
            typename MeshTraits<MeshType>::marker_element_const_iterator> type;
};

template<typename MeshType>
struct marked2elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
            typename MeshTraits<MeshType>::marker2_element_const_iterator,
            typename MeshTraits<MeshType>::marker2_element_const_iterator> type;
};

template<typename MeshType>
struct marked3elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
            typename MeshTraits<MeshType>::marker3_element_const_iterator,
            typename MeshTraits<MeshType>::marker3_element_const_iterator> type;
};

} // meta

/// \cond detail
namespace detail
{


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::element_const_iterator,
      typename MeshTraits<MeshType>::element_const_iterator>
      allelements( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              mesh.beginElement(),
                              mesh.endElement() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::element_const_iterator,
      typename MeshTraits<MeshType>::element_const_iterator>
      allelements( MeshType const& mesh,mpl::bool_<true> )
{
    return allelements( *mesh, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::element_const_iterator,
      typename MeshTraits<MeshType>::element_const_iterator>
      elements( MeshType const& mesh, rank_type pid, mpl::bool_<false> )
{

    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              mesh.beginElementWithProcessId( pid ),
                              mesh.endElementWithProcessId( pid ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::element_const_iterator,
      typename MeshTraits<MeshType>::element_const_iterator>
      elements( MeshType const& mesh, rank_type pid, mpl::bool_<true> )
{
    return elements( *mesh, pid, mpl::bool_<false>() );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::location_element_const_iterator,
      typename MeshTraits<MeshType>::location_element_const_iterator>
boundaryelements( MeshType const& mesh, uint16_type entity_min_dim, uint16_type entity_max_dim, rank_type pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::location_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.boundaryElements( entity_min_dim, entity_max_dim, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::location_element_const_iterator,
      typename MeshTraits<MeshType>::location_element_const_iterator>
boundaryelements( MeshType const& mesh, uint16_type entity_min_dim, uint16_type entity_max_dim, rank_type pid, mpl::bool_<true> )
{
    return boundaryelements( *mesh, entity_min_dim, entity_max_dim, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::location_element_const_iterator,
      typename MeshTraits<MeshType>::location_element_const_iterator>
      internalelements( MeshType const& mesh, rank_type pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::location_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.internalElements( pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::location_element_const_iterator,
      typename MeshTraits<MeshType>::location_element_const_iterator>
      internalelements( MeshType const& mesh, rank_type pid, mpl::bool_<true> )
{
    return internalelements( *mesh, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker_element_const_iterator,
      typename MeshTraits<MeshType>::marker_element_const_iterator>
      markedelements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<false>  )
{
    typedef typename MeshTraits<MeshType>::marker_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.elementsWithMarker( flag, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker_element_const_iterator,
      typename MeshTraits<MeshType>::marker_element_const_iterator>
      markedelements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<true> )
{
    return markedelements( *mesh, flag, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker2_element_const_iterator,
      typename MeshTraits<MeshType>::marker2_element_const_iterator>
      marked2elements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker2_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.elementsWithMarker2( flag, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker2_element_const_iterator,
      typename MeshTraits<MeshType>::marker2_element_const_iterator>
      marked2elements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<true> )
{
    return marked2elements( *mesh, flag, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker3_element_const_iterator,
      typename MeshTraits<MeshType>::marker3_element_const_iterator>
      marked3elements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker3_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.elementsWithMarker3( flag, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker3_element_const_iterator,
      typename MeshTraits<MeshType>::marker3_element_const_iterator>
      marked3elements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<true> )
{
    return marked3elements( *mesh, flag, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::element_const_iterator,
      typename MeshTraits<MeshType>::element_const_iterator>
      idedelements( MeshType const& mesh, size_type id, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              mesh.beginElementWithId( id ),
                              mesh.endElementWithId( id ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::element_const_iterator,
      typename MeshTraits<MeshType>::element_const_iterator>
      idedelements( MeshType const& mesh, size_type id, mpl::bool_<true> )
{
    return idedelements( *mesh, id, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::pid_face_const_iterator,
      typename MeshTraits<MeshType>::pid_face_const_iterator>
      faces( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::pid_face_const_iterator pid_face_const_iterator;
    pid_face_const_iterator it,en;
    boost::tie( it, en ) = mesh.facesWithProcessId( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(), it, en );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::pid_face_const_iterator,
      typename MeshTraits<MeshType>::pid_face_const_iterator>
      faces( MeshType const& mesh, rank_type __pid, mpl::bool_<true> )
{
    return faces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::face_const_iterator,
      typename MeshTraits<MeshType>::face_const_iterator>
      idedfaces( MeshType const& mesh, size_type id, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              mesh.beginFaceWithId( id ),
                              mesh.endFaceWithId( id ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::face_const_iterator,
      typename MeshTraits<MeshType>::face_const_iterator>
      idedfaces( MeshType const& mesh, size_type id, mpl::bool_<true> )
{
    return idedfaces( *mesh, id, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker_face_const_iterator,
      typename MeshTraits<MeshType>::marker_face_const_iterator>
      markedfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker_face_const_iterator iterator;
    auto beg = mesh.beginFaceWithMarker();
    auto end = mesh.endFaceWithMarker();
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              beg, end );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker_face_const_iterator,
      typename MeshTraits<MeshType>::marker_face_const_iterator>
      markedfaces( MeshType const& mesh, flag_type __marker, rank_type __pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.facesWithMarker( __marker, __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker2_face_const_iterator,
      typename MeshTraits<MeshType>::marker2_face_const_iterator>
      marked2faces( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker2_face_const_iterator iterator;
    auto beg = mesh.beginFaceWithMarker2();
    auto end = mesh.endFaceWithMarker2();
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              beg, end );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker2_face_const_iterator,
      typename MeshTraits<MeshType>::marker2_face_const_iterator>
      marked2faces( MeshType const& mesh, flag_type __marker, rank_type __pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker2_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.facesWithMarker2( __marker, __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker3_face_const_iterator,
      typename MeshTraits<MeshType>::marker3_face_const_iterator>
      marked3faces( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker3_face_const_iterator iterator;
    auto beg = mesh.beginFaceWithMarker3();
    auto end = mesh.endFaceWithMarker3();
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              beg, end );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker3_face_const_iterator,
      typename MeshTraits<MeshType>::marker3_face_const_iterator>
      marked3faces( MeshType const& mesh, flag_type __marker, rank_type __pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker3_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.facesWithMarker3( __marker, __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker_face_const_iterator,
      typename MeshTraits<MeshType>::marker_face_const_iterator>
      markedfaces( MeshType const& mesh, flag_type __marker, rank_type __pid, mpl::bool_<true> )
{
    return markedfaces( *mesh,  __marker, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker_face_const_iterator,
      typename MeshTraits<MeshType>::marker_face_const_iterator>
      markedfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<true> )
{
    return markedfaces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker2_face_const_iterator,
      typename MeshTraits<MeshType>::marker2_face_const_iterator>
      marked2faces( MeshType const& mesh, flag_type __marker, rank_type __pid, mpl::bool_<true> )
{
    return marked2faces( *mesh,  __marker, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker2_face_const_iterator,
      typename MeshTraits<MeshType>::marker2_face_const_iterator>
      marked2faces( MeshType const& mesh, rank_type __pid, mpl::bool_<true> )
{
    return marked2faces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker3_face_const_iterator,
      typename MeshTraits<MeshType>::marker3_face_const_iterator>
      marked3faces( MeshType const& mesh, flag_type __marker, rank_type __pid, mpl::bool_<true> )
{
    return marked3faces( *mesh,  __marker, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker3_face_const_iterator,
      typename MeshTraits<MeshType>::marker3_face_const_iterator>
      marked3faces( MeshType const& mesh, rank_type __pid, mpl::bool_<true> )
{
    return marked3faces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::location_face_const_iterator,
      typename MeshTraits<MeshType>::location_face_const_iterator>
      boundaryfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<false>  )
{
    typedef typename MeshTraits<MeshType>::location_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.facesOnBoundary( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::location_face_const_iterator,
      typename MeshTraits<MeshType>::location_face_const_iterator>
      boundaryfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<true>  )
{
    return boundaryfaces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::location_face_const_iterator,
      typename MeshTraits<MeshType>::location_face_const_iterator>
      internalfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{

    typedef typename MeshTraits<MeshType>::location_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.internalFaces( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::location_face_const_iterator,
      typename MeshTraits<MeshType>::location_face_const_iterator>
      internalfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<true> )
{
    return internalfaces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::interprocess_face_const_iterator,
      typename MeshTraits<MeshType>::interprocess_face_const_iterator>
      interprocessfaces( MeshType const& mesh, rank_type neighbor_pid, mpl::bool_<false> )
{

    typedef typename MeshTraits<MeshType>::interprocess_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.interProcessFaces( neighbor_pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::interprocess_face_const_iterator,
      typename MeshTraits<MeshType>::interprocess_face_const_iterator>
      interprocessfaces( MeshType const& mesh, rank_type neighbor_pid, mpl::bool_<true> )
{
    return interprocessfaces( *mesh, neighbor_pid, mpl::bool_<false>() );

}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::pid_edge_const_iterator,
      typename MeshTraits<MeshType>::pid_edge_const_iterator>
      edges( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::pid_edge_const_iterator pid_edge_const_iterator;
    pid_edge_const_iterator it,en;
    boost::tie( it, en ) = mesh.edgesWithProcessId( __pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(), it, en );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::pid_edge_const_iterator,
      typename MeshTraits<MeshType>::pid_edge_const_iterator>
      edges( MeshType const& mesh, rank_type __pid, mpl::bool_<true> )
{
    return edges( *mesh, __pid, mpl::bool_<false>() );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::marker_edge_const_iterator,
      typename MeshTraits<MeshType>::marker_edge_const_iterator>
      markededges( MeshType const& mesh, flag_type __marker, rank_type __pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker_edge_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.edgesWithMarker( __marker, __pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::marker_edge_const_iterator,
      typename MeshTraits<MeshType>::marker_edge_const_iterator>
      markededges( MeshType const& mesh, flag_type __marker, rank_type __pid, mpl::bool_<true> )
{
    return markededges( *mesh, __marker, __pid, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::location_edge_const_iterator,
      typename MeshTraits<MeshType>::location_edge_const_iterator>
      boundaryedges( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              mesh.beginEdgeOnBoundary(),
                              mesh.endEdgeOnBoundary() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::location_edge_const_iterator,
      typename MeshTraits<MeshType>::location_edge_const_iterator>
      boundaryedges( MeshType const& mesh, mpl::bool_<true> )
{
    return boundaryedges( *mesh, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::location_edge_const_iterator,
      typename MeshTraits<MeshType>::location_edge_const_iterator>
      internaledges( MeshType const& mesh, mpl::bool_<true> )
{
    return internaledges( mesh, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::location_edge_const_iterator,
      typename MeshTraits<MeshType>::location_edge_const_iterator>
      internaledges( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              mesh.beginInternalEdge(),
                              mesh.endInternalEdge() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::point_const_iterator,
      typename MeshTraits<MeshType>::point_const_iterator>
      points( MeshType const& mesh, mpl::bool_<true> )
{
    return points( *mesh, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::point_const_iterator,
      typename MeshTraits<MeshType>::point_const_iterator>
      points( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginPoint(),
                              mesh.endPoint() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::marker_point_const_iterator,
      typename MeshTraits<MeshType>::marker_point_const_iterator>
      markedpoints( MeshType const& mesh, size_type flag, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginPointWithMarker( flag ),
                              mesh.endPointWithMarker( flag ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::marker_point_const_iterator,
      typename MeshTraits<MeshType>::marker_point_const_iterator>
      markedpoints( MeshType const& mesh, size_type flag, mpl::bool_<true> )
{
    return markedpoints( *mesh, flag, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::location_point_const_iterator,
      typename MeshTraits<MeshType>::location_point_const_iterator>
      boundarypoints( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginPointOnBoundary( ),
                              mesh.endPointOnBoundary() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::location_point_const_iterator,
      typename MeshTraits<MeshType>::location_point_const_iterator>
      boundarypoints( MeshType const& mesh, mpl::bool_<true> )
{
    return boundarypoints( *mesh, mpl::bool_<false>() );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::location_point_const_iterator,
      typename MeshTraits<MeshType>::location_point_const_iterator>
      internalpoints( MeshType const& mesh, mpl::bool_<true> )
{
    return internalpoints( *mesh, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::location_point_const_iterator,
      typename MeshTraits<MeshType>::location_point_const_iterator>
      internalpoints( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginInternalPoint( ),
                              mesh.endInternalPoint() );
}

} // detail
/// \endcond

template<typename MeshType>
rank_type meshrank ( MeshType const& mesh, mpl::bool_<true> )
{
    return mesh->worldComm().localRank();
}

template<typename MeshType>
rank_type meshrank ( MeshType const& mesh, mpl::bool_<false> )
{
    return mesh.worldComm().localRank();
}
/**
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::element_const_iterator,
      typename MeshTraits<MeshType>::element_const_iterator>
      allelements( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::allelements( mesh, is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::element_const_iterator,
      typename MeshTraits<MeshType>::element_const_iterator>
      elements( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::elements( mesh, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
    //return elements( mesh, flag, is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the mesh
 * which share a face with the boundary
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::location_element_const_iterator,
      typename MeshTraits<MeshType>::location_element_const_iterator>
boundaryelements( MeshType const& mesh, uint16_type entity_min_dim = 0, uint16_type entity_max_dim = 2 )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::boundaryelements( mesh, entity_min_dim, entity_max_dim, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the mesh
 * which are stricly within the domain that is to say they do not
 * share a face with the boundary
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::location_element_const_iterator,
      typename MeshTraits<MeshType>::location_element_const_iterator>
      internalelements( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::internalelements( mesh, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::marker_element_const_iterator,
             typename MeshTraits<MeshType>::marker_element_const_iterator>
markedelements( MeshType const& mesh, std::string const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedelements( mesh, mesh->markerName( flag ), meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::marker_element_const_iterator,
             typename MeshTraits<MeshType>::marker_element_const_iterator>
markedelements( MeshType const& mesh, const char* flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedelements( mesh, mesh->markerName( flag ), meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::marker_element_const_iterator,
             typename MeshTraits<MeshType>::marker_element_const_iterator>
markedelements( MeshType const& mesh, boost::any const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;

    flag_type theflag = mesh->markerId( flag );
    VLOG(2) << "[markedelements] flag: " << theflag << "\n";
    return Feel::detail::markedelements( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                       typename MeshTraits<MeshType>::marker_element_const_iterator,
                       typename MeshTraits<MeshType>::marker_element_const_iterator> >
markedelements( MeshType const& mesh, std::initializer_list<boost::any> const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                           typename MeshTraits<MeshType>::marker_element_const_iterator,
                           typename MeshTraits<MeshType>::marker_element_const_iterator> > list_elements;
    for ( auto const& it : flag )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedelements] flag: " << theflag << "\n";
        list_elements.push_back( Feel::detail::markedelements( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() ) );
    }
    return list_elements;
}

template<typename MeshType>
std::list<boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                       typename MeshTraits<MeshType>::marker_element_const_iterator,
                       typename MeshTraits<MeshType>::marker_element_const_iterator> >
markedelements( MeshType const& mesh, std::list<std::string> const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                           typename MeshTraits<MeshType>::marker_element_const_iterator,
                           typename MeshTraits<MeshType>::marker_element_const_iterator> > list_elements;
    for ( auto const& it : flag )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedelements] flag: " << theflag << "\n";
        list_elements.push_back( Feel::detail::markedelements( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() ) );
    }
    return list_elements;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker2 \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker2_element_const_iterator,
      typename MeshTraits<MeshType>::marker2_element_const_iterator>
      marked2elements( MeshType const& mesh, flag_type flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked2elements( mesh, flag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker2 string
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker2_element_const_iterator,
      typename MeshTraits<MeshType>::marker2_element_const_iterator>
      marked2elements( MeshType const& mesh, std::string const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked2elements( mesh, mesh->markerName( flag ), meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}
template<typename MeshType>
std::list<boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                       typename MeshTraits<MeshType>::marker2_element_const_iterator,
                       typename MeshTraits<MeshType>::marker2_element_const_iterator> >
marked2elements( MeshType const& mesh, std::initializer_list<boost::any> const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                           typename MeshTraits<MeshType>::marker2_element_const_iterator,
                           typename MeshTraits<MeshType>::marker2_element_const_iterator> > list_elements;
    for ( auto const& it : flag )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedelements] flag: " << theflag << "\n";
        list_elements.push_back( Feel::detail::marked2elements( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() ) );
    }
    return list_elements;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker3 \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker3_element_const_iterator,
      typename MeshTraits<MeshType>::marker3_element_const_iterator>
      marked3elements( MeshType const& mesh, flag_type flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked3elements( mesh, flag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker3 string
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::marker3_element_const_iterator,
      typename MeshTraits<MeshType>::marker3_element_const_iterator>
      marked3elements( MeshType const& mesh, std::string const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked3elements( mesh, mesh->markerName( flag ), meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                       typename MeshTraits<MeshType>::marker3_element_const_iterator,
                       typename MeshTraits<MeshType>::marker3_element_const_iterator> >
marked3elements( MeshType const& mesh, std::initializer_list<boost::any> const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                           typename MeshTraits<MeshType>::marker3_element_const_iterator,
                           typename MeshTraits<MeshType>::marker3_element_const_iterator> > list_elements;
    for ( auto const& it : flag )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedelements] flag: " << theflag << "\n";
        list_elements.push_back( Feel::detail::marked3elements( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() ) );
    }
    return list_elements;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with id \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename MeshTraits<MeshType>::element_const_iterator,
      typename MeshTraits<MeshType>::element_const_iterator>
      idedelements( MeshType const& mesh, flag_type flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::idedelements( mesh, flag, is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over faces of the
 * mesh on processor \c __pid
 *
 * @param mesh a mesh data structure
 * @param __pid process id
 *
 * @return a pair of face iterators (begin,end)
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::pid_face_const_iterator,
      typename MeshTraits<MeshType>::pid_face_const_iterator>
      faces( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::faces( mesh, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with id \p id
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::face_const_iterator,
      typename MeshTraits<MeshType>::face_const_iterator>
      idedfaces( MeshType const& mesh, size_type id )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::idedfaces( mesh, id, is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over faces of the
 * mesh marked
 *
 * @param mesh a mesh data structure
 * @param __pid process id
 *
 * @return a pair of iterators (begin,end) for the set of marked faces
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker_face_const_iterator,
      typename MeshTraits<MeshType>::marker_face_const_iterator>
      markedfaces( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedfaces( mesh, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over faces of the
 * mesh marked with \c __marker
 *
 * @param mesh a mesh data structure
 * @param __marker a marker that identifies faces
 * @param __pid process id
 *
 * @return a pair of iterators (begin,end) for the set of marked faces
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker_face_const_iterator,
      typename MeshTraits<MeshType>::marker_face_const_iterator>
      markedfaces( MeshType const& mesh,
                   std::string const&__marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedfaces( mesh, mesh->markerName( __marker ), meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker_face_const_iterator,
      typename MeshTraits<MeshType>::marker_face_const_iterator>
      markedfaces( MeshType const& mesh,
                   const char*__marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedfaces( mesh, mesh->markerName( __marker ), meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}
template<typename MeshType>
std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                       typename MeshTraits<MeshType>::marker_face_const_iterator,
                       typename MeshTraits<MeshType>::marker_face_const_iterator> >
markedfaces( MeshType const& mesh,
             std::list<std::string> const& __markers )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                           typename MeshTraits<MeshType>::marker_face_const_iterator,
                           typename MeshTraits<MeshType>::marker_face_const_iterator> > list_faces;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_faces.push_back( Feel::detail::markedfaces( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() ) );
    }
    return list_faces;
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker_face_const_iterator,
      typename MeshTraits<MeshType>::marker_face_const_iterator>
      markedfaces( MeshType const& mesh,
                   boost::any const& __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    flag_type theflag = mesh->markerId( __marker );
    VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
    return Feel::detail::markedfaces( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );

}

template<typename MeshType>
std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                       typename MeshTraits<MeshType>::marker_face_const_iterator,
                       typename MeshTraits<MeshType>::marker_face_const_iterator> >
markedfaces( MeshType const& mesh,
             std::initializer_list<boost::any> __markers )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                           typename MeshTraits<MeshType>::marker_face_const_iterator,
                           typename MeshTraits<MeshType>::marker_face_const_iterator> > list_faces;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_faces.push_back( Feel::detail::markedfaces( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() ) );
    }
    return list_faces;
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker2_face_const_iterator,
      typename MeshTraits<MeshType>::marker2_face_const_iterator>
      marked2faces( MeshType const& mesh,
                    flag_type __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked2faces( mesh, __marker, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                       typename MeshTraits<MeshType>::marker2_face_const_iterator,
                       typename MeshTraits<MeshType>::marker2_face_const_iterator> >
marked2faces( MeshType const& mesh,
             std::initializer_list<boost::any> __markers )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                           typename MeshTraits<MeshType>::marker2_face_const_iterator,
                           typename MeshTraits<MeshType>::marker2_face_const_iterator> > list_faces;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_faces.push_back( Feel::detail::marked2faces( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() ) );
    }
    return list_faces;
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker3_face_const_iterator,
      typename MeshTraits<MeshType>::marker3_face_const_iterator>
      marked3faces( MeshType const& mesh,
                    flag_type __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked3faces( mesh, __marker, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker2_face_const_iterator,
      typename MeshTraits<MeshType>::marker2_face_const_iterator>
      marked2faces( MeshType const& mesh,
                    std::string const& __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked2faces( mesh, mesh->markerName( __marker ), meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker3_face_const_iterator,
      typename MeshTraits<MeshType>::marker3_face_const_iterator>
      marked3faces( MeshType const& mesh,
                    std::string const& __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked3faces( mesh, mesh->markerName( __marker ), meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                       typename MeshTraits<MeshType>::marker3_face_const_iterator,
                       typename MeshTraits<MeshType>::marker3_face_const_iterator> >
marked3faces( MeshType const& mesh,
             std::initializer_list<boost::any> __markers )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                           typename MeshTraits<MeshType>::marker3_face_const_iterator,
                           typename MeshTraits<MeshType>::marker3_face_const_iterator> > list_faces;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_faces.push_back( Feel::detail::marked3faces( mesh, theflag, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() ) );
    }
    return list_faces;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all boundary faces of the
 * mesh
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::location_face_const_iterator,
      typename MeshTraits<MeshType>::location_face_const_iterator>
      boundaryfaces( MeshType const& mesh  )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::boundaryfaces( mesh, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all internal faces of the
 * mesh belong to process domain \p __pid
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::location_face_const_iterator,
      typename MeshTraits<MeshType>::location_face_const_iterator>
      internalfaces( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::internalfaces( mesh, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all interprocess faces of the
 * mesh belonging to process \p __pid
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::interprocess_face_const_iterator,
      typename MeshTraits<MeshType>::interprocess_face_const_iterator>
      interprocessfaces( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::interprocessfaces( mesh, invalid_rank_type_value, is_ptr_or_shared_ptr() );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all interprocess faces of the
 * mesh belonging to process \p __pid
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::interprocess_face_const_iterator,
      typename MeshTraits<MeshType>::interprocess_face_const_iterator>
interprocessfaces( MeshType const& mesh, rank_type neighbor_pid )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::interprocessfaces( mesh, neighbor_pid, is_ptr_or_shared_ptr() );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over edges of the
 * mesh on processor \c __pid
 *
 * @param mesh a mesh data structure
 * @param __pid process id
 *
 * @return a pair of edge iterators (begin,end)
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::pid_edge_const_iterator,
      typename MeshTraits<MeshType>::pid_edge_const_iterator>
      edges( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::edges( mesh, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over edges of the
 * mesh marked with \c __marker
 *
 * @param mesh a mesh data structure
 * @param __marker a marker that identifies edges
 * @param __pid process id
 *
 * @return a pair of iterators (begin,end) for the set of marked edges
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::marker_edge_const_iterator,
      typename MeshTraits<MeshType>::marker_edge_const_iterator>
      markededges( MeshType const& mesh,
                   flag_type __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markededges( mesh, __marker, meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::marker_edge_const_iterator,
      typename MeshTraits<MeshType>::marker_edge_const_iterator>
      markededges( MeshType const& mesh,
                   std::string const& __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markededges( mesh, mesh->markerName(__marker), meshrank( mesh, is_ptr_or_shared_ptr() ), is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all boundary edges of the
 * mesh
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::location_edge_const_iterator,
      typename MeshTraits<MeshType>::location_edge_const_iterator>
      boundaryedges( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::boundaryedges( mesh, is_ptr_or_shared_ptr() );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all internal edges of the
 * mesh
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
      typename MeshTraits<MeshType>::location_edge_const_iterator,
      typename MeshTraits<MeshType>::location_edge_const_iterator>
      internaledges( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::internaledges( mesh, is_ptr_or_shared_ptr() );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::point_const_iterator,
      typename MeshTraits<MeshType>::point_const_iterator>
      points( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::points( mesh, is_ptr_or_shared_ptr() );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the marked points with \p flag of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::marker_point_const_iterator,
      typename MeshTraits<MeshType>::marker_point_const_iterator>
      markedpoints( MeshType const& mesh, size_type flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedpoints( mesh, flag, is_ptr_or_shared_ptr() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::marker_point_const_iterator,
             typename MeshTraits<MeshType>::marker_point_const_iterator>
markedpoints( MeshType const& mesh, std::string const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedpoints( mesh, mesh->markerName(flag), is_ptr_or_shared_ptr() );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the boundary points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::location_point_const_iterator,
      typename MeshTraits<MeshType>::location_point_const_iterator>
      boundarypoints( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::boundarypoints( mesh, is_ptr_or_shared_ptr() );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the internal(not on the boundary) points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::location_point_const_iterator,
      typename MeshTraits<MeshType>::location_point_const_iterator>
      internalpoints( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;

    return Feel::detail::internalpoints( mesh, is_ptr_or_shared_ptr() );
}

/**
 * \ingroup MeshIterators
 * \return the number of elements given element iterators constructed
 * using the mesh filters
 * \param its the mesh iterators
 * \param global all reduce number of elements if set to true
 *
 * The following code prints in the logfile the number of elements in
 * the mesh that are marked with marker1 equal to 1:
 *
 * \code
 * LOG(INFO) << "number of elements = " << nelements( markedelements(mesh,1) ) << "\n";
 * \endcode
 *
 */
template<typename MT, typename Iterator>
size_type
nelements( boost::tuple<MT,Iterator,Iterator> const& its, bool global = false )
{
    size_type d = std::distance( boost::get<1>( its ), boost::get<2>( its ) );
    size_type gd = d;
    if ( global )
        mpi::all_reduce(Environment::worldComm().globalComm(),
                        d,
                        gd,
                        std::plus<size_type>());
    return gd;

}
/**
 * \ingroup MeshIterators
 * \return the number of elements given element iterators constructed
 * using the mesh filters
 * \param a list of mesh iterators
 * \param global all reduce number of elements if set to true
 *
 * The following code prints in the logfile the number of elements in
 * the mesh that are marked with marker1 equal to 1:
 *
 * \code
 * LOG(INFO) << "number of elements = " << nelements( markedelements(mesh,{1,2,3}) ) << "\n";
 * \endcode
 *
 */
template<typename MT, typename Iterator>
size_type
nelements( std::list<boost::tuple<MT,Iterator,Iterator> > const& its, bool global = false )
{
    size_type d = 0;
    std::for_each( its.begin(), its.end(),
                   [&d]( boost::tuple<MT,Iterator,Iterator> const& t )
                   {
                       d+=std::distance( boost::get<1>( t ), boost::get<2>( t ) );
                   } );
    size_type gd = d;
    if ( global )
        mpi::all_reduce(Environment::worldComm().globalComm(),
                        d,
                        gd,
                        std::plus<size_type>());
    return gd;
}


/**
 * \ingroup MeshIterators
 * \return the number of elements given element iterators constructed
 * using custom range
 * \param its the mesh iterators
 * \param global all reduce number of elements if set to true
 *
 * The following code prints in the logfile the number of elements in
 * the mesh that are marked with marker1 equal to 1:
 *
 * \code
 * LOG(INFO) << "number of elements = " << nelements( myCustomRange ) << "\n";
 * \endcode
 *
 */
template<typename MT, typename Iterator,typename Container>
size_type
nelements( boost::tuple<MT,Iterator,Iterator,Container> const& its, bool global = false )
{
    size_type d = std::distance( boost::get<1>( its ), boost::get<2>( its ) );
    size_type gd = d;
    if ( global )
        mpi::all_reduce(Environment::worldComm().globalComm(),
                        d,
                        gd,
                        std::plus<size_type>());
    return gd;

}



template<typename ElementType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
      typename std::list<ElementType>::const_iterator,
      typename std::list<ElementType>::const_iterator>
      element( ElementType const& elt  )
{
    std::list<ElementType> lst;
    lst.push_back(  elt );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              lst.begin(),
                              lst.end() );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> >::const_iterator,
             typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> >::const_iterator,
             boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > >
             >
elements( MeshType const& mesh, EntityProcessType entity )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : elements(mesh) )
        {
            myelts->push_back(boost::cref(elt));
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto face_it = mesh->interProcessFaces().first;
        auto const face_en = mesh->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;

            // add elt in range
            myelts->push_back(boost::cref(eltOffProc));

            eltGhostDone.insert( eltOffProc.id() );
        }
    }

    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> >::const_iterator,
             typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> >::const_iterator,
             boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > >
             >
markedelements( MeshType const& mesh, boost::any const& flag, EntityProcessType entity )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    flag_type theflag = mesh->markerId( flag );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : markedelements(mesh,flag) )
        {
            myelts->push_back(boost::cref(elt));
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto face_it = mesh->interProcessFaces().first;
        auto const face_en = mesh->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;


            // add elt in range
            if ( eltOffProc.marker().value() == theflag )
                myelts->push_back(boost::cref(eltOffProc));

            eltGhostDone.insert( eltOffProc.id() );
        }
    }

    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> >::const_iterator,
             typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> >::const_iterator,
             boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> > >
             >
boundaryfaces( MeshType const& mesh, EntityProcessType entity )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& theface : boundaryfaces(mesh) )
        {
            myelts->push_back(boost::cref(theface));
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        //std::set<size_type> faceGhostDone;
        auto face_it = mesh->interProcessFaces().first;
        auto const face_en = mesh->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            for ( size_type f = 0; f < mesh->numLocalFaces(); f++ )
            {
                auto const& theface = eltOffProc.face(f);
                if ( theface.isOnBoundary() ) //&& faceGhostDone.find( theface.id() ) == faceGhostDone.end() )
                    myelts->push_back(boost::cref(theface));
            }
        }
    }

    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts );

}

/**
 * this function takes two sets of iterators and generate a data structure that
 * will hold the corresponding reference of mesh elenents. This data structure
 * can then be used to iterate over meshes like other mesh filters.
 *
 * \param it1 (begin,end) iterators to concatenate \param it2 (begin,end)
 * iterators to concatenate \return a data structure that holds the merge
 * between two sets of pair of iterators
 */
template<typename IteratorType>
boost::tuple<mpl::size_t<boost::tuples::element<0,IteratorType>::type::value>,
             typename std::vector<boost::reference_wrapper<typename boost::tuples::element<1,IteratorType>::type::value_type const> >::const_iterator,
             typename std::vector<boost::reference_wrapper<typename boost::tuples::element<1,IteratorType>::type::value_type const> >::const_iterator,
             boost::shared_ptr<std::vector<boost::reference_wrapper<typename boost::tuples::element<1,IteratorType>::type::value_type const> > >
             >
concatenate( IteratorType it1, IteratorType it2 )
{
    typedef std::vector<boost::reference_wrapper<typename boost::tuples::element<1,IteratorType>::type::value_type  const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    auto append = [&myelts]( typename boost::tuples::element<1,IteratorType>::type::value_type  const& e ) { myelts->push_back( boost::cref(e) ); };
    std::for_each( begin( it1 ), end( it1 ), append );
    std::for_each( begin( it2 ), end( it2 ), append );
    
    return boost::make_tuple( mpl::size_t<boost::tuples::element<0,IteratorType>::type::value>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts );

}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::edge_type const> >::const_iterator,
             typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::edge_type const> >::const_iterator,
             boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::edge_type const> > >
             >
interprocessedges( MeshType const& mesh, rank_type neighbor_pid = invalid_rank_type_value, EntityProcessType entity = EntityProcessType::ALL )
{
    typedef typename MeshTraits<MeshType>::face_type face_type;
    typedef typename MeshTraits<MeshType>::edge_type edge_type;
    typedef std::vector<boost::reference_wrapper<edge_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myedges( new cont_range_type );

    if ( entity == EntityProcessType::ALL )
    {
        if ( neighbor_pid == invalid_rank_type_value )
        {
            for ( auto const& theedge : edges(mesh) )
            {
                if ( theedge.numberOfProcGhost() > 0 )
                    myedges->push_back(boost::cref(theedge));
            }
        }
        else
        {
            for ( auto const& theedge : edges(mesh) )
            {
                if ( theedge.numberOfProcGhost() > 0 && theedge.isSubEntityOfGhostElement( neighbor_pid ) )
                    myedges->push_back(boost::cref(theedge));
            }
        }
    }
    else if (  entity == EntityProcessType::IGNORE_ENTITY_ON_INTERPROCESS_FACE )
    {
        std::set<size_type> saveIdEdges;
        for ( auto const& theface : interprocessfaces(mesh,neighbor_pid) )
        {
            for ( uint16_type e=0 ; e<face_type::numLocalEdges ; ++e )
                saveIdEdges.insert( theface.face(e).id() );
        }

        if ( neighbor_pid == invalid_rank_type_value )
        {
            for ( auto const& theedge : edges(mesh) )
            {
                if ( theedge.numberOfProcGhost() > 0 && saveIdEdges.find(theedge.id()) == saveIdEdges.end() )
                    myedges->push_back(boost::cref(theedge));
            }
        }
        else
        {
            for ( auto const& theedge : edges(mesh) )
            {
                if ( theedge.numberOfProcGhost() > 0 && theedge.isSubEntityOfGhostElement( neighbor_pid ) &&
                     saveIdEdges.find(theedge.id()) == saveIdEdges.end() )
                    myedges->push_back(boost::cref(theedge));
            }
        }
    }

    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              myedges->begin(),
                              myedges->end(),
                              myedges );

}

} // namespace Feel


#endif
