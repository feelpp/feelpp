/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-08-27

  Copyright (C) 2005,2006 EPFL

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
   \file filters.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-08-27
 */
#ifndef __LIFE_FILTERS_HPP
#define __LIFE_FILTERS_HPP 1

#include <utility>

#include <life/lifecore/application.hpp>
#include <life/lifemesh/traits.hpp>

namespace Life
{

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
                         typename mesh::Traits<MeshType>::element_const_iterator,
                         typename mesh::Traits<MeshType>::element_const_iterator> type;
};

template<typename MeshType>
struct markedelements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                         typename mesh::Traits<MeshType>::marker_element_const_iterator,
                         typename mesh::Traits<MeshType>::marker_element_const_iterator> type;
};

template<typename MeshType>
struct marked2elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                         typename mesh::Traits<MeshType>::marker2_element_const_iterator,
                         typename mesh::Traits<MeshType>::marker2_element_const_iterator> type;
};

template<typename MeshType>
struct marked3elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                         typename mesh::Traits<MeshType>::marker3_element_const_iterator,
                         typename mesh::Traits<MeshType>::marker3_element_const_iterator> type;
};

} // meta

/// \cond detail
namespace detail
{

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::element_const_iterator,
             typename mesh::Traits<MeshType>::element_const_iterator>
elements( MeshType const& mesh, flag_type flag, mpl::bool_<false> )
{

    Debug(4000) << "[filters] elements on proc " << Application::processId() << " flag = " << flag << "\n";
    Debug(4000) << "[filters] elements begin id : " << mesh.beginElementWithProcessId( flag )->id() << "\n";
    Debug(4000) << "[filters] elements n : " << std::distance( mesh.beginElementWithProcessId( flag ),
                                                               mesh.endElementWithProcessId( flag ) )
                << "\n";
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              mesh.beginElementWithProcessId( flag ),
                              mesh.endElementWithProcessId( flag ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::element_const_iterator,
             typename mesh::Traits<MeshType>::element_const_iterator>
elements( MeshType const& mesh, flag_type flag, mpl::bool_<true> )
{
    return elements( *mesh, flag, mpl::bool_<false>() );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::location_element_const_iterator,
             typename mesh::Traits<MeshType>::location_element_const_iterator>
boundaryelements( MeshType const& mesh, size_type pid, mpl::bool_<false> )
{
    typedef typename mesh::Traits<MeshType>::location_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.boundaryElements( pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::location_element_const_iterator,
             typename mesh::Traits<MeshType>::location_element_const_iterator>
boundaryelements( MeshType const& mesh, size_type pid, mpl::bool_<true> )
{
    return boundaryelements( *mesh, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::location_element_const_iterator,
             typename mesh::Traits<MeshType>::location_element_const_iterator>
internalelements( MeshType const& mesh, size_type pid, mpl::bool_<false> )
{
    typedef typename mesh::Traits<MeshType>::location_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.internalElements( pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::location_element_const_iterator,
             typename mesh::Traits<MeshType>::location_element_const_iterator>
internalelements( MeshType const& mesh, size_type pid, mpl::bool_<true> )
{
    return internalelements( *mesh, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::marker_element_const_iterator,
             typename mesh::Traits<MeshType>::marker_element_const_iterator>
markedelements( MeshType const& mesh, flag_type flag, size_type pid, mpl::bool_<false>  )
{
    typedef typename mesh::Traits<MeshType>::marker_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.elementsWithMarker( flag, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::marker_element_const_iterator,
             typename mesh::Traits<MeshType>::marker_element_const_iterator>
markedelements( MeshType const& mesh, flag_type flag, size_type pid, mpl::bool_<true> )
{
    return markedelements( *mesh, flag, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::marker2_element_const_iterator,
             typename mesh::Traits<MeshType>::marker2_element_const_iterator>
marked2elements( MeshType const& mesh, flag_type flag, size_type pid, mpl::bool_<false> )
{
    typedef typename mesh::Traits<MeshType>::marker2_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.elementsWithMarker2( flag, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::marker2_element_const_iterator,
             typename mesh::Traits<MeshType>::marker2_element_const_iterator>
marked2elements( MeshType const& mesh, flag_type flag, size_type pid, mpl::bool_<true> )
{
    return marked2elements( *mesh, flag, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::marker3_element_const_iterator,
             typename mesh::Traits<MeshType>::marker3_element_const_iterator>
marked3elements( MeshType const& mesh, flag_type flag, size_type pid, mpl::bool_<false> )
{
    typedef typename mesh::Traits<MeshType>::marker3_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.elementsWithMarker3( flag, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::marker3_element_const_iterator,
             typename mesh::Traits<MeshType>::marker3_element_const_iterator>
marked3elements( MeshType const& mesh, flag_type flag, size_type pid, mpl::bool_<true> )
{
    return marked3elements( mesh, flag, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::element_const_iterator,
             typename mesh::Traits<MeshType>::element_const_iterator>
idedelements( MeshType const& mesh, flag_type flag, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              mesh.beginElementWithId( flag ),
                              mesh.endElementWithId( flag ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::element_const_iterator,
             typename mesh::Traits<MeshType>::element_const_iterator>
idedelements( MeshType const& mesh, flag_type flag, mpl::bool_<true> )
{
    return idedelements( *mesh, flag, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::pid_face_const_iterator,
             typename mesh::Traits<MeshType>::pid_face_const_iterator>
faces( MeshType const& mesh, size_type __pid, mpl::bool_<false> )
{
    typedef typename mesh::Traits<MeshType>::pid_face_const_iterator pid_face_const_iterator;
    pid_face_const_iterator it,en;
    boost::tie( it, en ) = mesh.facesWithProcessId( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(), it, en );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::pid_face_const_iterator,
             typename mesh::Traits<MeshType>::pid_face_const_iterator>
faces( MeshType const& mesh, size_type __pid, mpl::bool_<true> )
{
    return faces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::face_const_iterator,
             typename mesh::Traits<MeshType>::face_const_iterator>
idedfaces( MeshType const& mesh, size_type id, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              mesh.beginFaceWithId( id ),
                              mesh.endFaceWithId( id ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::face_const_iterator,
             typename mesh::Traits<MeshType>::face_const_iterator>
idedfaces( MeshType const& mesh, size_type id, mpl::bool_<true> )
{
    return idedfaces( *mesh, id, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::marker_face_const_iterator,
             typename mesh::Traits<MeshType>::marker_face_const_iterator>
markedfaces( MeshType const& mesh, flag_type __marker, size_type __pid, mpl::bool_<false> )
{
    typedef typename mesh::Traits<MeshType>::marker_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.facesWithMarker( __marker, __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::marker_face_const_iterator,
             typename mesh::Traits<MeshType>::marker_face_const_iterator>
markedfaces( MeshType const& mesh, flag_type __marker, size_type __pid, mpl::bool_<true> )
{
    return markedfaces( *mesh,  __marker, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::location_face_const_iterator,
             typename mesh::Traits<MeshType>::location_face_const_iterator>
boundaryfaces( MeshType const& mesh, size_type __pid, mpl::bool_<false>  )
{
    typedef typename mesh::Traits<MeshType>::location_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.facesOnBoundary( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::location_face_const_iterator,
             typename mesh::Traits<MeshType>::location_face_const_iterator>
boundaryfaces( MeshType const& mesh, size_type __pid, mpl::bool_<true>  )
{
    return boundaryfaces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::location_face_const_iterator,
             typename mesh::Traits<MeshType>::location_face_const_iterator>
internalfaces( MeshType const& mesh, size_type __pid, mpl::bool_<false> )
{

    typedef typename mesh::Traits<MeshType>::location_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.internalFaces( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::location_face_const_iterator,
             typename mesh::Traits<MeshType>::location_face_const_iterator>
internalfaces( MeshType const& mesh, size_type __pid, mpl::bool_<true> )
{
    return internalfaces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::interprocess_face_const_iterator,
             typename mesh::Traits<MeshType>::interprocess_face_const_iterator>
interprocessfaces( MeshType const& mesh, size_type __pid, mpl::bool_<false> )
{

    typedef typename mesh::Traits<MeshType>::interprocess_face_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.interProcessFaces( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::interprocess_face_const_iterator,
             typename mesh::Traits<MeshType>::interprocess_face_const_iterator>
interprocessfaces( MeshType const& mesh, size_type __pid, mpl::bool_<true> )
{
    return interprocessfaces( *mesh, __pid, mpl::bool_<false>() );

}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::pid_edge_const_iterator,
             typename mesh::Traits<MeshType>::pid_edge_const_iterator>
edges( MeshType const& mesh, size_type __pid, mpl::bool_<false> )
{
    typedef typename mesh::Traits<MeshType>::pid_edge_const_iterator pid_edge_const_iterator;
    pid_edge_const_iterator it,en;
    boost::tie( it, en ) = mesh.edgesWithProcessId( __pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(), it, en );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::pid_edge_const_iterator,
             typename mesh::Traits<MeshType>::pid_edge_const_iterator>
edges( MeshType const& mesh, size_type __pid, mpl::bool_<true> )
{
    return edges( *mesh, __pid, mpl::bool_<false>() );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::marker_edge_const_iterator,
             typename mesh::Traits<MeshType>::marker_edge_const_iterator>
markededges( MeshType const& mesh, flag_type __marker, size_type __pid, mpl::bool_<false> )
{
    typedef typename mesh::Traits<MeshType>::marker_edge_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.edgesWithMarker( __marker, __pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::marker_edge_const_iterator,
             typename mesh::Traits<MeshType>::marker_edge_const_iterator>
markededges( MeshType const& mesh, flag_type __marker, size_type __pid, mpl::bool_<true> )
{
    return markededges( *mesh, __marker, __pid, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::location_edge_const_iterator,
             typename mesh::Traits<MeshType>::location_edge_const_iterator>
boundaryedges( MeshType const& mesh, mpl::bool_<true> )
{
    return boundaryedges( *mesh, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::location_edge_const_iterator,
             typename mesh::Traits<MeshType>::location_edge_const_iterator>
boundaryedges( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              mesh.beginEdgeOnBoundary(),
                              mesh.endEdgeOnBoundary() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::location_edge_const_iterator,
             typename mesh::Traits<MeshType>::location_edge_const_iterator>
internaledges( MeshType const& mesh, mpl::bool_<true> )
{
    return internaledges( mesh, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::location_edge_const_iterator,
             typename mesh::Traits<MeshType>::location_edge_const_iterator>
internaledges( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              mesh.beginInternalEdge(),
                              mesh.endInternalEdge() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::point_const_iterator,
             typename mesh::Traits<MeshType>::point_const_iterator>
points( MeshType const& mesh, mpl::bool_<true> )
{
    return points( *mesh, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::point_const_iterator,
             typename mesh::Traits<MeshType>::point_const_iterator>
points( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginPoint(),
                              mesh.endPoint() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::marked_point_const_iterator,
             typename mesh::Traits<MeshType>::marked_point_const_iterator>
markedpoints( MeshType const& mesh, size_type flag, mpl::int_<true> )
{
    return markedpoints( *mesh, flag, mpl::int_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::marked_point_const_iterator,
             typename mesh::Traits<MeshType>::marked_point_const_iterator>
markedpoints( MeshType const& mesh, size_type flag, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginPointWithMarker( flag ),
                              mesh.endPointWithMarker( flag) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::location_point_const_iterator,
             typename mesh::Traits<MeshType>::location_point_const_iterator>
boundarypoints( MeshType const& mesh, mpl::bool_<true> )
{
    return boundarypoints( *mesh, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::location_point_const_iterator,
             typename mesh::Traits<MeshType>::location_point_const_iterator>
boundarypoints( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginPointOnBoundary( ),
                              mesh.endPointOnBoundary() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::location_point_const_iterator,
             typename mesh::Traits<MeshType>::location_point_const_iterator>
internalpoints( MeshType const& mesh, mpl::bool_<true> )
{
    return internalpoints( *mesh, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::location_point_const_iterator,
             typename mesh::Traits<MeshType>::location_point_const_iterator>
internalpoints( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginInternalPoint( ),
                              mesh.endInternalPoint() );
}

} // detail
/// \endcond


/**
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::element_const_iterator,
             typename mesh::Traits<MeshType>::element_const_iterator>
allelements( MeshType const& mesh )
{
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              mesh.beginElement(),
                              mesh.endElement() );
}

/**
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::element_const_iterator,
             typename mesh::Traits<MeshType>::element_const_iterator>
elements( MeshType const& mesh, flag_type flag = Application::processId() )
{
    return detail::elements( mesh, flag, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
    //return elements( mesh, flag, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}
/**
 * \return a pair of iterators to iterate over elements of the mesh
 * which share a face with the boundary
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::location_element_const_iterator,
             typename mesh::Traits<MeshType>::location_element_const_iterator>
boundaryelements( MeshType const& mesh, size_type pid  = Application::processId() )
{
    return detail::boundaryelements( mesh, pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}


/**
 * \return a pair of iterators to iterate over elements of the mesh
 * which are stricly within the domain that is to say they do not
 * share a face with the boundary
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::location_element_const_iterator,
             typename mesh::Traits<MeshType>::location_element_const_iterator>
internalelements( MeshType const& mesh, size_type pid  = Application::processId() )
{
    return detail::internalelements( mesh, pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::marker_element_const_iterator,
             typename mesh::Traits<MeshType>::marker_element_const_iterator>
markedelements( MeshType const& mesh, flag_type flag, size_type pid  = Application::processId() )
{
    return detail::markedelements( mesh, flag, pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker2 \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::marker2_element_const_iterator,
             typename mesh::Traits<MeshType>::marker2_element_const_iterator>
marked2elements( MeshType const& mesh, flag_type flag, size_type pid  = Application::processId() )
{
    return detail::marked2elements( mesh, flag, pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker3 \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::marker3_element_const_iterator,
             typename mesh::Traits<MeshType>::marker3_element_const_iterator>
marked3elements( MeshType const& mesh, flag_type flag, size_type pid  = Application::processId() )
{
    return detail::marked3elements( mesh, flag, pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * \return a pair of iterators to iterate over elements with id \p flag
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename mesh::Traits<MeshType>::element_const_iterator,
             typename mesh::Traits<MeshType>::element_const_iterator>
idedelements( MeshType const& mesh, flag_type flag )
{
    return detail::idedelements( mesh, flag, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
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
             typename mesh::Traits<MeshType>::pid_face_const_iterator,
             typename mesh::Traits<MeshType>::pid_face_const_iterator>
faces( MeshType const& mesh, size_type __pid = Application::processId() )
{
    return detail::faces( mesh, __pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * \return a pair of iterators to iterate over elements with id \p id
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::face_const_iterator,
             typename mesh::Traits<MeshType>::face_const_iterator>
idedfaces( MeshType const& mesh, size_type id )
{
    return detail::idedfaces( mesh, id, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
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
             typename mesh::Traits<MeshType>::marker_face_const_iterator,
             typename mesh::Traits<MeshType>::marker_face_const_iterator>
markedfaces( MeshType const& mesh,
             flag_type __marker, size_type __pid = Application::processId() )
{
    Debug(4000) << "[markedfaces] marker = " << __marker << "\n";
    return detail::markedfaces( mesh, __marker, __pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * \return a pair of iterators to iterate over all boundary faces of the
 * mesh
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::location_face_const_iterator,
             typename mesh::Traits<MeshType>::location_face_const_iterator>
boundaryfaces( MeshType const& mesh, size_type __pid = Application::processId()  )
{
    return detail::boundaryfaces( mesh, __pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}


/**
 * \return a pair of iterators to iterate over all internal faces of the
 * mesh belong to process domain \p __pid
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::location_face_const_iterator,
             typename mesh::Traits<MeshType>::location_face_const_iterator>
internalfaces( MeshType const& mesh, size_type __pid = Application::processId() )
{
    return detail::internalfaces( mesh, __pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * \return a pair of iterators to iterate over all interprocess faces of the
 * mesh belonging to process \p __pid
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename mesh::Traits<MeshType>::interprocess_face_const_iterator,
             typename mesh::Traits<MeshType>::interprocess_face_const_iterator>
interprocessfaces( MeshType const& mesh, size_type __pid = Application::processId() )
{
    return detail::interprocessfaces( mesh, __pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}


/**
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
             typename mesh::Traits<MeshType>::pid_edge_const_iterator,
             typename mesh::Traits<MeshType>::pid_edge_const_iterator>
edges( MeshType const& mesh,
       size_type __pid = Application::processId() )
{
    return detail::edges( mesh, __pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}


/**
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
             typename mesh::Traits<MeshType>::marker_edge_const_iterator,
             typename mesh::Traits<MeshType>::marker_edge_const_iterator>
markededges( MeshType const& mesh,
             flag_type __marker, size_type __pid = Application::processId() )
{
    return detail::markededges( mesh, __marker, __pid, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * \return a pair of iterators to iterate over all boundary edges of the
 * mesh
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::location_edge_const_iterator,
             typename mesh::Traits<MeshType>::location_edge_const_iterator>
boundaryedges( MeshType const& mesh )
{
    return detail::boundaryedges( mesh, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}


/**
 * \return a pair of iterators to iterate over all internal edges of the
 * mesh
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename mesh::Traits<MeshType>::location_edge_const_iterator,
             typename mesh::Traits<MeshType>::location_edge_const_iterator>
internaledges( MeshType const& mesh )
{
    return detail::internaledges( mesh, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * return the range of iterators [begin,end] over the points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::point_const_iterator,
             typename mesh::Traits<MeshType>::point_const_iterator>
points( MeshType const& mesh )
{
    return detail::points( mesh, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * return the range of iterators [begin,end] over the marked points with \p flag of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::marked_point_const_iterator,
             typename mesh::Traits<MeshType>::marked_point_const_iterator>
markedpoints( MeshType const& mesh, size_type flag )
{
    return detail::markedpoints( mesh, flag, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * return the range of iterators [begin,end] over the boundary points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::location_point_const_iterator,
             typename mesh::Traits<MeshType>::location_point_const_iterator>
boundarypoints( MeshType const& mesh )
{
    return detail::boundarypoints( mesh, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

/**
 * return the range of iterators [begin,end] over the internal(not on the boundary) points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename mesh::Traits<MeshType>::location_point_const_iterator,
             typename mesh::Traits<MeshType>::location_point_const_iterator>
internalpoints( MeshType const& mesh )
{

    return detail::internalpoints( mesh, mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >() );
}

}


#endif
