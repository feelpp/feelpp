/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 29 juil. 2015

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
#ifndef FEELPP_FEELMESH_DETAIL_FILTERS_HPP
#define FEELPP_FEELMESH_DETAIL_FILTERS_HPP 1

namespace Feel {


/// \cond detail
namespace detail
{

template <typename RangeType>
struct submeshrangetype
{
    typedef typename mpl::if_< boost::is_std_list<RangeType>,
                               mpl::identity<RangeType>,
                               mpl::identity<std::list<RangeType> > >::type::type::value_type type;
};



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
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
markedelements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<false>  )
{
    auto rangeElementsWithMarker = mesh.elementsWithMarker( flag );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              std::get<0>( rangeElementsWithMarker ),
                              std::get<1>( rangeElementsWithMarker ),
                              std::get<2>( rangeElementsWithMarker ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
markedelements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<true> )
{
    return markedelements( *mesh, flag, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
marked2elements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker2_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.elementsWithMarker2( flag, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
marked2elements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<true> )
{
    return marked2elements( *mesh, flag, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
marked3elements( MeshType const& mesh, flag_type flag, rank_type pid, mpl::bool_<false> )
{
    typedef typename MeshTraits<MeshType>::marker3_element_const_iterator iterator;
    std::pair<iterator, iterator> p = mesh.elementsWithMarker3( flag, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              p.first, p.second );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
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
decltype(auto)
edges( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{
    auto r = mesh.edgesWithProcessId( __pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(), r.first, r.second );
}
template<typename MeshType>
decltype(auto)
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
points( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginPoint(),
                              mesh.endPoint() );
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
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
markedpoints( MeshType const& mesh, size_type flag, mpl::bool_<false> )
{
    auto rangePointsWithMarker = mesh.pointsWithMarker( flag );
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              std::get<0>( rangePointsWithMarker ),
                              std::get<1>( rangePointsWithMarker ),
                              std::get<2>( rangePointsWithMarker ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
markedpoints( MeshType const& mesh, size_type flag, mpl::bool_<true> )
{
    return markedpoints( *mesh, flag, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
boundarypoints( MeshType const& mesh, mpl::bool_<false> )
{
    auto rangeBoundaryPoints = mesh.boundaryPoints();
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              std::get<0>( rangeBoundaryPoints ),
                              std::get<1>( rangeBoundaryPoints ),
                              std::get<2>( rangeBoundaryPoints ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
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

}
#endif
