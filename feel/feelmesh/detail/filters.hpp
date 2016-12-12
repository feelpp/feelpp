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
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
elements( MeshType const& mesh, rank_type pid, mpl::bool_<false> )
{
    auto rangeElements = mesh.elementsWithProcessId( pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              std::get<0>( rangeElements ),
                              std::get<1>( rangeElements ),
                              std::get<2>( rangeElements ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
elements( MeshType const& mesh, rank_type pid, mpl::bool_<true> )
{
    return elements( *mesh, pid, mpl::bool_<false>() );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
boundaryelements( MeshType const& mesh, uint16_type entity_min_dim, uint16_type entity_max_dim, rank_type pid, mpl::bool_<false> )
{
    auto rangeElements = mesh.boundaryElements( entity_min_dim,entity_max_dim,pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              std::get<0>( rangeElements ),
                              std::get<1>( rangeElements ),
                              std::get<2>( rangeElements ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
boundaryelements( MeshType const& mesh, uint16_type entity_min_dim, uint16_type entity_max_dim, rank_type pid, mpl::bool_<true> )
{
    return boundaryelements( *mesh, entity_min_dim, entity_max_dim, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
internalelements( MeshType const& mesh, rank_type pid, mpl::bool_<false> )
{
    auto rangeElements = mesh.internalElements( pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              std::get<0>( rangeElements ),
                              std::get<1>( rangeElements ),
                              std::get<2>( rangeElements ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
internalelements( MeshType const& mesh, rank_type pid, mpl::bool_<true> )
{
    return internalelements( *mesh, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
markedelements( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, mpl::bool_<false>  )
{
    auto rangeElementsWithMarker = mesh.elementsWithMarkerByType( markerType, markersFlag, pid );
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
markedelements( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, mpl::bool_<true> )
{
    return markedelements( *mesh, markerType, markersFlag, pid, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
markedelements( MeshType const& mesh, uint16_type markerType, rank_type pid, mpl::bool_<false>  )
{
    auto rangeElementsWithMarker = mesh.elementsWithMarkerByType( markerType, pid );
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
markedelements( MeshType const& mesh, uint16_type markerType, rank_type pid, mpl::bool_<true> )
{
    return markedelements( *mesh, markerType, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
idedelements( MeshType const& mesh, size_type id, mpl::bool_<false> )
{
    auto rangeElementsWithId = mesh.elementsWithId( id );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              std::get<0>( rangeElementsWithId ),
                              std::get<1>( rangeElementsWithId ),
                              std::get<2>( rangeElementsWithId ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
idedelements( MeshType const& mesh, size_type id, mpl::bool_<true> )
{
    return idedelements( *mesh, id, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
faces( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{
    auto rangeFaces = mesh.facesWithProcessId( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              std::get<0>( rangeFaces ),
                              std::get<1>( rangeFaces ),
                              std::get<2>( rangeFaces ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
faces( MeshType const& mesh, rank_type __pid, mpl::bool_<true> )
{
    return faces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
idedfaces( MeshType const& mesh, size_type id, mpl::bool_<false> )
{
    auto rangeFaces = mesh.facesWithId( id );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              std::get<0>( rangeFaces ),
                              std::get<1>( rangeFaces ),
                              std::get<2>( rangeFaces ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
idedfaces( MeshType const& mesh, size_type id, mpl::bool_<true> )
{
    return idedfaces( *mesh, id, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
markedfaces( MeshType const& mesh, uint16_type markerType, rank_type pid, mpl::bool_<false> )
{
    auto rangeMarkedFaces = mesh.facesWithMarkerByType( markerType, pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              std::get<0>( rangeMarkedFaces ),
                              std::get<1>( rangeMarkedFaces ),
                              std::get<2>( rangeMarkedFaces ) );

}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
markedfaces( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, mpl::bool_<false> )
{
    auto rangeMarkedFaces = mesh.facesWithMarkerByType( markerType, markersFlag, pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              std::get<0>( rangeMarkedFaces ),
                              std::get<1>( rangeMarkedFaces ),
                              std::get<2>( rangeMarkedFaces ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
markedfaces( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type __pid, mpl::bool_<true> )
{
    return markedfaces( *mesh,   markerType, markersFlag, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
markedfaces( MeshType const& mesh, uint16_type markerType, rank_type __pid, mpl::bool_<true> )
{
    return markedfaces( *mesh,  markerType, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
boundaryfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<false>  )
{
    auto rangeBoundaryFaces = mesh.facesOnBoundary( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              std::get<0>( rangeBoundaryFaces ),
                              std::get<1>( rangeBoundaryFaces ),
                              std::get<2>( rangeBoundaryFaces ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
boundaryfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<true>  )
{
    return boundaryfaces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
internalfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{
    auto rangeInternalFaces = mesh.internalFaces( __pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              std::get<0>( rangeInternalFaces ),
                              std::get<1>( rangeInternalFaces ),
                              std::get<2>( rangeInternalFaces ) );

}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
internalfaces( MeshType const& mesh, rank_type __pid, mpl::bool_<true> )
{
    return internalfaces( *mesh, __pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
interprocessfaces( MeshType const& mesh, rank_type neighbor_pid, mpl::bool_<false> )
{
    auto rangeInterprocessFaces = mesh.interProcessFaces( neighbor_pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              std::get<0>( rangeInterprocessFaces ),
                              std::get<1>( rangeInterprocessFaces ),
                              std::get<2>( rangeInterprocessFaces ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
interprocessfaces( MeshType const& mesh, rank_type neighbor_pid, mpl::bool_<true> )
{
    return interprocessfaces( *mesh, neighbor_pid, mpl::bool_<false>() );

}

template<typename MeshType>
decltype(auto)
edges( MeshType const& mesh, rank_type __pid, mpl::bool_<false> )
{
    auto rangeMarkedEdges = mesh.edgesWithProcessId( __pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              std::get<0>( rangeMarkedEdges ),
                              std::get<1>( rangeMarkedEdges ),
                              std::get<2>( rangeMarkedEdges ) );
}
template<typename MeshType>
decltype(auto)
edges( MeshType const& mesh, rank_type __pid, mpl::bool_<true> )
{
    return edges( *mesh, __pid, mpl::bool_<false>() );
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
markededges( MeshType const& mesh, uint16_type markerType, rank_type pid, mpl::bool_<false> )
{
    auto rangeMarkedEdges = mesh.edgesWithMarkerByType( markerType, pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              std::get<0>( rangeMarkedEdges ),
                              std::get<1>( rangeMarkedEdges ),
                              std::get<2>( rangeMarkedEdges ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
markededges( MeshType const& mesh, uint16_type markerType, rank_type pid, mpl::bool_<true> )
{
    return markededges( *mesh, markerType, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
markededges( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, mpl::bool_<false> )
{
    auto rangeMarkedEdges = mesh.edgesWithMarkerByType( markerType, markersFlag, pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              std::get<0>( rangeMarkedEdges ),
                              std::get<1>( rangeMarkedEdges ),
                              std::get<2>( rangeMarkedEdges ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
markededges( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, mpl::bool_<true> )
{
    return markededges( *mesh, markerType, markersFlag, pid, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
boundaryedges( MeshType const& mesh, mpl::bool_<false> )
{
    auto rangeBoundaryEdges = mesh.edgesOnBoundary();
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              std::get<0>( rangeBoundaryEdges ),
                              std::get<1>( rangeBoundaryEdges ),
                              std::get<2>( rangeBoundaryEdges ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
boundaryedges( MeshType const& mesh, mpl::bool_<true> )
{
    return boundaryedges( *mesh, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
internaledges( MeshType const& mesh, mpl::bool_<false> )
{
    auto rangeInternalEdges = mesh.internalEdges();
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                              std::get<0>( rangeInternalEdges ),
                              std::get<1>( rangeInternalEdges ),
                              std::get<2>( rangeInternalEdges ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
internaledges( MeshType const& mesh, mpl::bool_<true> )
{
    return internaledges( *mesh, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_const_iterator,
             typename MeshTraits<MeshType>::point_const_iterator>
allpoints( MeshType const& mesh, mpl::bool_<false> )
{
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              mesh.beginPoint(),
                              mesh.endPoint() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
      typename MeshTraits<MeshType>::point_const_iterator,
      typename MeshTraits<MeshType>::point_const_iterator>
allpoints( MeshType const& mesh, mpl::bool_<true> )
{
    return allpoints( *mesh, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
points( MeshType const& mesh, mpl::bool_<false> )
{
    auto rangePoints = mesh.pointsWithProcessId();
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              std::get<0>( rangePoints ),
                              std::get<1>( rangePoints ),
                              std::get<2>( rangePoints ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
points( MeshType const& mesh, mpl::bool_<true> )
{
    return points( *mesh, mpl::bool_<false>() );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
markedpoints( MeshType const& mesh, uint16_type markerType, rank_type pid, mpl::bool_<false> )
{
    auto rangePointsWithMarker = mesh.pointsWithMarkerByType( markerType, pid );
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
markedpoints( MeshType const& mesh, uint16_type markerType, rank_type pid, mpl::bool_<true> )
{
    return markedpoints( *mesh, markerType, pid, mpl::bool_<false>() );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
markedpoints( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, mpl::bool_<false> )
{
    auto rangePointsWithMarker = mesh.pointsWithMarkerByType( markerType, markersFlag, pid );
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
markedpoints( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, mpl::bool_<true> )
{
    return markedpoints( *mesh, markerType, markersFlag, pid, mpl::bool_<false>() );
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
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
internalpoints( MeshType const& mesh, mpl::bool_<false> )
{
    auto rangeInternalPoints = mesh.internalPoints();
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(),
                              std::get<0>( rangeInternalPoints ),
                              std::get<1>( rangeInternalPoints ),
                              std::get<2>( rangeInternalPoints ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
internalpoints( MeshType const& mesh, mpl::bool_<true> )
{
    return internalpoints( *mesh, mpl::bool_<false>() );
}

} // detail
/// \endcond

}
#endif
