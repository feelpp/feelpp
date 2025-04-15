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

#include <feel/feelcore/unwrapptr.hpp>

namespace Feel {

#pragma GCC visibility push(hidden)
/// \cond detail
namespace detail
{

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
elements( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangeElements = Feel::unwrap_ptr( mesh ).template elementsFilter<entity_filter_t::PROCESS_ID>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), std::get<0>( rangeElements ), std::get<1>( rangeElements ), std::get<2>( rangeElements ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
boundaryelements( MeshType const& mesh, uint16_type entity_min_dim, uint16_type entity_max_dim, rank_type pid, entity_process_t ept )
{
    auto rangeElements = Feel::unwrap_ptr( mesh ).template elementsFilter<entity_filter_t::ON_BOUNDARY>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), std::get<0>( rangeElements ), std::get<1>( rangeElements ), std::get<2>( rangeElements ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
internalelements( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangeElements = Feel::unwrap_ptr( mesh ).template elementsFilter<entity_filter_t::INTERNAL>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), std::get<0>( rangeElements ), std::get<1>( rangeElements ), std::get<2>( rangeElements ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
markedelements( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, entity_process_t ept )
{
    auto rangeElements = Feel::unwrap_ptr( mesh ).template elementsFilter<entity_filter_t::MARKER>( ept, markerType, markersFlag, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), std::get<0>( rangeElements ), std::get<1>( rangeElements ), std::get<2>( rangeElements ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
markedelements( MeshType const& mesh, uint16_type markerType, rank_type pid, entity_process_t ept )
{
    auto rangeElements = Feel::unwrap_ptr( mesh ).template elementsFilter<entity_filter_t::MARKER>( ept, markerType, invalid_flag_type_value, pid );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), std::get<0>( rangeElements ), std::get<1>( rangeElements ), std::get<2>( rangeElements ) );
}

#if 0
template<int TheType,typename MeshType>
auto//std::map<int,markedelements_t<MeshType> >
collectionOfMarkedelements( MeshType const& mesh, std::any const& collectionOfMarkersFlag )
{
    std::map<int,std::set<flag_type>> collectionOfMarkerFlagSet;
    if ( auto argCasted = std::any_cast<std::map<int,int>>( &collectionOfMarkersFlag) )
    {
        for ( auto const& [part,markersFlag] : *argCasted )
            collectionOfMarkerFlagSet[part] = { Feel::unwrap_ptr( mesh ).markerId( markersFlag ) };
    }
    else if ( auto argCasted = std::any_cast<std::map<int,std::set<std::string>>>( &collectionOfMarkersFlag) )
    {
        for ( auto const& [part,markerNames] : *argCasted )
            collectionOfMarkerFlagSet[part] = Feel::unwrap_ptr( mesh ).markersId( markerNames );
    }
    else
        CHECK( false ) << "TODO : others cast";

    uint16_type markerType = 1;
    rank_type pid = rank( mesh );
    auto collectionOfRangeElement = Feel::unwrap_ptr( mesh ).template collectionOfElementsWithMarkerByType<TheType>( markerType, collectionOfMarkerFlagSet, pid );
    if constexpr ( TheType == 0 || TheType == 1 )
    {
        std::map<int,markedelements_t<MeshType> > res;
        for ( auto const& [part,rangeElementsWithMarker] : collectionOfRangeElement )
            res[part] = boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                           std::get<0>( rangeElementsWithMarker ),
                                           std::get<1>( rangeElementsWithMarker ),
                                           std::get<2>( rangeElementsWithMarker ) );
        return res;
    }
    else
    {
        std::map<int,std::tuple<markedelements_t<MeshType>,markedelements_t<MeshType>> > res;
        for ( auto const& [part,pairRangeElementsWithMarker] : collectionOfRangeElement )
        {
            auto rangeActiveElementsWithMarker = std::get<0>( pairRangeElementsWithMarker );
            auto rangeGhostElementsWithMarker = std::get<1>( pairRangeElementsWithMarker );
            res[part] = std::make_tuple( boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                                            std::get<0>( rangeActiveElementsWithMarker ),
                                                            std::get<1>( rangeActiveElementsWithMarker ),
                                                            std::get<2>( rangeActiveElementsWithMarker ) ),
                                         boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                                            std::get<0>( rangeGhostElementsWithMarker ),
                                                            std::get<1>( rangeGhostElementsWithMarker ),
                                                            std::get<2>( rangeGhostElementsWithMarker ) ) );
        }
        return res;
    }

}
#else
template<int TheType, typename MeshType>
using CollectionOfMarkedElementsResultType
                = std::conditional_t<TheType == 0 || TheType == 1,
                                     std::map<int, Range<MeshType,ElementsType::MESH_ELEMENTS>>,
                                     std::map<int, std::tuple<Range<MeshType,ElementsType::MESH_ELEMENTS>, Range<MeshType,ElementsType::MESH_ELEMENTS>>>>;

template<typename MeshType,typename TupleRange>
Range<decay_type<MeshType>,MESH_ELEMENTS> makeResultRange( std::shared_ptr<MeshType> const& mesh, TupleRange && r )
{
    auto&& meshRange = std::forward<TupleRange>(r);
    wc( mesh )->print( fmt::format( "makeResultRange: mesh count={}", mesh.use_count() ), FLAGS_v > 0, FLAGS_v > 0, FLAGS_v > 0 );

    //std::cout << fmt::format( "mesh: {} ptr:{}, count: {}\n", ( is_shared_ptr_v<decay_type<MeshType>>  ) ? "shared_ptr" : "ref_or_val", mesh, mesh.use_count()  ) << std::endl;
    Range<decay_type<MeshType>, MESH_ELEMENTS> res( mesh );
    for ( auto const& elt : *std::get<2>( meshRange ) )
        res.push_back( elt );
    res.shrink_to_fit();
    wc( mesh )->print( fmt::format( "makeResultRange: res size={}\n", res.size() ));
    return res;
}

template<int TheType, typename MeshType>
CollectionOfMarkedElementsResultType<TheType, MeshType>
collectionOfMarkedelements(MeshType const& mesh, std::any const& collectionOfMarkersFlag)
{
    std::map<int, std::set<flag_type>> collectionOfMarkerFlagSet;

    if (auto argCasted = std::any_cast<std::map<int, int>>(&collectionOfMarkersFlag))
    {
        for (auto const& [part, markersFlag] : *argCasted)
            collectionOfMarkerFlagSet[part] = { Feel::unwrap_ptr(mesh).markerId(markersFlag) };
    }
    else if (auto argCasted = std::any_cast<std::map<int, std::set<std::string>>>(&collectionOfMarkersFlag))
    {
        for (auto const& [part, markerNames] : *argCasted)
            collectionOfMarkerFlagSet[part] = Feel::unwrap_ptr(mesh).markersId(markerNames);
    }
    else
    {
        throw std::runtime_error("Unsupported data type in collectionOfMarkersFlag");
    }

    auto collectionOfRangeElement = Feel::unwrap_ptr(mesh).template collectionOfElementsWithMarkerByType<TheType>(1, collectionOfMarkerFlagSet, rank(mesh));

    CollectionOfMarkedElementsResultType<TheType, MeshType> res;
    for (auto const& [part, rangeElementsWithMarker] : collectionOfRangeElement)
    {
        if constexpr (TheType == 0 || TheType == 1)
        {
            res[part] = makeResultRange(mesh,rangeElementsWithMarker);
        }
        else
        {
            auto const& [rangeActive, rangeGhost] = rangeElementsWithMarker;
            res[part] = std::make_tuple( makeResultRange( mesh, rangeActive ), makeResultRange(mesh,rangeGhost) );
            CHECK( std::get<0>(res[part]).mesh() ) << "mesh with active elements is null";
            CHECK( std::get<1>(res[part]).mesh() ) << "mesh with ghost elements is null";
            CHECK( std::get<0>(res[part]).mesh() == std::get<1>(res[part]).mesh() );
        }
    }

    return res;
}
#endif



template<typename MeshType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype >
idedelements( MeshType const& mesh, size_type id )
{
    auto rangeElements = Feel::unwrap_ptr( mesh ).elementsWithId( id );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(), std::get<0>( rangeElements ), std::get<1>( rangeElements ), std::get<2>( rangeElements ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
faces( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangeFaces = Feel::unwrap_ptr( mesh ).template facesFilter<entity_filter_t::PROCESS_ID>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(), std::get<0>( rangeFaces ), std::get<1>( rangeFaces ), std::get<2>( rangeFaces ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
idedfaces( MeshType const& mesh, size_type id )
{
    auto rangeFaces = Feel::unwrap_ptr( mesh ).facesWithId( id );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(), std::get<0>( rangeFaces ), std::get<1>( rangeFaces ), std::get<2>( rangeFaces ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
markedfaces( MeshType const& mesh, uint16_type markerType, rank_type pid, entity_process_t ept )
{
    auto rangeFaces = Feel::unwrap_ptr( mesh ).template facesFilter<entity_filter_t::MARKER>( ept, markerType, invalid_flag_type_value, pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(), std::get<0>( rangeFaces ), std::get<1>( rangeFaces ), std::get<2>( rangeFaces ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
markedfaces( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, entity_process_t ept )
{
    auto rangeFaces = Feel::unwrap_ptr( mesh ).template facesFilter<entity_filter_t::MARKER>( ept, markerType, markersFlag, pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(), std::get<0>( rangeFaces ), std::get<1>( rangeFaces ), std::get<2>( rangeFaces ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
boundaryfaces( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangeFaces = Feel::unwrap_ptr( mesh ).template facesFilter<entity_filter_t::ON_BOUNDARY>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(), std::get<0>( rangeFaces ), std::get<1>( rangeFaces ), std::get<2>( rangeFaces ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
internalfaces( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangeFaces = Feel::unwrap_ptr( mesh ).template facesFilter<entity_filter_t::INTERNAL>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(), std::get<0>( rangeFaces ), std::get<1>( rangeFaces ), std::get<2>( rangeFaces ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype >
interprocessfaces( MeshType const& mesh, rank_type neighbor_pid )
{
    auto rangeInterprocessFaces = Feel::unwrap_ptr( mesh ).interProcessFaces( neighbor_pid );
    return boost::make_tuple( mpl::size_t<MESH_FACES>(), std::get<0>( rangeInterprocessFaces ), std::get<1>( rangeInterprocessFaces ), std::get<2>( rangeInterprocessFaces ) );
}

template<typename MeshType>
decltype(auto)
edges( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangeEdges = Feel::unwrap_ptr( mesh ).template edgesFilter<entity_filter_t::PROCESS_ID>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(), std::get<0>( rangeEdges ), std::get<1>( rangeEdges ), std::get<2>( rangeEdges ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
markededges( MeshType const& mesh, uint16_type markerType, rank_type pid, entity_process_t ept )
{
    auto rangeEdges = Feel::unwrap_ptr( mesh ).template edgesFilter<entity_filter_t::MARKER>( ept, markerType, invalid_flag_type_value, pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(), std::get<0>( rangeEdges ), std::get<1>( rangeEdges ), std::get<2>( rangeEdges ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
markededges( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, entity_process_t ept )
{
    auto rangeEdges = Feel::unwrap_ptr( mesh ).template edgesFilter<entity_filter_t::MARKER>( ept, markerType, markersFlag, pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(), std::get<0>( rangeEdges ), std::get<1>( rangeEdges ), std::get<2>( rangeEdges ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
boundaryedges( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangeEdges = Feel::unwrap_ptr( mesh ).template edgesFilter<entity_filter_t::ON_BOUNDARY>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(), std::get<0>( rangeEdges ), std::get<1>( rangeEdges ), std::get<2>( rangeEdges ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_EDGES>,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype >
internaledges( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangeEdges = Feel::unwrap_ptr( mesh ).template facesFilter<entity_filter_t::INTERNAL>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_EDGES>(), std::get<0>( rangeEdges ), std::get<1>( rangeEdges ), std::get<2>( rangeEdges ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
points( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangePoints = Feel::unwrap_ptr( mesh ).template pointsFilter<entity_filter_t::PROCESS_ID>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(), std::get<0>( rangePoints ), std::get<1>( rangePoints ), std::get<2>( rangePoints ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
markedpoints( MeshType const& mesh, uint16_type markerType, rank_type pid, entity_process_t ept )
{
    auto rangePoints = Feel::unwrap_ptr( mesh ).template pointsFilter<entity_filter_t::MARKER>( ept, markerType, invalid_flag_type_value, pid );
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(), std::get<0>( rangePoints ), std::get<1>( rangePoints ), std::get<2>( rangePoints ) );
}
template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
markedpoints( MeshType const& mesh, uint16_type markerType, std::set<flag_type> const& markersFlag, rank_type pid, entity_process_t ept )
{
    auto rangePoints = Feel::unwrap_ptr( mesh ).template pointsFilter<entity_filter_t::MARKER>( ept, markerType, markersFlag, pid );
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(), std::get<0>( rangePoints ), std::get<1>( rangePoints ), std::get<2>( rangePoints ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
boundarypoints( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangePoints = Feel::unwrap_ptr( mesh ).template pointsFilter<entity_filter_t::ON_BOUNDARY>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(), std::get<0>( rangePoints ), std::get<1>( rangePoints ), std::get<2>( rangePoints ) );
}

template<typename MeshType>
boost::tuple<mpl::size_t<MESH_POINTS>,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
             typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype >
internalpoints( MeshType const& mesh, rank_type pid, entity_process_t ept )
{
    auto rangePoints = Feel::unwrap_ptr( mesh ).template pointsFilter<entity_filter_t::INTERNAL>( ept, pid );
    return boost::make_tuple( mpl::size_t<MESH_POINTS>(), std::get<0>( rangePoints ), std::get<1>( rangePoints ), std::get<2>( rangePoints ) );
}

} // detail
/// \endcond
#pragma GCC visibility pop

}
#endif
