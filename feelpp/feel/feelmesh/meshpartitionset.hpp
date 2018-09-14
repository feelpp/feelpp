/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-10-16

  Copyright (C) 2013-2016 Feel++ Consortium

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
#ifndef FEELPP_MESHPARTITIONSET_HPP
#define FEELPP_MESHPARTITIONSET_HPP 1

namespace Feel
{


template<typename MeshType>
class MeshPartitionSet
{
public:
    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef std::vector<boost::reference_wrapper<typename MeshTraits<mesh_type>::point_type const> > point_container_type;
    typedef typename point_container_type::const_iterator point_const_iterator;
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<mesh_type>::edge_type const> > edge_container_type;
    typedef typename edge_container_type::const_iterator edge_const_iterator;
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<mesh_type>::face_type const> > face_container_type;
    typedef typename face_container_type::const_iterator face_const_iterator;
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<mesh_type>::element_type const> > element_container_type;
    typedef typename element_container_type::const_iterator element_const_iterator;

    MeshPartitionSet( mesh_ptrtype const& mesh )
        :
        M_mesh( mesh ),
        M_numGlobalPartition( mesh->worldComm().localSize() )
        {
            rank_type partId = M_mesh->worldComm().localRank();
            M_localPartitionIds.insert( partId );
            for ( rank_type p=0;p<M_numGlobalPartition;++p )
            {
                M_statistic[p].resize( 6 );
                M_statistic[p][0] = M_mesh->statNumPointsAll( p );
                M_statistic[p][1] = M_mesh->statNumElementsAll( p );
                M_statistic[p][2] = M_mesh->statNumElementsActive( p );
                M_statistic[p][3] = M_mesh->statNumFacesMarkedAll( p );
                M_statistic[p][4] = M_mesh->statNumEdgesMarkedAll( p );
                M_statistic[p][5] = M_mesh->statNumPointsMarkedAll( p );
            }

            M_containerPoints[partId].clear();
            M_containerActiveElements[partId].clear();
            M_containerGhostElements[partId].clear();
            M_containerMarkedFaces[partId].clear();
            M_containerMarkedEdges[partId].clear();
            M_containerMarkedPoints[partId].clear();

            auto rangeElements = M_mesh->elementsWithProcessId( partId );
            auto elt_it = std::get<0>( rangeElements );
            auto const elt_en = std::get<1>( rangeElements );
            for( ; elt_it != elt_en; ++ elt_it )
                M_containerActiveElements[partId].push_back( boost::cref( boost::unwrap_ref( *elt_it ) ) );
            CHECK( M_containerActiveElements[partId].size() == M_statistic[partId][2] ) << "something is wrong in active element : "
                                                                                   << M_containerActiveElements[partId].size() << " vs "
                                                                                   << M_statistic[partId][2];
            auto rangeGhostElement = M_mesh->ghostElements();
            auto ghostelt_it = std::get<0>( rangeGhostElement );
            auto ghostelt_en = std::get<1>( rangeGhostElement );
            for( ; ghostelt_it != ghostelt_en; ++ghostelt_it )
                M_containerGhostElements[partId].push_back( boost::cref( boost::unwrap_ref( *ghostelt_it ) ) );
            CHECK( M_containerGhostElements[partId].size() == (M_statistic[partId][1]-M_statistic[partId][2]) ) << "something is wrong in active element : "
                                                                                                      << M_containerGhostElements[partId].size() << " vs "
                                                                                                      << (M_statistic[partId][1]-M_statistic[partId][2]);

            auto pt_it = M_mesh->beginOrderedPoint();
            auto pt_en = M_mesh->endOrderedPoint();
            for( ; pt_it != pt_en; ++pt_it )
                M_containerPoints[partId].push_back(boost::cref(unwrap_ref( *pt_it ) ));

            this->updateMarkedSubEntitiesOnePartPerProcess<mesh_type>();
        }

    MeshPartitionSet( mesh_ptrtype const& mesh, rank_type nGlobalPart, std::set<rank_type> const& localPartitionIds )
        :
        M_mesh( mesh ),
        M_numGlobalPartition( nGlobalPart ),
        M_localPartitionIds( localPartitionIds )
        {
            CHECK( M_localPartitionIds.size() <= M_numGlobalPartition  ) << "number of local partition (in process) can not be greater than number of partition in full mesh";

            for ( rank_type partId : M_localPartitionIds )
            {
                M_containerPoints[partId].clear();
                M_containerActiveElements[partId].clear();
                M_containerGhostElements[partId].clear();
                M_containerMarkedFaces[partId].clear();
                M_containerMarkedEdges[partId].clear();
                M_containerMarkedPoints[partId].clear();
            }

            if ( M_localPartitionIds.size() == M_numGlobalPartition )
                this->buildAllPartInOneProcess();
            else
                CHECK( false ) << "TODO";
        }

    MeshPartitionSet( MeshPartitionSet const& ) = default;
    MeshPartitionSet( MeshPartitionSet && ) = default;
    ~MeshPartitionSet() {}

    MeshPartitionSet& operator=( MeshPartitionSet const& ) = default;
    MeshPartitionSet& operator=( MeshPartitionSet && ) noexcept = default;

    mesh_ptrtype const& mesh() const { return M_mesh; }

    rank_type numGlobalPartition() const { return M_numGlobalPartition; }
    rank_type numLocalPartition() const { return M_localPartitionIds.size(); }
    std::set<rank_type> const& localPartitionIds() const { return M_localPartitionIds; }
    bool hasLocalPartition( rank_type p ) const { return M_localPartitionIds.find(p) != M_localPartitionIds.end(); }
    void checkHasLocalPartition( rank_type p ) const { CHECK( this->hasLocalPartition(p) ) << "local partition " << p << " is not present"; }

    bool hasStatisticPartition( rank_type p ) const { return M_statistic.find(p) != M_statistic.end(); }
    void checkHasStatisticPartition( rank_type p ) const { CHECK( this->hasStatisticPartition(p) ) << "partition " << p << " is not present"; }
    size_type statNumPointsAll( rank_type p ) const { this->checkHasStatisticPartition(p); return M_statistic.find(p)->second[0]; }
    size_type statNumElementsAll( rank_type p ) const { this->checkHasStatisticPartition(p); return M_statistic.find(p)->second[1]; }
    size_type statNumElementsActive( rank_type p ) const { this->checkHasStatisticPartition(p); return M_statistic.find(p)->second[2]; }
    size_type statNumFacesMarkedAll( rank_type p ) const { this->checkHasStatisticPartition(p); return M_statistic.find(p)->second[3]; }
    size_type statNumEdgesMarkedAll( rank_type p ) const { this->checkHasStatisticPartition(p); return M_statistic.find(p)->second[4]; }
    size_type statNumPointsMarkedAll( rank_type p ) const { this->checkHasStatisticPartition(p); return M_statistic.find(p)->second[5]; }

    point_const_iterator beginPoint( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerPoints.find(p)->second.begin(); }
    point_const_iterator endPoint( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerPoints.find(p)->second.end(); }
    element_const_iterator beginActiveElement( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerActiveElements.find(p)->second.begin(); }
    element_const_iterator endActiveElement( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerActiveElements.find(p)->second.end(); }
    element_const_iterator beginGhostElement( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerGhostElements.find(p)->second.begin(); }
    element_const_iterator endGhostElement( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerGhostElements.find(p)->second.end(); }

    point_const_iterator beginMarkedPoint( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerMarkedPoints.find(p)->second.begin(); }
    point_const_iterator endMarkedPoint( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerMarkedPoints.find(p)->second.end(); }
    edge_const_iterator beginMarkedEdge( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerMarkedEdges.find(p)->second.begin(); }
    edge_const_iterator endMarkedEdge( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerMarkedEdges.find(p)->second.end(); }
    face_const_iterator beginMarkedFace( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerMarkedFaces.find(p)->second.begin(); }
    face_const_iterator endMarkedFace( rank_type p ) const { this->checkHasLocalPartition(p); return M_containerMarkedFaces.find(p)->second.end(); }

private:

    template<typename MT>
    void updateMarkedSubEntitiesOnePartPerProcess( typename std::enable_if<is_3d<MT>::value>::type* = nullptr );
    template<typename MT>
    void updateMarkedSubEntitiesOnePartPerProcess( typename std::enable_if<is_2d<MT>::value>::type* = nullptr );
    template<typename MT>
    void updateMarkedSubEntitiesOnePartPerProcess( typename std::enable_if<is_1d<MT>::value>::type* = nullptr );

    void buildAllPartInOneProcess();
    template<typename MT>
    void updateMarkedSubEntitiesAllPartInOneProcess( std::map<size_type,std::set<rank_type> > const& mapPointInterProcess,
                                                     typename std::enable_if<is_3d<MT>::value>::type* = nullptr );
    template<typename MT>
    void updateMarkedSubEntitiesAllPartInOneProcess( std::map<size_type,std::set<rank_type> > const& mapPointInterProcess,
                                                     typename std::enable_if<is_2d<MT>::value>::type* = nullptr );
    template<typename MT>
    void updateMarkedSubEntitiesAllPartInOneProcess( std::map<size_type,std::set<rank_type> > const& mapPointInterProcess,
                                                     typename std::enable_if<is_1d<MT>::value>::type* = nullptr );


private :
    mesh_ptrtype M_mesh;
    rank_type M_numGlobalPartition;
    std::set<rank_type> M_localPartitionIds;
    std::map<rank_type,std::vector<size_type> > M_statistic;

    std::map<rank_type,point_container_type> M_containerPoints;
    std::map<rank_type,element_container_type> M_containerActiveElements, M_containerGhostElements;

    std::map<rank_type,face_container_type> M_containerMarkedFaces;
    std::map<rank_type,edge_container_type> M_containerMarkedEdges;
    std::map<rank_type,point_container_type> M_containerMarkedPoints;
};


template<typename MeshType>
template<typename MT>
void
MeshPartitionSet<MeshType>::updateMarkedSubEntitiesOnePartPerProcess( typename std::enable_if<is_3d<MT>::value>::type* )
{
    rank_type partId = M_mesh->worldComm().localRank();
    auto face_it = M_mesh->beginOrderedFace();
    auto face_en = M_mesh->endOrderedFace();
    for( ; face_it != face_en; ++face_it )
    {
        auto const& face = unwrap_ref( *face_it );
        if ( !face.hasMarker() ) continue;
        M_containerMarkedFaces[partId].push_back(boost::cref(face));
    }
    auto edge_it = M_mesh->beginOrderedEdge();
    auto edge_en = M_mesh->endOrderedEdge();
    for( ; edge_it != edge_en; ++edge_it )
    {
        auto const& edge = unwrap_ref( *edge_it );
        if ( !edge.hasMarker() ) continue;
        M_containerMarkedEdges[partId].push_back(boost::cref(edge));
    }
    auto point_it = M_mesh->beginOrderedPoint();
    auto point_en = M_mesh->endOrderedPoint();
    for( ; point_it != point_en; ++point_it )
    {
        auto const& point = unwrap_ref( *point_it );
        if ( !point.hasMarker() ) continue;
        M_containerMarkedPoints[partId].push_back(boost::cref(point));
    }
}
template<typename MeshType>
template<typename MT>
void
MeshPartitionSet<MeshType>::updateMarkedSubEntitiesOnePartPerProcess( typename std::enable_if<is_2d<MT>::value>::type* )
{
    rank_type partId = M_mesh->worldComm().localRank();
    auto face_it = M_mesh->beginOrderedFace();
    auto face_en = M_mesh->endOrderedFace();
    for( ; face_it != face_en; ++face_it )
    {
        auto const& face = unwrap_ref( *face_it );
        if ( !face.hasMarker() ) continue;
        M_containerMarkedFaces[partId].push_back(boost::cref(face));
    }
    auto point_it = M_mesh->beginOrderedPoint();
    auto point_en = M_mesh->endOrderedPoint();
    for( ; point_it != point_en; ++point_it )
    {
        auto const& point = unwrap_ref( *point_it );
        if ( !point.hasMarker() ) continue;
        M_containerMarkedPoints[partId].push_back(boost::cref(point));
    }
}
template<typename MeshType>
template<typename MT>
void
MeshPartitionSet<MeshType>::updateMarkedSubEntitiesOnePartPerProcess( typename std::enable_if<is_1d<MT>::value>::type* )
{
    rank_type partId = M_mesh->worldComm().localRank();
    auto point_it = M_mesh->beginOrderedPoint();
    auto point_en = M_mesh->endOrderedPoint();
    for( ; point_it != point_en; ++point_it )
    {
        auto const& point = unwrap_ref( *point_it );
        if ( !point.hasMarker() ) continue;
        M_containerMarkedPoints[partId].push_back(boost::cref(point));
    }
}


template<typename MeshType>
void
MeshPartitionSet<MeshType>::buildAllPartInOneProcess()
{
    auto point_it = M_mesh->beginOrderedPoint();
    auto point_en = M_mesh->endOrderedPoint();
    for( ; point_it != point_en; ++point_it )
        unwrap_ref( *point_it ).setProcessId( invalid_rank_type_value );

    rank_type interprocessPointPidDetection = this->numGlobalPartition();

    std::map<size_type,std::set<rank_type> > mapPointInterProcess;

    // update active element containers + detect interprocess points
    for ( rank_type partId : M_localPartitionIds )
    {
        auto rangeElements = M_mesh->elementsWithProcessId( partId );
        auto elt_it = std::get<0>( rangeElements );
        auto elt_en = std::get<1>( rangeElements );
        // auto elt_it = M_mesh->beginElementWithProcessId( partId );
        // auto elt_en = M_mesh->endElementWithProcessId( partId );
        for ( ; elt_it != elt_en ; ++elt_it )
        {
            auto const& elt = boost::unwrap_ref( *elt_it );
            M_containerActiveElements[partId].push_back(boost::cref(elt));

            for ( uint16_type vLocId = 0 ; vLocId < mesh_type::element_type::numPoints; ++vLocId )
            {
                auto const& thepoint = elt.point( vLocId );
                size_type ptPid = thepoint.processId();
                if ( ptPid == invalid_rank_type_value ) // first time that we see this point
                    M_mesh->pointIterator( thepoint.id() )->second.setProcessId( partId );
                else if ( ptPid != partId ) // we found an interprocess point
                {
                    size_type ptId = thepoint.id();
                    mapPointInterProcess[ptId].insert( partId );
                    if ( ptPid != interprocessPointPidDetection ) // special treatment for the first interprocess detection
                    {
                        mapPointInterProcess[ptId].insert( ptPid );
                        M_mesh->pointIterator( ptId )->second.setProcessId( interprocessPointPidDetection );
                    }
                }
            }
        }
    }

#if 0
    for ( rank_type partId : M_localPartitionIds )
        std::cout << "nActiveElement in part " << partId << " : " << M_containerActiveElements[partId].size() << "\n";
#endif
    // update ghost element containers
    auto allelt_it = M_mesh->beginElement();
    auto allelt_en = M_mesh->endElement();
    for ( ; allelt_it != allelt_en ; ++allelt_it )
    {
        auto const& elt = allelt_it->second;
        rank_type eltPid = elt.processId();
        std::set<rank_type> neighborPids;
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::element_type::numPoints; ++vLocId )
        {
            auto const& thepoint = elt.point( vLocId );
            size_type ptPid = thepoint.processId();
            if ( ptPid != interprocessPointPidDetection )
                continue;

            auto pidsTouchPoint = mapPointInterProcess.find( thepoint.id() )->second;
            neighborPids.insert(pidsTouchPoint.begin(), pidsTouchPoint.end() );
        }
        for ( rank_type neighborPid : neighborPids )
        {
            if ( neighborPid == eltPid ) // not a gost elt
                continue;
            M_containerGhostElements[neighborPid].push_back(boost::cref(elt));
        }
    }
#if 0
    for ( rank_type partId : M_localPartitionIds )
        std::cout << "nGhostElement in part " << partId << " : " << M_containerGhostElements[partId].size() << "\n";
#endif
    // update point containers
    for ( rank_type partId : M_localPartitionIds )
    {
        // store point ids in ghost elements
        std::set<size_type> pointIdsInGhost;
        auto ghostelt_it = M_containerGhostElements[partId].begin();
        auto ghostelt_en = M_containerGhostElements[partId].end();
        for ( ; ghostelt_it != ghostelt_en ; ++ghostelt_it )
        {
            auto const& theghostelt = boost::unwrap_ref(*ghostelt_it);
            for ( uint16_type vLocId = 0 ; vLocId < mesh_type::element_type::numPoints; ++vLocId )
            {
                pointIdsInGhost.insert( theghostelt.point(vLocId).id() );
            }
        }
        // update containers
        auto rangePoints = M_mesh->pointsWithProcessId( partId );
        auto pt_it = std::get<0>( rangePoints );
        auto const pt_en = std::get<1>( rangePoints );
        for ( ; pt_it != pt_en ; ++pt_it )
        {
            M_containerPoints[partId].push_back( boost::cref(*pt_it) );
        }
        for ( size_type ptId : pointIdsInGhost )
        {
            M_containerPoints[partId].push_back( boost::cref(M_mesh->pointIterator( ptId )->second ) );
        }
        std::sort( M_containerPoints[partId].begin(), M_containerPoints[partId].end(),
                   []( auto const& a, auto const& b) -> bool
                   {
                       return unwrap_ref( a ).id() < unwrap_ref( b ).id();
                   });

    }

#if 0
    for ( rank_type partId : M_localPartitionIds )
        std::cout << "nPoints in part " << partId << " : " << M_containerPoints[partId].size() << "\n";
#endif


    this->updateMarkedSubEntitiesAllPartInOneProcess<mesh_type>( mapPointInterProcess );


    for ( rank_type p=0;p<M_numGlobalPartition;++p )
    {
        M_statistic[p].resize( 6 );
        M_statistic[p][0] = M_containerPoints[p].size();
        M_statistic[p][1] = M_containerActiveElements[p].size() + M_containerGhostElements[p].size();
        M_statistic[p][2] = M_containerActiveElements[p].size();
        M_statistic[p][3] = M_containerMarkedFaces[p].size();
        M_statistic[p][4] = M_containerMarkedEdges[p].size();
        M_statistic[p][5] = M_containerMarkedPoints[p].size();
    }

}



template<typename MeshType>
template<typename MT>
void
MeshPartitionSet<MeshType>::updateMarkedSubEntitiesAllPartInOneProcess( std::map<size_type,std::set<rank_type> > const& mapPointInterProcess,
                                                                        typename std::enable_if<is_3d<MT>::value>::type* )
{
    rank_type interprocessPointPidDetection = this->numGlobalPartition();

    auto face_it = M_mesh->beginOrderedFace();
    auto face_en = M_mesh->endOrderedFace();
    for ( ; face_it!=face_en ; ++face_it )
    {
        auto const& theface = unwrap_ref( *face_it );
        if ( !theface.hasMarker() ) continue;
        std::set<rank_type> facePids;
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::face_type::numVertices; ++vLocId )
        {
            auto const& thepoint = theface.point( vLocId );
            rank_type ptPid = thepoint.processId();
            if ( ptPid != interprocessPointPidDetection && ptPid != invalid_rank_type_value )
                facePids.insert( ptPid );
            else if ( ptPid == interprocessPointPidDetection )
            {
                auto itPidsTouchPoint = mapPointInterProcess.find( thepoint.id() );
                if ( itPidsTouchPoint ==  mapPointInterProcess.end() ) continue;
                for ( rank_type neighPid : itPidsTouchPoint->second )
                    facePids.insert( neighPid );
            }
        }
        for ( rank_type partId : facePids )
        {
            M_containerMarkedFaces[partId].push_back(boost::cref(theface));
        }
    }

    auto edge_it = M_mesh->beginOrderedEdge();
    auto edge_en = M_mesh->endOrderedEdge();
    for ( ; edge_it!=edge_en ; ++edge_it )
    {
        auto const& theedge = unwrap_ref( *edge_it );
        if ( !theedge.hasMarker() ) continue;
        std::set<rank_type> edgePids;
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::edge_type::numVertices; ++vLocId )
        {
            auto const& thepoint = theedge.point( vLocId );
            rank_type ptPid = thepoint.processId();
            if ( ptPid != interprocessPointPidDetection && ptPid != invalid_rank_type_value )
                edgePids.insert( ptPid );
            else if ( ptPid == interprocessPointPidDetection )
            {
                auto itPidsTouchPoint = mapPointInterProcess.find( thepoint.id() );
                if ( itPidsTouchPoint ==  mapPointInterProcess.end() ) continue;
                for ( rank_type neighPid : itPidsTouchPoint->second )
                    edgePids.insert( neighPid );
            }
        }
        for ( rank_type partId : edgePids )
        {
            M_containerMarkedEdges[partId].push_back(boost::cref(theedge));
        }
    }

    auto point_it = M_mesh->beginOrderedPoint();
    auto point_en = M_mesh->endOrderedPoint();
    for ( ; point_it!=point_en ; ++point_it )
    {
        auto const& thepoint = unwrap_ref( *point_it );
        if ( !thepoint.hasMarker() ) continue;
        std::set<rank_type> ptPids;
        rank_type ptPid = thepoint.processId();
        if ( ptPid != interprocessPointPidDetection && ptPid != invalid_rank_type_value )
            ptPids.insert( ptPid );
        else if ( ptPid == interprocessPointPidDetection )
        {
            auto itPidsTouchPoint = mapPointInterProcess.find( thepoint.id() );
            if ( itPidsTouchPoint ==  mapPointInterProcess.end() ) continue;
            for ( rank_type neighPid : itPidsTouchPoint->second )
                ptPids.insert( neighPid );
        }
        for ( rank_type partId : ptPids )
        {
            M_containerMarkedPoints[partId].push_back(boost::cref(thepoint));
        }
    }

}
template<typename MeshType>
template<typename MT>
void
MeshPartitionSet<MeshType>::updateMarkedSubEntitiesAllPartInOneProcess( std::map<size_type,std::set<rank_type> > const& mapPointInterProcess,
                                                                        typename std::enable_if<is_2d<MT>::value>::type* )
{
    rank_type interprocessPointPidDetection = this->numGlobalPartition();

    auto face_it = M_mesh->beginOrderedFace();
    auto face_en = M_mesh->endOrderedFace();
    for ( ; face_it!=face_en ; ++face_it )
    {
        auto const& theface = unwrap_ref( *face_it );
        if ( !theface.hasMarker() ) continue;
        std::set<rank_type> facePids;
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::face_type::numVertices; ++vLocId )
        {
            auto const& thepoint = theface.point( vLocId );
            rank_type ptPid = thepoint.processId();
            if ( ptPid != interprocessPointPidDetection && ptPid != invalid_rank_type_value )
                facePids.insert( ptPid );
            else if ( ptPid == interprocessPointPidDetection )
            {
                auto itPidsTouchPoint = mapPointInterProcess.find( thepoint.id() );
                if ( itPidsTouchPoint ==  mapPointInterProcess.end() ) continue;
                for ( rank_type neighPid : itPidsTouchPoint->second )
                    facePids.insert( neighPid );
            }
        }
        for ( rank_type partId : facePids )
        {
            M_containerMarkedFaces[partId].push_back(boost::cref(theface));
        }
    }

    auto point_it = M_mesh->beginOrderedPoint();
    auto point_en = M_mesh->endOrderedPoint();
    for ( ; point_it!=point_en ; ++point_it )
    {
        auto const& thepoint = unwrap_ref( *point_it );
        if ( !thepoint.hasMarker() ) continue;
        std::set<rank_type> ptPids;
        rank_type ptPid = thepoint.processId();
        if ( ptPid != interprocessPointPidDetection && ptPid != invalid_rank_type_value )
            ptPids.insert( ptPid );
        else if ( ptPid == interprocessPointPidDetection )
        {
            auto itPidsTouchPoint = mapPointInterProcess.find( thepoint.id() );
            if ( itPidsTouchPoint ==  mapPointInterProcess.end() ) continue;
            for ( rank_type neighPid : itPidsTouchPoint->second )
                ptPids.insert( neighPid );
        }
        for ( rank_type partId : ptPids )
        {
            M_containerMarkedPoints[partId].push_back(boost::cref(thepoint));
        }
    }

}
template<typename MeshType>
template<typename MT>
void
MeshPartitionSet<MeshType>::updateMarkedSubEntitiesAllPartInOneProcess( std::map<size_type,std::set<rank_type> > const& mapPointInterProcess,
                                                                        typename std::enable_if<is_1d<MT>::value>::type* )
{
    rank_type interprocessPointPidDetection = this->numGlobalPartition();

    auto point_it = M_mesh->beginOrderedPoint();
    auto point_en = M_mesh->endOrderedPoint();
    for ( ; point_it!=point_en ; ++point_it )
    {
        auto const& thepoint = unwrap_ref( *point_it );
        if ( !thepoint.hasMarker() ) continue;
        std::set<rank_type> ptPids;
        rank_type ptPid = thepoint.processId();
        if ( ptPid != interprocessPointPidDetection && ptPid != invalid_rank_type_value )
            ptPids.insert( ptPid );
        else if ( ptPid == interprocessPointPidDetection )
        {
            auto itPidsTouchPoint = mapPointInterProcess.find( thepoint.id() );
            if ( itPidsTouchPoint ==  mapPointInterProcess.end() ) continue;
            for ( rank_type neighPid : itPidsTouchPoint->second )
                ptPids.insert( neighPid );
        }
        for ( rank_type partId : ptPids )
        {
            M_containerMarkedPoints[partId].push_back(boost::cref(thepoint));
        }
    }
}

} // namespace Feel

#endif // FEELPP_MESHPARTITION_HPP
