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
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

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

            auto elt_it = M_mesh->beginElementWithProcessId( partId );
            auto elt_en = M_mesh->endElementWithProcessId( partId );
            for( ; elt_it != elt_en; ++ elt_it )
                M_containerActiveElements[partId].push_back(boost::cref(*elt_it));
            auto ghostelt_it = M_mesh->beginGhostElement();
            auto ghostelt_en = M_mesh->endGhostElement();
            for( ; ghostelt_it != ghostelt_en; ++ghostelt_it )
                M_containerGhostElements[partId].push_back(boost::cref(*ghostelt_it));

            auto pt_it = M_mesh->beginPoint();
            auto pt_en = M_mesh->endPoint();
            for( ; pt_it != pt_en; ++pt_it )
                M_containerPoints[partId].push_back(boost::cref(*pt_it));

            this->updateMarkedSubEntitiesOnePartPerProcess<mesh_type>();
        }

    MeshPartitionSet( mesh_ptrtype const& mesh, rank_type nGlobalPart, std::set<rank_type> const& localPartitionIds )
        :
        M_mesh( mesh ),
        M_numGlobalPartition( nGlobalPart ),
        M_localPartitionIds( localPartitionIds )
        {
            CHECK( M_localPartitionIds.size() <= M_numGlobalPartition  ) << "number of local partition (in process) can not be greater than number of partition in full mesh";
            CHECK( false ) << "TODO";
        }

    MeshPartitionSet( MeshPartitionSet const& ) = default;
    ~MeshPartitionSet() {}

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
    void updateMarkedSubEntitiesOnePartPerProcess( typename std::enable_if<is_3d<MT>::value>::type* = nullptr )
        {
            rank_type partId = M_mesh->worldComm().localRank();
            auto face_it = M_mesh->beginFace();
            auto face_en = M_mesh->endFace();
            for( ; face_it != face_en; ++face_it )
            {
                if ( face_it->marker().isOff() ) continue;
                M_containerMarkedFaces[partId].push_back(boost::cref(*face_it));
            }
            auto edge_it = M_mesh->beginEdge();
            auto edge_en = M_mesh->endEdge();
            for( ; edge_it != edge_en; ++edge_it )
            {
                if ( edge_it->marker().isOff() ) continue;
                M_containerMarkedEdges[partId].push_back(boost::cref(*edge_it));
            }
            auto point_it = M_mesh->beginPoint();
            auto point_en = M_mesh->endPoint();
            for( ; point_it != point_en; ++point_it )
            {
                if ( point_it->marker().isOff() ) continue;
                M_containerMarkedPoints[partId].push_back(boost::cref(*point_it));
            }
        }
    template<typename MT>
    void updateMarkedSubEntitiesOnePartPerProcess( typename std::enable_if<is_2d<MT>::value>::type* = nullptr )
        {
            rank_type partId = M_mesh->worldComm().localRank();
            auto face_it = M_mesh->beginFace();
            auto face_en = M_mesh->endFace();
            for( ; face_it != face_en; ++face_it )
            {
                if ( face_it->marker().isOff() ) continue;
                M_containerMarkedFaces[partId].push_back(boost::cref(*face_it));
            }
            auto point_it = M_mesh->beginPoint();
            auto point_en = M_mesh->endPoint();
            for( ; point_it != point_en; ++point_it )
            {
                if ( point_it->marker().isOff() ) continue;
                M_containerMarkedPoints[partId].push_back(boost::cref(*point_it));
            }
        }
    template<typename MT>
    void updateMarkedSubEntitiesOnePartPerProcess( typename std::enable_if<is_1d<MT>::value>::type* = nullptr )
        {
            rank_type partId = M_mesh->worldComm().localRank();
            auto point_it = M_mesh->beginPoint();
            auto point_en = M_mesh->endPoint();
            for( ; point_it != point_en; ++point_it )
            {
                if ( point_it->marker().isOff() ) continue;
                M_containerMarkedPoints[partId].push_back(boost::cref(*point_it));
            }
        }


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

} // namespace Feel

#endif // FEELPP_MESHPARTITION_HPP
