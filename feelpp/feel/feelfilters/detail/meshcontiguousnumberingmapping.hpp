//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 25 Oct 2017
//! @copyright 2019 Feel++ Consortium

#ifndef FEELPP_FILTERS_DETAILS_MESHCONTIGUOUSNUMBERINGMAPPING_HPP
#define FEELPP_FILTERS_DETAILS_MESHCONTIGUOUSNUMBERINGMAPPING_HPP 1

namespace Feel
{
namespace detail
{

template <typename MeshType,typename StorageNodeValueType>
struct MeshContiguousNumberingMapping
{
    using mesh_type = MeshType;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using index_type = typename mesh_type::index_type;
    using storage_node_value_type = StorageNodeValueType;
    using range_element_type = elements_reference_wrapper_t<mesh_type>;
    using point_ref_type = boost::reference_wrapper< typename mesh_type::point_type const>;

    explicit MeshContiguousNumberingMapping( mesh_type* mesh, bool interprocessPointAreDuplicated = false )
        :
        M_mesh( mesh ),
        M_interprocessPointAreDuplicated( interprocessPointAreDuplicated )
        {
            this->updateForUse();
        }

    void updateForUse()
        {
            mesh_type* mesh = M_mesh;
            rank_type currentPid = mesh->worldComm().localRank();
            rank_type worldSize = mesh->worldComm().localSize();

            if ( M_partIdToRangeElement.empty() )
            {
                std::map<int,int> collectionOfMarkersFlag;
                auto const en_part = mesh->endParts();
                for ( auto it_part = mesh->beginParts() ; it_part!=en_part;++it_part )
                    collectionOfMarkersFlag[it_part->first] = it_part->first;

                auto allRanges = collectionOfMarkedelements( mesh, collectionOfMarkersFlag );
                for ( auto const& [part,rangeElt] : allRanges )
                {
                    std::string markerName = mesh->markerName( part );
                    if ( markerName.empty() || !mesh->hasElementMarker( markerName ) )
                        markerName = "";
                    M_partIdToRangeElement[part] = std::make_tuple(markerName, rangeElt );
                }
            }

            // point id -> (  ( map of idsInOtherPart ), ( vector of ( marker, element id, id in elt) ) )
            std::unordered_map<index_type, std::tuple< std::map<rank_type, index_type>, std::vector< std::tuple<int,index_type,uint16_type>> >> dataPointsInterProcess;

            for ( auto const& [part,nameAndRangeElt] : M_partIdToRangeElement )
            {
                auto const& rangeElt = std::get<1>( nameAndRangeElt );
                index_type nEltInRange = nelements(rangeElt);
                auto & elementIdToContiguous = M_elementIdToContiguous[part];
                auto & pointIdsInElements = M_pointIdsInElements[part];
                auto & pointIdToContiguous = M_pointIdToContiguous[part];
                pointIdsInElements.resize( nEltInRange*mesh_type::element_type::numPoints, invalid_v<index_type> );
                index_type countPtId = 0, countEltId = 0;
                for ( auto const& eltWrap : rangeElt )
                {
                    auto const& elt = unwrap_ref( eltWrap );
                    index_type eltId = elt.id();
                    auto [itElt,eltIsInserted] = elementIdToContiguous.try_emplace( eltId, countEltId++ );
                    DCHECK( eltIsInserted ) << "something wrong, element already inserted";
                    index_type newEltId = itElt->second;
                    for ( uint16_type j = 0; j < mesh_type::element_type::numPoints; j++ )
                    {
                        auto const& pt = elt.point( j );
                        index_type ptid = pt.id();
                        auto const& ptIdnOthersPartitions = pt.idInOthersPartitions();
                        if ( M_interprocessPointAreDuplicated || ptIdnOthersPartitions.empty() ) // not a interprocess point
                        {
                            auto [itPt,isInserted] = pointIdToContiguous.try_emplace( ptid,std::make_pair(countPtId,boost::cref(pt)) );
                            if ( isInserted )
                                ++countPtId;
                            index_type newPtId = itPt->second.first;

                            DCHECK( (newEltId*mesh_type::element_type::numPoints+j) < pointIdsInElements.size() ) << "invalid size : " << (newEltId*mesh_type::element_type::numPoints+j) << " vs " <<  pointIdsInElements.size();
                            pointIdsInElements[newEltId*mesh_type::element_type::numPoints+j] = newPtId;
                        }
                        else
                        {
                            auto infoIpElt = std::make_tuple( part, newEltId/*eltId*/, j );
                            auto itFindPtIP = dataPointsInterProcess.find( ptid );
                            if ( itFindPtIP == dataPointsInterProcess.end() )
                            {
                                dataPointsInterProcess.emplace( ptid,  std::make_tuple( ptIdnOthersPartitions,  std::vector< std::tuple<int,index_type,uint16_type>>( { infoIpElt } ) ) );
                            }
                            else
                                std::get<1>( itFindPtIP->second ).push_back( infoIpElt );
                        }
                    }
                }
            }

            // --------------------------------------------------------------------------------------- //
            // treatment of interprocess point
            std::map<int,std::map<rank_type,std::map<index_type,index_type>>> dataPointsNotInProcess;
            if ( !M_interprocessPointAreDuplicated && worldSize > 1 )
            {
                std::map<rank_type, std::vector<std::pair<int,index_type>>> dataToSend;
                std::map<rank_type, std::vector<std::pair<int,index_type>>> dataToRecv;
                for ( auto const& [ ptId, dataIP ] : dataPointsInterProcess )
                {
                    auto const& dataIdInOtherPartition =  std::get<0>( dataIP );
                    auto const& dataEltInfo = std::get<1>( dataIP );
                    for ( auto const& [ procId, ptIdOtherPart ] : dataIdInOtherPartition )
                    {
                        for ( auto const& [marker,eltId,ptIdInElt ] : dataEltInfo )
                            dataToSend[procId].push_back( std::make_pair(marker, ptIdOtherPart) );
                    }
                }
                int neighborSubdomains = mesh->neighborSubdomains().size();
                int nbRequest = 2 * neighborSubdomains;
                mpi::request* reqs = new mpi::request[nbRequest];
                int cptRequest = 0;
                std::map<rank_type,size_type> sizeRecv;

                // get size of data to transfer
                for ( rank_type neighborRank : mesh->neighborSubdomains() )
                {
                    reqs[cptRequest++] = mesh->worldComm().localComm().isend( neighborRank , 0, (size_type)dataToSend[neighborRank].size() );
                    reqs[cptRequest++] = mesh->worldComm().localComm().irecv( neighborRank , 0, sizeRecv[neighborRank] );
                }
                // wait all requests
                mpi::wait_all(reqs, reqs + cptRequest);

                cptRequest = 0;
                for ( rank_type neighborRank : mesh->neighborSubdomains() )
                {
                    int nSendData = dataToSend[neighborRank].size();
                    if ( nSendData > 0 )
                        reqs[cptRequest++] = mesh->worldComm().localComm().isend( neighborRank, 0, &(dataToSend[neighborRank][0]), nSendData );
                    int nRecvData = sizeRecv[neighborRank];
                    dataToRecv[neighborRank].resize( nRecvData );
                    if ( nRecvData > 0 )
                        reqs[cptRequest++] = mesh->worldComm().localComm().irecv( neighborRank, 0, &(dataToRecv[neighborRank][0]), nRecvData );
                }
                // wait all requests
                mpi::wait_all( reqs, reqs + cptRequest );

                std::map<int,std::map<index_type, std::set<rank_type> > > treatRecv;
                // others process
                for ( auto const& [pid, ptData ] : dataToRecv )
                {
                    for ( auto const& [marker, ptId] : ptData )
                        treatRecv[marker][ptId].insert( pid );
                }
                // self process
                for ( auto const& [ ptId, dataIP ] : dataPointsInterProcess )
                {
                    auto const& dataEltInfo = std::get<1>( dataIP );
                    for ( auto const& [marker,eltId,ptIdInElt ] : dataEltInfo )
                        treatRecv[marker][ptId].insert( currentPid );
                }
                // determine which process have points
                std::map<int,std::set<index_type>> ipPointsOnCurrentProcess;
                for ( auto const& [marker, mapIdToPids] : treatRecv )
                {
                    for ( auto const& [ptId,allPids] : mapIdToPids )
                        if ( currentPid == *allPids.begin() )
                            ipPointsOnCurrentProcess[marker].insert( ptId );
                }
                // update M_pointIdToContiguous
                std::map<rank_type,std::vector<boost::tuple<int,index_type,index_type>>> dataToReSend;
                std::map<rank_type,std::vector<boost::tuple<int,index_type,index_type>>> dataToReRecv;
                for ( auto const& [marker, ptIds] : ipPointsOnCurrentProcess )
                {
                    auto & pointIdToContiguous = M_pointIdToContiguous[marker];
                    auto & pointIdsInElements = M_pointIdsInElements[marker];
                    index_type countPtId = pointIdToContiguous.size();
                    for ( index_type ptId : ptIds )
                    {
                        auto const& pt = mesh->point( ptId );
                        auto [itPt,isInserted] = pointIdToContiguous.try_emplace( ptId,std::make_pair(countPtId,boost::cref(pt)) );
                        index_type newPtId = itPt->second.first;
                        if ( isInserted )
                        {
                            ++countPtId;
                            for (auto const& [opid, optId] : pt.idInOthersPartitions() )
                                dataToReSend[opid].push_back( boost::make_tuple( marker,optId,newPtId ) );
                        }

                        auto itFindDataIP = dataPointsInterProcess.find( ptId );
                        CHECK( itFindDataIP != dataPointsInterProcess.end() ) << "point not register";
                        auto const& infoEltAssociated = std::get<1>( itFindDataIP->second );
                        for ( auto const& [ marker2, newEltId, j ] : infoEltAssociated )
                        {
                            if ( marker != marker2 )
                                continue;
                            pointIdsInElements[newEltId*mesh_type::element_type::numPoints+j] = newPtId;
                        }
                    }
                }


                // get size of data to transfer
                cptRequest = 0;
                for ( rank_type neighborRank : mesh->neighborSubdomains() )
                {
                    reqs[cptRequest++] = mesh->worldComm().localComm().isend( neighborRank , 0, (size_type)dataToReSend[neighborRank].size() );
                    reqs[cptRequest++] = mesh->worldComm().localComm().irecv( neighborRank , 0, sizeRecv[neighborRank] );
                }
                // wait all requests
                mpi::wait_all(reqs, reqs + cptRequest);

                cptRequest = 0;
                for ( rank_type neighborRank : mesh->neighborSubdomains() )
                {
                    int nSendData = dataToReSend[neighborRank].size();
                    if ( nSendData > 0 )
                        reqs[cptRequest++] = mesh->worldComm().localComm().isend( neighborRank, 0, &(dataToReSend[neighborRank][0]), nSendData );
                    int nRecvData = sizeRecv[neighborRank];
                    dataToReRecv[neighborRank].resize( nRecvData );
                    if ( nRecvData > 0 )
                        reqs[cptRequest++] = mesh->worldComm().localComm().irecv( neighborRank, 0, &(dataToReRecv[neighborRank][0]), nRecvData );
                }
                // wait all requests
                mpi::wait_all( reqs, reqs + cptRequest );
                // delete reqs because finish comm
                delete[] reqs;


                for ( auto const& [pid, dataByProc] : dataToReRecv )
                {
                    for ( auto const&  dataPt : dataByProc )
                    {
                        int marker = boost::get<0>( dataPt );
                        index_type ptId =  boost::get<1>( dataPt );
                        index_type newPtId =  boost::get<2>( dataPt );
                        dataPointsNotInProcess[marker][pid][ptId] = newPtId;
                    }
                }

            }

            // --------------------------------------------------------------------------------------- //
            // build nodes vector
            for ( auto const& [markerId,pointIdToContiguous] : M_pointIdToContiguous )
            {
                auto & nodes = M_nodes[markerId];
                nodes.resize( 3*pointIdToContiguous.size(),0 );
                for ( auto const& [ptId,ptData] : pointIdToContiguous )
                {
                    index_type newPtId = ptData.first;
                    auto const& pt = unwrap_ref( ptData.second );
                    for ( uint16_type d=0 ; d<mesh_type::nRealDim ;++d )
                    {
                        DCHECK( (3*newPtId+d) < nodes.size() ) << "invalid size : " << (3*newPtId+d) << " vs " << nodes.size() ;
                        nodes[3*newPtId+d] = pt.node()[d];
                    }
                }
            }

            // --------------------------------------------------------------------------------------- //
            // information in world comm
            int k=0;
            std::vector<boost::tuple<index_type,index_type>> nPointElementByMarker( M_pointIdToContiguous.size() );
            for ( auto const& [marker,pointIdToContiguous] : M_pointIdToContiguous )
            {
                nPointElementByMarker[k] = boost::make_tuple( pointIdToContiguous.size(), M_elementIdToContiguous.find( marker)->second.size() );
                ++k;
            }

            std::vector<std::vector<boost::tuple<index_type,index_type>>> recvInfos;
            mpi::all_gather( mesh->worldComm().comm(), nPointElementByMarker, recvInfos );

            k=0;
            for ( auto const& [marker,pointIdToContiguous] : M_pointIdToContiguous )
            {
                auto & numberOfPointElement = M_numberOfPointElement[marker];
                numberOfPointElement.resize( worldSize );
                index_type startPtId = 0, startEltId = 0;
                for ( rank_type p=0;p<worldSize;++p )
                {
                    index_type nPt = boost::get<0>( recvInfos[p][k] );
                    index_type nElt = boost::get<1>( recvInfos[p][k] );
                    numberOfPointElement[p] = std::make_tuple( startPtId, nPt, startEltId, nElt );
                    startPtId+=nPt;
                    startEltId+=nElt;
                }
                ++k;
                M_numberOfPointElementAllProcess[marker] = std::make_tuple( startPtId,startEltId );
            }

            // --------------------------------------------------------------------------------------- //
            // shift ids
            for ( auto const& [marker,pointIdToContiguous] : M_pointIdToContiguous )
            {
                index_type spi = this->startPointIds(marker,currentPid);
                index_type sei = this->startElementIds(marker,currentPid);
                for ( auto & [ptId,newPtData] :  M_pointIdToContiguous[marker] )
                    newPtData.first += spi;
                for ( auto & [eltId,newEltId] : M_elementIdToContiguous[marker] )
                    newEltId += sei;
                for ( index_type & ptId : M_pointIdsInElements[marker] )
                {
                    if ( ptId != invalid_v<index_type> )
                        ptId += spi;
                }
            }

            // --------------------------------------------------------------------------------------- //
            // update M_pointIdsInElements with interprocess points
            for ( auto const& [marker, dataByMarker ] : dataPointsNotInProcess )
            {
                auto & pointIdsInElements = M_pointIdsInElements[marker];
                for( auto const& [pid,dataByProc] : dataByMarker )
                {
                    index_type shiftPointId = this->startPointIds( marker, pid );
                    for ( auto const& [ ptId,newPtId ] : dataByProc )
                    {
                        auto itFindPtIp = dataPointsInterProcess.find( ptId );
                        CHECK( itFindPtIp != dataPointsInterProcess.end() ) << "invalid point ";
                        auto const& infosElt = std::get<1>( itFindPtIp->second );
                        // up pointIdsInElements
                        for (auto const& [ marker2, newEltId, j ] : infosElt )
                        {
                            if ( marker == marker2 )
                                pointIdsInElements[newEltId*mesh_type::element_type::numPoints+j] = newPtId + shiftPointId;
                        }
                    }
                }
            }

        }

    const mesh_type* mesh() const { return M_mesh; }

    std::map<int,std::tuple<std::string,range_element_type>> const& partIdToRangeElement() const { return M_partIdToRangeElement; }

    std::string const& name( int part ) const
        {
            auto itFindPart = M_partIdToRangeElement.find( part );
            CHECK( itFindPart !=  M_partIdToRangeElement.end() ) << "part not registerd";
            return std::get<0>( itFindPart->second );
        }
    range_element_type const& rangeElement( int part ) const
        {
            auto itFindPart = M_partIdToRangeElement.find( part );
            CHECK( itFindPart !=  M_partIdToRangeElement.end() ) << "part not registerd";
            return std::get<1>( itFindPart->second );
        }

    std::unordered_map<index_type,std::pair<index_type,point_ref_type>> const& pointIdToContiguous( int part ) const
        {
            auto itFindData = M_pointIdToContiguous.find( part );
            CHECK( itFindData != M_pointIdToContiguous.end() ) << "part not registerd";
            return itFindData->second;
        }
    index_type pointIdToContiguous( int part, index_type ptId ) const
        {
            auto itFindData = M_pointIdToContiguous.find( part );
            if ( itFindData == M_pointIdToContiguous.end() )
                return invalid_v<index_type>;
            auto const& data =  itFindData->second;
            auto itFindPt = data.find( ptId );
            if (  itFindPt == data.end() )
                return invalid_v<index_type>;
            return itFindPt->second.first;
        }
    std::unordered_map<index_type,index_type> const& elementIdToContiguous( int part ) const
        {
            auto itFindData = M_elementIdToContiguous.find( part );
            CHECK( itFindData == M_elementIdToContiguous.end() ) << "part not registerd";
            return itFindData->second;
        }
    index_type elementIdToContiguous( int part, index_type eltId ) const
        {
            auto itFindData = M_elementIdToContiguous.find( part );
            if ( itFindData == M_elementIdToContiguous.end() )
                return invalid_v<index_type>;
            auto const& data =  itFindData->second;
            auto itFindElt = data.find( eltId );
            if (  itFindElt == data.end() )
                return invalid_v<index_type>;
            return itFindElt->second;
        }
    std::vector<index_type> const& pointIdsInElements( int part ) const
        {
            auto itFindPointIdsInElements =  M_pointIdsInElements.find( part );
            CHECK( itFindPointIdsInElements != M_pointIdsInElements.end() ) << "part not registerd";
            return itFindPointIdsInElements->second;
        }
    std::vector<storage_node_value_type> const& nodes( int part ) const
        {
            auto itFindNodes =  M_nodes.find( part );
            CHECK( itFindNodes != M_nodes.end() ) << "part not registerd";
            return itFindNodes->second;
        }

    index_type startPointIds( int part, rank_type therank ) const { return genericInfo<0>( part, therank ); }
    index_type numberOfPoint( int part, rank_type therank ) const { return genericInfo<1>( part, therank ); }
    index_type startElementIds( int part, rank_type therank ) const { return genericInfo<2>( part, therank ); }
    index_type numberOfElement( int part, rank_type therank ) const { return genericInfo<3>( part, therank ); }

    index_type numberOfPointAllProcess( int part ) const { return genericInfoAllProcess<0>( part ); }
    index_type numberOfElementAllProcess( int part ) const { return genericInfoAllProcess<1>( part ); }

    void updateNodesCoordinates()
        {
            rank_type currentPid = M_mesh->worldComm().localRank();
            for ( auto const& [part,pointIdToContiguous] : M_pointIdToContiguous )
            {
                index_type spi = this->startPointIds(part,currentPid);
                auto & nodes = M_nodes[part];
                CHECK( nodes.size() == 3*pointIdToContiguous.size() ) << "wrong size";
                //nodes.resize( 3*pointIdToContiguous.size(),0 );
                for ( auto const& [ptId,ptData] : pointIdToContiguous )
                {
                    index_type newPtId = ptData.first - spi;
                    auto const& pt = unwrap_ref( ptData.second );
                    for ( uint16_type d=0 ; d<mesh_type::nRealDim ;++d )
                    {
                        CHECK( (3*newPtId+d) < nodes.size() ) << "invalid size : " << (3*newPtId+d) << " vs " << nodes.size() ;
                        nodes[3*newPtId+d] = pt.node()[d];
                    }
                }
            }
        }


    //! reorder the nodes ids in the element and put the new ordering in arg \newPointsIdsInElt
    template <typename TheNodeIndexType>
    void updateOrderingOfPointsIdsInElt( int part, rank_type therank, std::vector<TheNodeIndexType> & newPointsIdsInElt,
                                         std::vector<uint16_type> const& mappingWithThisKindOfElement, int shiftId = 0,
                                         int nPointsUsedInElt = mesh_type::element_type::numPoints ) const
        {
            CHECK( nPointsUsedInElt <= mappingWithThisKindOfElement.size() ) << "incomplete ordering";

            index_type _nElt = this->numberOfElement( part,therank );
            auto const& pointsIdsInElt_B = this->pointIdsInElements( part );
            newPointsIdsInElt.resize( _nElt*nPointsUsedInElt );

            for ( int k=0;k<_nElt;++k )
            {
                for ( uint16_type p=0;p<nPointsUsedInElt;++p )
                    newPointsIdsInElt[ k*nPointsUsedInElt + mappingWithThisKindOfElement[p] ] = pointsIdsInElt_B[ k*mesh_type::element_type::numPoints+p ] + shiftId;
            }
        }
private :

    template <int TupleId>
    index_type genericInfo( int part, rank_type therank ) const
        {
            auto itFindNumberOfPointElement = M_numberOfPointElement.find( part );
            CHECK( itFindNumberOfPointElement!= M_numberOfPointElement.end()) << "invalid part";
            CHECK( therank < itFindNumberOfPointElement->second.size() ) << "invalid rank";
            return std::get<TupleId>( itFindNumberOfPointElement->second[therank] );
        }
    template <int TupleId>
    index_type genericInfoAllProcess( int part ) const
        {
            auto itFindNumberOfPointElementAllProcess = M_numberOfPointElementAllProcess.find( part );
            CHECK( itFindNumberOfPointElementAllProcess != M_numberOfPointElementAllProcess.end()) << "invalid part";
            return std::get<TupleId>( itFindNumberOfPointElementAllProcess->second );
        }
private:
    mesh_type* M_mesh;
    bool M_interprocessPointAreDuplicated;
    std::map<int,std::tuple<std::string,range_element_type>> M_partIdToRangeElement;
    std::map<int,std::unordered_map<index_type,std::pair<index_type,point_ref_type> >> M_pointIdToContiguous;
    std::map<int,std::unordered_map<index_type,index_type>> M_elementIdToContiguous;
    std::map<int,std::vector<index_type>> M_pointIdsInElements;
    std::map<int,std::vector<storage_node_value_type>> M_nodes;
    std::map<int,std::vector<std::tuple<index_type,index_type,index_type,index_type>>> M_numberOfPointElement;
    std::map<int,std::tuple<index_type,index_type>> M_numberOfPointElementAllProcess;
};



template <typename T>
struct MeshPoints
{
    template <typename MeshType, typename IteratorType>
    MeshPoints( MeshType* mesh, const WorldComm&, IteratorType it, IteratorType en, const bool outer = false, const bool renumber = false, const bool fill = false, const int startIndex = 1 );

    int translatePointIds( std::vector<int32_t>& ids );
    int translateElementIds( std::vector<int32_t>& ids );

    int globalNumberOfPoints() const { return global_npts; }
    int globalNumberOfElements() const { return global_nelts; }

    std::vector<int> numberOfPoints, numberOfElements;
    int global_nelts{0}, global_npts{0};
    std::vector<int32_t> ids;
    std::unordered_map<int32_t, int32_t> new2old;
    std::unordered_map<int32_t, int32_t> old2new;
    std::unordered_map<int32_t, int32_t> nodemap;
    std::vector<T> coords;
    std::vector<int32_t> elemids;
    std::vector<int32_t> elem;
    size_type offsets_pts, global_offsets_pts;
    size_type offsets_elts, global_offsets_elts;
};

//!
//!  Builds information around faces/elements for exporting data
//!  @param mesh The mesh from which data is extracted
//!  @param it Starting iterator over the faces/elements
//!  @param en Endoing iterator over the faces/elements
//!  @param outer If false, the vertices are place in an x1 y1 z1 ... xn yn zn order, otherwise in the x1 ... xn y1 ... yn z1 ... zn
//!  @param renumber If true, the vertices will be renumbered with maps to keep the correspondance between the twoi, otherwise the original ids are kept
//!  @param fill It true, the method will generate points coordinates that are 3D, even if the point is specified with 1D or 2D coordinates (filled with 0)
//!  @param Specify the startIndex of the renumbered points (typically set to 0 or 1, but no restriction). This is only used when renumber is true, otherwise it is not used.
//!
template <typename T>
template <typename MeshType, typename IteratorType>
MeshPoints<T>::MeshPoints( MeshType* mesh, const WorldComm& worldComm, IteratorType it, IteratorType en, const bool outer, const bool renumber, const bool fill, const int startIndex )
{
    std::set<int> nodeset;
    size_type p = 0;
    auto elt_it = it;

    //!  Gather all the vertices of which the elements are made up with into a std::set */
    //!  build up correspondance arrays between index in nodeset and previous id */
    for ( auto eit = it; eit != en; ++eit )
    {
        auto const& elt = boost::unwrap_ref( *eit );
        for ( size_type j = 0; j < MeshType::element_type::numPoints; j++ )
        {
            int pid = elt.point( j ).id();
            auto ins = nodeset.insert( pid );
            if ( ins.second )
            {
                if ( renumber )
                {
                    ids.push_back( p + startIndex );
                }
                else
                {
                    ids.push_back( pid );
                }
                //!  old id -> new id */
                old2new[pid] = ids[p];
                //!  old id -> new id */
                new2old[ids[p]] = pid;
                //!  old id -> index of the new id */
                nodemap[pid] = p;
                ++p;
            }
        }
    }
    CHECK( p == ids.size() ) << "Invalid number of points " << ids.size() << "!=" << p;
    int nv = ids.size();

    coords.resize( 3 * nv, 0 );

    auto pit = ids.begin();
    auto pen = ids.end();
    //! for( auto i = 0; i < nv; ++i )

    //!  put coords of each point into the coords array */
    //!  if outer is true, the coords are placed like: x1 x2 ... xn y1 y2 ... yn z1 z2 ... zn */
    //!  otherwise, the coords are placed like: x1 y1 z1 x2 y2 z2 ... xn yn zn */
    for ( int i = 0; pit != pen; ++pit, ++i )
    {
        //! CHECK( *pit > 0 ) << "invalid id " << *pit;
        //! LOG(INFO) << "p " << i << "/" << nv << " =" << *pit;
        //! int pid = (renumber)?nodemap[*pit]+1:*pit;
        int pid = *pit;

        auto const& p = mesh->point( new2old[*pit] );
        if ( outer )
        {
            coords[i] = (T)p.node()[0];
        }
        else
        {
            coords[3 * i] = (T)p.node()[0];
        }

        if ( MeshType::nRealDim >= 2 )
        {
            if ( outer )
            {
                coords[nv + i] = ( T )( p.node()[1] );
            }
            else
            {
                coords[3 * i + 1] = ( T )( p.node()[1] );
            }
        }
        //!  Fill 2nd components with 0 if told to do so */
        else
        {
            if ( fill )
            {
                if ( outer )
                {
                    coords[nv + i] = (T)0;
                }
                else
                {
                    coords[3 * i + 1] = (T)0;
                }
            }
        }

        if ( MeshType::nRealDim >= 3 )
        {
            if ( outer )
            {
                coords[2 * nv + i] = ( T )( p.node()[2] );
            }
            else
            {
                coords[3 * i + 2] = ( T )( p.node()[2] );
            }
        }
        //!  Fill 3nd components with 0 if told to do so */
        else
        {
            if ( fill )
            {
                if ( outer )
                {
                    coords[2 * nv + i] = (T)0;
                }
                else
                {
                    coords[3 * i + 2] = (T)0;
                }
            }
        }
    }

    //!  number of local elements */
    int __ne = std::distance( it, en );

    //!  only do this resize if we have at least one element in the iterator */
    //!  otherwise it will segfault */
    if ( it != en )
    {
        elem.resize( __ne * MeshType::element_type::numPoints );
        //! elem.resize( __ne*mesh->numLocalVertices() );
        elemids.resize( __ne );
    }


    //!  build the array containing the id of each vertex for each element */
    elt_it = it;
    size_type e = 0;
    for ( ; elt_it != en; ++elt_it, ++e )
    {
        auto const& elt = boost::unwrap_ref( *elt_it );
        elemids[e] = elt.id() + 1;
        //! std::cout << "LocalV = " << elt.numLocalVertices << std::endl;
        //! for ( size_type j = 0; j < mesh->numLocalVertices(); j++ )
        for ( size_type j = 0; j < MeshType::element_type::numPoints; j++ )
        {
            //! std::cout << "LocalVId = " << j << " " << e*elt.numLocalVertices+j << std::endl;
            //! std::cout << elt.point( j ).id() << std::endl;
            //!  ensight id start at 1
            elem[e * MeshType::element_type::numPoints + j] = old2new[elt.point( j ).id()];
#if 0
            DCHECK( (elem[e*mesh->numLocalVertices()+j] > 0) && (elem[e*mesh->numLocalVertices()+j] <= nv ) )
                << "Invalid entry : " << elem[e*mesh->numLocalVertices()+j]
                << " at index : " << e*mesh->numLocalVertices()+j
                << " element :  " << e
                << " vertex :  " << j;
#endif
        }
    }
#if 0
    CHECK( e==__ne) << "Invalid number of elements, e= " << e << "  should be " << __ne;
    std::for_each( elem.begin(), elem.end(), [=]( int e )
                   { CHECK( ( e > 0) && e <= __nv ) << "invalid entry e = " << e << " nv = " << nv; } );
#endif


    //!  gather the number of points and elements fo each process */
    std::vector<int> ost{nv, __ne};
    std::vector<std::vector<int>> ospe;

    mpi::all_gather( worldComm.comm(), ost, ospe );

    //!  copy information about number of points/elements
    //!  per process in a local array */
    for ( size_type i = 0; i < ospe.size(); i++ )
    {
        numberOfPoints.push_back( ospe[i][0] );
        numberOfElements.push_back( ospe[i][1] );
    }

    //!  compute offsets to shift the point and element ids */
    //!  regarding to the processor rank */
    offsets_pts = 0;
    global_offsets_pts = 0;
    offsets_elts = 0;
    global_offsets_elts = 0;
    for ( size_type i = 0; i < ospe.size(); i++ )
    {
        if ( i < worldComm.localRank() )
        {
            offsets_pts += ospe[i][0];
            offsets_elts += ospe[i][1];
        }
        global_offsets_pts += ospe[i][0];
        global_offsets_elts += ospe[i][1];
    }
    global_npts = global_offsets_pts;
    global_nelts = global_offsets_elts;

    //!
    //! std::cout << "local offset pts : " << offsets_pts << std::endl;
    //! std::cout << "local offset elts : " << offsets_elts << std::endl;
    //! std::cout << "global offset pts : " << global_offsets_pts << std::endl;
    //! std::cout << "global offset elts : " << global_offsets_elts << std::endl;
    //! std::cout << "done with offsets" << std::endl;
}

//!
//!  Translate the list of points ids to the new global layout
//!  @param ids Array of local point ids to be translated
//!
template <typename T>
int MeshPoints<T>::translatePointIds( std::vector<int32_t>& ptids )
{
    for ( int i = 0; i < ptids.size(); i++ )
    {
        ptids[i] = offsets_pts + old2new[ptids[i]];
    }

    return 0;
}

//!
//!  Translate the list of element ids to the new global layout
//!  @param ids Array of local point ids to be translated
//!
template <typename T>
int MeshPoints<T>::translateElementIds( std::vector<int32_t>& elids )
{
    for ( int i = 0; i < elids.size(); i++ )
    {
        elids[i] = offsets_elts + elids[i];
    }

    return 0;
}



} // namespace detail

} // namespace Feel

#endif
