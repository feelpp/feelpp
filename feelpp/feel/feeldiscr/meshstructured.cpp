/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 04 Jan 2017

 Copyright (C) 2017 Feel++ Consortium

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
// #include <feel/feelconfig.h>
// #include <boost/geometry.hpp>
// #include <boost/geometry/geometries/point_xy.hpp>
// #include <boost/geometry/geometries/polygon.hpp>
#include <feel/feeldiscr/meshstructured.hpp>

namespace Feel
{

namespace detail
{
#if 0
class PolygonMeshStructured
{
public :
    typedef boost::geometry::model::d2::point_xy<double> poly_point_type;
    typedef boost::geometry::model::segment<poly_point_type> segment_type;
    typedef boost::geometry::model::polygon<poly_point_type> polygon_type;

    PolygonMeshStructured( std::string const& pathPoly, double pixelsize )
        {
            polygon_type p;
            std::ifstream flux( pathPoly );
            int nbPolygons;
            int nbVertices;
            double coordX, coordY;
            flux >> nbPolygons;
            //std::cout<< "nP " << nbPolygons <<  std::endl;
            for ( int j = 0; j < nbPolygons; j++ )
            {
                flux >> nbVertices;
                //std::cout<< "nV: " << nbVertices <<  std::endl;
                for ( int i = 0; i < nbVertices; i++ )
                {
                    flux >> coordX;
                    flux >> coordY;
                    p.outer().push_back( poly_point_type( coordX * pixelsize, coordY * pixelsize ) );
                }
                p.outer().push_back( poly_point_type( p.outer().at( 0 ).x(), p.outer().at( 0 ).y() ) );

                M_list_poly.push_back( p );
                p.clear();
            }
        }

    std::vector<polygon_type> const& polygons() const { return M_list_poly; }

private :
    std::vector<polygon_type> M_list_poly;
};
#endif

template <int Dim, typename IndexT>
class PartitioningMeshStructured
{
public:
    static constexpr int nDim = Dim;

    using index_type = IndexT;

    PartitioningMeshStructured(std::vector<index_type> const& nGlobalPtByAxis, rank_type nPartition ,rank_type partId )
        :
        M_nGlobalPtByAxis( nGlobalPtByAxis )
        {
            CHECK( M_nGlobalPtByAxis.size() == nDim ) << "invalid dim";

            M_nLocalPtByAxis.resize( nPartition, std::vector<index_type>(nDim,0) );
            M_startPtIdByAxis.resize( nPartition, std::vector<index_type>(nDim,0) );
            M_nLocalPtByAxisWithGhost.resize( nPartition, std::vector<index_type>(nDim,0) );
            M_startPtIdByAxisWithGhost.resize( nPartition, std::vector<index_type>(nDim,0) );


            std::string partitioningType = "striped";
            int partitioningStripedAxis = 0;

            rank_type nPartNonEmpty = std::min( M_nGlobalPtByAxis[partitioningStripedAxis]-1, (index_type)nPartition );
            bool currentPartIsEmpty = partId >= nPartNonEmpty;

            for ( int a=0;a<nDim;++a )
            {
                if ( a == partitioningStripedAxis )
                {
                    index_type nTotalPointsByCol = M_nGlobalPtByAxis[a] + ( nPartNonEmpty - 1 );
                    index_type nPtByColByProc = ( M_nGlobalPtByAxis[a] + ( nPartNonEmpty - 1 ) ) / nPartNonEmpty;
                    if ( nTotalPointsByCol <= nPtByColByProc * ( nPartNonEmpty - 1 ) )
                        --nPtByColByProc;
                    for ( rank_type p = 0 ; p<nPartNonEmpty ; ++p )
                        M_nLocalPtByAxis[p][a] = nPtByColByProc;
                    if ( nPtByColByProc * nPartNonEmpty < nTotalPointsByCol )
                    {
                        index_type nPointToDistribute = nTotalPointsByCol - nPtByColByProc * nPartNonEmpty;
                        for ( index_type k = 0; k < nPointToDistribute; ++k )
                        {
                            if ( k < nPartNonEmpty )
                                ++M_nLocalPtByAxis[nPartNonEmpty - 1 - k][a];
                        }
                    }

                    for ( rank_type p = 1 ; p<nPartNonEmpty ; ++p )
                    {
                        M_startPtIdByAxis[p][a] = M_nLocalPtByAxis[p-1][a] > 0? M_startPtIdByAxis[p-1][a] + (M_nLocalPtByAxis[p-1][a] -1 ) : M_startPtIdByAxis[p-1][a];
                        M_startPtIdByAxisWithGhost[p][a] = M_startPtIdByAxis[p][a] - 1;
                    }

                    for ( rank_type p = 0 ; p<nPartNonEmpty; ++p )
                    {
                        M_nLocalPtByAxisWithGhost[p][a] = M_nLocalPtByAxis[p][a];
                        if ( p < (nPartNonEmpty-1) )
                            ++M_nLocalPtByAxisWithGhost[p][a];
                        if ( p > 0 )
                            ++M_nLocalPtByAxisWithGhost[p][a];
                    }
                }
                else
                {
                    for ( rank_type p = 0 ; p<nPartNonEmpty ; ++p )
                    {
                        //startPtIdByAxis[p][a] = 0;
                        M_nLocalPtByAxis[p][a] = M_nGlobalPtByAxis[a];
                        M_nLocalPtByAxisWithGhost[p][a] = M_nGlobalPtByAxis[a];
                    }
                }
            }

            // ghost elements
            if ( partId > 0 && !currentPartIsEmpty )
            {
                std::vector<index_type> startEltGhost( nDim, 0 );
                std::vector<index_type> endEltGhost( nDim, 0 );
                for ( int a=0;a<nDim;++a )
                {
                    if ( a == partitioningStripedAxis )
                    {
                        startEltGhost[a] = M_startPtIdByAxisWithGhost[partId][a];
                        endEltGhost[a] =  startEltGhost[a] + 1;
                    }
                    else
                    {
                        startEltGhost[a] = 0;
                        endEltGhost[a] = M_nGlobalPtByAxis[a]-1;
                    }
                }
                M_eltGhosts.emplace( std::make_pair( partId-1, std::make_tuple( std::move(startEltGhost),std::move(endEltGhost) ) ) );
            }
            if ( partId < (nPartNonEmpty-1) && !currentPartIsEmpty )
            {
                std::vector<index_type> startEltGhost( nDim, 0 );
                std::vector<index_type> endEltGhost( nDim, 0 );
                for ( int a=0;a<nDim;++a )
                {
                    if ( a == partitioningStripedAxis )
                    {
                        startEltGhost[a] = M_startPtIdByAxis[partId][a] + (M_nLocalPtByAxis[partId][a] - 1);
                        endEltGhost[a] = startEltGhost[a] +1;
                    }
                    else
                    {
                        startEltGhost[a] = 0;
                        endEltGhost[a] = M_nGlobalPtByAxis[a]-1;
                    }
                }
                M_eltGhosts.emplace( std::make_pair( partId+1, std::make_tuple( std::move(startEltGhost),std::move(endEltGhost) ) ) );
            }

        }

    bool pointIsGhost( rank_type partId, std::array<index_type,nDim> const& indexes ) const
        {
            if constexpr ( nDim == 1 )
            {
                index_type i = indexes[0];
                return
                    ( i < M_startPtIdByAxis[partId][0] ) ||
                    ( i >= (M_startPtIdByAxis[partId][0]+M_nLocalPtByAxis[partId][0]) );
            }
            else if constexpr ( nDim == 2 )
            {
                index_type i = indexes[0];
                index_type j = indexes[1];
                return
                    ( i < M_startPtIdByAxis[partId][0] ) ||
                    ( i >= (M_startPtIdByAxis[partId][0]+M_nLocalPtByAxis[partId][0]) ) ||
                    ( j < M_startPtIdByAxis[partId][1] ) ||
                    ( j >= (M_startPtIdByAxis[partId][1]+M_nLocalPtByAxis[partId][1]) );
            }
            else
            {
                index_type i = indexes[0];
                index_type j = indexes[1];
                index_type k = indexes[2];
                return
                    ( i < M_startPtIdByAxis[partId][0] ) ||
                    ( i >= (M_startPtIdByAxis[partId][0]+M_nLocalPtByAxis[partId][0]) ) ||
                    ( j < M_startPtIdByAxis[partId][1] ) ||
                    ( j >= (M_startPtIdByAxis[partId][1]+M_nLocalPtByAxis[partId][1]) );
                    ( k < M_startPtIdByAxis[partId][2] ) ||
                    ( k >= (M_startPtIdByAxis[partId][2]+M_nLocalPtByAxis[partId][2]) );
            }
        }

    std::vector<index_type> const& nLocalPtByAxis(rank_type p) const { return M_nLocalPtByAxis.at( p ); }
    std::vector<index_type> const& startPtIdByAxis(rank_type p) const { return M_startPtIdByAxis.at( p ); }
    std::vector<index_type> const& nLocalPtByAxisWithGhost(rank_type p) const { return M_nLocalPtByAxisWithGhost.at( p ); }
    std::vector<index_type> const& startPtIdByAxisWithGhost(rank_type p) const { return M_startPtIdByAxisWithGhost.at( p ); }

    std::map<rank_type,std::tuple<std::vector<index_type>,std::vector<index_type>>> const& eltGhost() const { return M_eltGhosts; }
private:
    std::vector<index_type> M_nGlobalPtByAxis;
    std::vector<std::vector<index_type>> M_nLocalPtByAxis;
    std::vector<std::vector<index_type>> M_startPtIdByAxis;
    std::vector<std::vector<index_type>> M_nLocalPtByAxisWithGhost;
    std::vector<std::vector<index_type>> M_startPtIdByAxisWithGhost;

    std::map<rank_type,std::tuple<std::vector<index_type>,std::vector<index_type>>> M_eltGhosts;

};
} // namespace detail


template <typename GeoShape, typename T, typename IndexT>
MeshStructured<GeoShape,T,IndexT>::MeshStructured( std::shared_ptr<setup_type> const& setup, worldcomm_ptr_t const& wc )
    :
    super( wc ),
    M_setup( setup )
{
    if ( M_setup )
        this->generateStructuredMesh();
}

template <typename GeoShape, typename T, typename IndexT>
void
MeshStructured<GeoShape,T,IndexT>::generateStructuredMesh()
{

    if ( M_setup->isCartesian() )
        this->setStructureProperty( "00110" );
    else
        this->setStructureProperty( "00010" );

    rank_type nProc = this->worldComm().localSize();
    rank_type partId = this->worldComm().localRank();

    std::unordered_map<size_type, std::tuple<size_type, rank_type>> mapGhostElt;
    std::unordered_map<size_type, size_type> idStructuredMeshToFeelMesh;
    //node_type coords( 2 );

    std::vector<index_type> M_nGlobalPtByAxis = M_setup->nPointByAxis();

    Feel::detail::PartitioningMeshStructured<nDim,index_type> partitionTool(M_nGlobalPtByAxis, nProc, partId);
    auto const& startPtIdByAxisOnProc = partitionTool.startPtIdByAxis(partId);
    auto const& nLocalPtByAxisOnProc = partitionTool.nLocalPtByAxis(partId);
    auto const& startPtIdByAxisWithGhostOnProc = partitionTool.startPtIdByAxisWithGhost(partId);
    auto const& nLocalPtByAxisWithGhostOnProc = partitionTool.nLocalPtByAxisWithGhost(partId);

    // points
    if constexpr ( nDim == 1 )
    {
        for ( index_type i = startPtIdByAxisWithGhostOnProc[0]; i < (startPtIdByAxisWithGhostOnProc[0] + nLocalPtByAxisWithGhostOnProc[0]); ++i )
        {
            bool currentPtIsGhost = partitionTool.pointIsGhost( partId, {i} );
            this->addStructuredPoint( {i}, partId, currentPtIsGhost );
        }
    }
    else if constexpr ( nDim == 2 )
    {
        for ( index_type i = startPtIdByAxisWithGhostOnProc[0]; i < (startPtIdByAxisWithGhostOnProc[0] + nLocalPtByAxisWithGhostOnProc[0]); ++i )
        {
            for ( index_type j = startPtIdByAxisWithGhostOnProc[1]; j < (startPtIdByAxisWithGhostOnProc[1] + nLocalPtByAxisWithGhostOnProc[1]); ++j )
            {
                bool currentPtIsGhost = partitionTool.pointIsGhost( partId, {i,j} );
                this->addStructuredPoint( {i,j}, partId, currentPtIsGhost );
            }
        }
    }
    else
    {
        for ( index_type i = startPtIdByAxisWithGhostOnProc[0]; i < (startPtIdByAxisWithGhostOnProc[0] + nLocalPtByAxisWithGhostOnProc[0]); ++i )
        {
            for ( index_type j = startPtIdByAxisWithGhostOnProc[1]; j < (startPtIdByAxisWithGhostOnProc[1] + nLocalPtByAxisWithGhostOnProc[1]); ++j )
            {
                for ( index_type k = startPtIdByAxisWithGhostOnProc[2]; k < (startPtIdByAxisWithGhostOnProc[2] + nLocalPtByAxisWithGhostOnProc[2]); ++k )
                {
                    bool currentPtIsGhost = partitionTool.pointIsGhost( partId, {i,j,k} );
                    this->addStructuredPoint( {i,j,k}, partId, currentPtIsGhost );
                }
            }
        }
    }

    // active elements
    std::vector<size_type> endEltIdByAxisOnProc(nDim,0);
    for ( int a=0;a<nDim;++a )
        endEltIdByAxisOnProc[a] = nLocalPtByAxisOnProc[a]>0? startPtIdByAxisOnProc[a] + (nLocalPtByAxisOnProc[a] - 1) : startPtIdByAxisOnProc[a];
    if constexpr ( nDim == 1 )
    {
        for ( index_type i = startPtIdByAxisOnProc[0]; i < endEltIdByAxisOnProc[0]; ++i )
        {
            auto [eid,eidFeel] = this->addStructuredElement({i},partId,partId,{} );
            idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eidFeel ) );
        }
    }
    else if constexpr ( nDim == 2 )
    {
        for ( index_type i = startPtIdByAxisOnProc[0]; i < endEltIdByAxisOnProc[0]; ++i )
        {
            for ( index_type j = startPtIdByAxisOnProc[1]; j < endEltIdByAxisOnProc[1]; ++j )
            {
                auto [eid,eidFeel] = this->addStructuredElement({i,j},partId,partId,{} );
                idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eidFeel ) );
            }
        }
    }
    else
    {
        for ( index_type i = startPtIdByAxisOnProc[0]; i < endEltIdByAxisOnProc[0]; ++i )
        {
            for ( index_type j = startPtIdByAxisOnProc[1]; j < endEltIdByAxisOnProc[1]; ++j )
            {
                for ( index_type k = startPtIdByAxisOnProc[2]; k < endEltIdByAxisOnProc[2]; ++k )
                {
                    auto [eid,eidFeel] = this->addStructuredElement({i,j,k},partId,partId,{} );
                    idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eidFeel ) );
                }
            }
        }
    }

    // ghost elements
    for ( auto const& [partIdGhost, ghostEltsDesc] : partitionTool.eltGhost() )
    {
        if constexpr ( nDim == 1 )
        {
            index_type start_i = std::get<0>( ghostEltsDesc )[0];
            index_type end_i = std::get<1>( ghostEltsDesc )[0];
            for ( index_type i = start_i; i < end_i; ++i )
            {
                auto [eid,eidFeel] = this->addStructuredElement({i},partId,partIdGhost,{} );
                idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eidFeel ) );
                mapGhostElt.insert( std::make_pair( eid, std::make_tuple( eidFeel, partIdGhost ) ) );
            }
        }
        else if constexpr ( nDim == 2 )
        {
            index_type start_i = std::get<0>( ghostEltsDesc )[0];
            index_type start_j = std::get<0>( ghostEltsDesc )[1];
            index_type end_i = std::get<1>( ghostEltsDesc )[0];
            index_type end_j = std::get<1>( ghostEltsDesc )[1];
            for ( index_type i = start_i; i < end_i; ++i )
            {
                for ( index_type j = start_j; j < end_j; ++j )
                {
                    auto [eid,eidFeel] = this->addStructuredElement({i,j},partId,partIdGhost,{} );
                    idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eidFeel ) );
                    mapGhostElt.insert( std::make_pair( eid, std::make_tuple( eidFeel, partIdGhost ) ) );
                }
            }
        }
        else
        {
            index_type start_i = std::get<0>( ghostEltsDesc )[0];
            index_type start_j = std::get<0>( ghostEltsDesc )[1];
            index_type start_k = std::get<0>( ghostEltsDesc )[2];
            index_type end_i = std::get<1>( ghostEltsDesc )[0];
            index_type end_j = std::get<1>( ghostEltsDesc )[1];
            index_type end_k = std::get<1>( ghostEltsDesc )[2];
            for ( index_type i = start_i; i < end_i; ++i )
            {
                for ( index_type j = start_j; j < end_j; ++j )
                {
                    for ( index_type k = start_k; k < end_k; ++k )
                    {
                        auto [eid,eidFeel] = this->addStructuredElement({i,j,k},partId,partIdGhost,{} );
                        idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eidFeel ) );
                        mapGhostElt.insert( std::make_pair( eid, std::make_tuple( eidFeel, partIdGhost ) ) );
                    }
                }
            }
        }
    }

    this->updateGhostCellInfo( idStructuredMeshToFeelMesh, mapGhostElt );
}

template <typename GeoShape, typename T, typename IndexT>
void
MeshStructured<GeoShape,T,IndexT>::addStructuredPoint( std::array<index_type,nDim> const& indexes, rank_type partId, bool isGhost )
{
    index_type ptid = M_setup->pointId( indexes );
    node_type coords = M_setup->pointCoordinates( indexes );
    point_type pt( ptid, coords );
    pt.setProcessId( partId );
    pt.setProcessIdInPartition( partId );
    this->addPoint( pt );
}

template <typename GeoShape, typename T, typename IndexT>
std::pair<typename MeshStructured<GeoShape,T,IndexT>::size_type,typename MeshStructured<GeoShape,T,IndexT>::size_type>
MeshStructured<GeoShape,T,IndexT>::addStructuredElement( std::array<index_type,nDim> const& indexes, rank_type processId, rank_type partId,
                                                         std::vector<rank_type> const& neighborPartitionIds )
{
    index_type eid = M_setup->cellId( indexes ); // StructuredMesh Id
    element_type e;
    e.setMarker( 1, 1 );
    e.setProcessIdInPartition( processId );
    e.setProcessId( partId );

    auto ptids = M_setup->pointIdsFromCellId( indexes );
    for ( uint16_type k = 0; k < element_type::numPoints; ++k )
        e.setPoint( k, this->point( ptids[k] ) );

    if ( neighborPartitionIds.empty() )
        e.setNeighborPartitionIds( neighborPartitionIds );

    auto [eit,inserted] = this->addElement( e, true ); // e.id() is defined by Feel++
    return std::make_pair( eid, eit->second.id() );
}

template <typename GeoShape, typename T, typename IndexT>
void
MeshStructured<GeoShape,T,IndexT>::updateGhostCellInfo( std::unordered_map<size_type, size_type> const& idStructuredMeshToFeelMesh,
                                                        std::unordered_map<size_type, std::tuple<size_type, rank_type>> const& mapGhostElt )
{
    const rank_type myRank = this->worldComm().localRank();
    const rank_type nProc = this->worldComm().localSize();

    std::map<rank_type, std::vector<size_type> > dataToSend;
    std::map<rank_type, std::vector<size_type> > dataToRecv;
    std::map<rank_type, std::vector<std::reference_wrapper<element_type>> > memory;

    for ( auto const& [eidMeshStructured,ghostData] : mapGhostElt )
    {
        size_type ghostEltIdFeel = std::get<0>( ghostData );
        rank_type ghostPid = std::get<1>( ghostData );
        dataToSend[ghostPid].push_back( eidMeshStructured );
        memory[ghostPid].push_back( std::ref( this->elementIterator( ghostEltIdFeel )->second ) );
    }

    int nbMaxRequest = 2*dataToSend.size();
    std::vector<mpi::request> reqs( nbMaxRequest );
    int countRequest = 0;

    // get size of data to transfer
    std::map<rank_type,std::size_t> sizeRecv;
    std::map<rank_type,std::size_t> sizeSend;
    for ( auto const& [neighborRank,data]: dataToSend )
    {
        sizeSend[neighborRank] = data.size();
        reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, sizeSend[neighborRank] );
        reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
    }
    // wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    // step 1 :send/recv of data
    for ( auto const& [neighborRank,data]: dataToSend )
    {
        std::size_t nSendData = dataToSend[neighborRank].size();
        if ( nSendData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, dataToSend[neighborRank].data(), nSendData );
        std::size_t nRecvData = sizeRecv[neighborRank];
        dataToRecv[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, dataToRecv[neighborRank].data(), nRecvData );
    }
    // step 1 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    std::map<rank_type, std::vector<size_type> > dataToSendStep2;
    std::map<rank_type, std::vector<size_type> > dataToRecvStep2;
    for ( auto const& [neighborRank,eids] : dataToRecv )
    {
        dataToSendStep2[neighborRank].resize( eids.size() );
        for ( int k = 0; k < eids.size(); ++k )
        {
            size_type eidFeel = idStructuredMeshToFeelMesh.at(eids[k]);
            dataToSendStep2[neighborRank][k] = eidFeel;
            this->elementIterator( eidFeel )->second.addNeighborPartitionId( neighborRank );
        }
    }

    // step 2 :send/recv of data
    for ( auto const& [neighborRank,data]: dataToSend )
    {
        std::size_t nSendData = dataToSendStep2[neighborRank].size();
        if ( nSendData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().isend( neighborRank, 0, dataToSendStep2[neighborRank].data(), nSendData );

        std::size_t nRecvData = sizeSend[neighborRank]; // use inverse size
        dataToRecvStep2[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[countRequest++] = this->worldComm().localComm().irecv( neighborRank, 0, dataToRecvStep2[neighborRank].data(), nRecvData );
    }
    // step 2 :wait all requests
    mpi::wait_all( std::begin(reqs), std::begin(reqs) + countRequest );
    countRequest = 0;

    // uodate element multiprocess mapping
    for ( auto const& [neighborRank,dataRecv] : dataToRecvStep2 )
    {
        std::size_t nDataRecv = dataRecv.size();
        for ( std::size_t k = 0; k < nDataRecv; ++k )
            memory[neighborRank][k].get().setIdInOtherPartitions( neighborRank, dataRecv[k] );
    }
}



template class MeshStructured<Hypercube<1>>;
template class MeshStructured<Hypercube<2>>;
template class MeshStructured<Hypercube<3>>;

} // namespace Feel
