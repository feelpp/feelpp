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
    static const int nDim = Dim;

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

    bool pointIsGhost( rank_type partId, int i, int j ) const
        {
            return
                ( i < M_startPtIdByAxis[partId][0] ) ||
                ( i >= (M_startPtIdByAxis[partId][0]+M_nLocalPtByAxis[partId][0]) ) ||
                ( j < M_startPtIdByAxis[partId][1] ) ||
                ( j >= (M_startPtIdByAxis[partId][1]+M_nLocalPtByAxis[partId][1]) );
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
MeshStructured<GeoShape,T,IndexT>::MeshStructured( int nx, int ny, double pixelsize,
                                                   std::optional<holo3_image<float>> const& cx,
                                                   std::optional<holo3_image<float>> const& cy,
                                                   worldcomm_ptr_t const& wc,
                                                   bool withCoord )
    : super( wc ),
      M_nx( nx ),
      M_ny( ny ),
      M_cx( cx ),
      M_cy( cy ),
      M_pixelsize( pixelsize )
{
    if ( !M_cx && !M_cy )
        this->setStructureProperty( "00110" );
    else
        this->setStructureProperty( "00010" );
    LOG(INFO) << "nx x ny = " << nx << " x " << ny << "\t" << nx * ny << std::endl;
    CHECK( nx > 1 ) << "number of point in x dir should be greater than 1";
    CHECK( ny > 1 ) << "number of point in y dir should be greater than 1";


    rank_type nProc = wc->localSize();
    rank_type partId = wc->localRank();

    std::unordered_map<size_type, boost::tuple<size_type, rank_type>> mapGhostElt;
    std::unordered_map<size_type, size_type> idStructuredMeshToFeelMesh;
    node_type coords( 2 );


    // std::shared_ptr<Feel::detail::PolygonMeshStructured> polygonTool;
    // if ( withPoly )
    //     polygonTool = std::make_shared<Feel::detail::PolygonMeshStructured>( pathPoly, M_pixelsize );
    // bool inPoly = false;

    std::vector<index_type> M_nGlobalPtByAxis;
    M_nGlobalPtByAxis = { M_nx, M_ny };

    Feel::detail::PartitioningMeshStructured<nDim,index_type> partitionTool(M_nGlobalPtByAxis, nProc, partId);
    auto const& startPtIdByAxisOnProc = partitionTool.startPtIdByAxis(partId);
    auto const& nLocalPtByAxisOnProc = partitionTool.nLocalPtByAxis(partId);
    auto const& startPtIdByAxisWithGhostOnProc = partitionTool.startPtIdByAxisWithGhost(partId);
    auto const& nLocalPtByAxisWithGhostOnProc = partitionTool.nLocalPtByAxisWithGhost(partId);

    // points
    for ( int i = startPtIdByAxisWithGhostOnProc[0]; i < (startPtIdByAxisWithGhostOnProc[0] + nLocalPtByAxisWithGhostOnProc[0]); ++i )
    {
        for ( int j = startPtIdByAxisWithGhostOnProc[1]; j < (startPtIdByAxisWithGhostOnProc[1] + nLocalPtByAxisWithGhostOnProc[1]); ++j )
        {
            bool currentPtIsGhost = partitionTool.pointIsGhost( partId, i, j );
            this->addStructuredPoint( i,j, partId, currentPtIsGhost, withCoord );
        }
    }

    // active elements
    std::vector<size_type> endEltIdByAxisOnProc(nDim,0);
    for ( int a=0;a<nDim;++a )
        endEltIdByAxisOnProc[a] = nLocalPtByAxisOnProc[a]>0? startPtIdByAxisOnProc[a] + (nLocalPtByAxisOnProc[a] - 1) : startPtIdByAxisOnProc[a];

    for ( int i = startPtIdByAxisOnProc[0]; i < endEltIdByAxisOnProc[0]; ++i )
    {
        for ( int j = startPtIdByAxisOnProc[1]; j < endEltIdByAxisOnProc[1]; ++j )
        {
            auto [eid,eidFeel] = this->addStructuredElement(i,j,partId,partId,{},withCoord );
            idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eidFeel ) );
        }
    }

    // ghost elements
    for ( auto const& [partIdGhost, ghostEltsDesc] : partitionTool.eltGhost() )
    {
        size_type start_i = std::get<0>( ghostEltsDesc )[0];
        size_type start_j = std::get<0>( ghostEltsDesc )[1];
        size_type end_i = std::get<1>( ghostEltsDesc )[0];
        size_type end_j = std::get<1>( ghostEltsDesc )[1];
        for ( size_type i = start_i; i < end_i; ++i )
        {
            for ( size_type j = start_j; j < end_j; ++j )
            {
                auto [eid,eidFeel] = this->addStructuredElement(i,j,partId,partIdGhost,{},withCoord );
                idStructuredMeshToFeelMesh.insert( std::make_pair( eid, eidFeel ) );
                mapGhostElt.insert( std::make_pair( eid, boost::make_tuple( eidFeel, partIdGhost ) ) );
            }
        }
    }

    this->updateGhostCellInfoByUsingNonBlockingComm( idStructuredMeshToFeelMesh, mapGhostElt );
}

template <typename GeoShape, typename T, typename IndexT>
void
MeshStructured<GeoShape,T,IndexT>::updateGhostCellInfoByUsingNonBlockingComm( std::unordered_map<size_type, size_type> const& idStructuredMeshToFeelMesh,
                                                           std::unordered_map<size_type, boost::tuple<size_type, rank_type>> const& mapGhostElt )
{
    DVLOG( 1 ) << "updateGhostCellInfoNonBlockingComm : start on rank " << this->worldComm().localRank() << "\n";

    const int nProc = this->worldComm().localSize();
    DVLOG( 1 ) << "updateGhostCellInfoNonBlockingComm : nProc = " << nProc;
    //std::cout << nProc << std::endl;

    //-----------------------------------------------------------//
    // compute size of container to send
    std::unordered_map<rank_type, int> nDataInVecToSend;
    auto it_map = mapGhostElt.begin();
    auto const en_map = mapGhostElt.end();
    for ( ; it_map != en_map; ++it_map )
    {
        const rank_type idProc = it_map->second.template get<1>();
        if ( nDataInVecToSend.find( idProc ) == nDataInVecToSend.end() )
            nDataInVecToSend[idProc] = 0;
        nDataInVecToSend[idProc]++;
    }
    //-----------------------------------------------------------//
    // init and resize the container to send
    std::unordered_map<rank_type, std::vector<int>> dataToSend;
    auto itNDataInVecToSend = nDataInVecToSend.begin();
    auto const enNDataInVecToSend = nDataInVecToSend.end();
    for ( ; itNDataInVecToSend != enNDataInVecToSend; ++itNDataInVecToSend )
    {
        const rank_type idProc = itNDataInVecToSend->first;
        const int nData = itNDataInVecToSend->second;
        dataToSend[idProc].resize( nData );
    }
    //-----------------------------------------------------------//
    // prepare container to send
    std::unordered_map<rank_type, std::unordered_map<int, int>> memoryMsgToSend;
    std::unordered_map<rank_type, int> nDataInVecToSendBis;
    it_map = mapGhostElt.begin();
    for ( ; it_map != en_map; ++it_map )
    {
        const int idGmsh = it_map->first;
        const int idFeel = it_map->second.template get<0>();
        const rank_type idProc = it_map->second.template get<1>();

        if ( nDataInVecToSendBis.find( idProc ) == nDataInVecToSendBis.end() )
            nDataInVecToSendBis[idProc] = 0;
        // save request
        memoryMsgToSend[idProc][nDataInVecToSendBis[idProc]] = idFeel;
        // update container
        dataToSend[idProc][nDataInVecToSendBis[idProc]] = idGmsh;
        // update counter
        nDataInVecToSendBis[idProc]++;
        // std::cout << idProc << std::endl;
    }
    //-----------------------------------------------------------//
    // counter of request
    int nbRequest = 0;
    for ( rank_type proc = 0; proc < nProc; ++proc )
    {
        if ( dataToSend.find( proc ) != dataToSend.end() )
            nbRequest +=2;
    }
    if ( nbRequest == 0 ) return;

    mpi::request* reqs = new mpi::request[nbRequest];
    int cptRequest = 0;
    //-----------------------------------------------------------//
    // first send
    std::unordered_map<rank_type, std::vector<int>> dataToRecv;
    for ( auto const& [procComm,dataToSendOnProc] : dataToSend )
    {
        reqs[cptRequest++] = this->worldComm().localComm().isend( procComm, 0, dataToSendOnProc.data(), dataToSendOnProc.size() );
        auto & dataToRecvOnProc = dataToRecv[procComm];
        dataToRecvOnProc.resize( dataToSendOnProc.size() );
        reqs[cptRequest++] = this->worldComm().localComm().irecv( procComm, 0, dataToRecvOnProc.data(), dataToRecvOnProc.size() );
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all( reqs, reqs + nbRequest );
    //-----------------------------------------------------------//
    // build the container to ReSend
    std::unordered_map<rank_type, std::vector<int>> dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv != enDataRecv; ++itDataRecv )
    {
        const rank_type idProc = itDataRecv->first;
        const int nDataRecv = itDataRecv->second.size();
        dataToReSend[idProc].resize( nDataRecv );
        //store the idFeel corresponding
        for ( int k = 0; k < nDataRecv; ++k )
        {
            dataToReSend[idProc][k] = idStructuredMeshToFeelMesh.find( itDataRecv->second[k] )->second;
        }
    }
    //-----------------------------------------------------------//
    // send respond to the request
    cptRequest = 0;
    std::unordered_map<rank_type, std::vector<int>> finalDataToRecv;
    for ( auto const& [procComm,dataToSendOnProc] : dataToReSend )
    {
        reqs[cptRequest++] = this->worldComm().localComm().isend( procComm, 0, dataToSendOnProc.data(), dataToSendOnProc.size() );
        auto & dataToRecvOnProc = finalDataToRecv[procComm];
        dataToRecvOnProc.resize( dataToSendOnProc.size() );
        reqs[cptRequest++] = this->worldComm().localComm().irecv( procComm, 0, dataToRecvOnProc.data(), dataToRecvOnProc.size() );
    }
    //-----------------------------------------------------------//
    // wait all requests
    mpi::wait_all( reqs, reqs + nbRequest );
    // delete reqs because finish comm
    delete[] reqs;
    //-----------------------------------------------------------//
    // update mesh : id in other partitions for the ghost cells
    auto itFinalDataToRecv = finalDataToRecv.begin();
    auto const enFinalDataToRecv = finalDataToRecv.end();

    for ( ; itFinalDataToRecv != enFinalDataToRecv; ++itFinalDataToRecv )
    {
        const rank_type idProc = itFinalDataToRecv->first;
        const int nDataRecv = itFinalDataToRecv->second.size();
        //std::cout << idProc << ":" << nDataRecv << std::endl;
        for ( int k = 0; k < nDataRecv; ++k )
        {
            /* std::cout << "I want element " << memoryMsgToSend[idProc][k] << ": " << idProc << std::endl;*/
            auto& eltToUpdate = this->elementIterator( memoryMsgToSend[idProc][k] /*,idProc*/ )->second;
#if 0
            std::cout << "k = " << k << std::endl;
            std::cout << "itFinalDataToRecv->second[k]  " << itFinalDataToRecv->second[k]  << std::endl;
            std::cout << "eltToUpdate->id()             " << eltToUpdate->id()             << std::endl;
            std::cout << "eltToUpdate->processId()      " << eltToUpdate->processId()      << std::endl;
            std::cout << "eltToUpdate->pidInPartition() " << eltToUpdate->pidInPartition() << std::endl;
            std::cout << "eltToUpdate->refDim()         " << eltToUpdate->refDim()         << std::endl;
            std::cout << "eltToUpdate->nPoints()        " << eltToUpdate->nPoints()        << std::endl;
#endif
            eltToUpdate.setIdInOtherPartitions( idProc, itFinalDataToRecv->second[k] );
        }
    }
    //-----------------------------------------------------------//
    DVLOG( 1 ) << "updateGhostCellInfoNonBlockingComm : finish on rank " << this->worldComm().localRank() << "\n";
}


template <typename GeoShape, typename T, typename IndexT>
void
MeshStructured<GeoShape,T,IndexT>::addStructuredPoint( size_type i, size_type j, rank_type partId, bool isGhost,
                                                       bool withCoord )
{
    int ptid = (M_ny)*i + j;
    node_type coords( 2 );
    if ( withCoord && M_cx && M_cy )
    {
        coords[0] = (*M_cy)( j, i );
        coords[1] = (*M_cx)( j, i );
    }
    else
    {
        coords[0] = M_pixelsize * i;
        coords[1] = M_pixelsize *(M_ny - 1 - j);
    }

    point_type pt( ptid, coords );
    if ( !isGhost )
        pt.setProcessId( partId );
    pt.setProcessIdInPartition( partId );
    this->addPoint( pt );
}

template <typename GeoShape, typename T, typename IndexT>
std::pair<typename MeshStructured<GeoShape,T,IndexT>::size_type,typename MeshStructured<GeoShape,T,IndexT>::size_type>
MeshStructured<GeoShape,T,IndexT>::addStructuredElement( size_type i, size_type j, rank_type processId, rank_type partId,
                                                         std::vector<rank_type> const& neighborPartitionIds,
                                                         bool withCoord )
{
    size_type eid = ( M_ny - 1 ) * i + j; // StructuredMesh Id
    element_type e;
    e.setMarker( 1, 1 );
    e.setProcessIdInPartition( processId );
    e.setProcessId( partId );

    std::vector<size_type> ptid = {( M_ny ) * ( i + 1 ) + j,     // 0
                                   ( M_ny ) * ( i + 1 ) + j + 1, // 1
                                   (M_ny)*i + j + 1,             // 2
                                   (M_ny)*i + j};                // 3

    for ( uint16_type k = 0; k < 4; ++k )
        e.setPoint( k, this->point( ptid[k] ) );

    if ( neighborPartitionIds.empty() )
        e.setNeighborPartitionIds( neighborPartitionIds );

    auto [eit,inserted] = this->addElement( e, true ); // e.id() is defined by Feel++
    return std::make_pair( eid, eit->second.id() );
}



template class MeshStructured<Hypercube<2>>;

} // namespace Feel
