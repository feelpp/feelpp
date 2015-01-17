/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-28

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
   \file datamap.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-28
 */
#include <feel/feelalg/datamap.hpp>
//#include <boost/thread/thread.hpp>

namespace Feel
{
DataMap::DataMap( WorldComm const& _worldComm )
    :
    M_closed( false ),
    M_n_dofs( 0 ),
    M_n_localWithoutGhost_df( _worldComm.globalSize(),0 ),
    M_n_localWithGhost_df( _worldComm.globalSize(),0 ),
    M_first_df( _worldComm.globalSize(),0 ),
    M_last_df( _worldComm.globalSize(),0 ),
    M_first_df_globalcluster( _worldComm.globalSize(),0 ),
    M_last_df_globalcluster( _worldComm.globalSize(),0 ),
    M_myglobalelements(),
    M_mapGlobalProcessToGlobalCluster(),
    M_mapGlobalClusterToGlobalProcess(),
    M_worldComm( _worldComm ),
    M_indexSplit()
{}

DataMap::DataMap( size_type n, size_type n_local, WorldComm const& _worldComm )
    :
    M_closed( false ),
    M_n_dofs( n ),
    M_n_localWithoutGhost_df( _worldComm.globalSize(),0 ),
    M_n_localWithGhost_df( _worldComm.globalSize(),0 ),
    M_first_df( _worldComm.globalSize(),0 ),
    M_last_df( _worldComm.globalSize(),0 ),
    M_first_df_globalcluster( _worldComm.globalSize(),0 ),
    M_last_df_globalcluster( _worldComm.globalSize(),0 ),
    M_myglobalelements(),
    M_mapGlobalProcessToGlobalCluster(),
    M_mapGlobalClusterToGlobalProcess(),
    M_worldComm( _worldComm ),
    M_indexSplit()
{

    FEELPP_ASSERT ( n_local <= n )
    ( n_local )( n )
    ( this->worldComm().globalRank() )
    ( this->worldComm().globalSize() ).error( "Invalid local vector size" );

#ifdef FEELPP_HAS_MPI
    std::vector<int> local_sizes     ( this->worldComm().size(), 0 );
    std::vector<int> local_sizes_send( this->worldComm().size(), 0 );

    if ( this->worldComm().size() > 1 )
    {
        LOG(WARNING) << "Not imlemented!";
    }
    else // sequential
    {
        local_sizes[this->worldComm().rank()] = n_local;

        M_first_df[this->worldComm().rank()] = 0;
        M_last_df[this->worldComm().rank()] = n_local-1;
        // mpi
        M_n_localWithoutGhost_df[this->worldComm().rank()]=n_local;
        M_n_localWithGhost_df[this->worldComm().rank()]=n_local;
        M_first_df_globalcluster[this->worldComm().rank()]=0;
        M_last_df_globalcluster[this->worldComm().rank()]=n_local-1;
        // todo! : les map
        M_mapGlobalProcessToGlobalCluster.resize( n_local );
        M_mapGlobalClusterToGlobalProcess.resize( n_local );
    }

    if ( n == invalid_size_type_value )
        M_n_dofs = M_last_df[this->worldComm().rank()]+1;

#  ifdef DEBUG
    // Make sure all the local sizes sum up to the global
    // size, otherwise there is big trouble!
    int sum=0;

    for ( int p=0; p< this->worldComm().size(); p++ )
        sum += local_sizes[p];

    if ( n != invalid_size_type_value )
        FEELPP_ASSERT ( sum == static_cast<int>( n ) )
        ( sum )( n )
        ( this->worldComm().rank() )
        ( this->worldComm().size() ).warn( "invalid distributed vector construction" );

#  endif

#else // FEELPP_HAS_MPI

    // No other options without MPI!
    if ( n != n_local )
    {
        std::cerr << "ERROR:  MPI is required for n != n_local!"
                  << std::endl;
        //error();
    }

#endif // FEELPP_HAS_MPI

    /*
    DVLOG(2) << "        global size = " << this->size() << "\n";
    DVLOG(2) << "        local  size = " << this->localSize() << "\n";
    DVLOG(2) << "  first local index = " << this->firstLocalIndex() << "\n";
    DVLOG(2) << "   last local index = " << this->lastLocalIndex() << "\n";
    */

}

DataMap::~DataMap()
{}


void
DataMap::close() const
{
    // we assume here that the data is contiguous
    M_myglobalelements.resize( nMyElements() );

    for ( size_type i = 0; i < nMyElements(); ++i )
    {
        M_myglobalelements[i] = minMyGID() + i;
    }

    M_closed = true;
}
std::vector<size_type> const&
DataMap::myGlobalElements() const
{
    if ( !this->closed() )
        this->close();

    return M_myglobalelements;
}


void
DataMap::setNDof( size_type ndof )
{
    M_n_dofs=ndof;
}

void
DataMap::setNLocalDofWithoutGhost( const rank_type proc, const size_type n, bool inWorld )
{
    M_n_localWithoutGhost_df[proc]=n;
}
void
DataMap::setNLocalDofWithGhost( const rank_type proc, const size_type n, bool inWorld )
{
    M_n_localWithGhost_df[proc]=n;
}
void
DataMap::setFirstDof( const rank_type proc, const size_type df, bool inWorld )
{
    FEELPP_ASSERT( proc < M_first_df.size() )( proc )( M_first_df.size() ).error( "invalid proc id or dof table" );
    M_first_df[proc]=df;
}
void
DataMap::setLastDof( const rank_type proc, const size_type df, bool inWorld )
{
    FEELPP_ASSERT( proc < M_first_df.size() )( proc )( M_first_df.size() ).error( "invalid proc id or dof table" );
    M_last_df[proc]=df;
}
void
DataMap::setFirstDofGlobalCluster( const rank_type proc, const size_type df, bool inWorld )
{
    FEELPP_ASSERT( proc < M_first_df_globalcluster.size() )( proc )( M_first_df_globalcluster.size() ).error( "invalid proc id or dof table" );
    M_first_df_globalcluster[proc]=df;
}
void
DataMap::setLastDofGlobalCluster( const rank_type proc, const size_type df, bool inWorld )
{
    FEELPP_ASSERT( proc < M_first_df_globalcluster.size() )( proc )( M_first_df_globalcluster.size() ).error( "invalid proc id or dof table" );
    M_last_df_globalcluster[proc]=df;
}


void
DataMap::setMapGlobalProcessToGlobalCluster( std::vector<size_type> const& map )
{
    M_mapGlobalProcessToGlobalCluster=map;
}
void
DataMap::setMapGlobalClusterToGlobalProcess( std::vector<size_type> const& map )
{
    M_mapGlobalClusterToGlobalProcess=map;
}
void
DataMap::setMapGlobalProcessToGlobalCluster( size_type i, size_type j )
{
    M_mapGlobalProcessToGlobalCluster[i]=j;
}
void
DataMap::setMapGlobalClusterToGlobalProcess( size_type i, size_type j )
{
    M_mapGlobalClusterToGlobalProcess[i]=j;
}
void
DataMap::resizeMapGlobalProcessToGlobalCluster( size_type n )
{
    M_mapGlobalProcessToGlobalCluster.resize( n );
}
void
DataMap::resizeMapGlobalClusterToGlobalProcess( size_type n )
{
    M_mapGlobalClusterToGlobalProcess.resize( n );
}

void
DataMap::updateDataInWorld()
{
    mpi::all_gather( this->worldComm().globalComm(),
                     this->M_n_localWithoutGhost_df[this->worldComm().globalRank()],
                     this->M_n_localWithoutGhost_df );
    mpi::all_gather( this->worldComm().globalComm(),
                     this->M_n_localWithGhost_df[this->worldComm().globalRank()],
                     this->M_n_localWithGhost_df );

    mpi::all_gather( this->worldComm().globalComm(),
                     this->M_first_df[this->worldComm().globalRank()],
                     this->M_first_df );
    mpi::all_gather( this->worldComm().globalComm(),
                     this->M_last_df[this->worldComm().globalRank()],
                     this->M_last_df );
    mpi::all_gather( this->worldComm().globalComm(),
                     this->M_first_df_globalcluster[this->worldComm().globalRank()],
                     this->M_first_df_globalcluster );
    mpi::all_gather( this->worldComm().globalComm(),
                     this->M_last_df_globalcluster[this->worldComm().globalRank()],
                     this->M_last_df_globalcluster );

    //_globalcluster
}


rank_type
DataMap::procOnGlobalCluster( size_type globDof ) const
{
    rank_type proc=0,res=invalid_rank_type_value;
    bool find=false;

    while ( ( !find ) && ( proc<this->nProcessors() ) )
    {
        if ( ( this->nLocalDofWithoutGhost(proc) > 0 ) && ( globDof <= M_last_df_globalcluster[proc] ) && ( globDof >= M_first_df_globalcluster[proc] ) )
        {
            res = proc;
            find=true;
        }
        else
            ++proc;
    }
    return res;
}

boost::tuple<bool,size_type>
DataMap::searchGlobalProcessDof( size_type gpdof ) const
{
    bool find=false;
    size_type gDofProcess = 0;
    const size_type startLoc = this->firstDof();
    const size_type endLoc = startLoc+this->nLocalDofWithGhost();
    for ( size_type k=startLoc ; k < endLoc && !find ; ++k )
        if ( this->mapGlobalProcessToGlobalCluster(k) == gpdof )
        {
            gDofProcess=k;
            find =true;
        }

    return boost::make_tuple( find,gDofProcess );
}

void
DataMap::updateIndexSetWithParallelMissingDof( std::vector<size_type> & _indexSet ) const
{
    auto res = this->buildIndexSetWithParallelMissingDof( _indexSet );
    _indexSet.assign( res.begin(), res.end() );
}

std::vector<size_type>
DataMap::buildIndexSetWithParallelMissingDof( std::vector<size_type> const& _indexSet ) const
{
    std::vector<size_type> indexSet;
    indexSet.insert( indexSet.begin(), _indexSet.begin(), _indexSet.end() );

    // if sequential return identical index set
    if ( this->worldComm().localSize() == 1 )
        return indexSet;

    // init data used in mpi comm
    std::map< rank_type, std::vector< size_type > > dataToSend, dataToRecv;
    for ( rank_type p : this->neighborSubdomains() )
        dataToSend[p].clear();

    // up data used in mpi comm
    for ( auto const& id : indexSet )
    {
        if ( !this->dofGlobalProcessIsGhost(id) )
        {
            if ( this->activeDofSharedOnCluster().find( id ) != this->activeDofSharedOnCluster().end() )
                for ( rank_type pNeighborId : this->activeDofSharedOnCluster().find( id )->second )
                    dataToSend[pNeighborId].push_back( this->mapGlobalProcessToGlobalCluster( id ) );
        }
        else
        {
            size_type gcdof = this->mapGlobalProcessToGlobalCluster( id );
            rank_type procIdFinded = procOnGlobalCluster( gcdof );
            CHECK ( procIdFinded != invalid_rank_type_value ) << " proc not find for gcdof : " << gcdof;
            dataToSend[procIdFinded].push_back( gcdof );
        }
    }

    // prepare mpi com
    int nbRequest=2*this->neighborSubdomains().size();
    mpi::request * reqs = new mpi::request[nbRequest];
    // apply isend/irecv
    int cptRequest=0;
    for ( rank_type p : this->neighborSubdomains() )
    {
        CHECK( dataToSend.find(p) != dataToSend.end() ) << " no data to send to proc " << p << "\n";
        reqs[cptRequest] = this->worldComm().localComm().isend( p , 0, dataToSend.find(p)->second );
        ++cptRequest;
        reqs[cptRequest] = this->worldComm().localComm().irecv( p , 0, dataToRecv[p] );
        ++cptRequest;
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    delete [] reqs;

    // update indexSet : can be usefull for fix missing dof entries in input at interprocess zone!
    for ( auto const& dataR : dataToRecv )
    {
        rank_type theproc = dataR.first;
        for ( size_type dataRfromproc : dataR.second )
        {
            //size_type thelocdof = this->mapGlobalClusterToGlobalProcess( dataRfromproc );
            auto thelocdof = this->searchGlobalProcessDof( dataRfromproc );
            CHECK( thelocdof.get<0>() ) << "local dof not find with cluster id : " << dataRfromproc;
            //indexSet.insert( thelocdof.get<1>() );
            indexSet.push_back( thelocdof.get<1>() );
        }
    }

    return indexSet;
}

boost::shared_ptr<DataMap>
DataMap::createSubDataMap( std::vector<size_type> const& _idExtract, bool _checkAndFixInputRange ) const
{
    rank_type curProcId = this->worldComm().localRank();

    bool checkAndFixInputRange = ( this->worldComm().localSize() > 1 &&  _checkAndFixInputRange );
    std::vector<size_type> idExtractVec = ( checkAndFixInputRange )?
        this->buildIndexSetWithParallelMissingDof( _idExtract ) : _idExtract;

    // convert input idExtract into std::set ( preserve ordering )
    std::set<size_type> idExtract;
    idExtract.insert( idExtractVec.begin(), idExtractVec.end() );

    // init new dataMap pointer
    boost::shared_ptr<DataMap> dataMapRes( new DataMap(this->worldComm()) );

    // define neighbor process as identical
    dataMapRes->setNeighborSubdomains( this->neighborSubdomains() );

    // number of local dof define from input
    size_type nLocalDofWithGhost = idExtract.size();
    // compute number of active dof on current proc
    size_type nLocalDofWithoutGhost=0;
    for ( auto const& id : idExtract )
    {
        if ( !this->dofGlobalProcessIsGhost(id) )
            ++nLocalDofWithoutGhost;
    }

    // collective mpi comm about nLocalDofWithGhost and nLocalDofWithoutGhost
    std::vector<boost::tuple<size_type,size_type> > dataRecvFromGather;
    auto dataSendToGather = boost::make_tuple(nLocalDofWithGhost,nLocalDofWithoutGhost);
    mpi::all_gather( this->worldComm(),
                     dataSendToGather,
                     dataRecvFromGather );

    // define first, last, nlocalDof for each process and get nDof in global cluster
    size_type nDof = 0;
    for (rank_type p=0;p<this->worldComm().localSize();++p)
    {
        const size_type nLocalDofWithGhostOnProc = dataRecvFromGather[p].template get<0>();
        dataMapRes->setFirstDof( p, 0 );
        dataMapRes->setLastDof( p, (nLocalDofWithGhostOnProc > 0)? nLocalDofWithGhostOnProc-1 :0  );
        dataMapRes->setNLocalDofWithGhost( p, nLocalDofWithGhostOnProc );

        const size_type nLocalDofWithoutGhostOnProc = dataRecvFromGather[p].template get<1>();
        if ( p == 0 )
        {
            dataMapRes->setFirstDofGlobalCluster( p, 0 );
            dataMapRes->setLastDofGlobalCluster( p, (nLocalDofWithoutGhostOnProc>0)? nLocalDofWithoutGhostOnProc-1 : 0 );
        }
        else
        {
            if ( dataMapRes->nLocalDofWithoutGhost(p-1) > 0 )
                dataMapRes->setFirstDofGlobalCluster( p, dataMapRes->lastDofGlobalCluster( p-1 ) + 1 );
            else
                dataMapRes->setFirstDofGlobalCluster( p, dataMapRes->lastDofGlobalCluster( p-1 ) );

            if ( nLocalDofWithoutGhostOnProc > 0 )
                dataMapRes->setLastDofGlobalCluster( p, dataMapRes->firstDofGlobalCluster(p) + nLocalDofWithoutGhostOnProc - 1 );
            else
                dataMapRes->setLastDofGlobalCluster( p, dataMapRes->firstDofGlobalCluster(p) );
        }
        dataMapRes->setNLocalDofWithoutGhost( p, nLocalDofWithoutGhostOnProc );
        nDof+=nLocalDofWithoutGhostOnProc;
    }
    dataMapRes->setNDof( nDof );


    // define mapping between global process and global cluster and store missing ghost dof to send
    dataMapRes->resizeMapGlobalProcessToGlobalCluster( dataMapRes->nLocalDofWithGhost(curProcId) );
    dataMapRes->resizeMapGlobalClusterToGlobalProcess( dataMapRes->nLocalDofWithoutGhost(curProcId) );
    size_type firstDofGC = dataMapRes->firstDofGlobalCluster(curProcId);
    size_type idLocalDof=0, idGlobalDof=firstDofGC;
    std::map< rank_type, std::vector< size_type > > dataToSend, memoryDataToSend;
    for ( rank_type p : this->neighborSubdomains() )
    {
        dataToSend[p].clear();
        memoryDataToSend[p].clear();
    }
    std::map< size_type, size_type> dofTableRelationGP;
    for ( auto const& id : idExtract )
    {
        if ( !this->dofGlobalProcessIsGhost(id) )
        {
            dataMapRes->M_mapGlobalProcessToGlobalCluster[idLocalDof]=idGlobalDof;
            dataMapRes->M_mapGlobalClusterToGlobalProcess[idGlobalDof-firstDofGC]=idLocalDof;
            if ( this->activeDofSharedOnCluster().find( id ) != this->activeDofSharedOnCluster().end() )
                dofTableRelationGP[id]=idLocalDof;
            ++idGlobalDof;
        }
        else
        {
            size_type gcdof = this->mapGlobalProcessToGlobalCluster( id );
            rank_type procIdFinded = procOnGlobalCluster( gcdof );
            CHECK ( procIdFinded != invalid_rank_type_value ) << " proc not find for gcdof : " << gcdof;
            dataToSend[procIdFinded].push_back( gcdof );
            memoryDataToSend[procIdFinded].push_back( idLocalDof );
            dataMapRes->M_mapGlobalProcessToGlobalCluster[idLocalDof]=invalid_size_type_value;
        }
        ++idLocalDof;
    }

    if ( this->worldComm().localSize() > 1 )
    {
        // init data used in mpi comm
        std::map< rank_type, std::vector< size_type > > dataToRecv;

        // prepare mpi com
        int nbRequest=2*this->neighborSubdomains().size();
        mpi::request * reqs = new mpi::request[nbRequest];
        // apply isend/irecv
        int cptRequest=0;
        for ( rank_type p : this->neighborSubdomains() )
        {
            CHECK( dataToSend.find(p) != dataToSend.end() ) << " no data to send to proc " << p << "\n";
            reqs[cptRequest] = this->worldComm().localComm().isend( p , 0, dataToSend.find(p)->second );
            ++cptRequest;
            reqs[cptRequest] = this->worldComm().localComm().irecv( p , 0, dataToRecv[p] );
            ++cptRequest;
        }
        // wait all requests
        mpi::wait_all(reqs, reqs + nbRequest);

        // treat recv and store request to re-send
        std::map< rank_type, std::vector< size_type > > dataToReSend,dataToReRecv;
        for ( rank_type p : this->neighborSubdomains() )
            dataToReSend[p].clear();
        size_type firstDofGC = this->firstDofGlobalCluster(curProcId);
        for ( auto const& dataR : dataToRecv )
        {
            rank_type theproc = dataR.first;
            for ( auto const& gcdof : dataR.second )
            {
                size_type gpdof = this->mapGlobalClusterToGlobalProcess(gcdof-firstDofGC);
                dataMapRes->M_activeDofSharedOnCluster[gpdof].insert( theproc );

                CHECK( dofTableRelationGP.find( gpdof ) != dofTableRelationGP.end() ) << "active dof gpdof not register";
                size_type gcdoffinded = dataMapRes->mapGlobalProcessToGlobalCluster( dofTableRelationGP.find( gpdof )->second );
                CHECK( gcdoffinded != invalid_size_type_value ) << " active dof not register in map";
                dataToReSend[theproc].push_back( gcdoffinded );
            }
        }

        // apply isend/irecv
        cptRequest=0;
        for ( rank_type p : this->neighborSubdomains() )
        {
            CHECK( dataToReSend.find(p) != dataToReSend.end() ) << " no data to send to proc " << p << "\n";
            reqs[cptRequest] = this->worldComm().localComm().isend( p , 0, dataToReSend.find(p)->second );
            ++cptRequest;
            reqs[cptRequest] = this->worldComm().localComm().irecv( p , 0, dataToReRecv[p] );
            ++cptRequest;
        }
        // wait all requests
        mpi::wait_all(reqs, reqs + nbRequest);
        delete [] reqs;

        // treat last recv to fix missing entries in mapGlobalProcessToGlobalCluster
        for ( auto const& dataRR : dataToReRecv )
        {
            rank_type theproc = dataRR.first;
            size_type k=0;
            for ( auto const& gcdof : dataRR.second )
            {
                size_type locdof = memoryDataToSend.find( theproc )->second[k];
                dataMapRes->M_mapGlobalProcessToGlobalCluster[locdof] = gcdof;
                ++k;
            }

        }
    } // if ( this->worldComm().localSize() > 1 )



    if ( true )
    {
        // build one split for all dof
        dataMapRes->buildIndexSplit();
    }
    else
    {
        // TODO : reuse parent indexsplit structure to build the new one
        boost::shared_ptr<IndexSplit> indexSplit( new IndexSplit( this->indexSplit()->size() ) );
        size_type startIS = dataMapRes->firstDofGlobalCluster( curProcId );
        for ( int k = 0 ; k < this->indexSplit()->size() ; ++k )
        {

        }
    }



    return dataMapRes;
}


void
DataMap::showMeMapGlobalProcessToGlobalCluster( bool showAll, std::ostream& __out2 ) const
{
    //__out << std::endl;
    this->comm().globalComm().barrier();

    for ( rank_type proc = 0; proc<this->comm().globalSize(); ++proc )
    {
        this->comm().globalComm().barrier();
        if ( proc==this->worldComm().globalRank() )//this->worldComm().masterRank() )
        {
            std::ostringstream __out;
            this->comm().globalComm().barrier();
            __out << "\n";
            __out << "-----------------------------------------------------------------------\n"
                  << "------------------showMeMapGlobalProcessToGlobalCluster----------------\n"
                  << "-----------------------------------------------------------------------\n"
                  << "god rank : " << this->comm().godRank()  << "\n"
                  << "global rank : " << this->comm().globalRank()  << "\n"
                  << "local rank : " << this->comm().localRank()  << "\n"
                  << "rank : " << proc  << "\n"
                  << "nDof : " << this->nDof() << "\n"
                  << "nLocalDof : " << this->nLocalDof() << "\n"
                  << "nLocalDofWithoutGhost : " << this->nLocalDofWithoutGhost() << "\n"
                  << "nLocalDofWithGhost : " << this->nLocalDofWithGhost() << "\n"
                  << "neighbor processors : (" << M_neighbor_processors.size() << ")";
            for ( rank_type pNeighborId : M_neighbor_processors )
                __out << " " << pNeighborId;
            __out << "\n";
            __out << "mapGlobalProcessToGlobalCluster().size() " << this->mapGlobalProcessToGlobalCluster().size() << "\n"
                  << "-----------------------------------------------------------------------\n";

            if (showAll)
            {
#if 1
                __out << "mapGlobalProcessToGlobalCluster : \n";
                for ( size_type i=0 ; i<this->mapGlobalProcessToGlobalCluster().size() ; ++i )
                {
                    __out << i << " " << this->mapGlobalProcessToGlobalCluster()[i]
                          << " real proc " << procOnGlobalCluster( /*this->*/mapGlobalProcessToGlobalCluster()[i] ) <<"\n";
                }
                __out << "-----------------------------------------------------------------------\n";
#endif
#if 0
                __out << "mapGlobalClusterToGlobalProcess : \n";
                for ( size_type i=0 ; i<this->mapGlobalClusterToGlobalProcess().size() ; ++i )
                {
                    __out << i << " " << this->mapGlobalClusterToGlobalProcess()[i]
                          <<"\n";
                }
                __out << "-----------------------------------------------------------------------\n";
#endif
            }
#if 1
            __out << " M_first_df : ";

            for ( rank_type i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this-> M_first_df[i] << " ";
            }

            __out << "\n";
            __out << " M_last_df : ";

            for ( rank_type i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this-> M_last_df[i] << " ";
            }

            __out << "\n";
            __out << " M_first_df_globalcluster : ";

            for ( rank_type i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this-> M_first_df_globalcluster[i] << " ";
            }

            __out << "\n";
            __out << " M_last_df_globalcluster : ";

            for ( rank_type i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this-> M_last_df_globalcluster[i] << " ";
            }

            __out << "\n";
            __out << " M_n_localWithGhost_df : ";

            for ( rank_type i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this->M_n_localWithGhost_df[i] << " ";
            }

            __out << "\n";
            __out << " M_n_localWithoutGhost_df : ";

            for ( rank_type i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this->M_n_localWithoutGhost_df[i] << " ";
            }

            __out << "\n";

#endif
            __out << "-----------------------------------------------------------------------\n";

            __out << "\n";// << std::endl;

            __out2 << __out.str() << std::endl;
        }
        //else  sleep(1);


        //this->comm().barrier();
        this->comm().globalComm().barrier();
        double mydelay=0.3;
        sleep(mydelay);

        //boost::this_thread::sleep( boost::posix_time::seconds(1) );

    }

    this->comm().globalComm().barrier();
    //sleep(1);

}


void
DataMap::buildIndexSplit()
{
    DVLOG(1) << "buildIndexSplit() start\n";
    M_indexSplit.reset( new indexsplit_type( 1 ) );
    M_indexSplit->operator[](0).resize( this->nLocalDofWithoutGhost() );

    M_indexSplit->setFirstIndex( 0, this->firstDofGlobalCluster() );
    M_indexSplit->setLastIndex( 0, this->lastDofGlobalCluster() );
    M_indexSplit->setNIndex( 0, this->nLocalDofWithoutGhost() );

    const size_type firstDof = this->firstDofGlobalCluster();

    for ( size_type index = 0; index < this->nLocalDofWithGhost() ; ++index )
    {
        if ( this->dofGlobalProcessIsGhost(index) ) continue;
        const size_type globalDof = this->mapGlobalProcessToGlobalCluster(index);
        M_indexSplit->operator[](0)[globalDof - firstDof] = globalDof;
    }

    size_type nDofForSmallerRankId=0;
    for ( rank_type proc=0;proc<this->worldComm().globalRank();++proc )
        nDofForSmallerRankId+=this->nLocalDofWithoutGhost( proc );
    M_indexSplit->setNIndexForSmallerRankId( 0, nDofForSmallerRankId );

    DVLOG(1) << "buildIndexSplit() done\n";
}

void
IndexSplit::FieldsDef::showMe() const
{
    std::cout << "FieldsDef showMe\n ";
    for (auto it=this->begin(),en=this->end();it!=en;++it)
    {
        std::cout << "field " << it->first << " -> ";
        for ( int val : it->second )
            std::cout << " " << val;
        std::cout << "\n";
    }
}


IndexSplit::FieldsDef
IndexSplit::parseFieldsDef( std::string s )
{
    // TODO : remove white space before

    IndexSplit::FieldsDef res;

    int index=0;
    while ( index < s.size() )
    {
        bool findArrow=false;
        int indexLenght=0;
        int fieldId=0;

        while (!findArrow && index < s.size())
        {
            if ( s.substr(index,2) == "->" )
            {
                fieldId=boost::lexical_cast<int>( s.substr(index-indexLenght,indexLenght).c_str() );
                findArrow=true;
                index+=2;
            }
            else
            {
                ++indexLenght;
                ++index;
            }
        }

        CHECK( findArrow ) << "fields def expression has not symbol -> : " << s << "\n";

        if ( s.substr(index,1) == "(" )
        {
            ++index;
            bool find=false;
            std::vector<int> allIndexLenght;
            /*int*/ indexLenght=0;
            int nSplit=0;
            while (!find && index < s.size() )
            {
                if ( s.substr(index,1) == ")" )
                {
                    ++nSplit;
                    int splitId = boost::lexical_cast<int>( s.substr(index-indexLenght,indexLenght).c_str() );
                    res[fieldId].insert( splitId );
                    find=true;
                    index+=indexLenght+1;
                }
                else
                {
                    if ( s.substr(index,1) == "," )
                    {
                        ++nSplit;
                        int splitId = boost::lexical_cast<int>( s.substr(index-indexLenght,indexLenght).c_str() );
                        res[fieldId].insert( splitId );
                        indexLenght=0;
                    }
                    else
                    {
                        ++indexLenght;
                    }
                    ++index;
                }
            }
            CHECK( find ) << "fields def expression has an invalid format " << s << "\n";
        } // (

    }

    //res.showMe();
    return res;
}


void
IndexSplit::resize( int size )
{
    super_type::resize( size );
    M_firstIndex.resize( size );
    M_lastIndex.resize( size );
    M_nIndex.resize( size );
    M_nIndexForSmallerRankId.resize( size );
}


void
IndexSplit::addSplit( size_type startSplit, self_ptrtype const& addedIndexSplit )
{
    const int sizeIS1 = this->size();
    const int sizeIS2 = addedIndexSplit->size();
    const int newSize = sizeIS1+sizeIS2;
    this->resize( newSize );

    size_type startIS = startSplit;
    for ( int k = 0 ; k < newSize ; ++k )
    {
        if ( k >= sizeIS1 )
        {
            int sizeSplitAdded = addedIndexSplit->split(k-sizeIS1).size();

            this->setFirstIndex( k, startIS );
            this->setLastIndex( k, (sizeSplitAdded>0)? startIS+sizeSplitAdded-1 : startIS  );
            this->setNIndex( k, sizeSplitAdded );

            this->operator[](k).resize( sizeSplitAdded );
            const size_type firstIndexAdded = addedIndexSplit->firstIndex(k-sizeIS1);
            for ( int l=0 ; l<sizeSplitAdded ; ++l )
                this->operator[](k)[l] = startIS + addedIndexSplit->split(k-sizeIS1)[l] - firstIndexAdded;
            this->setNIndexForSmallerRankId( k, addedIndexSplit->nIndexForSmallerRankId( k-sizeIS1 ) );
        }

        startIS += this->operator[](k).size();
    }

}


IndexSplit::self_ptrtype
IndexSplit::applyFieldsDef( IndexSplit::FieldsDef const& fieldsDef ) const
{
#if 0
    int nField = fieldsDef.size();
    self_ptrtype newIS( new self_type( nField ) );
    auto it = fieldsDef.begin();
    auto const en = fieldsDef.end();
    for ( ; it!=en ; ++it)
    {
        int k = it->first;
        int sizeNewField = 0;
        for ( int field : it->second )
        {
            sizeNewField += this->operator[]( field ).size();
        }
        newIS->operator[](k).resize( sizeNewField );

        int startField=0;
        for ( int field : it->second )
        {
            int sizeField = this->operator[]( field ).size();
            for ( int i = 0 ; i < sizeField; ++i )
                newISoperator[](k)[startField+i] = this->operator[]( field )[i];
            startField += sizeField;
        }
    }

    return newIS;
#else
    int nField = fieldsDef.size();

    self_ptrtype newIS( new self_type( nField ) );

    // get split id present in def
    std::set<int> splitSetInDef;
    auto it = fieldsDef.begin();
    auto const en = fieldsDef.end();
    for ( ; it != en ; ++it)
        for ( int fieldId : it->second )
            splitSetInDef.insert( fieldId );

    //std::cout << "splitSetInDef.size() " << splitSetInDef.size() << "\n";

    // init start index
    std::vector<int> startSplit( this->size() );
    for ( int i = 0 ; i < this->size(); ++i )
        startSplit[i] = this->firstIndex( i );
    // fix start index if not all field are present in def
    for ( int splitId = 0 ; splitId < this->size(); ++splitId )
    {
        auto itFindSplit = std::find_if( splitSetInDef.begin(), splitSetInDef.end(),
                                         [&splitId]( int s ) { return s == splitId; } );
        if ( itFindSplit == splitSetInDef.end() )
        {
            for ( int splitId2 : splitSetInDef )
            {
                if ( splitId2 > splitId )
                    startSplit[splitId2] -= this->nIndex(splitId);
                startSplit[splitId2] -= this->nIndexForSmallerRankId( splitId );
            }
        }
    }

    // update new index split
    for ( it = fieldsDef.begin() ; it != en ; ++it)
    {
        int fieldId = it->first;
        int sizeNewSplit = 0;
        for ( int splitId : it->second )
            sizeNewSplit += this->operator[]( splitId ).size();

        newIS->operator[](fieldId).resize( sizeNewSplit );

        bool hasInitFirstIndex = false;
        size_type startIndexSplit=0, nIndexForSmallerRank=0;
        for ( int splitId : it->second )
        {
            if ( !hasInitFirstIndex )
            {
                newIS->setFirstIndex( fieldId, startSplit[splitId] );
                hasInitFirstIndex=true;
            }

            int sizeSplit = this->operator[]( splitId ).size();
            for ( int i = 0 ; i < sizeSplit; ++i )
                newIS->operator[](fieldId)[startIndexSplit+i] = startSplit[splitId] + this->operator[]( splitId )[i] - this->firstIndex( splitId );
            //newIS[fieldId][startIndexSplit+i] =  this->operator[]( splitId )[i];
            startIndexSplit += sizeSplit;

            nIndexForSmallerRank += this->nIndexForSmallerRankId( splitId );
        }

        newIS->setLastIndex( fieldId, (sizeNewSplit>0)? newIS->firstIndex(fieldId)+sizeNewSplit-1 : newIS->firstIndex(fieldId) );
        newIS->setNIndex( fieldId, sizeNewSplit );
        newIS->setNIndexForSmallerRankId( fieldId, nIndexForSmallerRank );
    }

    return newIS;
#endif

}



void
IndexSplit::showMe() const
{
    int nSplit = this->size();
    std::ostringstream ostr;
    ostr << "------------------IndexSplit::showMe()----------------\n";
    ostr << "-proc " << Environment::worldComm().globalRank() << "\n"
         << "-number of split : " << nSplit << "\n";
    for (size_type k=0;k<nSplit;++k)
    {
        size_type nDofInSplit = this->operator[](k).size();
        ostr << "-split : " << k << "\n"
             << "-firstIndex : " << this->firstIndex(k) << "\n"
             << "-lastIndex : " << this->lastIndex(k) << "\n"
             << "-nIndex : " << this->nIndex(k) << "\n"
             << "-nIndexForSmallerRankId : " << this->nIndexForSmallerRankId(k) << "\n";
        if ( true )
        {
            for ( int i=0;i<nDofInSplit;++i)
            {
                ostr << this->operator[](k)[i] << " ";
            }
            ostr << "\n";
        }
        ostr << "------------------------------------------------------\n";
    }

    Environment::worldComm().barrier();
    for ( rank_type proc = 0; proc<Environment::worldComm().globalSize(); ++proc )
    {
        Environment::worldComm().barrier();
        if ( proc==Environment::worldComm().globalRank() )
            std::cout<< ostr.str();
    }
    Environment::worldComm().barrier();

}


} // namespace Feel
