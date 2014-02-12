/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
    M_worldComm( _worldComm )
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
    M_worldComm( _worldComm )
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

DataMap::DataMap( DataMap const & dm )
    :
    M_closed( dm.M_closed ),
    M_n_dofs( dm.M_n_dofs ),
    M_n_localWithoutGhost_df( dm.M_n_localWithoutGhost_df ),
    M_n_localWithGhost_df( dm.M_n_localWithGhost_df ),
    M_first_df( dm.M_first_df ),
    M_last_df( dm.M_last_df ),
    M_first_df_globalcluster( dm.M_first_df_globalcluster ),
    M_last_df_globalcluster( dm.M_last_df_globalcluster ),
    M_myglobalelements(),
    M_mapGlobalProcessToGlobalCluster( dm.M_mapGlobalProcessToGlobalCluster ),
    M_mapGlobalClusterToGlobalProcess( dm.M_mapGlobalClusterToGlobalProcess ),
    M_neighbor_processors( dm.M_neighbor_processors ),
    M_worldComm( dm.M_worldComm )
{}
DataMap::~DataMap()
{}

DataMap&
DataMap::operator=( DataMap const& dm )
{
    if ( this != &dm )
    {
        M_worldComm = dm.M_worldComm;
        M_closed = dm.M_closed;
        M_n_dofs = dm.M_n_dofs;
        M_n_localWithoutGhost_df = dm.M_n_localWithoutGhost_df;
        M_n_localWithGhost_df = dm.M_n_localWithGhost_df;
        M_first_df = dm.M_first_df;
        M_last_df = dm.M_last_df;
        M_first_df_globalcluster = dm.M_first_df_globalcluster;
        M_last_df_globalcluster = dm.M_last_df_globalcluster;
        M_myglobalelements = dm.M_myglobalelements;
        M_mapGlobalProcessToGlobalCluster = dm.M_mapGlobalProcessToGlobalCluster;
        M_mapGlobalClusterToGlobalProcess = dm.M_mapGlobalClusterToGlobalProcess;
        M_neighbor_processors = dm.M_neighbor_processors;
    }

    return *this;
}
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
DataMap::setNLocalDofWithoutGhost( const size_type proc, const size_type n, bool inWorld )
{
    M_n_localWithoutGhost_df[proc]=n;
}
void
DataMap::setNLocalDofWithGhost( const size_type proc, const size_type n, bool inWorld )
{
    M_n_localWithGhost_df[proc]=n;
}
void
DataMap::setFirstDof( const size_type proc, const size_type df, bool inWorld )
{
    FEELPP_ASSERT( proc < M_first_df.size() )( proc )( M_first_df.size() ).error( "invalid proc id or dof table" );
    M_first_df[proc]=df;
}
void
DataMap::setLastDof( const size_type proc, const size_type df, bool inWorld )
{
    FEELPP_ASSERT( proc < M_first_df.size() )( proc )( M_first_df.size() ).error( "invalid proc id or dof table" );
    M_last_df[proc]=df;
}
void
DataMap::setFirstDofGlobalCluster( const size_type proc, const size_type df, bool inWorld )
{
    FEELPP_ASSERT( proc < M_first_df_globalcluster.size() )( proc )( M_first_df_globalcluster.size() ).error( "invalid proc id or dof table" );
    M_first_df_globalcluster[proc]=df;
}
void
DataMap::setLastDofGlobalCluster( const size_type proc, const size_type df, bool inWorld )
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


uint16_type
DataMap::procOnGlobalCluster( size_type globDof ) const
{
    uint16_type proc=0,res=0;
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
DataMap::showMeMapGlobalProcessToGlobalCluster( bool showAll, std::ostream& __out2 ) const
{
    //__out << std::endl;
    this->comm().globalComm().barrier();

    for ( int proc = 0; proc<this->comm().globalSize(); ++proc )
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
                  << "mapGlobalProcessToGlobalCluster().size() " << this->mapGlobalProcessToGlobalCluster().size() << "\n"
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

            for ( int i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this-> M_first_df[i] << " ";
            }

            __out << "\n";
            __out << " M_last_df : ";

            for ( int i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this-> M_last_df[i] << " ";
            }

            __out << "\n";
            __out << " M_first_df_globalcluster : ";

            for ( int i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this-> M_first_df_globalcluster[i] << " ";
            }

            __out << "\n";
            __out << " M_last_df_globalcluster : ";

            for ( int i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this-> M_last_df_globalcluster[i] << " ";
            }

            __out << "\n";
            __out << " M_n_localWithGhost_df : ";

            for ( int i=0; i<this->worldComm().globalSize(); ++i )
            {
                __out << this->M_n_localWithGhost_df[i] << " ";
            }

            __out << "\n";
            __out << " M_n_localWithoutGhost_df : ";

            for ( int i=0; i<this->worldComm().globalSize(); ++i )
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

} // namespace Feel
