/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-28
 */
#include <feel/feelalg/datamap.hpp>

namespace Feel
{
DataMap::DataMap(WorldComm const& _worldComm)
    :
    M_closed( false ),
    _M_n_dofs( 0 ),
    _M_n_localWithoutGhost_df( ),
    _M_n_localWithGhost_df( ),
    _M_first_df( ),
    _M_last_df(),
    _M_first_df_globalcluster(),
    _M_last_df_globalcluster(),
    M_myglobalelements(),
    M_mapGlobalProcessToGlobalCluster(),
    M_mapGlobalClusterToGlobalProcess(),
    M_worldComm(_worldComm)
{
    //std::cout << "\nDataMap : build empty! on godrank " << this->worldComm().godRank() << std::endl;
    _M_n_localWithoutGhost_df.resize( this->worldComm().globalSize() );
    _M_n_localWithGhost_df.resize( this->worldComm().globalSize() );
    _M_first_df.resize( this->worldComm().globalSize() );
    _M_last_df.resize( this->worldComm().globalSize() );
    _M_first_df_globalcluster.resize( this->worldComm().globalSize() );
    _M_last_df_globalcluster.resize( this->worldComm().globalSize() );
}
    DataMap::DataMap( size_type n, size_type n_local, WorldComm const& _worldComm )
    :
    M_closed( false ),
    _M_n_dofs( n ),
    _M_n_localWithoutGhost_df( ),
    _M_n_localWithGhost_df( ),
    _M_first_df(),
    _M_last_df(),
    _M_first_df_globalcluster(),
    _M_last_df_globalcluster(),
    M_myglobalelements(),
    M_mapGlobalProcessToGlobalCluster(),
    M_mapGlobalClusterToGlobalProcess(),
    M_worldComm(_worldComm)
{
    _M_n_localWithoutGhost_df.resize( this->worldComm().globalSize() );
    _M_n_localWithGhost_df.resize( this->worldComm().globalSize() );
    _M_first_df.resize( this->worldComm().globalSize() );
    _M_last_df.resize( this->worldComm().globalSize() );
    _M_first_df_globalcluster.resize( this->worldComm().globalSize() );
    _M_last_df_globalcluster.resize( this->worldComm().globalSize() );

    FEELPP_ASSERT (n_local <= n)
        ( n_local )( n )
        ( this->worldComm().globalRank() )
        ( this->worldComm().globalSize() ).error( "Invalid local vector size" );

    _M_n_dofs = n;

#ifdef HAVE_MPI
    std::vector<int> local_sizes     ( this->worldComm().size(), 0);
    std::vector<int> local_sizes_send( this->worldComm().size(), 0);

    if ( this->worldComm().size() > 1 )
        {
#if !defined(FEELPP_ENABLE_MPI_MODE)
            local_sizes_send[this->worldComm().rank()] = n_local;
#if 1
            MPI_Allreduce (&local_sizes_send[0],
                           &local_sizes[0],
                           local_sizes.size(),
                           MPI_INT,
                           MPI_SUM,
                           this->worldComm() );
#else
            mpi::all_reduce( this->worldComm(), local_sizes_send, local_sizes, std::plus<std::vector<int> >() );
#endif
            std::vector<int> M_recvcounts( this->worldComm().size() );
            std::vector<int> M_displs( this->worldComm().size() );

            int _local_index = 0;
            for ( int p=0; p<this->worldComm().size(); p++ )
                {
                    // size of data per processor
                    M_recvcounts[p] = local_sizes[p];
                    Debug( 5600 ) << "VectorUblas::init M_recvcounts[" <<p<< "]=" << M_recvcounts[p] << "\n";

                    // first index on the p-th processor
                    M_displs[p] = _local_index;
                    Debug( 5600 ) << "VectorUblas::init M_displs[" <<p<< "]=" << M_displs[p] << "\n";

                    _local_index += local_sizes[p];

                    _M_first_df[p] = _local_index;
                    _M_last_df[p] = _local_index+n_local-1;
                }
#else // defined(FEELPP_ENABLE_MPI_MODE)

            // Vincent TODO!

#endif

        }
    else // sequential
        {
            local_sizes[this->worldComm().rank()] = n_local;

            _M_first_df[this->worldComm().rank()] = 0;
            _M_last_df[this->worldComm().rank()] = n_local-1;
            // mpi
            _M_n_localWithoutGhost_df[this->worldComm().rank()]=n_local;
            _M_n_localWithGhost_df[this->worldComm().rank()]=n_local;
            _M_first_df_globalcluster[this->worldComm().rank()]=0;
            _M_last_df_globalcluster[this->worldComm().rank()]=n_local-1;
            // todo! : les map
            M_mapGlobalProcessToGlobalCluster.resize(n_local);
            M_mapGlobalClusterToGlobalProcess.resize(n_local);
        }
    if ( n == invalid_size_type_value )
        _M_n_dofs = _M_last_df[this->worldComm().rank()]+1;
#  ifdef DEBUG
    // Make sure all the local sizes sum up to the global
    // size, otherwise there is big trouble!
    int sum=0;

    for (int p=0; p< this->worldComm().size(); p++)
        sum += local_sizes[p];

    if ( n != invalid_size_type_value )
        FEELPP_ASSERT (sum == static_cast<int>(n))
            ( sum )( n )
            ( this->worldComm().rank() )
            ( this->worldComm().size() ).warn( "invalid distributed vector construction" );

#  endif

#else

    // No other options without MPI!
    if (n != n_local)
        {
            std::cerr << "ERROR:  MPI is required for n != n_local!"
                      << std::endl;
            //error();
        }

#endif

    /*
    Debug( 5600 ) << "        global size = " << this->size() << "\n";
    Debug( 5600 ) << "        local  size = " << this->localSize() << "\n";
    Debug( 5600 ) << "  first local index = " << this->firstLocalIndex() << "\n";
    Debug( 5600 ) << "   last local index = " << this->lastLocalIndex() << "\n";
    */

}

DataMap::DataMap( DataMap const & dm )
    :
    M_closed( dm.M_closed ),
    _M_n_dofs( dm._M_n_dofs ),
    _M_n_localWithoutGhost_df( dm._M_n_localWithoutGhost_df ),
    _M_n_localWithGhost_df( dm._M_n_localWithGhost_df ),
    _M_first_df( dm._M_first_df ),
    _M_last_df( dm._M_last_df ),
    _M_first_df_globalcluster( dm._M_first_df_globalcluster),
    _M_last_df_globalcluster( dm._M_last_df_globalcluster),
    M_myglobalelements(),
    M_mapGlobalProcessToGlobalCluster(dm.M_mapGlobalProcessToGlobalCluster),
    M_mapGlobalClusterToGlobalProcess(dm.M_mapGlobalClusterToGlobalProcess),
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
            _M_n_dofs = dm._M_n_dofs;
            _M_n_localWithoutGhost_df = dm._M_n_localWithoutGhost_df;
            _M_n_localWithGhost_df = dm._M_n_localWithGhost_df;
            _M_first_df = dm._M_first_df;
            _M_last_df = dm._M_last_df;
            _M_first_df_globalcluster = dm._M_first_df_globalcluster;
            _M_last_df_globalcluster = dm._M_last_df_globalcluster;
            M_myglobalelements = dm.M_myglobalelements;
            M_mapGlobalProcessToGlobalCluster = dm.M_mapGlobalProcessToGlobalCluster;
            M_mapGlobalClusterToGlobalProcess = dm.M_mapGlobalClusterToGlobalProcess;
        }
    return *this;
}
void
DataMap::close() const
{
    // we assume here that the data is contiguous
    M_myglobalelements.resize( nMyElements() );
    for( size_type i = 0; i < nMyElements(); ++i )
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
DataMap::setNDof(size_type ndof) { _M_n_dofs=ndof; }

void
DataMap::setNLocalDofWithoutGhost(const size_type proc, const size_type n, bool inWorld)
{
    _M_n_localWithoutGhost_df[this->worldComm().globalRank()]=n;
}
void
DataMap::setNLocalDofWithGhost(const size_type proc, const size_type n, bool inWorld)
{
    _M_n_localWithGhost_df[this->worldComm().globalRank()]=n;
}
void
DataMap::setFirstDof(const size_type proc, const size_type df, bool inWorld)
{
    FEELPP_ASSERT(proc < _M_first_df.size())( proc )( _M_first_df.size() ).error( "invalid proc id or dof table" );
    _M_first_df[this->worldComm().globalRank()]=df;
}
void
DataMap::setLastDof(const size_type proc, const size_type df, bool inWorld)
{
    FEELPP_ASSERT(proc < _M_first_df.size())( proc )( _M_first_df.size() ).error( "invalid proc id or dof table" );
    _M_last_df[this->worldComm().globalRank()]=df;
}
void
DataMap::setFirstDofGlobalCluster(const size_type proc, const size_type df, bool inWorld)
{
    FEELPP_ASSERT(proc < _M_first_df_globalcluster.size())( proc )( _M_first_df_globalcluster.size() ).error( "invalid proc id or dof table" );
    _M_first_df_globalcluster[this->worldComm().globalRank()]=df;
}
void
DataMap::setLastDofGlobalCluster(const size_type proc, const size_type df, bool inWorld)
{
    FEELPP_ASSERT(proc < _M_first_df_globalcluster.size())( proc )( _M_first_df_globalcluster.size() ).error( "invalid proc id or dof table" );
    _M_last_df_globalcluster[this->worldComm().globalRank()]=df;
}


void
DataMap::setMapGlobalProcessToGlobalCluster( std::vector<size_type> const& map) { M_mapGlobalProcessToGlobalCluster=map; };
void
DataMap::setMapGlobalClusterToGlobalProcess( std::vector<size_type> const& map) { M_mapGlobalClusterToGlobalProcess=map; };
void
DataMap::setMapGlobalProcessToGlobalCluster( size_type i, size_type j) { M_mapGlobalProcessToGlobalCluster[i]=j; };
void
DataMap::setMapGlobalClusterToGlobalProcess( size_type i, size_type j) { M_mapGlobalClusterToGlobalProcess[i]=j; };
void
DataMap::resizeMapGlobalProcessToGlobalCluster( size_type n) { M_mapGlobalProcessToGlobalCluster.resize(n); };
void
DataMap::resizeMapGlobalClusterToGlobalProcess( size_type n) { M_mapGlobalClusterToGlobalProcess.resize(n); };

void
DataMap::updateDataInWorld()
{
    mpi::all_gather( this->worldComm().globalComm(),
                     this->_M_n_localWithoutGhost_df[this->worldComm().globalRank()],
                     this->_M_n_localWithoutGhost_df );
    mpi::all_gather( this->worldComm().globalComm(),
                     this->_M_n_localWithGhost_df[this->worldComm().globalRank()],
                     this->_M_n_localWithGhost_df );

    mpi::all_gather( this->worldComm().globalComm(),
                     this->_M_first_df[this->worldComm().globalRank()],
                     this->_M_first_df );
    mpi::all_gather( this->worldComm().globalComm(),
                     this->_M_last_df[this->worldComm().globalRank()],
                     this->_M_last_df );
    mpi::all_gather( this->worldComm().globalComm(),
                     this->_M_first_df_globalcluster[this->worldComm().globalRank()],
                     this->_M_first_df_globalcluster );
    mpi::all_gather( this->worldComm().globalComm(),
                     this->_M_last_df_globalcluster[this->worldComm().globalRank()],
                     this->_M_last_df_globalcluster );

    //_globalcluster
}


void
DataMap::showMeMapGlobalProcessToGlobalCluster( std::ostream& __out  ) const
{
    __out << std::endl;
    //this->comm().barrier();
    this->comm().godComm().barrier();

    //for (int proc = 0;proc<this->comm().size();++proc
    for (int proc = 0;proc<this->comm().godSize();++proc)
        {
            //if (proc==this->comm().rank())
            if (proc==this->comm().godRank())
                {
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
#if 0
                    for (size_type i=0 ; i<this->mapGlobalProcessToGlobalCluster().size() ;++i)
                        {
                            __out << i << " " << this->mapGlobalProcessToGlobalCluster()[i]
                                  << " real proc " << procOnGlobalCluster(/*this->*/mapGlobalProcessToGlobalCluster()[i]) <<"\n";
                        }
                    __out << "-----------------------------------------------------------------------\n";
#endif
#if 1
                    __out << " _M_first_df : ";
                    for (int i=0;i<this->worldComm().globalSize();++i)
                        {
                            __out << this-> _M_first_df[i] << " ";
                        }
                    __out << "\n";
                    __out << " _M_last_df : ";
                    for (int i=0;i<this->worldComm().globalSize();++i)
                        {
                            __out << this-> _M_last_df[i] << " ";
                        }
                    __out << "\n";
                    __out << " _M_first_df_globalcluster : ";
                    for (int i=0;i<this->worldComm().globalSize();++i)
                        {
                            __out << this-> _M_first_df_globalcluster[i] << " ";
                        }
                    __out << "\n";
                    __out << " _M_last_df_globalcluster : ";
                    for (int i=0;i<this->worldComm().globalSize();++i)
                        {
                            __out << this-> _M_last_df_globalcluster[i] << " ";
                        }
                    __out << "\n";

#endif
                    __out << "-----------------------------------------------------------------------\n";

                    __out << std::endl;
                }
            //this->comm().barrier();
            this->comm().godComm().barrier();
        }


}

} // namespace Feel

