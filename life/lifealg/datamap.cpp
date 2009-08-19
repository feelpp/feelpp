/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-28

  Copyright (C) 2008, 2009 Université Joseph Fourier (Grenoble I)

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
#include <life/lifealg/datamap.hpp>

namespace Life
{
DataMap::DataMap()
    :
    M_closed( false ),
    _M_n_dofs( 0 ),
    _M_first_df( Application::nProcess()),
    _M_last_df( Application::nProcess() ),
    M_myglobalelements()
{}
DataMap::DataMap( size_type n, size_type n_local )
    :
    M_closed( false ),
    _M_n_dofs( n ),
    _M_first_df( Application::nProcess()),
    _M_last_df( Application::nProcess() ),
    M_myglobalelements()
{
    LIFE_ASSERT (n_local <= n)
        ( n_local )( n )
        ( Application::processId() )
        ( Application::nProcess() ).error( "Invalid local vector size" );

    _M_n_dofs = n;


    size_type M_local_size  = n_local;
    size_type M_first_local_index = 0;

#ifdef HAVE_MPI
    std::vector<int> local_sizes     ( Application::nProcess(), 0);
    std::vector<int> local_sizes_send( Application::nProcess(), 0);

    if ( Application::nProcess() > 1 )
        {
            local_sizes_send[Application::processId()] = n_local;
#if 1
            MPI_Allreduce (&local_sizes_send[0],
                           &local_sizes[0],
                           local_sizes.size(),
                           MPI_INT,
                           MPI_SUM,
                           Application::COMM_WORLD );
#else
            mpi::all_reduce( Application::comm(), local_sizes_send, local_sizes, std::plus<std::vector<int> >() );
#endif
            std::vector<int> M_recvcounts( Application::nProcess() );
            std::vector<int> M_displs( Application::nProcess() );

            // _first_local_index is the sum of _local_size
            // for all processor ids less than ours
            for ( int p=0; p<Application::processId(); p++ )
                {
                    M_first_local_index += local_sizes[p];
                }

            int _local_index = 0;
            for ( int p=0; p<Application::nProcess(); p++ )
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

        }
    else
        {
            M_first_local_index = 0;
            local_sizes[Application::processId()] = n_local;

            _M_first_df[Application::processId()] = 0;
            _M_last_df[Application::processId()] = n_local-1;
        }
    if ( n == invalid_size_type_value )
        _M_n_dofs = _M_last_df[Application::processId()]+1;
#  ifdef DEBUG
    // Make sure all the local sizes sum up to the global
    // size, otherwise there is big trouble!
    int sum=0;

    for (int p=0; p< Application::nProcess(); p++)
        sum += local_sizes[p];

    if ( n != invalid_size_type_value )
        LIFE_ASSERT (sum == static_cast<int>(n))
            ( sum )( n )
            ( Application::processId() )
            ( Application::nProcess() ).warn( "invalid distributed vector construction" );

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

    size_type M_last_local_index = M_first_local_index + M_local_size;


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
    _M_first_df( dm._M_first_df ),
    _M_last_df( dm._M_last_df ),
    M_myglobalelements()
{}
DataMap::~DataMap()
{}

DataMap&
DataMap::operator=( DataMap const& dm )
{
    if ( this != &dm )
        {
            M_closed = dm.M_closed;
            _M_n_dofs = dm._M_n_dofs;
            _M_first_df = dm._M_first_df;
            _M_last_df = dm._M_last_df;
            M_myglobalelements = dm.M_myglobalelements;

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

}
