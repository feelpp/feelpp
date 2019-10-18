/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-03-10

  Copyright (C) 2014-2016 Feel++ Consortium

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
#include <feel/feelcore/feel.hpp>
#include <feel/feelfilters/detail/fileindex.hpp>

/* handle cases where we are not using the 2.2 MPI Standard */
#if !defined(MPI_INT32_T) && defined(FEELPP_MPI_INT32)
    #define MPI_INT32_T FEELPP_MPI_INT32
#endif
#if !defined(MPI_INT64_T) && defined(FEELPP_MPI_INT64)
    #define MPI_INT64_T FEELPP_MPI_INT64
#endif

namespace Feel { namespace detail {

FileIndex::FileIndex( worldcomm_ptr_t const& w )
    :
    CommObject( w ),
    M_nextFreePosFile( -1 )
{
    MPI_Type_size( MPI_INT32_T , &M_sizeOfInt32_t );
    MPI_Type_size( MPI_INT64_T , &M_sizeOfInt64_t );
}

void FileIndex::read( MPI_File fh )
{
    char buffer[80];

    MPI_Offset offset;
    MPI_Offset prevOffset;

    MPI_Status status;

    LOG(INFO) << "Start reading FILE_INDEX (if any)";
    /* Stored current position */
    // MPI_File_get_position_shared(fh, &prevOffset);
    /* Check file size */
    MPI_File_get_size(fh, &offset);
    LOG(INFO) << "file length: " << offset;

    if ( offset <= 0 )
    {
        this->worldComm().barrier();
        LOG(INFO) << "no FILE_INDEX (stop trying to read it)";
        return;
    }


    if ( this->worldComm().isMasterRank() )
    {
        // read last line FILE_INDEX
        //MPI_File_seek(fh, -80, MPI_SEEK_END);
        MPI_File_read_at(fh, offset-80, buffer, 80, MPI_CHAR, &status);

        LOG(INFO) <<"buffer:" << buffer;
        if (strncmp(buffer, "FILE_INDEX", 10) == 0)
        {
            LOG(INFO) << "found FILE_INDEX";
            // right before the FILE_INDEX entry we find the address of the index start
            //MPI_File_seek(fh, -80-sizeof(int64_type), MPI_SEEK_END);

            int64_type addr;
            MPI_File_read_at(fh, offset-80-M_sizeOfInt64_t, &addr, 1, MPI_INT64_T, &status);
            this->M_nextFreePosFile = addr;

            // MPI_File_seek(fh, addr, MPI_SEEK_SET);
            offset = addr;

            int32_type nBlock;
            MPI_File_read_at(fh, offset, &nBlock, 1, MPI_INT32_T, &status);
            LOG(INFO) << "read in FILE_INDEX number of steps: " << nBlock;
            // need some check here regarding the number of time steps probably
            offset += M_sizeOfInt32_t;
            // now we can read the fileblocks
            M_fileblocks.clear();
            for( int i = 0; i < nBlock; ++i )
            {
                int64_type fb;
                MPI_File_read_at(fh, offset,  &fb, 1, MPI_INT64_T, &status);
                M_fileblocks.push_back( fb );
                offset+=M_sizeOfInt64_t;
            }

            int32_type flag;
            MPI_File_read_at(fh, offset, &flag, 1, MPI_INT32_T, &status);
            offset+=M_sizeOfInt32_t;

            CHECK( flag == 0 ) << "invalid FILE_INDEX, flag must be equal to 0 but we have " << flag ;
            LOG(INFO) << "Done reading FILE_INDEX";
        }
    }

    mpi::broadcast( this->worldComm(), M_fileblocks, this->worldComm().masterRank() );
    mpi::broadcast( this->worldComm(), M_nextFreePosFile, this->worldComm().masterRank() );
}

void FileIndex::write( MPI_File fh,  MPI_Offset & offset )
{
    int size;
    char buffer[80];

    //MPI_Offset offset;
    MPI_Status status;

    LOG(INFO) << "Start writing FILE_INDEX";

    /* only process 0 writes data */
    if ( this->worldComm().isMasterRank() )
    {
        // go to end of file
        // MPI_File_seek_shared(fh, 0, MPI_SEEK_END);

        //MPI_File_get_size(fh, &offset);
        // get position of fileblock for number of steps
        // MPI_File_get_position_shared(fh, &offset);
        int64_type fb_n_step = offset;
        int32_type n = this->numberOfBlock();

        // write number of steps
        MPI_File_write_at(fh, offset, &n, 1, MPI_INT32_T, &status);
        offset+=M_sizeOfInt32_t;
        //MPI_File_write_ordered(fh, &n, size, MPI_INT32_T, &status);
        LOG(INFO) << "Writing "<< this->numberOfBlock() << " fileblocks in FILE_INDEX";

        // write fileblocks stored
        for( int64_type fb : M_fileblocks )
        {
            MPI_File_write_at(fh, offset, &fb, 1, MPI_INT64_T, &status);
            offset+=M_sizeOfInt64_t;
            //MPI_File_write_ordered(fh, &fb, size, MPI_INT64_T, &status);
        }

        // write 32bit integer flag (set to 0)
        int32_type flag = 0;
        MPI_File_write_at(fh, offset, &flag, 1, MPI_INT32_T, &status);
        offset+=M_sizeOfInt32_t;
        //MPI_File_write_ordered(fh, &flag, size, MPI_INT32_T, &status);

        // write position of fileblock for number of steps
        MPI_File_write_at(fh, offset, &fb_n_step, 1, MPI_INT64_T, &status);
        offset+=M_sizeOfInt64_t;
        //MPI_File_write_ordered(fh, &fb_n_step, size, MPI_INT64_T, &status);

        // write string FILE_INDEX
        memset(buffer, '\0', sizeof(buffer));
        strcpy( buffer, "FILE_INDEX" );

        MPI_File_write_at(fh, offset, buffer, sizeof(buffer), MPI_CHAR, &status);
        offset+=80;
        //MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        LOG(INFO) << "Done writing FILE_INDEX";
    }
}

} // detail

} // Feel
