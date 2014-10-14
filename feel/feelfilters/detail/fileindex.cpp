/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-03-10

  Copyright (C) 2014 Feel++ Consortium

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

void FileIndex::read( MPI_File fh )
{
    char buffer[80];

    MPI_Offset offset;
    MPI_Offset prevOffset;

    MPI_Status status;

    LOG(INFO) << "Start reading FILE_INDEX (if any)";
    /* Stored current position */
    MPI_File_get_position_shared(fh, &prevOffset);
    
    /* Check file size */
    MPI_File_get_size(fh, &offset);

    LOG(INFO) << "file length: " << offset;

    if ( offset <= 0 )
    {
        LOG(INFO) << "no FILE_INDEX (stop trying to read it)";
        return;
    }

    // read last line FILE_INDEX
    MPI_File_seek(fh, -80, MPI_SEEK_END);
    MPI_File_read_all(fh, buffer, 80, MPI_CHAR, &status);

    LOG(INFO) <<"buffer:" << buffer;
    if (strncmp(buffer, "FILE_INDEX", 10) == 0)
    {
        LOG(INFO) << "found FILE_INDEX";
        // right before the FILE_INDEX entry we find the address of the index start
        MPI_File_seek(fh, -80-sizeof(int64_type), MPI_SEEK_END);

        int64_type addr;
        MPI_File_read_all(fh, &addr, 1, MPI_INT64_T, &status);
        this->fileblock_n_steps = addr;

        MPI_File_seek(fh, addr, MPI_SEEK_SET);

        int32_type n;
        MPI_File_read_all(fh, &n, 1, MPI_INT32_T, &status);
        this->n_steps = n;
        LOG(INFO) << "read in FILE_INDEX number of steps: " << this->n_steps;
        // need some check here regarding the number of time steps probably

        // now we can read the fileblocks
        this->fileblocks.clear();
        for( int i = 0; i < n_steps; ++i )
        {
            int64_type fb;
            MPI_File_read_all(fh, &fb, 1, MPI_INT64_T, &status);
            this->fileblocks.push_back( fb );
        }

        int32_type flag;
        MPI_File_read_all(fh, &flag, 1, MPI_INT32_T, &status);
        CHECK( flag == 0 ) << "invalid FILE_INDEX, flag must be equal to 0";
        LOG(INFO) << "Done reading FILE_INDEX";
    }

    /* restore previous offset */
    MPI_File_seek_shared(fh, prevOffset, MPI_SEEK_SET);
}

void FileIndex::write( MPI_File fh )
{
    int size;
    char buffer[80];

    MPI_Offset offset;
    MPI_Status status;

    LOG(INFO) << "Start writing FILE_INDEX";

    /* only process 0 writes data */
    if( Environment::isMasterRank() ) 
    { size = 1; }
    else 
    { size = 0; }

    // go to end of file
    MPI_File_seek_shared(fh, 0, MPI_SEEK_END);

    // get position of fileblock for number of steps
    MPI_File_get_position_shared(fh, &offset);
    int64_type fb_n_step = offset;
    int32_type n = this->n_steps;

    // write number of steps
    MPI_File_write_ordered(fh, &n, size, MPI_INT32_T, &status);
    LOG(INFO) << "Writing "<< this->n_steps << " fileblocks in FILE_INDEX";
    LOG(INFO) << "Writing "<< this->fileblocks.size() << " fileblocks in FILE_INDEX";

    // write fileblocks stored
    for( int64_type fb : this->fileblocks )
    {
        MPI_File_write_ordered(fh, &fb, size, MPI_INT64_T, &status);
    }

    LOG(INFO) << "Writing flag==0 in FILE_INDEX";
    // write 32bit integer flag (set to 0)
    int32_type flag = 0;
    MPI_File_write_ordered(fh, &flag, size, MPI_INT32_T, &status);
    
    // write position of fileblock for number of steps
    MPI_File_write_ordered(fh, &fb_n_step, size, MPI_INT64_T, &status);

    // write string FILE_INDEX
    if( Environment::isMasterRank() ) 
    { size = sizeof(buffer); }
    else 
    { size = 0; }
    memset(buffer, '\0', sizeof(buffer));
    strcpy( buffer, "FILE_INDEX" );
    MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
    LOG(INFO) << "Done writing FILE_INDEX";
}

} // detail

} // Feel
