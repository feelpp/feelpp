/*

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
#include <limits>

/* handle cases where we are not using the 2.2 MPI Standard */
#if !defined( MPI_INT32_T ) && defined( FEELPP_MPI_INT32 )
#define MPI_INT32_T FEELPP_MPI_INT32
#endif
#if !defined( MPI_INT64_T ) && defined( FEELPP_MPI_INT64 )
#define MPI_INT64_T FEELPP_MPI_INT64
#endif

namespace Feel
{
namespace detail
{

namespace
{
//! \brief Fail collectively in effect, without adding a collective per I/O call.
void checkIndexIo( int error, WorldComm const &comm )
{
    if ( error == MPI_SUCCESS )
        return;
    char message[MPI_MAX_ERROR_STRING];
    int length = 0;
    MPI_Error_string( error, message, &length );
    std::cerr << "EnSight FILE_INDEX: " << std::string( message, length ) << std::endl;
    MPI_Abort( comm.comm(), error );
    std::terminate();
}

//! \brief Read exactly count entries; short reads are corrupt indexes.
void readIndexAt( MPI_File file, MPI_Offset offset, void *data, int count, MPI_Datatype type,
                  MPI_Status *status, WorldComm const &comm )
{
    checkIndexIo( MPI_File_read_at( file, offset, data, count, type, status ), comm );
    int actual = 0;
    checkIndexIo( MPI_Get_count( status, type, &actual ), comm );
    if ( actual != count )
        checkIndexIo( MPI_ERR_IO, comm );
}

//! \brief Write exactly count entries without allowing rank-local exceptions.
void writeIndexAt( MPI_File file, MPI_Offset offset, void const *data, int count, MPI_Datatype type,
                   MPI_Status *status, WorldComm const &comm )
{
    checkIndexIo( MPI_File_write_at( file, offset, data, count, type, status ), comm );
    int actual = 0;
    checkIndexIo( MPI_Get_count( status, type, &actual ), comm );
    if ( actual != count )
        checkIndexIo( MPI_ERR_IO, comm );
}
} // namespace

FileIndex::FileIndex( worldcomm_ptr_t const &w ) : CommObject( w ), M_nextFreePosFile( -1 )
{
    MPI_Type_size( MPI_INT32_T, &M_sizeOfInt32_t );
    MPI_Type_size( MPI_INT64_T, &M_sizeOfInt64_t );
}

void FileIndex::read( MPI_File fh )
{
    char buffer[80] = {};
    clear();

    MPI_Offset offset;
    MPI_Offset prevOffset;

    MPI_Status status;

    LOG( INFO ) << "Start reading FILE_INDEX (if any)";
    /* Stored current position */
    // MPI_File_get_position_shared(fh, &prevOffset);
    /* Check file size */
    if ( this->worldComm().comm().rank() == this->worldComm().masterRank() )
        checkIndexIo( MPI_File_get_size( fh, &offset ), this->worldComm() );
    mpi::broadcast( this->worldComm().comm(), offset, this->worldComm().masterRank() );
    LOG( INFO ) << "file length: " << offset;

    if ( offset < 80 )
    {
        this->worldComm().comm().barrier();
        LOG( INFO ) << "no FILE_INDEX (stop trying to read it)";
        return;
    }

    if ( this->worldComm().comm().rank() == this->worldComm().masterRank() )
    {
        // read last line FILE_INDEX
        // MPI_File_seek(fh, -80, MPI_SEEK_END);
        readIndexAt( fh, offset - 80, buffer, 80, MPI_CHAR, &status, this->worldComm() );

        if ( strncmp( buffer, "FILE_INDEX", 10 ) == 0 )
        {
            LOG( INFO ) << "found FILE_INDEX";
            // right before the FILE_INDEX entry we find the address of the index start
            // MPI_File_seek(fh, -80-sizeof(int64_type), MPI_SEEK_END);

            MPI_Offset const fileSize = offset;
            if ( fileSize < 80 + M_sizeOfInt64_t + 2 * M_sizeOfInt32_t )
                checkIndexIo( MPI_ERR_IO, this->worldComm() );
            int64_type addr;
            readIndexAt( fh, offset - 80 - M_sizeOfInt64_t, &addr, 1, MPI_INT64_T, &status,
                         this->worldComm() );
            if ( addr < 0 || addr > fileSize - 80 - M_sizeOfInt64_t - 2 * M_sizeOfInt32_t )
                checkIndexIo( MPI_ERR_IO, this->worldComm() );
            this->M_nextFreePosFile = addr;

            // MPI_File_seek(fh, addr, MPI_SEEK_SET);
            offset = addr;

            int32_type nBlock;
            readIndexAt( fh, offset, &nBlock, 1, MPI_INT32_T, &status, this->worldComm() );
            if ( nBlock < 0 || addr + MPI_Offset( nBlock ) * M_sizeOfInt64_t + 2 * M_sizeOfInt32_t +
                                       M_sizeOfInt64_t + 80 !=
                                   fileSize )
                checkIndexIo( MPI_ERR_IO, this->worldComm() );
            LOG( INFO ) << "read in FILE_INDEX number of steps: " << nBlock;
            // need some check here regarding the number of time steps probably
            offset += M_sizeOfInt32_t;
            // now we can read the fileblocks
            M_fileblocks.clear();
            M_fileblocks.resize( nBlock );
            readIndexAt( fh, offset, M_fileblocks.data(), nBlock, MPI_INT64_T, &status,
                         this->worldComm() );
            offset += MPI_Offset( nBlock ) * M_sizeOfInt64_t;
            int64_type previous = -1;
            for ( auto block : M_fileblocks )
            {
                if ( block <= previous || block >= addr )
                    checkIndexIo( MPI_ERR_IO, this->worldComm() );
                previous = block;
            }

            int32_type flag;
            readIndexAt( fh, offset, &flag, 1, MPI_INT32_T, &status, this->worldComm() );
            offset += M_sizeOfInt32_t;

            if ( flag != 0 )
                checkIndexIo( MPI_ERR_IO, this->worldComm() );
            LOG( INFO ) << "Done reading FILE_INDEX";
        }
    }

    mpi::broadcast( this->worldComm().comm(), M_fileblocks, this->worldComm().masterRank() );
    mpi::broadcast( this->worldComm().comm(), M_nextFreePosFile, this->worldComm().masterRank() );
}

void FileIndex::write( MPI_File fh, MPI_Offset &offset )
{
    int size;
    char buffer[80];

    // MPI_Offset offset;
    MPI_Status status;

    LOG( INFO ) << "Start writing FILE_INDEX";
    // The next append overwrites this trailer, not the end of the whole file.
    M_nextFreePosFile = offset;

    /* only process 0 writes data */
    if ( this->worldComm().comm().rank() == this->worldComm().masterRank() )
    {
        // go to end of file
        // MPI_File_seek_shared(fh, 0, MPI_SEEK_END);

        // MPI_File_get_size(fh, &offset);
        //  get position of fileblock for number of steps
        //  MPI_File_get_position_shared(fh, &offset);
        int64_type fb_n_step = offset;
        if ( M_fileblocks.size() > std::size_t( std::numeric_limits<int32_type>::max() ) ||
             offset < 0 )
            checkIndexIo( MPI_ERR_COUNT, this->worldComm() );
        int32_type n = this->numberOfBlock();

        // write number of steps
        writeIndexAt( fh, offset, &n, 1, MPI_INT32_T, &status, this->worldComm() );
        offset += M_sizeOfInt32_t;
        // MPI_File_write_ordered(fh, &n, size, MPI_INT32_T, &status);
        LOG( INFO ) << "Writing " << this->numberOfBlock() << " fileblocks in FILE_INDEX";

        // write fileblocks stored
        writeIndexAt( fh, offset, M_fileblocks.data(), n, MPI_INT64_T, &status, this->worldComm() );
        offset += MPI_Offset( n ) * M_sizeOfInt64_t;

        // write 32bit integer flag (set to 0)
        int32_type flag = 0;
        writeIndexAt( fh, offset, &flag, 1, MPI_INT32_T, &status, this->worldComm() );
        offset += M_sizeOfInt32_t;
        // MPI_File_write_ordered(fh, &flag, size, MPI_INT32_T, &status);

        // write position of fileblock for number of steps
        writeIndexAt( fh, offset, &fb_n_step, 1, MPI_INT64_T, &status, this->worldComm() );
        offset += M_sizeOfInt64_t;
        // MPI_File_write_ordered(fh, &fb_n_step, size, MPI_INT64_T, &status);

        // write string FILE_INDEX
        memset( buffer, '\0', sizeof( buffer ) );
        strcpy( buffer, "FILE_INDEX" );

        writeIndexAt( fh, offset, buffer, sizeof( buffer ), MPI_CHAR, &status, this->worldComm() );
        offset += 80;
        // MPI_File_write_ordered(fh, buffer, size, MPI_CHAR, &status);
        LOG( INFO ) << "Done writing FILE_INDEX";
    }
}

} // namespace detail

} // namespace Feel
