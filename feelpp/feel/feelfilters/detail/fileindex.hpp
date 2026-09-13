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
#ifndef FEELPP_INDEX_HPP
#define FEELPP_INDEX_HPP 1

#include <feel/feelcore/commobject.hpp>

namespace Feel
{

namespace detail
{
//! \brief EnSight Gold transient-file index with 64-bit block offsets.
//! Reads/writes the standard trailing FILE_INDEX without changing its layout.
//! The offset array is transferred in one MPI request, not one per time step.
//! @author Christophe Prud'homme
//! @see
class FileIndex : public CommObject
{
  public:
    //! \brief Bind index operations to the file's communicator.
    explicit FileIndex( worldcomm_ptr_t const &w );

    //! \brief Collectively read and broadcast a validated trailing index.
    //! Empty/non-indexed files yield an undefined index. Invalid indexed files
    //! or MPI-I/O failures abort the communicator rather than stranding peers.
    void read( MPI_File fh );

    //! \brief Write the index on the master and advance its byte offset.
    void write( MPI_File fh, MPI_Offset &offset );

    //! \brief Append the byte offset immediately after BEGIN TIME STEP.
    void add( int64_type tellp ) { M_fileblocks.push_back( tellp ); }

    //! \return Whether one or more temporal blocks are recorded.
    bool defined() const { return !M_fileblocks.empty(); }

    //! \return Ordered block offsets in the binary file.
    std::vector<int64_type> const &fileBlocks() const { return M_fileblocks; }

    //! \return Number of temporal blocks.
    int64_type numberOfBlock() const { return M_fileblocks.size(); }

    //! \return Position at which the old index can be replaced by a new step.
    int64_type nextFreePosFile() const { return M_nextFreePosFile; }

    //! \brief Discard a previous pack's offsets and append position, retaining capacity.
    void clear()
    {
        M_fileblocks.clear();
        M_nextFreePosFile = -1;
    }

  private:
    std::vector<int64_type> M_fileblocks;
    int64_type M_nextFreePosFile;
    int M_sizeOfInt32_t, M_sizeOfInt64_t;
};
} // namespace detail
} // namespace Feel
#endif /* FEELPP_INDEX_HPP */
