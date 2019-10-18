/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

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

namespace Feel {

namespace detail {
/**
 * \class FileIndex
 *
 * @author Christophe Prud'homme
 * @see
 */
class FileIndex : public CommObject
{
public:
    FileIndex( worldcomm_ptr_t const& w );
    void read( MPI_File fh );
    void write( MPI_File fh,  MPI_Offset & offset );
    void add( int64_type tellp ) { M_fileblocks.push_back( tellp ); }
    bool defined() const { return !M_fileblocks.empty(); }

    std::vector<int64_type> const& fileBlocks() const { return M_fileblocks; }
    int64_type numberOfBlock() const { return M_fileblocks.size(); }
    int64_type nextFreePosFile() const { return M_nextFreePosFile; }
private :
    std::vector<int64_type> M_fileblocks;
    int64_type M_nextFreePosFile;
    int M_sizeOfInt32_t, M_sizeOfInt64_t;
};
} // detail
} // Feel
#endif /* FEELPP_INDEX_HPP */
