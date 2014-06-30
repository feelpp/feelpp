/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

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
#ifndef FEELPP_INDEX_HPP
#define FEELPP_INDEX_HPP 1

#include <feel/feelcore/environment.hpp>

namespace Feel {

namespace detail {
/**
 * \class FileIndex
 *
 * @author Christophe Prud'homme
 * @see
 */
class FileIndex
{
public:
    FileIndex() : n_steps( 0 ), fileblock_n_steps( -1 ) {}
    void read( MPI_File fh );
    void write( MPI_File fh );
    void add( int64_type tellp ) { ++n_steps ; fileblocks.push_back( tellp ); }
    bool defined() const { return n_steps != 0; }
    int64_type n_steps;
    std::vector<int64_type> fileblocks;
    int64_type fileblock_n_steps;
};
} // detail
} // Feel
#endif /* FEELPP_INDEX_HPP */
