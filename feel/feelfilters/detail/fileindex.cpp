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

namespace Feel { namespace detail {

void FileIndex::read( std::istream& is )
{
    LOG(INFO) << "Start reading FILE_INDEX (if any)";
    is.seekg (0, is.end);
    int length = is.tellg();
    is.seekg (0, is.beg);
    LOG(INFO) << "file length: " << length;
    if ( length <= 0 )
    {
        LOG(INFO) << "no FILE_INDEX (stop trying to read it)";
        return;
    }

    // read last line FILE_INDEX
    is.seekg(-80*sizeof(char), std::ios::end );
    char buffer[80];
    is.read( (char*)&buffer, sizeof(buffer) );
    LOG(INFO) <<"buffer:" << buffer;
    if (strncmp(buffer, "FILE_INDEX", 10) == 0)
    {
        LOG(INFO) << "found FILE_INDEX";
        // right before the FILE_INDEX entry we find the address of the index start
        is.seekg(-80*sizeof(char)-sizeof(int64_type), std::ios::end);
        int64_type addr;
        is.read( (char*)&addr, sizeof(int64_type) );
        this->fileblock_n_steps = addr;

        is.seekg(addr, std::ios::beg);
        int32_type n;
        is.read( (char*)&n, sizeof(int32_type) );
        this->n_steps = n;
        LOG(INFO) << "read in FILE_INDEX number of steps: " << this->n_steps;
        // need some check here regarding the number of time steps probably

        // now we can read the fileblocks
        this->fileblocks.clear();
        for( int i = 0; i < n_steps; ++i )
        {
            int64_type fb;
            is.read( (char*)&fb, sizeof(int64_type) );
            this->fileblocks.push_back( fb );
        }
        int32_type flag;
        is.read( (char*)&flag, sizeof(int32_type) );
        CHECK( flag == 0 ) << "invalid FILE_INDEX, flag must be equal to 0";
        LOG(INFO) << "Done reading FILE_INDEX";
    }

}

void FileIndex::write( std::ostream& os )
{
    LOG(INFO) << "Start writing FILE_INDEX";
    // go to end of file
    os.seekp( 0, std::ios::end );
    // get position of fileblock for number of steps
    int64_type fb_n_step = os.tellp();
    int32_type n = this->n_steps;
    // write number of steps
    os.write( (char*)&n, sizeof(int32_type) );
    LOG(INFO) << "Writing "<< this->n_steps << " fileblocks in FILE_INDEX";
    LOG(INFO) << "Writing "<< this->fileblocks.size() << " fileblocks in FILE_INDEX";
    // write fileblocks stored
    for( int64_type fb : this->fileblocks )
    {
        os.write( (char*)&fb, sizeof(int64_type) );
    }
    LOG(INFO) << "Writing flag==0 in FILE_INDEX";
    // write 32bit integer flag (set to 0)
    int32_type flag = 0;
    os.write( (char*)&flag, sizeof(int32_type) );
    // write position of fileblock for number of steps
    os.write( (char*)&fb_n_step, sizeof(int64_type) );

    // write string FILE_INDEX
    char buffer[80];
    strcpy( buffer, "FILE_INDEX" );
    os.write( (char*)&buffer, sizeof(buffer) );
    LOG(INFO) << "Done writing FILE_INDEX";
}

} // detail

} // Feel
