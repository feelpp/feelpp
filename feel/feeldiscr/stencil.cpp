/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2012-04-07

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
/**
   \file stencil.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-04-07
 */
#include <sstream>
#include <feel/feeldiscr/stencil.hpp>

namespace Feel
{
std::vector<size_type> default_block_pattern( 1,size_type( Pattern::HAS_NO_BLOCK_PATTERN ) );

void
stencilManagerGarbageCollect()
{
    BOOST_FOREACH( StencilManagerImpl::value_type& entry, StencilManager::instance() )
    {
        // each entry is a pair of tuple and graph

        if ( entry.second.unique() )
        {
            std::ostringstream ostr;
            std::for_each( entry.first.get<3>().begin(),
                           entry.first.get<3>().end(),
                           [&ostr]( size_type i ) { ostr << i << ","; } );
            std::cout << "[stencilManagerGarbageCollect] deleting entry space: ( "
                      << entry.first.get<0>()
                      << "," << entry.first.get<1>()
                      << "," << int(entry.first.get<2>())
                      << ",( " << ostr.str()
                      << "),"
                      << "," << int(entry.first.get<4>())
                      << " ): "
                      << entry.second << "\n";

            StencilManager::instance().erase( entry.first );
        }
    }
}

void
stencilManagerPrint()
{
    std::cout << "********************************************************************************\n";
    if ( StencilManager::instance().empty() )
    {
        std::cout << "[stencilManagerPrint] no entries in StencilManager\n";
    }
    BOOST_FOREACH( StencilManagerImpl::value_type& entry, StencilManager::instance() )
    {
        std::ostringstream ostr;
        std::for_each( entry.first.get<3>().begin(),
                       entry.first.get<3>().end(),
                       [&ostr]( size_type i ) { ostr << i << ","; } );
            std::cout << "[stencilManagerPrint] ("
                      << entry.first.get<0>() << "[" << entry.first.get<0>().use_count() << "]"
                      << "," << entry.first.get<1>() << "[" << entry.first.get<1>().use_count() << "]"
                      << "," << int(entry.first.get<2>())
                      << ", (" << ostr.str()
                      << "),"
                      << "," << int(entry.first.get<4>())
                      << "): "
                      << " " << entry.second << "[" << entry.second.use_count() << "]"<< "\n";
    }
    std::cout << "********************************************************************************\n";
}

}
