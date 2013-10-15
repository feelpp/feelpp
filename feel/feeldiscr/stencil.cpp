/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-04-07
 */
#include <sstream>
#include <feel/feeldiscr/stencil.hpp>

namespace Feel
{

BlocksStencilPattern default_block_pattern( 1,1,size_type( Pattern::HAS_NO_BLOCK_PATTERN ) );

void
stencilManagerGarbageCollect()
{

    std::list<StencilManagerImpl::key_type> eltToDelete;

    BOOST_FOREACH( StencilManagerImpl::value_type & entry, StencilManager::instance() )
    {
        auto fspace1 = entry.first.get<0>().lock();
        auto fspace2 = entry.first.get<1>().lock();
        // each entry is a pair of tuple and graph
        if ( entry.second.unique() || entry.first.get<0>().expired() || entry.first.get<1>().expired() )
        {
#if !defined ( NDEBUG )
            std::ostringstream ostr;
            std::for_each( entry.first.get<3>().begin(),
                           entry.first.get<3>().end(),
                           [&ostr]( size_type i ) { ostr << i << ","; } );
            LOG(INFO) << "[stencilManagerGarbageCollect] deleting entry space: ( "
                      << fspace1
                      << "," << fspace2
                      << "," << int(entry.first.get<2>())
                      << ",( " << ostr.str()
                      << "),"
                      << "," << int(entry.first.get<4>())
                      << " ): "
                      << entry.second << "\n";
#endif
            eltToDelete.push_back(entry.first);
        }
    }
    VLOG(1) << "Deleting " << eltToDelete.size() << " stencils...";
    for ( auto it=eltToDelete.begin(),en=eltToDelete.end();it!=en;++it)
    {
        auto ki=StencilManager::instance().erase( *it );
    }
}

void
stencilManagerGarbage(StencilManagerImpl::key_type const& key)
{
    auto git = StencilManager::instance().find( key );
    if (  git != StencilManager::instance().end() )
    {
        if ( git->second.unique() )
            {
                StencilManager::instance().erase( git->first );
            }
    }
}

void
stencilManagerAdd(StencilManagerImpl::key_type const& key,StencilManagerImpl::graph_ptrtype graph)
{
    auto git = StencilManager::instance().find( key );
    if (  git == StencilManager::instance().end() )
        {
            StencilManager::instance().operator[]( key ) = graph;
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
        auto fspace1 = entry.first.get<0>().lock();
        auto fspace2 = entry.first.get<1>().lock();
        std::ostringstream ostr;
        std::for_each( entry.first.get<3>().begin(),
                       entry.first.get<3>().end(),
                       [&ostr]( size_type i ) { ostr << i << ","; } );
        LOG(INFO) << "[stencilManagerPrint] ("
                  << fspace1 << "[" << fspace2.use_count() << "]"
                  << "," << fspace2 << "[" << fspace2.use_count() << "]"
                  << "," << int(entry.first.get<2>())
                  << ", (" << ostr.str()
                  << "),"
                  << "," << int(entry.first.get<4>())
                  << "): "
                  << " " << entry.second << "[" << entry.second.use_count() << "]"<< "\n";
    }
    std::cout << "********************************************************************************\n";
}

stencilRangeMap0Type
stencilRangeMap( )
{
    return stencilRangeMap0Type();
}

}
