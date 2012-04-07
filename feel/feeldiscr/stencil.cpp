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
#include <feel/feeldiscr/stencil.hpp>

namespace Feel
{
namespace detail
{
void
runGarbageCollector()
{
    BOOST_FOREACH( StencilManagerImpl::value_type& entry, StencilManager::instance() )
    {
        // each entry is a pair of tuple and graph

        if ( entry.first.get<0>().unique() && entry.first.get<1>().unique() && entry.second.unique() )
        {
            std::cout << "[runGarbageCollector] deleting entry space:"
                      << entry.first.get<0>()
                      << "," << entry.first.get<1>()
                      << " and graph " << entry.second << "\n";

            StencilManager::instance().erase( entry.first );
        }
    }
}
} // detail
}
