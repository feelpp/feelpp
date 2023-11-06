/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 21 Apr 2016

 Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_DETAIL_CONCATENATE_HPP
#define FEELPP_DETAIL_CONCATENATE_HPP 1

namespace Feel
{

#pragma GCC visibility push(hidden)

namespace detail
{


/**
 * implementation for faces
 */
template<typename RangeResult, typename RangeType>
void
concatenate_entities( RangeResult& elts, std::unordered_set<size_type> & idsIn, RangeType&& it )
{
    auto append = [&elts,&idsIn]( element_t<RangeType> const& e )
        {
            size_type eid = e.id();
            if ( idsIn.find( eid ) == idsIn.end() )
            {
                elts.push_back( e );
                idsIn.insert( eid );
            }
        };
    std::for_each( std::forward<RangeType>(it).begin(), std::forward<RangeType>(it).end(), append );
}

template<typename RangeResult, typename RangeType, typename ...Args>
void
concatenate_entities( RangeResult& elts, std::unordered_set<size_type> & idsIn, RangeType&& it, Args&&... args )
{
    concatenate_entities( elts, idsIn, std::forward<RangeType>(it) );
    concatenate_entities( elts, idsIn, std::forward<Args>(args)... );
}


template<typename RangeType, typename ...Args>
auto
concatenate_impl( RangeType&& it, Args&&... args )
{
    decay_type<RangeType> myelts{ std::forward<RangeType>(it).mesh() };
    std::unordered_set<size_type> idsIn;
    auto append = [&myelts,&idsIn]( element_t<RangeType> const& e ) { myelts.push_back( e ); idsIn.insert( e.id() ); };
    std::for_each( begin( it ), end( it ), append );

    concatenate_entities( myelts, idsIn, std::forward<Args>(args)... );
    return myelts;
}

} // detail
#pragma GCC visibility pop
}
#endif
