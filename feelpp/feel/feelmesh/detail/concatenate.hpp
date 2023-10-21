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
template<typename ContainerType, typename RangeType>
void
concatenate_entities( std::shared_ptr<ContainerType>& elts, std::unordered_set<size_type> & idsIn, RangeType it )
{
    using face_t = filter_entity_t<RangeType>;
    auto append = [&elts,&idsIn]( face_t const& e )
        {
            size_type eid = e.id();
            if ( idsIn.find( eid ) == idsIn.end() )
            {
                elts->push_back( boost::cref(e) );
                idsIn.insert( eid );
            }
        };
    std::for_each( begin( it ), end( it ), append );
}

template<typename ContainerType, typename RangeType, typename ...Args>
void
concatenate_entities( std::shared_ptr<ContainerType>& elts, std::unordered_set<size_type> & idsIn, RangeType it, Args... args )
{
    using face_t = filter_entity_t<RangeType>;
    auto append = [&elts,&idsIn]( face_t const& e )
        {
            size_type eid = e.id();
            if ( idsIn.find( eid ) == idsIn.end() )
            {
                elts->push_back( boost::cref(e) );
                idsIn.insert( eid );
            }
        };
    std::for_each( begin( it ), end( it ), append );
    concatenate_entities( elts, idsIn, args... );
}


template<typename RangeType, typename ...Args>
auto
concatenate_impl( RangeType&& it, Args... args )
{
    using face_t = filter_entity_t<RangeType>;
    typedef std::vector<boost::reference_wrapper<face_t const> > cont_range_type;
    std::shared_ptr<cont_range_type> myelts( new cont_range_type );
    std::unordered_set<size_type> idsIn;
    auto append = [&myelts,&idsIn]( face_t const& e ) { myelts->push_back( boost::cref(e) ); idsIn.insert( e.id() ); };
    std::for_each( begin( it ), end( it ), append );
    concatenate_entities( myelts, idsIn, args... );
    return range( _range = boost::make_tuple( filter_enum_t<RangeType>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts ), _mesh = it.mesh() );
}

} // detail
#pragma GCC visibility pop
}
#endif
