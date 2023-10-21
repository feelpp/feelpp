/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 25 Apr 2017

 Copyright (C) 2017 Feel++ Consortium

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
#ifndef FEELPP_DETAIL_INTERSECT_HPP
#define FEELPP_DETAIL_INTERSECT_HPP 1

namespace Feel
{

#pragma GCC visibility push(hidden)

namespace detail
{


/**
 * implementation for faces
 */
template<typename ContainerType, typename IteratorType>
std::shared_ptr<ContainerType>
intersect_entities( std::shared_ptr<ContainerType> const& elts, IteratorType it )
{
    using entity_t = filter_entity_t<IteratorType>;
    std::shared_ptr<ContainerType> myelts( new ContainerType );

    std::unordered_set<size_type> idsIn;
    for ( auto const& eltWrap1 : *elts )
        idsIn.insert( unwrap_ref( eltWrap1 ).id() );
    for ( auto const& eltWrap2 : it )
    {
        auto const& elt2 = unwrap_ref( eltWrap2 );
        auto itFindElt = idsIn.find( elt2.id() );
        if ( itFindElt != idsIn.end() )
            myelts->push_back( boost::cref(elt2) );
    }
    return myelts;
}

template<typename ContainerType, typename IteratorType, typename ...Args>
std::shared_ptr<ContainerType>
intersect_entities( std::shared_ptr<ContainerType> const& elts, IteratorType it, Args... args )
{
    return intersect_entities( intersect_entities( elts,it ), args... );
}

template<typename IteratorType, typename ...Args>
auto
intersect_impl( IteratorType it, Args... args )
{
    using entity_t = filter_entity_t<IteratorType>;
    typedef std::vector<boost::reference_wrapper<entity_t const> > cont_range_type;
    std::shared_ptr<cont_range_type> firstelts( new cont_range_type );

    auto append = [&firstelts]( entity_t const& e ) { firstelts->push_back( boost::cref(e) ); };
    std::for_each( begin( it ), end( it ), append );
    auto myelts = intersect_entities( firstelts, args... );
    return range( _range=boost::make_tuple( filter_enum_t<IteratorType>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts ),
                  _mesh=it.mesh() );
}

} // detail
#pragma GCC visibility pop
}
#endif
