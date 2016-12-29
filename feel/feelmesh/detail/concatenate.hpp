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
template<typename ContainerType, typename IteratorType>
void
concatenate_entities( boost::shared_ptr<ContainerType>& elts, IteratorType it )
{
    using face_t = filter_entity_t<IteratorType>;
    auto append = [&elts]( face_t const& e ) { elts->push_back( boost::cref(e) ); };
    std::for_each( begin( it ), end( it ), append );
}

template<typename ContainerType, typename IteratorType, typename ...Args>
void
concatenate_entities( boost::shared_ptr<ContainerType>& elts, IteratorType it, Args... args )
{
    using face_t = filter_entity_t<IteratorType>;
    auto append = [&elts]( face_t const& e ) { elts->push_back( boost::cref(e) ); };
    std::for_each( begin( it ), end( it ), append );
    concatenate_entities( elts, args... );
}


template<typename IteratorType, typename ...Args>
ext_entities_from_iterator_t<IteratorType>
concatenate_impl( IteratorType it, Args... args )
{
    using face_t = filter_entity_t<IteratorType>;
    typedef std::vector<boost::reference_wrapper<face_t const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    auto append = [&myelts]( face_t const& e ) { myelts->push_back( boost::cref(e) ); };
    std::for_each( begin( it ), end( it ), append );
    concatenate_entities( myelts, args... );
    return boost::make_tuple( filter_enum_t<IteratorType>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts );
}

} // detail
#pragma GCC visibility pop
}
#endif
