/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 24 Apr 2016

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
#ifndef FEELPP_COMPLEMENT_HPP
#define FEELPP_COMPLEMENT_HPP 1

#include <utility>
#include <feel/feelmesh/filters.hpp>


namespace Feel {

/**
 * Provides the complement of the set of entities given by \c
 * [begin(it),end(it)] with respect to the negation of predicate \p pred
 *
 * \code
 * // get all faces not marked by marker1 "R"
 * auto set = complement( faces(mesh), [&mesh]( auto const& e ) { return e.marker() == mesh->markerName( "R" ); });
 * \endcode
 */
template<typename IteratorType, typename Predicate>
auto
complement( IteratorType&& it, Predicate pred )
{
    using entity_t = filter_entity_t<IteratorType>;
    typedef std::vector<boost::reference_wrapper<entity_t const> > cont_range_type;
    std::shared_ptr<cont_range_type> myelts( new cont_range_type );

    auto append = [&myelts,&pred]( entity_t const& e ) { if ( !pred( e ) ) myelts->push_back( boost::cref(e) ); };
    std::for_each( begin( it ), end( it ), append );
    return range( _range=boost::make_tuple( filter_enum_t<IteratorType>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts ), _mesh=it.mesh() );
}

} // Feel
#endif
