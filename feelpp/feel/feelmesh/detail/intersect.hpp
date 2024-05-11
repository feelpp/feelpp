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
template<typename RangeType1, typename RangeType2>
auto intersect_entities(RangeType1&& elts, RangeType2&& it)
{
    // If RangeType has a non-trivial constructor, use std::decay_t to get the underlying type
    using DecayedRangeType = std::remove_reference_t<decay_type<RangeType1>>;
    using mesh_ptr_const_t = typename DecayedRangeType::mesh_ptr_const_t;
    //mesh_ptr_const_t p{std::forward<RangeType1>(elts).mesh()};
    //mesh_ptr_const_t p = elts.mesh();//std::const_pointer_cast<const mesh_ptr_non_const_t>(std::forward<RangeType1>(elts).mesh());

    DecayedRangeType myelts{elts.mesh()};

    std::unordered_set<size_type> idsIn;
    std::for_each(std::forward<RangeType1>(elts).begin(), std::forward<RangeType1>(elts).end(), 
                  [&idsIn](element_t<DecayedRangeType> const& elt) { idsIn.insert(elt.id()); });

    std::for_each(std::forward<RangeType2>(it).begin(), std::forward<RangeType2>(it).end(), 
                  [&myelts, &idsIn](element_t<DecayedRangeType> const& elt) {
                        auto itFindElt = idsIn.find(elt.id());
                        if (itFindElt != idsIn.end())
                            myelts.push_back(elt);
                 });
    
    return myelts;
}


template<typename RangeType1, typename RangeType2, typename ...Args>
auto
intersect_entities(RangeType1&& elts, RangeType2&& it, Args&&... args)
{
    return intersect_entities(intersect_entities(std::forward<RangeType1>(elts), std::forward<RangeType2>(it)), std::forward<Args>(args)...);
}


template<typename RangeType, typename ...Args>
auto
intersect_impl( RangeType&& it, Args&&... args )
{
    return intersect_entities( std::forward<RangeType>(it), std::forward<Args>(args)... );
}

} // detail
#pragma GCC visibility pop
}
#endif
