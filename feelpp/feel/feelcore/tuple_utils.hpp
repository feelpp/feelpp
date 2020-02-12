/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Thibaut Metivet <thibaut.metivet@inria.fr>
       Date: 2020-01-23

  Copyright (C) 2020 Inria

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file tuple_utils.hpp
   \author Thibaut Metivet <thibaut.metivet@inria.fr>
   \date 2020-01-23
 */

#ifndef FEELPP_TUPLE_UTILS_HPP
#define FEELPP_TUPLE_UTILS_HPP 1

#include <feel/feelcore/traits.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/for_each.hpp>
#include <algorithm>
#include <boost/iterator/transform_iterator.hpp>

namespace Feel {

/* Dispatching for_each to hana or STL depending on the iterated sequence */
template<typename Xs, typename F,
    std::enable_if_t<boost::hana::Foldable<Xs>::value, int> = 0 >
void for_each( Xs && xs, F && f )
{
    boost::hana::for_each( xs, f );
}

template<typename Xs, typename F,
    std::enable_if_t<Feel::is_iterable<Xs>::value, int> = 0 >
void for_each( Xs && xs, F && f )
{
    std::for_each( begin(xs), end(xs), f );
}

/* Zip an arbitrary -- runtime -- sequence of tuples to a tuple of runtime sequences.
 * The zip function creates a tuple of vectors by default.
 * The zip_with function uses the provided functional to create the runtime sequences in the returned tuple
 */
namespace detail {
template <typename S>
struct zip_with_impl
{
    template <std::size_t N, typename F, typename ItT1, typename ItT2>
    static decltype(auto) transverse( F&& f, ItT1&& itBegin, ItT2&& itEnd ) {
        auto const at = []( auto const t ) { 
            return hana::at_c<N>( t );
        };
        return static_cast<F&&>(f)( iter::transform_iterator( static_cast<ItT1&&>(itBegin), at ), iter::transform_iterator( static_cast<ItT2&&>(itEnd), at ) );
    }

    template <std::size_t ...N, typename F, typename ItT1, typename ItT2>
    static auto
    zip_helper(std::index_sequence<N...>, F&& f, ItT1&& itBegin, ItT2&& itEnd ) {
        return hana::make<S>(transverse<N>( static_cast<F&&>(f), static_cast<ItT1&&>(itBegin), static_cast<ItT2&&>(itEnd) )...);
    }

    template <typename F, typename ItT1, typename ItT2>
    static auto
    apply( F&& f, ItT1&& itBegin, ItT2&& itEnd ) {
        constexpr std::size_t N = decltype(hana::length(*itBegin))::value;
        return zip_helper( std::make_index_sequence<N>{}, static_cast<F&&>(f), static_cast<ItT1&&>(itBegin), static_cast<ItT2&&>(itEnd) );
    }
};
} // namespace detail

template <typename F, typename ItT1, typename ItT2>
static auto
zip_with( F&& f, ItT1&& itBegin, ItT2&& itEnd ) {
    return detail::zip_with_impl<typename hana::tag_of<decltype(*itBegin)>::type>::apply(
            static_cast<F&&>(f), static_cast<ItT1&&>(itBegin), static_cast<ItT2&&>(itEnd)
            );
}

template <typename ItT1, typename ItT2>
static auto
zip( ItT1&& itBegin, ItT2&& itEnd ) {
    return zip_with( []( auto it1, auto it2 ) { return std::vector( it1, it2 ); }, itBegin, itEnd );
}

} // namespace Feel
#endif
