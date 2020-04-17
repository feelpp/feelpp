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

namespace detail {

/* Zip an arbitrary -- runtime -- sequence of tuples to a tuple of runtime sequences.
 * The zip function creates a tuple of vectors by default.
 * The zip_with function uses the provided functional to create the runtime sequences in the returned tuple
 */

template <typename S>
struct zip_with_impl
{
    template <std::size_t N, typename F, typename ItT1, typename ItT2>
    static decltype(auto) transverse( F&& f, ItT1&& itBegin, ItT2&& itEnd ) {
        //auto const at = []( auto const t ) { 
            //return hana::at_c<N>( t );
        //};
        struct at {
            typedef std::decay_t<decltype( *itBegin )> it_value_type;
            typedef std::decay_t<decltype( hana::at_c<N>( *itBegin ) )> result_type;
            result_type operator()( it_value_type const& t ) const { return hana::at_c<N>( t ); }
        };
        return static_cast<F&&>(f)( 
                boost::iterators::transform_iterator( static_cast<ItT1&&>(itBegin), at{} ), 
                boost::iterators::transform_iterator( static_cast<ItT2&&>(itEnd), at{} ) 
                );
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


/**
 * Concat hana tuple by storing inside an iterable container (typically a vector) the data of same type
 * and avoid to increase the tuple size.
 */
template <typename FeelppTagOfTupleType,typename FeelppTagOfContainerType>
struct AdvancedConcatOfTupleContainerType
{
    template <typename Tag, typename T>
    struct is_a_t :
        hana::integral_constant<bool, std::is_same<Tag, typename T::feelpp_tag >::value >
    {};

    template<typename... SomeType>
    static constexpr auto apply( const SomeType&... v )
        {
            return applyImpl( hana::tuple<>{}, v... );
        }
private :

    template<typename ResType >
    static constexpr auto applyImpl( ResType && res )
        {
            return std::move( res );
        }

    template<int Index,typename ResType, typename... SomeType >
    static constexpr auto applyImplFromTupleType( ResType && res, hana::tuple<SomeType...> const& t )
        {
            constexpr int nTupleElt = std::decay_t<decltype(hana::size( t ))>::value;
            if constexpr ( Index < nTupleElt )
                return applyImplFromTupleType<Index+1>( applyImpl( std::forward<ResType>( res ), hana::at( t, hana::int_c<Index> ) ), t );
            else
                return std::move( res );
        }

    template < typename T1, typename... SomeType >
    static constexpr auto applyFromContainerType( T1 const& t1, hana::tuple<SomeType...> && res )
        {
            if constexpr ( hana::find( hana::to_tuple(hana::tuple_t<SomeType...> ), hana::type_c<T1>) == hana::nothing )
                         {
                             return hana::append( std::forward<hana::tuple<SomeType...>>( res ), t1 );
                         }
            else
            {
                hana::for_each( res, [&t1]( auto & e )
                                {
                                    if constexpr ( std::is_same_v<std::decay_t<decltype(e)>, T1> )
                                        {
                                            for ( auto const& se : t1 )
                                                e.push_back( se );
                                        }
                                });
                return std::move( res );
            }
        }
    template<typename ResType, typename T1, typename... SomeType2 >
    static constexpr auto applyImpl( ResType && res, T1 const& t1, const SomeType2&... tothers )
        {
            if constexpr ( is_a_t<FeelppTagOfContainerType, T1 >::value )
                {
                    return applyImpl( applyFromContainerType( t1,std::forward<ResType>(res) ), tothers... );
                }
            else if constexpr ( is_a_t<FeelppTagOfTupleType, T1 >::value )
                {
                    return applyImpl( applyImplFromTupleType<0>( std::forward<ResType>( res ), t1.tuple() ), tothers... );
                }
            else
                return std::move( res );
        }
};


} // namespace detail

template <typename F, typename ItT1, typename ItT2>
static auto
zip_with( F&& f, ItT1&& itBegin, ItT2&& itEnd ) {
    return Feel::detail::zip_with_impl<typename hana::tag_of<decltype(*itBegin)>::type>::apply(
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
