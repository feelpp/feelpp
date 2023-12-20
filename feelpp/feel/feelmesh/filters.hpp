/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-27

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2010 Universite de Grenoble 1 (Joseph Fourier)
  Copyright (C) 2010-2016 Feel++ Consortium

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
   \file filters.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-27
 */
#ifndef __FEELPP_FILTERS_HPP
#define __FEELPP_FILTERS_HPP 1

#include <iterator>
#include <type_traits>
#include <utility>
#include <unordered_set>
#include <any>
#include <tuple>
#include <boost/mp11/function.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/rank.hpp>
#include <feel/feelcore/enums.hpp>
#include <feel/feelmesh/meshbase.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/iterator.hpp>
#include <feel/feelmesh/ranges.hpp>
#include <feel/feelmesh/detail/filters.hpp>
#include <feel/feelvf/expr.hpp>

namespace std
{
template< std::size_t I, class T >
struct tuple_element;
template <size_t I, typename Tuple>
using tuple_element_t = typename tuple_element<I,Tuple>::type;
}
namespace Feel
{
#if 0   
template <int I, typename Tuple>
struct tuple_element;
template <int I, typename Head, typename... Tail>
struct tuple_element<I, boost::tuple<Head, Tail...>>
{
    typedef typename tuple_element<I - 1, boost::tuple<Tail...>>::type type;
};
template <typename Head, typename... Tail>
struct tuple_element<0, boost::tuple<Head, Tail...>>
{
    typedef Head type;
};
template <int I, typename Tuple>
using tuple_element_t = typename std::tuple_element<I,Tuple>::type;

template <typename> struct is_tuple: std::false_type {};

template <typename ...T> struct is_tuple<boost::tuple<T...>>: std::true_type {};

template <typename T> 
inline constexpr bool is_tuple_v = is_tuple<T>::value;
#else // 0
template <typename> struct is_tuple: std::false_type {};
template <typename ...T> struct is_tuple<std::tuple<T...>>: std::true_type {};
template <typename T>
inline constexpr bool is_tuple_v = is_tuple<T>::value;  
#endif // 0

template<typename MeshType>
decltype(auto)
make_elements_wrapper()
{
    using elt_wrapper_t = typename MeshTraits<MeshType>::elements_reference_wrapper_type;
    return std::make_shared<elt_wrapper_t>();
}

template <typename RangeT>
using range_elementtype_t = std::tuple_element_t<0, RangeT>;

template <typename RangeT>
constexpr size_type range_elementtype_v = range_elementtype_t<RangeT>::value;

template <typename RangeT>//, std::enable_if_t < std::is_same_v < std::tuple_element<1, RangeT>,std::tuple_element<2, RangeT>>,int> = 0 >
using range_iterators_t = std::tuple_element_t<1, RangeT>;

template <typename RangeT>
using range_container_t = std::tuple_element_t<3, RangeT>;

//!
//! @return the worldcomm associated to the mesh stored in the range
//! @ingroup Mesh
//!
template <typename RangeT, std::enable_if_t<is_range_v<RangeT>, int> = 0>
WorldComm const& worldComm( RangeT const& range )
{
    return boost::unwrap_ref( *begin(range) ).mesh()->worldComm();
}

/**
 * @brief get the rank of a range of entities
 * 
 * @tparam S 
 * @tparam ITERATOR 
 * @tparam CONTAINER 
 * @param range 
 * @return int the rank of the entity range
 */
template <typename RangeT,std::enable_if_t<is_range_v<RangeT>,int> = 0>
rank_type rank( RangeT const& range )
{
    return worldComm( range ).globalRank();
}

#if 0
/**
 * namespace for meta mesh computation data structure
 */
namespace meta
{

template<typename MeshType>
struct elements
{
    static inline const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<MeshType>::element_const_iterator,
                         typename MeshTraits<MeshType>::element_const_iterator> type;
};

template<typename MeshType>
struct markedelements
{
    static inline const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<MeshType>::marker_element_const_iterator,
                         typename MeshTraits<MeshType>::marker_element_const_iterator> type;
};

template<typename MeshType>
struct marked2elements
{
    static inline const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<MeshType>::marker2_element_const_iterator,
                         typename MeshTraits<MeshType>::marker2_element_const_iterator> type;
};

template<typename MeshType>
struct marked3elements
{
    static inline const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                         typename MeshTraits<MeshType>::marker3_element_const_iterator,
                         typename MeshTraits<MeshType>::marker3_element_const_iterator> type;
};

} // meta
#endif

/**
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
allelements( MeshType const& mesh   )
{
    return range( _range=Feel::detail::allelements( mesh ), _mesh = mesh );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,decay_type<std::remove_pointer_t<MeshType>>>,int> = 0>
auto
elements( MeshType const& mesh )
{
    return range( _range=Feel::detail::elements( mesh, rank( mesh ) ), _mesh=mesh );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,decay_type<std::remove_pointer_t<MeshType>>>,int> = 0>
auto
elements( MeshType const& mesh, int pid )
{
    return range( _range=Feel::detail::elements( mesh, pid ), _mesh=mesh );
}
/**
 * @brief EXTRACT_ELEM
 */
enum class select_elements_from_expression {
    with_positive_values = 1,
    with_negative_values = 2,
    with_changing_sign = 3,
    with_value = 4,
    with_value_range = 5
};

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with pid \p flag and
 *  verified expr > 0 for all nodes of the element
 */
template <typename MeshType, typename ExprType, typename... Ts, std::enable_if_t<std::is_base_of_v<MeshBase<>, decay_type<std::remove_pointer_t<MeshType>>>, int> = 0>
Range<MeshType,MESH_ELEMENTS>
elements( MeshType const& mesh, vf::Expr<ExprType> const& expr, Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    select_elements_from_expression selector = args.get_else(_selector, select_elements_from_expression::with_positive_values);
    auto && threshold = args.get_else(_threshold, 1e-9);
    auto && value = args.get_else(_value, 0 );
    auto && min = args.get_else(_min, 0 );
    auto && max = args.get_else(_max, 1 );
    using expr_t = decay_type<decltype( expr )>;
    rank_type pid = rank( mesh );
    Range<MeshType,MESH_ELEMENTS> myelts(mesh);
    auto const& imesh = Feel::unwrap_ptr( mesh );
    typedef typename MeshTraits<MeshType>::mesh_type mesh_type;
    auto it = imesh.beginOrderedElement();
    auto en = imesh.endOrderedElement();
    if ( it != en )
    {
        auto gm = imesh.gm();
        auto const& initElt = unwrap_ref( *it );
        typename mesh_type::reference_convex_type refConvex;
        auto geopc = gm->preCompute( refConvex.points() );
        const size_type context = expr_t::context|vm::POINT;
        auto ctx = gm->template context<context>( initElt, geopc );
        auto expr_evaluator = expr.evaluator( vf::mapgmc(ctx) );
        for ( ; it!=en;++it )
        {
            auto const& elt = unwrap_ref( *it );
            if ( elt.processId() != pid )
                continue;
            ctx->template update<context>( elt );
            expr_evaluator.update( vf::mapgmc( ctx ) );
            int elt_count_positive = 0, elt_count_negative = 0;
            for ( uint16_type q=0;q<ctx->nPoints();++q )
            {
                double val = expr_evaluator.evalq( 0,0,q );
                if ( selector == select_elements_from_expression::with_positive_values && val > -threshold )
                {
                    elt_count_positive++;
                    break;
                }
                else if ( selector == select_elements_from_expression::with_negative_values && val < threshold )
                {
                    elt_count_negative++;
                    break;
                }
                else if ( selector == select_elements_from_expression::with_changing_sign ) 
                {
                    if ( val < 0 )
                        elt_count_negative++;
                    else if ( val >= 0 )
                        elt_count_positive++;
                }
                else if ( selector == select_elements_from_expression::with_value ) 
                {
                    if ( std::abs(val - value) < threshold )
                    {
                        elt_count_positive++;
                        break;
                    }
                }
                else if ( selector == select_elements_from_expression::with_value_range ) 
                {
                    if ( ( val > min-threshold ) && ( val < max+threshold ) )
                        elt_count_positive++;
                }
            }
            if ( ( selector == select_elements_from_expression::with_positive_values ) && elt_count_positive > 0 )
            {
                myelts.push_back(elt);
            }
            else if ( ( selector == select_elements_from_expression::with_negative_values ) && elt_count_negative > 0 )
            {
                myelts.push_back(elt);
            }
            else if ( ( selector == select_elements_from_expression::with_changing_sign ) && ( elt_count_positive > 0 && elt_count_negative > 0 ) )
            {
                myelts.push_back(elt);
            }
            else if ( ( selector == select_elements_from_expression::with_value ) && elt_count_positive > 0 )
            {
                myelts.push_back(elt);
            }
            else if ( ( selector == select_elements_from_expression::with_value_range ) && elt_count_positive > 0 )
            {
                myelts.push_back(elt);
            }
        }
    }
    myelts.shrink_to_fit();

    return myelts;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the mesh
 * which share a face with the boundary
 */
template<typename MeshType>
auto
boundaryelements( MeshType const& mesh, uint16_type entity_min_dim = 0, uint16_type entity_max_dim = 2 )
{
    return range(_range=Feel::detail::boundaryelements( mesh, entity_min_dim, entity_max_dim, rank( mesh ) ), _mesh=mesh );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the mesh
 * which are strictly within the domain that is to say they do not
 * share a face with the boundary
 */
template<typename MeshType>
auto
internalelements( MeshType const& mesh )
{
    return range(_range=Feel::detail::internalelements( mesh, rank( mesh ) ), _mesh=mesh );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
auto
markedelementsByType( MeshType const& mesh, uint16_type markerType )
{
    return range(_range=Feel::detail::markedelements( mesh, markerType, rank( mesh ) ), _mesh=mesh );
}
/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
auto
markedelementsByType( MeshType const& mesh, uint16_type markerType,
                      boost::any const& markersFlag )
{
    std::set<flag_type> markerFlagSet = Feel::unwrap_ptr( mesh ).markersId( markersFlag );
    return range( _range=Feel::detail::markedelements( mesh, markerType, markerFlagSet, rank( mesh ) ), _mesh=mesh );
}
/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
auto
markedelementsByType( MeshType const& mesh, uint16_type markerType,
                      std::initializer_list<boost::any> const& markersFlag )
{
    std::set<flag_type> markerFlagSet;
    for ( auto const& it : markersFlag )
    {
        flag_type theflag = Feel::unwrap_ptr( mesh ).markerId( it );
        VLOG(2) << "[markedfacesByType] flag: " << theflag << "\n";
        markerFlagSet.insert( theflag );
    }
    return range( _range=Feel::detail::markedelements( mesh, markerType, markerFlagSet, rank( mesh ) ), _mesh=mesh );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markedelements( MeshType const& mesh )
{
    return markedelementsByType( mesh, 1 );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markedelements( MeshType const& mesh, boost::any const& markersFlag )
{
    return markedelementsByType( mesh, 1, markersFlag );
}
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markedelements( MeshType const& mesh, std::initializer_list<boost::any> const& markersFlag )
{
    return markedelementsByType( mesh, 1, markersFlag );
}
/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker2 string
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
marked2elements( MeshType const& mesh, boost::any const& markersFlag )
{
    return markedelementsByType( mesh, 2, markersFlag );
}
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
marked2elements( MeshType const& mesh, std::initializer_list<boost::any> const& markersFlag )
{
    return markedelementsByType( mesh, 2, markersFlag );
}
/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker3 string
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
marked3elements( MeshType const& mesh, boost::any const& markersFlag )
{
    return markedelementsByType( mesh, 3, markersFlag );
}
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
marked3elements( MeshType const& mesh, std::initializer_list<boost::any> const& markersFlag )
{
    return markedelementsByType( mesh, 3, markersFlag );
}



/**
 *
 * \ingroup MeshIterators
 * \return a mapping between a id and pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
std::map<int, Range<MeshType, MESH_ELEMENTS>>
collectionOfMarkedelements( MeshType const& mesh, std::any const& collectionOfMarkersFlag )
{
    return Feel::detail::collectionOfMarkedelements<0>( mesh,collectionOfMarkersFlag );
}
/**
 *
 * \ingroup MeshIterators
 * \return a mapping between a id and pair of (range active element, range ghost elements)
 */
template<typename MeshType>
auto
collectionOfMarkedelementsWithGhostsSplitted( MeshType const& mesh, std::any const& collectionOfMarkersFlag )
{
    return Feel::detail::collectionOfMarkedelements<2>( mesh,collectionOfMarkersFlag );
}



/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with id \p flag
 */
template<typename MeshType>
auto
idedelements( MeshType const& mesh, flag_type flag )
{
    return range( _range=Feel::detail::idedelements( mesh, flag ), _mesh=mesh );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over faces of the
 * mesh on processor \c __pid
 *
 * @param mesh a mesh data structure
 * @param __pid process id
 *
 * @return a pair of face iterators (begin,end)
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
faces( MeshType const& mesh )
{
    return range(_range=Feel::detail::faces( mesh, rank( mesh ) ), _mesh=mesh, _pid=rank(mesh) );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with id \p id
 */
template<typename MeshType>
auto
idedfaces( MeshType const& mesh, size_type id )
{
    return range(_range=Feel::detail::idedfaces( mesh, id ), _mesh=mesh );
}



template<typename MeshType>
auto
markedfacesByType( MeshType const& mesh, uint16_type markerType )
{
    return range( _range=Feel::detail::markedfaces( mesh, markerType, rank( mesh ) ), _mesh=mesh );
}

template<typename MeshType>
auto
markedfacesByType( MeshType const& mesh, uint16_type markerType,
                   boost::any const& __marker )
{
    std::set<flag_type> markerFlagSet = Feel::unwrap_ptr( mesh ).markersId( __marker );
    return range( _range=Feel::detail::markedfaces( mesh, markerType, markerFlagSet, rank( mesh ) ), _mesh=mesh );
}

template<typename MeshType>
auto
markedfacesByType( MeshType const& mesh, uint16_type markerType,
                   std::initializer_list<boost::any> const& markersFlag )
{
    std::set<flag_type> markerFlagSet;
    for ( auto const& it : markersFlag )
    {
        flag_type theflag = Feel::unwrap_ptr( mesh ).markerId( it );
        VLOG(2) << "[markedfacesByType] flag: " << theflag << "\n";
        markerFlagSet.insert( theflag );
    }
    return range( _range=Feel::detail::markedfaces( mesh, markerType, markerFlagSet, rank( mesh ) ), _mesh=mesh );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over faces of the
 * mesh marked
 *
 * @param mesh a mesh data structure
 * @param __pid process id
 *
 * @return a pair of iterators (begin,end) for the set of marked faces
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markedfaces( MeshType const& mesh )
{
    return markedfacesByType( mesh, 1 );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over faces of the
 * mesh marked with \c __marker
 *
 * @param mesh a mesh data structure
 * @param markersFlag a marker that identifies faces
 *
 * @return a pair of iterators (begin,end) for the set of marked faces
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markedfaces( MeshType const& mesh, boost::any const& markersFlag )
{
    return markedfacesByType( mesh, 1, markersFlag );
}
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markedfaces( MeshType const& mesh, std::initializer_list<boost::any> const& markersFlag )
{
    return markedfacesByType( mesh, 1, markersFlag );
}

template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
marked2faces( MeshType const& mesh, boost::any const& markersFlag )
{
    return markedfacesByType( mesh, 2, markersFlag );
}
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
marked2faces( MeshType const& mesh,
              std::initializer_list<boost::any> const& markersFlag )
{
    return markedfacesByType( mesh, 2, markersFlag );
}

template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
marked3faces( MeshType const& mesh, boost::any const& markersFlag )
{
    return markedfacesByType( mesh, 3, markersFlag );
}
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
marked3faces( MeshType const& mesh,
              std::initializer_list<boost::any> markersFlag )
{
    return markedfacesByType( mesh, 3, markersFlag );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all boundary faces of the
 * mesh
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
boundaryfaces( MeshType const& mesh  )
{
    return range( _range=Feel::detail::boundaryfaces( mesh, rank( mesh ) ), _mesh=mesh, _pid=rank(mesh) );
}
template <typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>, unwrap_ptr_t<MeshType>>, int> = 0>
Range<MeshType,MESH_FACES>
boundaryfaces( MeshType const& mesh, Range<MeshType,MESH_ELEMENTS> const& r )
{
    Range<MeshType,MESH_FACES> ret(mesh);

    for ( auto const& e : r )
    {
        auto const& elt = boost::unwrap_ref( e );
        auto const& eltfaces = elt.faces();
        for ( auto it = eltfaces.first, en = eltfaces.second; it != en; ++it )
        {
            auto const& face = *it;
            if ( face->isConnectedTo0() && !face->isConnectedTo1() )
            {
                ret.push_back( *face );
            }
        }
    }
    ret.shrink_to_fit();
    return ret;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all internal faces of the
 * mesh belong to process domain \p __pid
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
internalfaces( MeshType const& mesh )
{
    return range( _range=Feel::detail::internalfaces( mesh, rank( mesh ) ), _mesh=mesh, _pid=rank(mesh) );
}

template <typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>, unwrap_ptr_t<MeshType>>, int> = 0>
Range<MeshType,MESH_FACES>
internalfaces( MeshType const& mesh, Range<MeshType,MESH_ELEMENTS> const& r )
{
    typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::faces_reference_wrapper_type );
    // store face ids to know if a face is already in the list
    std::set<int> fids;
    for(auto const&e : r )
    {
        auto const& elt = boost::unwrap_ref( e );
        auto const& eltfaces = elt.faces();
        for(auto it = eltfaces.first, en = eltfaces.second; it != en; ++it)
        {
            auto const& face = *it;
            if ( face->isConnectedTo0() && face->isConnectedTo1()  )
            {
                if ( fids.find( face->id() ) == fids.end() )
                {
                    fids.insert( face->id() );
                    myelts->push_back( boost::cref( *face ) );
                }
            }
        }
    }
    myelts->shrink_to_fit();
    return range( _range=boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts ), _mesh=mesh );
}
/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all interprocess faces of the
 * mesh belonging to process \p __pid
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
interprocessfaces( MeshType const& mesh )
{
    return range( _range=Feel::detail::interprocessfaces( mesh, invalid_rank_type_value ),  _mesh=mesh );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all interprocess faces of the
 * mesh belonging to process \p __pid
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
interprocessfaces( MeshType const& mesh, rank_type neighbor_pid )
{
    return range( _range=Feel::detail::interprocessfaces( mesh, neighbor_pid ), _mesh=mesh );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over edges of the
 * mesh on processor \c __pid
 *
 * @param mesh a mesh data structure
 * @param __pid process id
 *
 * @return a pair of edge iterators (begin,end)
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
edges( MeshType const& mesh )
{
    return range( _range = Feel::detail::edges( mesh, rank( mesh ) ), _mesh=mesh, _pid=rank(mesh) );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over edges of the
 * mesh marked with \c __marker
 *
 * @param mesh a mesh data structure
 * @param __marker a marker that identifies edges
 * @param __pid process id
 *
 * @return a pair of iterators (begin,end) for the set of marked edges
 */
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markededgesByType( MeshType const& mesh, uint16_type markerType,
                   typename std::enable_if<is_3d<MeshType>::value>::type* = nullptr )
{
    return range( _range = Feel::detail::markededges( mesh, markerType, rank( mesh ) ), _mesh=mesh );
}
template<typename MeshType>
auto
markededgesByType( MeshType const& mesh, uint16_type markerType,
                   typename std::enable_if<is_2d<MeshType>::value>::type* = nullptr )
{
    return markedfacesByType( mesh,markerType );
}
template<typename MeshType>
auto
markededgesByType( MeshType const& mesh, uint16_type markerType,
                   boost::any const& __marker,
                   typename std::enable_if<is_3d<MeshType>::value>::type* = nullptr )
{
    std::set<flag_type> markerFlagSet = Feel::unwrap_ptr( mesh ).markersId( __marker );
    return range( _range = Feel::detail::markededges( mesh, markerType, markerFlagSet, rank( mesh ) ), _mesh=mesh );
}
template<typename MeshType>
auto
markededgesByType( MeshType const& mesh, uint16_type markerType,
                   boost::any const& markersFlag,
                   typename std::enable_if<is_2d<MeshType>::value>::type* = nullptr )
{
    return markedfacesByType( mesh,markerType,markersFlag );
}
template<typename MeshType>
auto
markededgesByType( MeshType const& mesh, uint16_type markerType,
                   std::initializer_list<boost::any> const& markersFlag,
                   typename std::enable_if<is_3d<MeshType>::value>::type* = nullptr )
{
    std::set<flag_type> markerFlagSet;
    for ( auto const& it : markersFlag )
    {
        flag_type theflag = Feel::unwrap_ptr( mesh ).markerId( it );
        VLOG(2) << "[markededgesByType] flag: " << theflag << "\n";
        markerFlagSet.insert( theflag );
    }
    return range( _range = Feel::detail::markededges( mesh, markerType, markerFlagSet, rank( mesh ) ), _mesh=mesh );
}
template<typename MeshType>
auto
markededgesByType( MeshType const& mesh, uint16_type markerType,
                   std::initializer_list<boost::any> const& markersFlag,
                   typename std::enable_if<is_2d<MeshType>::value>::type* = nullptr )
{
    return markedfacesByType( mesh,markerType,markersFlag );
}

template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markededges( MeshType const& mesh )
{
    return markededgesByType( mesh, 1 );
}

template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markededges( MeshType const& mesh, boost::any const& markersFlag )
{
    return markededgesByType( mesh, 1, markersFlag );
}

template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0>
auto
markededges( MeshType const& mesh, std::initializer_list<boost::any> const& markersFlag )
{
    return markededgesByType( mesh, 1, markersFlag );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all boundary edges of the
 * mesh
 */
template<typename MeshType>
auto
boundaryedges( MeshType const& mesh )
{
    return range( _range = Feel::detail::boundaryedges( mesh ), _mesh = mesh );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all internal edges of the
 * mesh
 */
template<typename MeshType>
auto
internaledges( MeshType const& mesh )
{
    return range( _range=Feel::detail::internaledges( mesh ), _mesh=mesh );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
auto
points( MeshType const& mesh )
{
    return range( _range = Feel::detail::points( mesh ), _mesh=mesh );
}


template<typename MeshType>
auto
markedpointsByType( MeshType const& mesh, uint16_type markerType )
{
    return range( _range = Feel::detail::markedpoints( mesh, markerType, rank( mesh ) ), _mesh = mesh );
}

template<typename MeshType>
auto
markedpointsByType( MeshType const& mesh, uint16_type markerType,
                   boost::any const& __marker )
{
    std::set<flag_type> markerFlagSet = Feel::unwrap_ptr( mesh ).markersId( __marker );
    return range( _range=Feel::detail::markedpoints( mesh, markerType, markerFlagSet, rank( mesh ) ), _mesh=mesh);
}

template<typename MeshType>
auto
markedpointsByType( MeshType const& mesh, uint16_type markerType,
                   std::initializer_list<boost::any> const& markersFlag )
{
    std::set<flag_type> markerFlagSet;
    for ( auto const& it : markersFlag )
    {
        flag_type theflag = Feel::unwrap_ptr( mesh ).markerId( it );
        VLOG(2) << "[markedpointsByType] flag: " << theflag << "\n";
        markerFlagSet.insert( theflag );
    }
    return range( _range=Feel::detail::markedpoints( mesh, markerType, markerFlagSet, rank( mesh ) ), _mesh=mesh );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over points of the
 * mesh marked
 *
 * @param mesh a mesh data structure
 * @param __pid process id
 *
 * @return a pair of iterators (begin,end) for the set of marked points
 */
template<typename MeshType>
auto
markedpoints( MeshType const& mesh )
{
    return markedpointsByType( mesh, 1 );
}

/**
 *
 * \ingroup MeshIterators
 *
 * @param mesh a mesh data structure
 * @param markersFlag a marker that identifies points
 * \return the range of iterators [begin,end,container] over the marked points with \p flag of the mesh
 *
 */
template<typename MeshType>
auto
markedpoints( MeshType const& mesh, boost::any const& markersFlag )
{
    return markedpointsByType( mesh, 1, markersFlag );
}
template<typename MeshType>
auto
markedpoints( MeshType const& mesh, std::initializer_list<boost::any> const& markersFlag )
{
    return markedpointsByType( mesh, 1, markersFlag );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the boundary points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
auto
boundarypoints( MeshType const& mesh )
{
    return range( _range=Feel::detail::boundarypoints( mesh ), _mesh=mesh );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the internal(not on the boundary) points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
auto
internalpoints( MeshType const& mesh )
{
    return range( _range=Feel::detail::internalpoints( mesh ), _mesh=mesh );
}



template<typename ElementType>
boost::tuple<mpl::size_t<MESH_ELEMENTS>,
             typename std::list<ElementType>::const_iterator,
             typename std::list<ElementType>::const_iterator>
element( ElementType const& elt  )
{
    std::list<ElementType> lst;
    lst.push_back(  elt );
    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              lst.begin(),
                              lst.end() );
}


template<typename MeshType>
auto
elements( MeshType const& imesh, EntityProcessType entity )
{
    typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::elements_reference_wrapper_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );
    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : elements(mesh.shared_from_this()) )
        {
            myelts->push_back( boost::cref( boost::unwrap_ref( elt ) ) );
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto rangeInterProcessFaces = mesh.interProcessFaces();
        auto face_it = std::get<0>( rangeInterProcessFaces );
        auto const face_en = std::get<1>( rangeInterProcessFaces );
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& faceip = boost::unwrap_ref( *face_it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;

            // add elt in range
            myelts->push_back(boost::cref(eltOffProc));

            eltGhostDone.insert( eltOffProc.id() );
        }
    }
    myelts->shrink_to_fit();
    return range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                            myelts->begin(),
                                            myelts->end(),
                                            myelts ), _mesh=mesh.shared_from_this() );
}

template<typename MeshType>
auto
boundaryelements( MeshType const& imesh, EntityProcessType entity )
{
    typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::elements_reference_wrapper_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : boundaryelements(mesh.shared_from_this()) )
        {
            myelts->push_back( boost::cref( boost::unwrap_ref( elt ) ) );
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto rangeInterProcessFaces = mesh.interProcessFaces();
        auto face_it = std::get<0>( rangeInterProcessFaces );
        auto const face_en = std::get<1>( rangeInterProcessFaces );
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& faceip = boost::unwrap_ref( *face_it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;

            // add elt in range
            if ( eltOffProc.isOnBoundary() )
            {
                myelts->push_back(boost::cref(eltOffProc));
            }

            eltGhostDone.insert( eltOffProc.id() );
        }
    }
    myelts->shrink_to_fit();
    return range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                            myelts->begin(),
                                            myelts->end(),
                                            myelts ), _mesh=mesh.shared_from_this() );
}

template<typename MeshType>
auto
elementsWithMarkedFaces( MeshType const& imesh, boost::any const& flag, EntityProcessType entity = EntityProcessType::LOCAL_ONLY )
{
    typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::elements_reference_wrapper_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );

    flag_type theflag = mesh.markerId( flag );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        for ( auto const& faceWrap : markedfaces(mesh,flag) )
        {
            auto const& face = boost::unwrap_ref( faceWrap );
            if ( face.isConnectedTo0() )
                myelts->push_back( boost::cref( face.element0() ) );
            if ( face.isConnectedTo1() )
                myelts->push_back( boost::cref( face.element1() ) );
        }
    }
    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto rangeInterProcessFaces = mesh.interProcessFaces();
        auto face_it = std::get<0>( rangeInterProcessFaces );
        auto const face_en = std::get<1>( rangeInterProcessFaces );
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& faceip = boost::unwrap_ref( *face_it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;


            // add elt in range
            if ( eltOffProc.hasMarker() && eltOffProc.marker().value() == theflag )
            {
                myelts->push_back(boost::cref(eltOffProc));
            }

            eltGhostDone.insert( eltOffProc.id() );
        }
    }
    myelts->shrink_to_fit();
    return range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                            myelts->begin(),
                                            myelts->end(),
                                            myelts ), _mesh=mesh.shared_from_this() );
}

template<typename MeshType>
auto
markedelements( MeshType const& imesh, boost::any const& flag, EntityProcessType entity )
{
    typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::elements_reference_wrapper_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );

    flag_type theflag = mesh.markerId( flag );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : markedelements(mesh,flag) )
        {
            // myelts->push_back(boost::cref(elt));
            myelts->push_back( boost::cref( boost::unwrap_ref( elt ) ) );
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto rangeInterProcessFaces = mesh.interProcessFaces();
        auto face_it = std::get<0>( rangeInterProcessFaces );
        auto const face_en = std::get<1>( rangeInterProcessFaces );
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& faceip = boost::unwrap_ref( *face_it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;


            // add elt in range
            if ( eltOffProc.hasMarker() && eltOffProc.marker().value() == theflag )
                myelts->push_back(boost::cref(eltOffProc));

            eltGhostDone.insert( eltOffProc.id() );
        }
    }
    myelts->shrink_to_fit();
    return range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                            myelts->begin(),
                                            myelts->end(),
                                            myelts ), _mesh=mesh.shared_from_this() );
}

template<typename MeshType>
auto
marked2elements( MeshType const& imesh, boost::any const& flag, EntityProcessType entity )
{
    typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::elements_reference_wrapper_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );

    flag_type theflag = mesh.markerId( flag );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : marked2elements(mesh,flag) )
        {
            myelts->push_back( boost::cref( boost::unwrap_ref( elt ) ) );
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto rangeInterProcessFaces = mesh.interProcessFaces();
        auto face_it = std::get<0>( rangeInterProcessFaces );
        auto const face_en = std::get<1>( rangeInterProcessFaces );
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& faceip = boost::unwrap_ref( *face_it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;


            // add elt in range
            if ( eltOffProc.hasMarker2() && eltOffProc.marker2().value() == theflag )
                myelts->push_back(boost::cref(eltOffProc));

            eltGhostDone.insert( eltOffProc.id() );
        }
    }
    myelts->shrink_to_fit();
    return range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                            myelts->begin(),
                                            myelts->end(),
                                            myelts ), _mesh=mesh.shared_from_this() );
}

template<typename MeshType>
auto
faces( MeshType const& imesh, EntityProcessType entity )
{
    typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::faces_reference_wrapper_type );
    //typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> > cont_range_type;
    //std::shared_ptr<cont_range_type> myelts( new cont_range_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );
    
    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& theface : faces(mesh) )
        {
            myelts->push_back( boost::cref( boost::unwrap_ref( theface ) ) );
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        auto rangeInterProcessFaces = mesh.interProcessFaces();
        auto face_it = std::get<0>( rangeInterProcessFaces );
        auto const face_en = std::get<1>( rangeInterProcessFaces );
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& faceip = boost::unwrap_ref( *face_it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            for ( size_type f = 0; f < mesh.numLocalFaces(); f++ )
            {
                auto const& theface = eltOffProc.face(f);
                myelts->push_back(boost::cref(theface));
            }
        }
    }
    myelts->shrink_to_fit();
    return range( _range=boost::make_tuple( mpl::size_t<MESH_FACES>(),
                                            myelts->begin(),
                                            myelts->end(),
                                            myelts ), _mesh=mesh.shared_from_this() );

}

template<typename MeshType>
//std::enable_if_t<std::is_base_of_v<MeshBase,unwrap_ptr_t<MeshType>>,ext_faces_t<MeshType>>
auto
boundaryfaces( MeshType const& imesh, EntityProcessType entity )
{
    typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::faces_reference_wrapper_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& theface : boundaryfaces(mesh) )
        {
            myelts->push_back( boost::cref( boost::unwrap_ref( theface ) ) );
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        //std::set<size_type> faceGhostDone;
        auto rangeInterProcessFaces = mesh.interProcessFaces();
        auto face_it = std::get<0>( rangeInterProcessFaces );
        auto const face_en = std::get<1>( rangeInterProcessFaces );
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& faceip = boost::unwrap_ref( *face_it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            for ( size_type f = 0; f < mesh.numLocalFaces(); f++ )
            {
                auto const& theface = eltOffProc.face(f);
                if ( theface.isOnBoundary() ) //&& faceGhostDone.find( theface.id() ) == faceGhostDone.end() )
                    myelts->push_back(boost::cref(theface));
            }
        }
    }
    myelts->shrink_to_fit();
    return range( _range=boost::make_tuple( mpl::size_t<MESH_FACES>(),
                                            myelts->begin(),
                                            myelts->end(),
                                            myelts ), _mesh=mesh.shared_from_this() );

}


template<typename MeshType>
auto
marked2faces( MeshType const& imesh, boost::any flag, EntityProcessType entity )
{
    typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::faces_reference_wrapper_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );

    flag_type theflag = mesh.markerId( flag );
    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& theface : marked2faces(mesh, theflag) )
        {
            myelts->push_back( boost::cref( boost::unwrap_ref( theface ) ) );
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        //std::set<size_type> faceGhostDone;
        auto rangeInterProcessFaces = mesh.interProcessFaces();
        auto face_it = std::get<0>( rangeInterProcessFaces );
        auto const face_en = std::get<1>( rangeInterProcessFaces );
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& faceip = boost::unwrap_ref( *face_it );
            auto const& elt0 = faceip.element0();
            auto const& elt1 = faceip.element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            for ( size_type f = 0; f < mesh.numLocalFaces(); f++ )
            {
                auto const& theface = eltOffProc.face(f);
                if ( theface.hasMarker2() && theface.marker2().value() == theflag ) //&& faceGhostDone.find( theface.id() ) == faceGhostDone.end() )
                    myelts->push_back(boost::cref(theface));
            }
        }
    }
    myelts->shrink_to_fit();
    return range( _range = boost::make_tuple( mpl::size_t<MESH_FACES>(),
                                              myelts->begin(),
                                              myelts->end(),
                                              myelts ), 
                  _mesh=mesh.shared_from_this(), _marker1=theflag );

}



template<typename MeshType>
auto
interprocessedges( MeshType const& imesh, rank_type neighbor_pid = invalid_rank_type_value, EntityProcessType entity = EntityProcessType::ALL )
{
    typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype myedges( new typename MeshTraits<MeshType>::edges_reference_wrapper_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );

    if ( entity == EntityProcessType::ALL )
    {
        if ( neighbor_pid == invalid_rank_type_value )
        {
            for ( auto const& theedge : edges(mesh) )
            {
                if ( theedge.numberOfProcGhost() > 0 )
                    myedges->push_back(boost::cref(theedge));
            }
        }
        else
        {
            for ( auto const& theedge : edges(mesh) )
            {
                if ( theedge.numberOfProcGhost() > 0 && theedge.isSubEntityOfGhostElement( neighbor_pid ) )
                    myedges->push_back(boost::cref(theedge));
            }
        }
    }
    else if (  entity == EntityProcessType::IGNORE_ENTITY_ON_INTERPROCESS_FACE )
    {
        std::set<size_type> saveIdEdges;
        for ( auto const& theface : interprocessfaces(mesh,neighbor_pid) )
        {
            for ( uint16_type e=0 ; e<MeshTraits<MeshType>::face_type::numLocalEdges ; ++e )
                saveIdEdges.insert( theface.face(e).id() );
        }

        if ( neighbor_pid == invalid_rank_type_value )
        {
            for ( auto const& theedge : edges(mesh) )
            {
                if ( theedge.numberOfProcGhost() > 0 && saveIdEdges.find(theedge.id()) == saveIdEdges.end() )
                    myedges->push_back(boost::cref(theedge));
            }
        }
        else
        {
            for ( auto const& theedge : edges(mesh) )
            {
                if ( theedge.numberOfProcGhost() > 0 && theedge.isSubEntityOfGhostElement( neighbor_pid ) &&
                     saveIdEdges.find(theedge.id()) == saveIdEdges.end() )
                    myedges->push_back(boost::cref(theedge));
            }
        }
    }
    myedges->shrink_to_fit();
    return range( _range=boost::make_tuple( mpl::size_t<MESH_EDGES>(),
                                            myedges->begin(),
                                            myedges->end(),
                                            myedges ), 
                  _mesh=mesh.shared_from_this() );

}

//!
//! @return build a list of elements based on iterators  [begin,end) from a mesh \p imesh
//! @ingroup Mesh
//!
template<typename MeshType, typename IteratorType>
auto
idelements( MeshType const& imesh, IteratorType begin, IteratorType end )
{
    //auto myelts = make_elements_wrapper<MeshType>();
    typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype myelts( new typename MeshTraits<MeshType>::elements_reference_wrapper_type );
    auto const& mesh = Feel::unwrap_ptr( imesh );
    std::unordered_set<size_type> eltIdsDone;
    for( auto it = begin; it != end; ++ it )
    {
        size_type eltId = *it;
        if ( eltIdsDone.find( eltId ) != eltIdsDone.end() )
            continue;
        eltIdsDone.insert( eltId );
        if ( mesh.hasElement( eltId ) )
            myelts->push_back( boost::cref( mesh.element( eltId ) ) );
    }
    myelts->shrink_to_fit();
    return range( _range= boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                                             myelts->begin(),
                                             myelts->end(),
                                             myelts ),
                  _mesh=mesh.shared_from_this() );
}
//!
//! build a list of elements based on a list of element ids \p l from a mesh \p imesh
//!
template<typename MeshType, typename T>
auto
idelements( MeshType const& imesh, std::vector<T> const& l )
{
    return idelements( imesh, l.begin(), l.end() );
}


//! \return a pair of iterators to iterate over a range of elements 'rangeElt' which
//!  touch the range of faces or edges \rangeEntity by a point/edge/faces (with respect to type arg)
template<typename MeshType,typename RangeType, typename EntityRangeType>
Range<MeshType,MESH_ELEMENTS>
elements( MeshType const& mesh, RangeType const& rangeElt, EntityRangeType const& rangeEntity, ElementsType type = ElementsType::MESH_POINTS )
{
    std::unordered_set<size_type> entityIds;
    if constexpr ( std::is_same_v<EntityRangeType,Range<MeshType,MESH_FACES> > )
    {
        for( auto const& faceWrap : rangeEntity )
        {
            auto const& face = unwrap_ref( faceWrap );
            if ( type == ElementsType::MESH_POINTS )
            {
                for ( uint16_type p=0;p<face.nVertices();++p )
                    entityIds.insert( face.point(p).id() );
            }
            else if ( type == ElementsType::MESH_EDGES && is_3d<MeshType>::value )
            {
                for ( uint16_type p=0;p<face.nEdges();++p )
                    entityIds.insert( face.edge(p).id() );
            }
            else if ( type == ElementsType::MESH_FACES || ( is_2d<MeshType>::value && type == ElementsType::MESH_EDGES ) )
            {
                entityIds.insert( face.id() );
            }
            else
                CHECK( false ) << "invalid type " << type;
        }
    }
    else if constexpr ( std::is_same_v<EntityRangeType,Range<MeshType,MESH_EDGES> > )
        {
            for( auto const& edgeWrap : rangeEntity )
            {
                auto const& edge = unwrap_ref( edgeWrap );
                if ( type == ElementsType::MESH_POINTS )
                {
                    for ( uint16_type p=0;p<edge.nVertices();++p )
                        entityIds.insert( edge.point(p).id() );
                }
                else if ( type == ElementsType::MESH_EDGES )
                {
                    entityIds.insert( edge.id() );
                }
                else
                    CHECK( false ) << "invalid type :" << type << " (should be MESH_POINTS or MESH_EDGES)";
            }
        }
    else
    {
        CHECK( false ) << "invalid range";
    }

    Range<MeshType, MESH_ELEMENTS> myelts( mesh );

    for ( auto const& eltWrap : rangeElt )
    {
        auto const& elt = unwrap_ref( eltWrap );
        bool addElt = false;
        if ( type == ElementsType::MESH_POINTS )
        {
            for ( uint16_type p=0;p<elt.nVertices();++p )
            {
                if  ( entityIds.find( elt.point(p).id() ) != entityIds.end() )
                {
                    addElt = true;
                    break;
                }
            }
        }
        else if ( type == ElementsType::MESH_EDGES && is_3d<MeshType>::value )
        {
            for ( uint16_type p=0;p<elt.nEdges();++p )
            {
                if  ( entityIds.find( elt.edge(p).id() ) != entityIds.end() )
                {
                    addElt = true;
                    break;
                }
            }
        }
        else if ( type == ElementsType::MESH_FACES || ( is_2d<MeshType>::value && type == ElementsType::MESH_EDGES ) )
        {
            for ( uint16_type p=0;p<elt.nTopologicalFaces();++p )
            {
                if  ( entityIds.find( elt.face(p).id() ) != entityIds.end() )
                {
                    addElt = true;
                    break;
                }
            }
        }
        else
            CHECK( false ) << "invalid type " << type;
        if ( addElt )
            myelts.push_back(elt);
    }
    myelts.shrink_to_fit();
    return myelts;
}

//! \return a pair of iterators to iterate over all elements (not ghost) which
//!  touch the range of faces or edges rangeEntity by a point/edge/faces (with respect to type arg)
template<typename MeshType,typename EntityRangeType>
auto
elements( MeshType const& mesh, EntityRangeType const& rangeEntity, ElementsType type = ElementsType::MESH_POINTS )
{
    return elements( mesh, elements(mesh), rangeEntity, type );
}


//! return  mesh fragmentatin by marked elements  : a mapping between  fragment id to to the tuple (range, marker ids, fragment name)
template<typename MeshType, std::enable_if_t<std::is_base_of_v<MeshBase<>,unwrap_ptr_t<MeshType>>,int> = 0  >
auto
fragmentationMarkedElements( MeshType const& mesh, EntityProcessType entity = EntityProcessType::LOCAL_ONLY, bool includeNonMarkedElement = true )
{
    CHECK( entity == EntityProcessType::LOCAL_ONLY ) << "TODO";
    rank_type pid = rank( mesh );
    using elements_reference_wrapper_ptrtype = typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype;
    using elements_reference_wrapper_type = typename MeshTraits<MeshType>::elements_reference_wrapper_type;
    using marker_type = typename MeshTraits<MeshType>::element_type::marker_type;

    std::map<int,std::tuple<Range<std::remove_pointer_t<std::remove_const_t<decay_type<MeshType>>>,MESH_ELEMENTS>,marker_type,std::string>> res;

    auto const& imesh = Feel::unwrap_ptr( mesh );
    using mesh_type = typename MeshTraits<MeshType>::mesh_type;
    auto it = imesh.beginOrderedElement();
    auto en = imesh.endOrderedElement();

    marker_type emptyMarker;
    std::map<marker_type, std::tuple<elements_reference_wrapper_ptrtype,int>> mapMarkersToElements;

    auto const& mappingFragments = imesh.meshFragmentationByMarker();
    for ( auto const& [fragmentId,fragmentMarkerIds] : mappingFragments )
    {
        if ( fragmentMarkerIds.empty() && !includeNonMarkedElement )
            continue;
        mapMarkersToElements.emplace( fragmentMarkerIds, std::make_tuple( std::make_shared<elements_reference_wrapper_type>(),fragmentId ) );
    }

    for ( ; it != en; ++it )
    {
        auto const& elt = unwrap_ref( *it );
        if ( elt.processId() != pid )
            continue;

        if ( !includeNonMarkedElement && !elt.hasMarker() )
            continue;

        bool hasMarker = elt.hasMarker();
        if ( !includeNonMarkedElement && !hasMarker )
            continue;
        marker_type const& eltMarkers = hasMarker? elt.marker() : emptyMarker;

        auto itFind = mapMarkersToElements.find( eltMarkers );
        DCHECK( itFind != mapMarkersToElements.end() ) << "missing fragment";
        std::get<0>(itFind->second)->push_back( boost::cref(elt) );
    }

    for ( auto & [mIds,fragementData] : mapMarkersToElements )
    {
        auto & [myelts,fragmentId] = fragementData;
        myelts->shrink_to_fit();
        auto therange = range( _range=boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),myelts->begin(), myelts->end(), myelts ), _mesh=mesh->shared_from_this() );

        std::string fragmentName;
        for ( auto mId : mIds )
        {
            std::string mName = imesh.markerName( mId, mesh_type::element_type::nDim );
            if ( !mName.empty() && !fragmentName.empty() )
                fragmentName += "_";
            fragmentName += mName;
        }
        res.emplace( fragmentId, std::make_tuple( std::move( therange ), mIds, fragmentName ) );
    }
    return res;
}

} // namespace Feel


#endif
