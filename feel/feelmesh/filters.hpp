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

#include <utility>
#if BOOST_VERSION >= 105600
#include <boost/phoenix/stl/algorithm/detail/is_std_list.hpp>
#else
#include <boost/spirit/home/phoenix/stl/algorithm/detail/is_std_list.hpp>
#endif

#include <feel/feelcore/environment.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/iterator.hpp>
#include <feel/feelmesh/detail/filters.hpp>

namespace Feel
{
using elements_c = mpl::size_t<MESH_ELEMENTS>;
using faces_c = mpl::size_t<MESH_FACES>;
using edges_c = mpl::size_t<MESH_EDGES>;
using points_c = mpl::size_t<MESH_POINTS>;
/**
 * a RangeType can be one or more filter/range objects of the same type, we
 * extract the underlying type by first casting everything to a list and then
 * retrieving consistently the type.
 */
template <typename RangeType>
using range_t = typename mpl::if_< boost::is_std_list<RangeType>,
                                   mpl::identity<RangeType>,
                                   mpl::identity<std::list<RangeType> > >::type::type::value_type;

template<typename MeshType>
using elements_t =  boost::tuple<elements_c,
                                 typename MeshTraits<MeshType>::element_const_iterator,
                                 typename MeshTraits<MeshType>::element_const_iterator>;

template<typename MeshType>
using ext_elements_t = boost::tuple<elements_c,
                                    typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> >::const_iterator,
                                    typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> >::const_iterator,
                                    boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > >
                                    >;

template<typename MeshType>
using boundaryelements_t =  boost::tuple<elements_c,
                                         typename MeshTraits<MeshType>::location_element_const_iterator,
                                         typename MeshTraits<MeshType>::location_element_const_iterator>;

template<typename MeshType>
using internalelements_t = boundaryelements_t<MeshType>;

template<typename MeshType>
using markedelements_t = boost::tuple<elements_c,
                                      typename MeshTraits<MeshType>::marker_element_const_iterator,
                                      typename MeshTraits<MeshType>::marker_element_const_iterator>;

template<typename MeshType>
using marked2elements_t = boost::tuple<elements_c,
                                       typename MeshTraits<MeshType>::marker2_element_const_iterator,
                                       typename MeshTraits<MeshType>::marker2_element_const_iterator>;

template<typename MeshType>
using marked3elements_t =  boost::tuple<elements_c,
                                        typename MeshTraits<MeshType>::marker3_element_const_iterator,
                                        typename MeshTraits<MeshType>::marker3_element_const_iterator> ;


template<typename MeshType>
using faces_t =  boost::tuple<faces_c,
                              typename MeshTraits<MeshType>::pid_face_const_iterator,
                              typename MeshTraits<MeshType>::pid_face_const_iterator>;

template<typename MeshType>
using ext_faces_t = boost::tuple<faces_c,
                                 typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> >::const_iterator,
                                 typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> >::const_iterator,
                                 boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> > >
                                 >;


template<typename IteratorType>
using filter_enum_t = typename boost::tuples::element<0,IteratorType>::type;
template<typename IteratorType>
using filter_iterator_t = typename boost::tuples::element<1,IteratorType>::type;
template<typename IteratorType>
using filter_entity_t = typename filter_iterator_t<IteratorType>::value_type;

template<typename IteratorType>
using ext_entities_from_iterator_t = boost::tuple<filter_enum_t<IteratorType>,
                                                  typename std::vector<boost::reference_wrapper<filter_entity_t<IteratorType> const> >::const_iterator,
                                                  typename std::vector<boost::reference_wrapper<filter_entity_t<IteratorType> const> >::const_iterator,
                                                  boost::shared_ptr<std::vector<boost::reference_wrapper<filter_entity_t<IteratorType> const> > >
                                                  >;

template<typename MeshType>
using idfaces_t =  boost::tuple<faces_c,
                                typename MeshTraits<MeshType>::face_const_iterator,
                                typename MeshTraits<MeshType>::face_const_iterator>;

template<typename MeshType>
using boundaryfaces_t =  boost::tuple<faces_c,
                                      typename MeshTraits<MeshType>::location_face_const_iterator,
                                      typename MeshTraits<MeshType>::location_face_const_iterator>;

template<typename MeshType>
using internalfaces_t = boundaryfaces_t<MeshType>;

template<typename MeshType>
using interprocessfaces_t =  boost::tuple<faces_c,
                                          typename MeshTraits<MeshType>::interprocess_face_const_iterator,
                                          typename MeshTraits<MeshType>::interprocess_face_const_iterator>;


template<typename MeshType>
using markedfaces_t = boost::tuple<faces_c,
                                   typename MeshTraits<MeshType>::marker_face_const_iterator,
                                   typename MeshTraits<MeshType>::marker_face_const_iterator>;

template<typename MeshType>
using marked2faces_t = boost::tuple<faces_c,
                                    typename MeshTraits<MeshType>::marker2_face_const_iterator,
                                    typename MeshTraits<MeshType>::marker2_face_const_iterator>;

template<typename MeshType>
using marked3faces_t =  boost::tuple<faces_c,
                                     typename MeshTraits<MeshType>::marker3_face_const_iterator,
                                     typename MeshTraits<MeshType>::marker3_face_const_iterator> ;

template<typename MeshType>
using edges_t =  boost::tuple<edges_c,
                              typename MeshTraits<MeshType>::edge_const_iterator,
                              typename MeshTraits<MeshType>::edge_const_iterator>;
template<typename MeshType>
using pid_edges_t =  boost::tuple<edges_c,
                                  typename MeshTraits<MeshType>::pid_edge_const_iterator,
                                  typename MeshTraits<MeshType>::pid_edge_const_iterator>;

template<typename MeshType>
using ext_edges_t =  boost::tuple<edges_c,
                                  typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::edge_type const> >::const_iterator,
                                  typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::edge_type const> >::const_iterator,
                                  boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::edge_type const> > >
                                  >;

template<typename MeshType>
using markededges_t = boost::tuple<edges_c,
                                   typename MeshTraits<MeshType>::marker_edge_const_iterator,
                                   typename MeshTraits<MeshType>::marker_edge_const_iterator>;

template<typename MeshType>
using boundaryedges_t = boost::tuple<edges_c,
                                     typename MeshTraits<MeshType>::location_edge_const_iterator,
                                     typename MeshTraits<MeshType>::location_edge_const_iterator>;

template<typename MeshType>
using internaledges_t = boundaryedges_t<MeshType>;

template<typename MeshType>
using points_t =  boost::tuple<points_c,
                               typename MeshTraits<MeshType>::point_const_iterator,
                               typename MeshTraits<MeshType>::point_const_iterator>;

template<typename MeshType>
using markedpoints_t = boost::tuple<points_c,
                                    typename MeshTraits<MeshType>::marker_point_const_iterator,
                                    typename MeshTraits<MeshType>::marker_point_const_iterator>;

template<typename MeshType>
using boundarypoints_t = boost::tuple<points_c,
                                      typename MeshTraits<MeshType>::location_point_const_iterator,
                                      typename MeshTraits<MeshType>::location_point_const_iterator>;
template<typename MeshType>
using internalpoints_t = boundarypoints_t<MeshType>;

template<typename IteratorRangeT>
using submeshrange_t = typename Feel::detail::submeshrangetype<IteratorRangeT>::type;
/**
 * namespace for meta mesh computation data structure
 */
namespace meta
{

template<typename MeshType>
struct elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<elements_c,
                         typename MeshTraits<MeshType>::element_const_iterator,
                         typename MeshTraits<MeshType>::element_const_iterator> type;
};

template<typename MeshType>
struct markedelements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<elements_c,
                         typename MeshTraits<MeshType>::marker_element_const_iterator,
                         typename MeshTraits<MeshType>::marker_element_const_iterator> type;
};

template<typename MeshType>
struct marked2elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<elements_c,
                         typename MeshTraits<MeshType>::marker2_element_const_iterator,
                         typename MeshTraits<MeshType>::marker2_element_const_iterator> type;
};

template<typename MeshType>
struct marked3elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<elements_c,
                         typename MeshTraits<MeshType>::marker3_element_const_iterator,
                         typename MeshTraits<MeshType>::marker3_element_const_iterator> type;
};

} // meta

/**
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType>
elements_t<MeshType>
allelements( MeshType const& mesh )
{
    return boost::make_tuple( elements_c{},
                              remove_ptr( mesh ).beginElement(),
                              remove_ptr( mesh ).endElement() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType>
elements_t<MeshType>
elements( MeshType const& mesh )
{
    return elements_t<MeshType>( elements_c{},
                                 remove_ptr( mesh ).beginElementWithProcessId(rank(mesh)),
                                 remove_ptr( mesh ).endElementWithProcessId(rank(mesh)) );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the mesh
 * which share a face with the boundary
 */
template<typename MeshType>
boundaryelements_t<MeshType>
boundaryelements( MeshType const& mesh, uint16_type entity_min_dim = 0, uint16_type entity_max_dim = 2 )
{
    // pair of iterators
    auto p = remove_ptr(mesh).boundaryElements( entity_min_dim, entity_max_dim, rank(mesh) );
    return boundaryelements_t<MeshType>( elements_c{}, p.first, p.second );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the mesh
 * which are stricly within the domain that is to say they do not
 * share a face with the boundary
 */
template<typename MeshType>
internalelements_t<MeshType>
internalelements( MeshType const& mesh )
{
    // pair of iterators
    auto p = remove_ptr(mesh).internalElements( rank(mesh) );
    return internalelements_t<MeshType>( elements_c{}, p.first, p.second );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
markedelements_t<MeshType>
markedelements( MeshType const& mesh, std::string const& flag )
{
    auto p = remove_ptr(mesh).elementsWithMarker( remove_ptr(mesh).markerName(flag), 
                                                  rank(mesh) );
    return markedelements_t<MeshType>( elements_c{}, p.first, p.second );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
markedelements_t<MeshType>
markedelements( MeshType const& mesh, const char* flag )
{
    auto p = remove_ptr(mesh).elementsWithMarker( remove_ptr(mesh).markerName(flag), 
                                                  rank(mesh) );
    return markedelements_t<MeshType>( elements_c{}, p.first, p.second );
}

template<typename MeshType>
markedelements_t<MeshType>
markedelements( MeshType const& mesh, boost::any const& flag )
{
    auto p = remove_ptr(mesh).elementsWithMarker( remove_ptr(mesh).markerId(flag), 
                                                  rank(mesh) );
    return markedelements_t<MeshType>( elements_c{}, p.first, p.second );
}

template<typename MeshType>
std::list<markedelements_t<MeshType>>
markedelements( MeshType const& mesh, std::initializer_list<boost::any> const& flag )
{
    std::list<markedelements_t<MeshType>> list_elements;
    for ( auto const& it : flag )
    {
        auto p = remove_ptr(mesh).elementsWithMarker( remove_ptr(mesh).markerId(it), 
                                                      rank(mesh) );
        list_elements.push_back( markedelements_t<MeshType>( elements_c{}, p.first, p.second ) );
    }
    return list_elements;
}

template<typename MeshType>
std::list<markedelements_t<MeshType>>
markedelements( MeshType const& mesh, std::list<std::string> const& flag )
{
    std::list<markedelements_t<MeshType>> list_elements;
    for ( auto const& it : flag )
    {
        auto p = remove_ptr(mesh).elementsWithMarker( remove_ptr(mesh).markerId(it), 
                                                      rank(mesh) );
        list_elements.push_back( markedelements_t<MeshType>( elements_c{}, p.first, p.second ) );
    }
    return list_elements;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker2 \p flag
 */
template<typename MeshType>
marked2elements_t<MeshType>
marked2elements( MeshType const& mesh, flag_type flag )
{
    auto p = remove_ptr(mesh).elementsWithMarker2( remove_ptr(mesh).markerId(flag), 
                                                   rank(mesh) );
    return marked2elements_t<MeshType>( elements_c{}, p.first, p.second );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker2 string
 */
template<typename MeshType>
marked2elements_t<MeshType>
marked2elements( MeshType const& mesh, std::string const& flag )
{
    auto p = remove_ptr(mesh).elementsWithMarker2( remove_ptr(mesh).markerName(flag), 
                                                   rank(mesh) );
    return marked2elements_t<MeshType>( elements_c{}, p.first, p.second );
}
/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with marker \p flag
 */
template<typename MeshType>
marked2elements_t<MeshType>
marked2elements( MeshType const& mesh, const char* flag )
{
    auto p = remove_ptr(mesh).elementsWithMarker2( remove_ptr(mesh).markerName(flag), 
                                                   rank(mesh) );
    return marked2elements_t<MeshType>( elements_c{}, p.first, p.second );
}

template<typename MeshType>
marked2elements_t<MeshType>
marked2elements( MeshType const& mesh, boost::any const& flag )
{
    auto p = remove_ptr(mesh).elementsWithMarker2( remove_ptr(mesh).markerId(flag), 
                                                   rank(mesh) );
    return marked2elements_t<MeshType>( elements_c{}, p.first, p.second );
}

template<typename MeshType>
std::list<marked2elements_t<MeshType>>
marked2elements( MeshType const& mesh, std::initializer_list<boost::any> const& flag )
{
    std::list<marked2elements_t<MeshType>> list_elements;
    for ( auto const& it : flag )
    {
        auto p = remove_ptr(mesh).elementsWithMarker2( remove_ptr(mesh).markerId(it), 
                                                       rank(mesh) );
        list_elements.push_back(  marked2elements_t<MeshType>( elements_c{}, p.first, p.second ) );
    }
    return list_elements;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker3 \p flag
 */
template<typename MeshType>
marked3elements_t<MeshType>
marked3elements( MeshType const& mesh, flag_type flag )
{
    auto p = remove_ptr(mesh).elementsWithMarker3( remove_ptr(mesh).markerId(flag), 
                                                   rank(mesh) );
    return marked3elements_t<MeshType>( elements_c{}, p.first, p.second );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements of the
 * mesh with \c Marker3 string
 */
template<typename MeshType>
marked3elements_t<MeshType>
marked3elements( MeshType const& mesh, std::string const& flag )
{
    auto p = remove_ptr(mesh).elementsWithMarker3( remove_ptr(mesh).markerName(flag), 
                                                   rank(mesh) );
    return marked3elements_t<MeshType>( elements_c{}, p.first, p.second );
}

template<typename MeshType>
std::list<marked3elements_t<MeshType>>
marked3elements( MeshType const& mesh, std::initializer_list<boost::any> const& flag )
{
    std::list<marked3elements_t<MeshType>> list_elements;
    for ( auto const& it : flag )
    {
        auto p = remove_ptr(mesh).elementsWithMarker3( remove_ptr(mesh).markerId(it), 
                                                       rank(mesh) );
        list_elements.push_back(  marked3elements_t<MeshType>( elements_c{}, p.first, p.second ) );
    }
    return list_elements;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with id \p flag
 */
template<typename MeshType>
elements_t<MeshType>
idedelements( MeshType const& mesh, flag_type flag )
{
    return boost::make_tuple( elements_c{},
                              remove_ptr(mesh).beginElementWithId( flag ),
                              remove_ptr(mesh).endElementWithId( flag ) );
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
template<typename MeshType>
faces_t<MeshType>
faces( MeshType const& mesh )
{
    auto p = remove_ptr(mesh).facesWithProcessId( rank(mesh) );
    return boost::make_tuple( faces_c{}, p.first, p.second );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with id \p id
 */
template<typename MeshType>
idfaces_t<MeshType>
idedfaces( MeshType const& mesh, size_type id )
{
    return idfaces_t<MeshType>( faces_c{}, 
                                remove_ptr(mesh).beginFaceWithId( id ),
                                remove_ptr(mesh).endFaceWithId( id ) );
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
template<typename MeshType>
markedfaces_t<MeshType>
markedfaces( MeshType const& mesh )
{
    return markedfaces_t<MeshType>( faces_c{}, 
                                    remove_ptr(mesh).beginFaceWithMarker(),
                                    remove_ptr(mesh).endFaceWithMarker() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over faces of the
 * mesh marked with \c __marker
 *
 * @param mesh a mesh data structure
 * @param __marker a marker that identifies faces
 * @param __pid process id
 *
 * @return a pair of iterators (begin,end) for the set of marked faces
 */
template<typename MeshType>
markedfaces_t<MeshType>
markedfaces( MeshType const& mesh,
             flag_type m )
{
    return markedfaces_t<MeshType>( faces_c{}, 
                                    remove_ptr(mesh).beginFaceWithMarker(m,rank(mesh)),
                                    remove_ptr(mesh).endFaceWithMarker(m,rank(mesh)) );
}

template<typename MeshType>
markedfaces_t<MeshType>
markedfaces( MeshType const& mesh,
             std::string const&__marker )
{
    return markedfaces( mesh, remove_ptr(mesh).markerName(__marker));
}

template<typename MeshType>
markedfaces_t<MeshType>
markedfaces( MeshType const& mesh,
             const char*__marker )
{
    return markedfaces( mesh, std::string( __marker ) );
}
template<typename MeshType>
std::list<markedfaces_t<MeshType>>
markedfaces( MeshType const& mesh,
             std::list<std::string> const& __markers )
{
    std::list<markedfaces_t<MeshType>> list_faces;
    for ( auto const& it : __markers )
    {
        list_faces.push_back( markedfaces(mesh,it));
    }
    return list_faces;
}

template<typename MeshType>
markedfaces_t<MeshType>
markedfaces( MeshType const& mesh,
             boost::any const& __marker )
{
    return markedfaces( mesh, mesh->markerId( __marker ) );
}

template<typename MeshType>
std::list<markedfaces_t<MeshType>>
markedfaces( MeshType const& mesh,
             std::initializer_list<boost::any> __markers )
{
    std::list<markedfaces_t<MeshType>> list_faces;
    for ( auto const& it : __markers )
    {
        list_faces.push_back( markedfaces( mesh, it ) );
    }
    return list_faces;
}


template<typename MeshType>
boost::tuple<faces_c,
             typename MeshTraits<MeshType>::marker2_face_const_iterator,
             typename MeshTraits<MeshType>::marker2_face_const_iterator>
marked2faces( MeshType const& mesh,
              flag_type __marker )
{
    return marked2faces_t<MeshType>( faces_c{}, 
                                     remove_ptr(mesh).beginFaceWithMarker2(__marker,rank(mesh)),
                                     remove_ptr(mesh).endFaceWithMarker2(__marker,rank(mesh)) );
}

template<typename MeshType>
std::list<boost::tuple<faces_c,
                       typename MeshTraits<MeshType>::marker2_face_const_iterator,
                       typename MeshTraits<MeshType>::marker2_face_const_iterator> >
marked2faces( MeshType const& mesh,
              std::initializer_list<boost::any> __markers )
{
    std::list<boost::tuple<faces_c,
                           typename MeshTraits<MeshType>::marker2_face_const_iterator,
                           typename MeshTraits<MeshType>::marker2_face_const_iterator> > list_faces;
    for ( auto const& it : __markers )
    {
        list_faces.push_back( marked2faces( mesh, remove_ptr(mesh).markerId(it) ) );
    }
    return list_faces;
}

template<typename MeshType>
decltype(auto)
marked2faces( MeshType const& mesh,
              std::string const& __marker )
{
    return marked2faces( mesh, mesh->markerName( __marker ) );
}

template<typename MeshType>
marked3faces_t<MeshType>
marked3faces( MeshType const& mesh,
              flag_type __marker )
{
    return marked3faces_t<MeshType>( faces_c{}, 
                                     remove_ptr(mesh).beginFaceWithMarker3(__marker,rank(mesh)),
                                     remove_ptr(mesh).endFaceWithMarker3(__marker,rank(mesh)) );
}


template<typename MeshType>
marked3faces_t<MeshType>
marked3faces( MeshType const& mesh,
              std::string const& __marker )
{
    return marked3faces( mesh, mesh->markerName( __marker ) );
}

template<typename MeshType>
std::list<marked3faces_t<MeshType>>
marked3faces( MeshType const& mesh,
              std::initializer_list<boost::any> __markers )
{
    std::list<marked3faces_t<MeshType>> list_faces;
    for ( auto const& it : __markers )
    {
        list_faces.push_back( marked3faces( mesh, it ) );
    }
    return list_faces;
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all boundary faces of the
 * mesh
 */
template<typename MeshType>
boundaryfaces_t<MeshType>
boundaryfaces( MeshType const& mesh  )
{
    auto p = remove_ptr(mesh).facesOnBoundary( rank(mesh) );
    return boundaryfaces_t<MeshType>( faces_c{}, p.first, p.second );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all internal faces of the
 * mesh belong to process domain \p __pid
 */
template<typename MeshType>
internalfaces_t<MeshType>
internalfaces( MeshType const& mesh )
{
    auto p = remove_ptr(mesh).internalFaces( rank(mesh) );
    return internalfaces_t<MeshType>( faces_c{}, p.first, p.second );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all interprocess faces of the
 * mesh belonging to process \p __pid
 */
template<typename MeshType>
interprocessfaces_t<MeshType>
interprocessfaces( MeshType const& mesh, rank_type neighbor_pid = invalid_rank_type_value )
{
    auto p =  remove_ptr(mesh).interProcessFaces( neighbor_pid );
    return boost::make_tuple( faces_c{}, p.first, p.second );
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
template<typename MeshType>
pid_edges_t<MeshType>
edges( MeshType const& mesh )
{
    auto p =  remove_ptr(mesh).edgesWithProcessId( rank(mesh) );
    return pid_edges_t<MeshType>( edges_c{}, p.first, p.second );
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
template<typename MeshType>
markededges_t<MeshType>
markededges( MeshType const& mesh,
             flag_type __marker )
{
    auto p = remove_ptr(mesh).edgesWithMarker( __marker, rank(mesh) );
    return markededges_t<MeshType>( edges_c{}, p.first, p.second );
}

template<typename MeshType>
markededges_t<MeshType>
markededges( MeshType const& mesh,
             std::string const& __marker,
             typename std::enable_if<is_3d<MeshType>::value>::type* = nullptr )
{
    
    return markededges( mesh, remove_ptr(mesh).markerName( __marker ) );
}

/**
 * @return the range of faces of the \p mesh associated to \p __marker
 */
template<typename MeshType>
decltype(auto)
markededges( MeshType const& mesh,
             std::string const& __marker,
             typename std::enable_if<is_2d<MeshType>::value>::type* = nullptr )
{
    return markedfaces( mesh, __marker );
}


template<typename MeshType>
std::list<markededges_t<MeshType>>
markededges( MeshType const& mesh,
             std::list<std::string> const& __markers,
             typename std::enable_if<is_3d<MeshType>::value>::type* = nullptr )
{
    std::list<markededges_t<MeshType>> list_edges;
    for ( auto const& it : __markers )
    {
        list_edges.push_back( markededges( mesh, it ) );
    }
    return list_edges;
}
/**
 * @return the range of faces of the \p mesh associated to \p __marker
 */
template<typename MeshType>
decltype(auto)
markededges( MeshType const& mesh,
             std::list<std::string> const& __markers,
             typename std::enable_if<is_2d<MeshType>::value>::type* = nullptr )
{
    return markedfaces(mesh, __markers);
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all boundary edges of the
 * mesh
 */
template<typename MeshType>
boundaryedges_t<MeshType>
boundaryedges( MeshType const& mesh )
{
    return boundaryedges_t<MeshType>( edges_c{},
                                      remove_ptr(mesh).beginEdgeOnBoundary(),
                                      remove_ptr(mesh).endEdgeOnBoundary() );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all internal edges of the
 * mesh
 */
template<typename MeshType>
internaledges_t<MeshType>
internaledges( MeshType const& mesh )
{
    return internaledges_t<MeshType>( edges_c{},
                                      remove_ptr(mesh).beginInternalEdge(),
                                      remove_ptr(mesh).endInternalEdge() );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
points_t<MeshType>
points( MeshType const& mesh )
{
    return points_t<MeshType>( points_c{},
                               remove_ptr(mesh).beginPoint(),
                               remove_ptr(mesh).endPoint() );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the marked points with \p flag of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
markedpoints_t<MeshType>
markedpoints( MeshType const& mesh, size_type flag )
{
    return markedpoints_t<MeshType>( points_c{},
                                     remove_ptr(mesh).beginPointWithMarker( flag ),
                                     remove_ptr(mesh).endPointWithMarker( flag ) );
}

template<typename MeshType>
markedpoints_t<MeshType>
markedpoints( MeshType const& mesh, std::string const& flag )
{
    auto n = remove_ptr(mesh).markerName(flag);
    return markedpoints_t<MeshType>( points_c{},
                                     remove_ptr(mesh).beginPointWithMarker( n ),
                                     remove_ptr(mesh).endPointWithMarker( n ) );
}

template<typename MeshType>
std::list<markedpoints_t<MeshType>>
markedpoints( MeshType const& mesh, std::list<std::string> const& __markers )
{
    std::list<markedpoints_t<MeshType>> list_points;
    for ( auto const& it : __markers )
    {
        list_points.push_back( markedpoints(mesh, it ) );
    }
    return list_points;
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the boundary points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
boundarypoints_t<MeshType>
boundarypoints( MeshType&& mesh )
{
    return boundarypoints_t<MeshType>( points_c{},
                                       remove_ptr(mesh).beginPointOnBoundary( ),
                                       remove_ptr(mesh).endPointOnBoundary() );
}

/**
 * \ingroup MeshIterators
 * \return the range of iterators [begin,end] over the internal(not on the boundary) points of the mesh
 * \warning this filter is not parallelized
 */
template<typename MeshType>
internalpoints_t<MeshType>
internalpoints( MeshType const& mesh )
{
    return boundarypoints_t<MeshType>( points_c{},
                                       remove_ptr(mesh).beginInternalPoint( ),
                                       remove_ptr(mesh).endInternalPoint() );
}

/**
 * \ingroup MeshIterators
 * \return the number of elements given element iterators constructed
 * using the mesh filters
 * \param its the mesh iterators
 * \param global all reduce number of elements if set to true
 *
 * The following code prints in the logfile the number of elements in
 * the mesh that are marked with marker1 equal to 1:
 *
 * \code
 * LOG(INFO) << "number of elements = " << nelements( markedelements(mesh,1) ) << "\n";
 * \endcode
 *
 */
template<typename MT, typename Iterator>
size_type
nelements( boost::tuple<MT,Iterator,Iterator> const& its, bool global = false )
{
    size_type d = std::distance( boost::get<1>( its ), boost::get<2>( its ) );
    size_type gd = d;
    if ( global )
        mpi::all_reduce(Environment::worldComm().globalComm(),
                        d,
                        gd,
                        std::plus<size_type>());
    return gd;

}
/**
 * \ingroup MeshIterators
 * \return the number of elements given element iterators constructed
 * using the mesh filters
 * \param a list of mesh iterators
 * \param global all reduce number of elements if set to true
 *
 * The following code prints in the logfile the number of elements in
 * the mesh that are marked with marker1 equal to 1:
 *
 * \code
 * LOG(INFO) << "number of elements = " << nelements( markedelements(mesh,{1,2,3}) ) << "\n";
 * \endcode
 *
 */
template<typename MT, typename Iterator>
size_type
nelements( std::list<boost::tuple<MT,Iterator,Iterator> > const& its, bool global = false )
{
    size_type d = std::accumulate( its.begin(), its.end(), size_type(0),
                                   [&d]( boost::tuple<MT,Iterator,Iterator> const& t )
                                   {
                                       d+=std::distance( boost::get<1>( t ), boost::get<2>( t ) );
                                   } );
    size_type gd = d;
    if ( global )
        mpi::all_reduce(Environment::worldComm().globalComm(),
                        d,
                        gd,
                        std::plus<size_type>());
    return gd;
}


/**
 * \ingroup MeshIterators
 * \return the number of elements given element iterators constructed
 * using custom range
 * \param its the mesh iterators
 * \param global all reduce number of elements if set to true
 *
 * The following code prints in the logfile the number of elements in
 * the mesh that are marked with marker1 equal to 1:
 *
 * \code
 * LOG(INFO) << "number of elements = " << nelements( myCustomRange ) << "\n";
 * \endcode
 *
 */
template<typename MT, typename Iterator,typename Container>
size_type
nelements( boost::tuple<MT,Iterator,Iterator,Container> const& its, bool global = false )
{
    size_type d = std::distance( boost::get<1>( its ), boost::get<2>( its ) );
    size_type gd = d;
    if ( global )
        mpi::all_reduce(Environment::worldComm().globalComm(),
                        d,
                        gd,
                        std::plus<size_type>());
    return gd;

}



template<typename ElementType>
boost::tuple<elements_c,
             typename std::list<ElementType>::const_iterator,
             typename std::list<ElementType>::const_iterator>
element( ElementType const& elt  )
{
    std::list<ElementType> lst;
    lst.push_back(  elt );
    return boost::make_tuple( elements_c{},
                              lst.begin(),
                              lst.end() );
}


template<typename MeshType>
ext_elements_t<MeshType>
elements( MeshType&&  mesh, EntityProcessType entity )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : elements(mesh) )
        {
            myelts->push_back(boost::cref(elt));
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto face_it = mesh.interProcessFaces().first;
        auto const face_en = mesh.interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;

            // add elt in range
            myelts->push_back(boost::cref(eltOffProc));

            eltGhostDone.insert( eltOffProc.id() );
        }
    }

    return boost::make_tuple( elements_c{},
                              myelts->begin(),
                              myelts->end(),
                              myelts );
}


template<typename MeshType>
ext_elements_t<MeshType>
elements( MeshType const& mesh, EntityProcessType entity )
{
    return elements( remove_ptr(mesh), entity );
}


template<typename MeshType>
ext_elements_t<MeshType>
boundaryelements( MeshType const& mesh, EntityProcessType entity )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : boundaryelements(mesh) )
        {
            myelts->push_back(boost::cref(elt));
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto face_it = mesh->interProcessFaces().first;
        auto const face_en = mesh->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
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

    return boost::make_tuple( elements_c{},
                              myelts->begin(),
                              myelts->end(),
                              myelts );
}

template<typename MeshType>
ext_elements_t<MeshType>
elementsWithMarkedFaces( MeshType const& mesh, boost::any const& flag, EntityProcessType entity = EntityProcessType::LOCAL_ONLY )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    flag_type theflag = mesh->markerId( flag );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        for ( auto const& elt : markedfaces(mesh,flag) )
        {
            if ( elt.isConnectedTo0() )
                myelts->push_back(boost::cref(elt.element0()));
            if ( elt.isConnectedTo1() )
                myelts->push_back(boost::cref(elt.element1()));
        }
    }
    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto face_it = mesh->interProcessFaces().first;
        auto const face_en = mesh->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;


            // add elt in range
            if ( eltOffProc.marker().value() == theflag )
            {
                myelts->push_back(boost::cref(eltOffProc));
            }

            eltGhostDone.insert( eltOffProc.id() );
        }
    }

    return boost::make_tuple( elements_c{},
                              myelts->begin(),
                              myelts->end(),
                              myelts );
}

template<typename MeshType>
ext_elements_t<MeshType>
markedelements( MeshType const& mesh, boost::any const& flag, EntityProcessType entity )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    flag_type theflag = mesh->markerId( flag );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : markedelements(mesh,flag) )
        {
            myelts->push_back(boost::cref(elt));
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto face_it = mesh->interProcessFaces().first;
        auto const face_en = mesh->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;


            // add elt in range
            if ( eltOffProc.marker().value() == theflag )
                myelts->push_back(boost::cref(eltOffProc));

            eltGhostDone.insert( eltOffProc.id() );
        }
    }

    return boost::make_tuple( elements_c{},
                              myelts->begin(),
                              myelts->end(),
                              myelts );
}

template<typename MeshType>
ext_elements_t<MeshType>
marked2elements( MeshType const& mesh, boost::any const& flag, EntityProcessType entity )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    flag_type theflag = mesh->markerId( flag );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& elt : marked2elements(mesh,flag) )
        {
            myelts->push_back(boost::cref(elt));
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        std::set<size_type> eltGhostDone;

        auto face_it = mesh->interProcessFaces().first;
        auto const face_en = mesh->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            if ( eltGhostDone.find( eltOffProc.id() ) != eltGhostDone.end() ) continue;


            // add elt in range
            if ( eltOffProc.marker2().value() == theflag )
                myelts->push_back(boost::cref(eltOffProc));

            eltGhostDone.insert( eltOffProc.id() );
        }
    }

    return boost::make_tuple( elements_c{},
                              myelts->begin(),
                              myelts->end(),
                              myelts );
}


template<typename MeshType>
ext_faces_t<MeshType>
boundaryfaces( MeshType const& mesh, EntityProcessType entity )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& theface : boundaryfaces(mesh) )
        {
            myelts->push_back(boost::cref(theface));
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        //std::set<size_type> faceGhostDone;
        auto face_it = mesh->interProcessFaces().first;
        auto const face_en = mesh->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            for ( size_type f = 0; f < mesh->numLocalFaces(); f++ )
            {
                auto const& theface = eltOffProc.face(f);
                if ( theface.isOnBoundary() ) //&& faceGhostDone.find( theface.id() ) == faceGhostDone.end() )
                    myelts->push_back(boost::cref(theface));
            }
        }
    }

    return boost::make_tuple( faces_c{},
                              myelts->begin(),
                              myelts->end(),
                              myelts );

}

template<typename MeshType>
ext_faces_t<MeshType>
marked2faces( MeshType const& mesh, boost::any flag, EntityProcessType entity )
{
    typedef std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );
    flag_type theflag = mesh->markerId( flag );
    if ( ( entity == EntityProcessType::LOCAL_ONLY ) || ( entity == EntityProcessType::ALL ) )
        for ( auto const& theface : marked2faces(mesh, theflag) )
        {
            myelts->push_back(boost::cref(theface));
        }

    if ( ( entity == EntityProcessType::GHOST_ONLY ) || ( entity == EntityProcessType::ALL ) )
    {
        //std::set<size_type> faceGhostDone;
        auto face_it = mesh->interProcessFaces().first;
        auto const face_en = mesh->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;

            for ( size_type f = 0; f < mesh->numLocalFaces(); f++ )
            {
                auto const& theface = eltOffProc.face(f);
                if ( theface.marker2().value() == theflag ) //&& faceGhostDone.find( theface.id() ) == faceGhostDone.end() )
                    myelts->push_back(boost::cref(theface));
            }
        }
    }

    return boost::make_tuple( faces_c{},
                              myelts->begin(),
                              myelts->end(),
                              myelts );

}



template<typename MeshType>
ext_edges_t<MeshType>
interprocessedges( MeshType const& mesh, rank_type neighbor_pid = invalid_rank_type_value, EntityProcessType entity = EntityProcessType::ALL )
{
    typedef typename MeshTraits<MeshType>::face_type face_type;
    typedef typename MeshTraits<MeshType>::edge_type edge_type;
    typedef std::vector<boost::reference_wrapper<edge_type const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myedges( new cont_range_type );

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
            for ( uint16_type e=0 ; e<face_type::numLocalEdges ; ++e )
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

    return boost::make_tuple( edges_c{},
                              myedges->begin(),
                              myedges->end(),
                              myedges );

}

template<size_t S, class ITERATOR, class CONTAINER>
WorldComm const&
worldComm( boost::tuple<mpl::size_t<S>,ITERATOR,ITERATOR,CONTAINER> const& range )
{
    return range.template get<1>()->mesh()->worldComm();
}


} // namespace Feel


#endif
