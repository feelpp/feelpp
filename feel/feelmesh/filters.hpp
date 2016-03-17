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
using elements_t =  boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                                 typename MeshTraits<MeshType>::element_const_iterator,
                                 typename MeshTraits<MeshType>::element_const_iterator>;

template<typename MeshType>
using pid_elements_t =  boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                                     typename MeshTraits<MeshType>::pid_element_const_iterator,
                                     typename MeshTraits<MeshType>::pid_element_const_iterator>;

template<typename MeshType>
using ext_elements_t = boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                                    typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> >::const_iterator,
                                    typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> >::const_iterator,
                                    boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::element_type const> > >
                                    >;

template<typename MeshType>
using boundaryelements_t =  boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                                         typename MeshTraits<MeshType>::location_element_const_iterator,
                                         typename MeshTraits<MeshType>::location_element_const_iterator>;

template<typename MeshType>
using internalelements_t = boundaryelements_t<MeshType>;

template<typename MeshType>
using markedelements_t = boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                                      typename MeshTraits<MeshType>::marker_element_const_iterator,
                                      typename MeshTraits<MeshType>::marker_element_const_iterator>;

template<typename MeshType>
using marked2elements_t = boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                                       typename MeshTraits<MeshType>::marker2_element_const_iterator,
                                       typename MeshTraits<MeshType>::marker2_element_const_iterator>;

template<typename MeshType>
using marked3elements_t =  boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                                        typename MeshTraits<MeshType>::marker3_element_const_iterator,
                                        typename MeshTraits<MeshType>::marker3_element_const_iterator> ;


template<typename MeshType>
using faces_t =  boost::tuple<mpl::size_t<MESH_FACES>,
                              typename MeshTraits<MeshType>::pid_face_const_iterator,
                              typename MeshTraits<MeshType>::pid_face_const_iterator>;

template<typename MeshType>
using ext_faces_t = boost::tuple<mpl::size_t<MESH_FACES>,
                                 typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> >::const_iterator,
                                 typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> >::const_iterator,
                                 boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> > >
                                 >;

template<typename MeshType>
using idfaces_t =  boost::tuple<mpl::size_t<MESH_FACES>,
                                typename MeshTraits<MeshType>::face_const_iterator,
                                typename MeshTraits<MeshType>::face_const_iterator>;

template<typename MeshType>
using boundaryfaces_t =  boost::tuple<mpl::size_t<MESH_FACES>,
                                      typename MeshTraits<MeshType>::location_face_const_iterator,
                                      typename MeshTraits<MeshType>::location_face_const_iterator>;

template<typename MeshType>
using internalfaces_t = boost::tuple<mpl::size_t<MESH_FACES>,
                                     typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> >::const_iterator,
                                     typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> >::const_iterator,
                                     boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::face_type const> > >
                                     >;

template<typename MeshType>
using interprocessfaces_t =  boost::tuple<mpl::size_t<MESH_FACES>,
                                          typename MeshTraits<MeshType>::interprocess_face_const_iterator,
                                          typename MeshTraits<MeshType>::interprocess_face_const_iterator>;


template<typename MeshType>
using markedfaces_t = boost::tuple<mpl::size_t<MESH_FACES>,
                                   typename MeshTraits<MeshType>::marker_face_const_iterator,
                                   typename MeshTraits<MeshType>::marker_face_const_iterator>;

template<typename MeshType>
using marked2faces_t = boost::tuple<mpl::size_t<MESH_FACES>,
                                    typename MeshTraits<MeshType>::marker2_face_const_iterator,
                                    typename MeshTraits<MeshType>::marker2_face_const_iterator>;

template<typename MeshType>
using marked3faces_t =  boost::tuple<mpl::size_t<MESH_FACES>,
                                     typename MeshTraits<MeshType>::marker3_face_const_iterator,
                                     typename MeshTraits<MeshType>::marker3_face_const_iterator> ;

template<typename MeshType>
using edges_t =  boost::tuple<mpl::size_t<MESH_EDGES>,
                               typename MeshTraits<MeshType>::edge_const_iterator,
                               typename MeshTraits<MeshType>::edge_const_iterator>;
template<typename MeshType>
using pid_edges_t =  boost::tuple<mpl::size_t<MESH_EDGES>,
                                typename MeshTraits<MeshType>::pid_edge_const_iterator,
                                typename MeshTraits<MeshType>::pid_edge_const_iterator>;

template<typename MeshType>
using ext_edges_t =  boost::tuple<mpl::size_t<MESH_EDGES>,
                              typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::edge_type const> >::const_iterator,
                              typename std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::edge_type const> >::const_iterator,
                              boost::shared_ptr<std::vector<boost::reference_wrapper<typename MeshTraits<MeshType>::edge_type const> > >
                              >;

template<typename MeshType>
using markededges_t = boost::tuple<mpl::size_t<MESH_EDGES>,
                                   typename MeshTraits<MeshType>::marker_edge_const_iterator,
                                   typename MeshTraits<MeshType>::marker_edge_const_iterator>;

template<typename MeshType>
using boundaryedges_t = boost::tuple<mpl::size_t<MESH_EDGES>,
                                     typename MeshTraits<MeshType>::location_edge_const_iterator,
                                     typename MeshTraits<MeshType>::location_edge_const_iterator>;

template<typename MeshType>
using internaledges_t = boundaryedges_t<MeshType>;

template<typename MeshType>
using points_t =  boost::tuple<mpl::size_t<MESH_POINTS>,
                               typename MeshTraits<MeshType>::point_const_iterator,
                               typename MeshTraits<MeshType>::point_const_iterator>;

template<typename MeshType>
using markedpoints_t = boost::tuple<mpl::size_t<MESH_POINTS>,
                                    typename MeshTraits<MeshType>::marker_point_const_iterator,
                                    typename MeshTraits<MeshType>::marker_point_const_iterator>;

template<typename MeshType>
using boundarypoints_t = boost::tuple<mpl::size_t<MESH_POINTS>,
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
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
            typename MeshTraits<MeshType>::element_const_iterator,
            typename MeshTraits<MeshType>::element_const_iterator> type;
};

template<typename MeshType>
struct markedelements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
            typename MeshTraits<MeshType>::marker_element_const_iterator,
            typename MeshTraits<MeshType>::marker_element_const_iterator> type;
};

template<typename MeshType>
struct marked2elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
            typename MeshTraits<MeshType>::marker2_element_const_iterator,
            typename MeshTraits<MeshType>::marker2_element_const_iterator> type;
};

template<typename MeshType>
struct marked3elements
{
    static const uint16_type nDim = MeshType::nDim;
    typedef boost::tuple<mpl::size_t<MESH_ELEMENTS>,
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::allelements( mesh, is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over elements with pid \p flag
 */
template<typename MeshType>
pid_elements_t<MeshType>
elements( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::elements( mesh, rank( mesh ), is_ptr_or_shared_ptr() );
    //return elements( mesh, flag, is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::boundaryelements( mesh, entity_min_dim, entity_max_dim, rank( mesh ), is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::internalelements( mesh, rank( mesh ), is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedelements( mesh, mesh->markerName( flag ), rank( mesh ), is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedelements( mesh, mesh->markerName( flag ), rank( mesh ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
markedelements_t<MeshType>
markedelements( MeshType const& mesh, boost::any const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;

    flag_type theflag = mesh->markerId( flag );
    VLOG(2) << "[markedelements] flag: " << theflag << "\n";
    return Feel::detail::markedelements( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<markedelements_t<MeshType>>
markedelements( MeshType const& mesh, std::initializer_list<boost::any> const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<markedelements_t<MeshType>> list_elements;
    for ( auto const& it : flag )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedelements] flag: " << theflag << "\n";
        list_elements.push_back( Feel::detail::markedelements( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() ) );
    }
    return list_elements;
}

template<typename MeshType>
std::list<markedelements_t<MeshType>>
markedelements( MeshType const& mesh, std::list<std::string> const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<markedelements_t<MeshType>> list_elements;
    for ( auto const& it : flag )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedelements] flag: " << theflag << "\n";
        list_elements.push_back( Feel::detail::markedelements( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() ) );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked2elements( mesh, flag, rank( mesh ), is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked2elements( mesh, mesh->markerName( flag ), rank( mesh ), is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked2elements( mesh, mesh->markerName( flag ), rank( mesh ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
marked2elements_t<MeshType>
marked2elements( MeshType const& mesh, boost::any const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;

    flag_type theflag = mesh->markerId( flag );
    VLOG(2) << "[markedelements] flag: " << theflag << "\n";
    return Feel::detail::marked2elements( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<marked2elements_t<MeshType>>
marked2elements( MeshType const& mesh, std::initializer_list<boost::any> const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<marked2elements_t<MeshType>> list_elements;
    for ( auto const& it : flag )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedelements] flag: " << theflag << "\n";
        list_elements.push_back( Feel::detail::marked2elements( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() ) );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked3elements( mesh, flag, rank( mesh ), is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked3elements( mesh, mesh->markerName( flag ), rank( mesh ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<marked3elements_t<MeshType>>
marked3elements( MeshType const& mesh, std::initializer_list<boost::any> const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<marked3elements_t<MeshType>> list_elements;
    for ( auto const& it : flag )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedelements] flag: " << theflag << "\n";
        list_elements.push_back( Feel::detail::marked3elements( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() ) );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::idedelements( mesh, flag, is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::faces( mesh, rank( mesh ), is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::idedfaces( mesh, id, is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedfaces( mesh, rank( mesh ), is_ptr_or_shared_ptr() );
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
             std::string const&__marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedfaces( mesh, mesh->markerName( __marker ), rank( mesh ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
markedfaces_t<MeshType>
markedfaces( MeshType const& mesh,
             const char*__marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedfaces( mesh, mesh->markerName( __marker ), rank( mesh ), is_ptr_or_shared_ptr() );
}
template<typename MeshType>
std::list<markedfaces_t<MeshType>>
markedfaces( MeshType const& mesh,
             std::list<std::string> const& __markers )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<markedfaces_t<MeshType>> list_faces;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_faces.push_back( Feel::detail::markedfaces( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() ) );
    }
    return list_faces;
}

template<typename MeshType>
markedfaces_t<MeshType>
markedfaces( MeshType const& mesh,
             boost::any const& __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    flag_type theflag = mesh->markerId( __marker );
    VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
    return Feel::detail::markedfaces( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() );

}

template<typename MeshType>
std::list<markedfaces_t<MeshType>>
markedfaces( MeshType const& mesh,
             std::initializer_list<boost::any> __markers )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<markedfaces_t<MeshType>> list_faces;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_faces.push_back( Feel::detail::markedfaces( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() ) );
    }
    return list_faces;
}


template<typename MeshType>
boost::tuple<mpl::size_t<MESH_FACES>,
      typename MeshTraits<MeshType>::marker2_face_const_iterator,
      typename MeshTraits<MeshType>::marker2_face_const_iterator>
marked2faces( MeshType const& mesh,
              flag_type __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked2faces( mesh, __marker, rank( mesh ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                       typename MeshTraits<MeshType>::marker2_face_const_iterator,
                       typename MeshTraits<MeshType>::marker2_face_const_iterator> >
marked2faces( MeshType const& mesh,
             std::initializer_list<boost::any> __markers )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<boost::tuple<mpl::size_t<MESH_FACES>,
                           typename MeshTraits<MeshType>::marker2_face_const_iterator,
                           typename MeshTraits<MeshType>::marker2_face_const_iterator> > list_faces;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_faces.push_back( Feel::detail::marked2faces( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() ) );
    }
    return list_faces;
}

template<typename MeshType>
marked3faces_t<MeshType>
marked3faces( MeshType const& mesh,
              flag_type __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked3faces( mesh, __marker, rank( mesh ), is_ptr_or_shared_ptr() );
}


template<typename MeshType>
marked2faces_t<MeshType>
marked2faces( MeshType const& mesh,
              std::string const& __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked2faces( mesh, mesh->markerName( __marker ), rank( mesh ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
marked3faces_t<MeshType>
marked3faces( MeshType const& mesh,
              std::string const& __marker )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::marked3faces( mesh, mesh->markerName( __marker ), rank( mesh ), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<marked3faces_t<MeshType>>
marked3faces( MeshType const& mesh,
             std::initializer_list<boost::any> __markers )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<marked3faces_t<MeshType>> list_faces;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_faces.push_back( Feel::detail::marked3faces( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() ) );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::boundaryfaces( mesh, rank( mesh ), is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::internalfaces( mesh, rank( mesh ), is_ptr_or_shared_ptr() );
}

/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all interprocess faces of the
 * mesh belonging to process \p __pid
 */
template<typename MeshType>
interprocessfaces_t<MeshType>
interprocessfaces( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::interprocessfaces( mesh, invalid_rank_type_value, is_ptr_or_shared_ptr() );
}


/**
 *
 * \ingroup MeshIterators
 * \return a pair of iterators to iterate over all interprocess faces of the
 * mesh belonging to process \p __pid
 */
template<typename MeshType>
interprocessfaces_t<MeshType>
interprocessfaces( MeshType const& mesh, rank_type neighbor_pid )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::interprocessfaces( mesh, neighbor_pid, is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::edges( mesh, rank( mesh ), is_ptr_or_shared_ptr() );
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
    return Feel::detail::markededges( mesh, __marker, rank( mesh ), is_ptr_or_shared_ptr<MeshType>() );
}

template<typename MeshType>
markededges_t<MeshType>
markededges( MeshType const& mesh,
             std::string const& __marker,
             typename std::enable_if<is_3d<MeshType>::value>::type* = nullptr )
{
    return Feel::detail::markededges( mesh,
                                      mesh->markerName( __marker ),
                                      rank( mesh ),
                                      is_ptr_or_shared_ptr<MeshType>() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<markededges_t<MeshType>> list_edges;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_edges.push_back( Feel::detail::markededges( mesh, theflag, rank( mesh ), is_ptr_or_shared_ptr() ) );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::boundaryedges( mesh, is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::internaledges( mesh, is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::points( mesh, is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedpoints( mesh, flag, is_ptr_or_shared_ptr() );
}

template<typename MeshType>
markedpoints_t<MeshType>
markedpoints( MeshType const& mesh, std::string const& flag )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::markedpoints( mesh, mesh->markerName(flag), is_ptr_or_shared_ptr() );
}

template<typename MeshType>
std::list<markedpoints_t<MeshType>>
markedpoints( MeshType const& mesh, std::list<std::string> const& __markers )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    std::list<markedpoints_t<MeshType>> list_points;
    for ( auto const& it : __markers )
    {
        flag_type theflag = mesh->markerId( it );
        VLOG(2) << "[markedfaces] flag: " << theflag << "\n";
        list_points.push_back( Feel::detail::markedpoints( mesh, theflag/*, rank( mesh )*/, is_ptr_or_shared_ptr() ) );
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
boundarypoints( MeshType const& mesh )
{
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;
    return Feel::detail::boundarypoints( mesh, is_ptr_or_shared_ptr() );
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
    typedef typename mpl::or_<is_shared_ptr<MeshType>, boost::is_pointer<MeshType> >::type is_ptr_or_shared_ptr;

    return Feel::detail::internalpoints( mesh, is_ptr_or_shared_ptr() );
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
    size_type d = 0;
    std::for_each( its.begin(), its.end(),
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
ext_elements_t<MeshType>
elements( MeshType const& mesh, EntityProcessType entity )
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
            myelts->push_back(boost::cref(eltOffProc));

            eltGhostDone.insert( eltOffProc.id() );
        }
    }

    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts );
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

    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
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

    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
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

    return boost::make_tuple( mpl::size_t<MESH_ELEMENTS>(),
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

    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
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

    return boost::make_tuple( mpl::size_t<MESH_FACES>(),
                              myelts->begin(),
                              myelts->end(),
                              myelts );

}

/**
 * this function takes two sets of iterators and generate a data structure that
 * will hold the corresponding reference of mesh elenents. This data structure
 * can then be used to iterate over meshes like other mesh filters.
 *
 * \param it1 (begin,end) iterators to concatenate \param it2 (begin,end)
 * iterators to concatenate \return a data structure that holds the merge
 * between two sets of pair of iterators
 */
template<typename IteratorType>
boost::tuple<mpl::size_t<boost::tuples::element<0,IteratorType>::type::value>,
             typename std::vector<boost::reference_wrapper<typename boost::tuples::element<1,IteratorType>::type::value_type const> >::const_iterator,
             typename std::vector<boost::reference_wrapper<typename boost::tuples::element<1,IteratorType>::type::value_type const> >::const_iterator,
             boost::shared_ptr<std::vector<boost::reference_wrapper<typename boost::tuples::element<1,IteratorType>::type::value_type const> > >
             >
concatenate( IteratorType it1, IteratorType it2 )
{
    typedef std::vector<boost::reference_wrapper<typename boost::tuples::element<1,IteratorType>::type::value_type  const> > cont_range_type;
    boost::shared_ptr<cont_range_type> myelts( new cont_range_type );

    auto append = [&myelts]( typename boost::tuples::element<1,IteratorType>::type::value_type  const& e ) { myelts->push_back( boost::cref(e) ); };
    std::for_each( begin( it1 ), end( it1 ), append );
    std::for_each( begin( it2 ), end( it2 ), append );

    return boost::make_tuple( mpl::size_t<boost::tuples::element<0,IteratorType>::type::value>(),
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

    return boost::make_tuple( mpl::size_t<MESH_EDGES>(),
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
