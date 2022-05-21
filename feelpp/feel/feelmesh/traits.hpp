/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomm@zion>
       Date: 2005-08-27

  Copyright (C) 2011-2022 Feel++ Consortium
  Copyright (C) 2011-2022 University of Strasbourg
  Copyright (C) 2005,2006 EPFL

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
   \file traits.hpp
   \author Christophe Prud'homme <prudhomm@zion>
   \date 2005-08-27
 */
#ifndef __FEELPP_MESH_TRAITS_HPP
#define __FEELPP_MESH_TRAITS_HPP 1

#include <feel/feelcore/traits.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <boost/phoenix/stl/algorithm/detail/is_std_list.hpp>
#include <feel/feelmesh/enums.hpp>

namespace Feel
{
/**
 * \class MeshTraits
 * \ingroup Mesh
 * \brief Traits for meshes
 *
 * @author Christophe Prud'homme
 */
template<typename MeshType>
struct MeshTraits
{
    /** @name Typedefs
     */
    //@{

    typedef MeshTraits<MeshType> self_type;
    typedef typename boost::remove_pointer<typename remove_shared_ptr<MeshType>::type >::type mesh_type;

    typedef typename mesh_type::shape_type element_shape_type;

    typedef typename mesh_type::element_type element_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::point_type point_type;

    // element iterators
    typedef typename mesh_type::element_iterator element_iterator;
    typedef typename mesh_type::element_const_iterator element_const_iterator;

    typedef typename mesh_type::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename mesh_type::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef typename mesh_type::element_reference_wrapper_iterator element_reference_wrapper_iterator;
    typedef typename mesh_type::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;

    // face iterators
    typedef typename mesh_type::face_iterator face_iterator;
    typedef typename mesh_type::face_const_iterator face_const_iterator;

    typedef typename mesh_type::faces_reference_wrapper_type faces_reference_wrapper_type;
    typedef typename mesh_type::faces_reference_wrapper_ptrtype faces_reference_wrapper_ptrtype;
    typedef typename mesh_type::face_reference_wrapper_iterator face_reference_wrapper_iterator;
    typedef typename mesh_type::face_reference_wrapper_const_iterator face_reference_wrapper_const_iterator;

    // edge iterators
    typedef typename mesh_type::edge_iterator edge_iterator;
    typedef typename mesh_type::edge_const_iterator edge_const_iterator;

    typedef typename mesh_type::edges_reference_wrapper_type edges_reference_wrapper_type;
    typedef typename mesh_type::edges_reference_wrapper_ptrtype edges_reference_wrapper_ptrtype;
    typedef typename mesh_type::edge_reference_wrapper_iterator edge_reference_wrapper_iterator;
    typedef typename mesh_type::edge_reference_wrapper_const_iterator edge_reference_wrapper_const_iterator;

    // point iterators
    typedef typename mesh_type::point_iterator point_iterator;
    typedef typename mesh_type::point_const_iterator point_const_iterator;

    typedef typename mesh_type::points_reference_wrapper_type points_reference_wrapper_type;
    typedef typename mesh_type::points_reference_wrapper_ptrtype points_reference_wrapper_ptrtype;
    typedef typename mesh_type::point_reference_wrapper_iterator point_reference_wrapper_iterator;
    typedef typename mesh_type::point_reference_wrapper_const_iterator point_reference_wrapper_const_iterator;

    //@}
};

template<typename T>
struct dimension_t : mpl::int_<decay_type<T>::nDim> {};
template<typename T>
constexpr uint16_type dimension_v = dimension_t<T>::value;
/**
 * @return topological dimension of a type T
 */
template<typename T>
inline constexpr int dimension( T const& t )
{
  return dimension_v<T>;
}

template<typename T>
struct real_dimension_t : mpl::int_<decay_type<T>::nRealDim> {};
template<typename T>
constexpr uint16_type real_dimension_v = real_dimension_t<T>::value;

/**
 * @return real dimension of a type T
 */
template<typename T>
inline constexpr int real_dimension( T const& t )
{
  return real_dimension_v<T>;
}



template<typename T>
struct is_3d : mpl::bool_<decay_type<T>::nDim == 3 /*|| decay_type<T>::nRealDim ==3*/> {};
template<typename T>
struct is_2d : mpl::bool_<decay_type<T>::nDim == 2 /*|| decay_type<T>::nRealDim ==2*/> {};
template<typename T>
struct is_1d : mpl::bool_<decay_type<T>::nDim == 1 /*|| decay_type<T>::nRealDim ==1*/> {};
template<typename T>
struct is_0d : mpl::bool_<decay_type<T>::nDim == 0 /*|| decay_type<T>::nRealDim ==0*/> {};

template<typename T>
struct is_3d_real : mpl::bool_< decay_type<T>::nRealDim ==3 > {};
template<typename T>
struct is_2d_real : mpl::bool_< decay_type<T>::nRealDim ==2 > {};
template<typename T>
struct is_1d_real : mpl::bool_< decay_type<T>::nRealDim ==1 > {};
template<typename T>
struct is_0d_real : mpl::bool_< decay_type<T>::nRealDim ==0 > {};


template<typename T>
struct is_topological_face : mpl::bool_<(decay_type<T>::nDim==decay_type<T>::nRealDim-1)> {};
template<typename T>
struct is_face : mpl::bool_<(decay_type<T>::nDim == 2 && decay_type<T>::nRealDim == 3)> {};
template<typename T>
struct is_edge : mpl::bool_<(decay_type<T>::nDim==1 && decay_type<T>::nRealDim == 3)||(decay_type<T>::nDim==1 && decay_type<T>::nRealDim == 2)> {};
template<typename T>
struct is_point : mpl::bool_<(decay_type<T>::nDim == 0)> {};

template<typename T>
struct is_convex : std::is_convertible<T,ConvexBase>::type {};

template<typename T>
struct is_simplex : std::is_base_of<SimplexBase, T>::type {};
template<typename T>
inline constexpr bool is_simplex_v = is_simplex<T>::value;

template<typename T>
struct is_triangle : mpl::and_<is_simplex<T>,is_2d<T>> {};
template<typename T>
struct is_tetrahedron : mpl::and_<is_simplex<T>,is_3d<T>> {};

template<typename T>
struct is_hypercube : std::is_base_of<HypercubeBase, T>::type {};
template <typename T>
inline constexpr bool is_hypercube_v = is_hypercube<T>::value;


template<typename T>
struct is_square : mpl::and_<is_hypercube<T>,is_2d<T>> {};
template<typename T>
struct is_cube : mpl::and_<is_hypercube<T>,is_3d<T>> {};

template<typename T>
struct is_segment : mpl::and_<is_convex<T>,is_1d<T>> {};

/**
 * Checks whether T is a GeoElement<n>D type. 
 * Provides the member constant value that is equal to true, if T is the type GeoElement<n>D. 
 * Otherwise, value is equal to false.
 */
template <typename T> struct is_geoelement: std::false_type {};
/**
 * Helper variable template
 * \return true is T is a GeoElement<n>D at compilation time, false otherwise
 */
template <typename... T>
inline constexpr bool is_geoelement_v = is_geoelement<T...>::value;





/**
 * Filters
 */

template <typename T> struct is_filter: std::false_type {};

template <typename... T> struct is_filter<boost::tuple<T...>>: std::true_type {};
template <typename... T>
constexpr bool is_filter_v = is_filter<T...>::value;

/**
 * a RangeType can be one or more filter/range objects of the same type, we
 * extract the underlying type by first casting everything to a list and then
 * retrieving consistently the type.
 */
template <typename RangeType>
using range_t = typename mpl::if_< boost::is_std_list<RangeType>,
                                   mpl::identity<RangeType>,
                                   mpl::identity<std::list<RangeType> > >::type::type::value_type;

template <typename RangeType>
using entity_range_t = typename  boost::unwrap_reference<typename boost::tuples::template element<1,range_t<RangeType> >::type::value_type>::type;

template<typename MeshType>
using elements_reference_wrapper_t = boost::tuple<mpl::size_t<MESH_ELEMENTS>,
                                                  typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
                                                  typename MeshTraits<MeshType>::element_reference_wrapper_const_iterator,
                                                  typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype>;

template<typename MeshType>
using allelements_t = elements_reference_wrapper_t<MeshType>;

template<typename MeshType>
using elements_pid_t = elements_reference_wrapper_t<MeshType>;

template<typename MeshType>
using idelements_t = elements_reference_wrapper_t<MeshType>;

template<typename MeshType>
using ext_elements_t = elements_reference_wrapper_t<MeshType>;
template <typename MeshType>
using range_element_t = typename MeshTraits<MeshType>::elements_reference_wrapper_type;
template <typename MeshType>
using range_element_ptr_t = typename MeshTraits<MeshType>::elements_reference_wrapper_ptrtype;

template <typename MeshType>
using boundaryelements_t = elements_reference_wrapper_t<MeshType>;
template<typename MeshType>
using internalelements_t = elements_reference_wrapper_t<MeshType>;

template<typename MeshType>
using markedelements_t = elements_reference_wrapper_t<MeshType>;
template<typename MeshType>
using marked2elements_t = elements_reference_wrapper_t<MeshType>;
template<typename MeshType>
using marked3elements_t = elements_reference_wrapper_t<MeshType>;


template<typename MeshType>
using faces_reference_wrapper_t = boost::tuple<mpl::size_t<MESH_FACES>,
                                                  typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
                                                  typename MeshTraits<MeshType>::face_reference_wrapper_const_iterator,
                                                  typename MeshTraits<MeshType>::faces_reference_wrapper_ptrtype>;

template<typename MeshType>
using allfaces_t =  faces_reference_wrapper_t<MeshType>;

template<typename MeshType>
using idfaces_t =  faces_reference_wrapper_t<MeshType>;

template<typename MeshType>
using faces_pid_t = faces_reference_wrapper_t<MeshType>;

template<typename MeshType>
using ext_faces_t = faces_reference_wrapper_t<MeshType>;

template<typename MeshType>
using boundaryfaces_t = faces_reference_wrapper_t<MeshType>;
template<typename MeshType>
using internalfaces_t = faces_reference_wrapper_t<MeshType>;
template<typename MeshType>
using interprocessfaces_t =  faces_reference_wrapper_t<MeshType>;

template<typename MeshType>
using markedfaces_t = faces_reference_wrapper_t<MeshType>;
template<typename MeshType>
using marked2faces_t = faces_reference_wrapper_t<MeshType>;
template<typename MeshType>
using marked3faces_t = faces_reference_wrapper_t<MeshType>;


template<typename MeshType>
using edges_reference_wrapper_t = boost::tuple<mpl::size_t<MESH_EDGES>,
                                                  typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
                                                  typename MeshTraits<MeshType>::edge_reference_wrapper_const_iterator,
                                                  typename MeshTraits<MeshType>::edges_reference_wrapper_ptrtype>;

template<typename MeshType>
using alledges_t = edges_reference_wrapper_t<MeshType>;

template<typename MeshType>
using pid_edges_t = edges_reference_wrapper_t<MeshType>;

template<typename MeshType>
using ext_edges_t =  edges_reference_wrapper_t<MeshType>;

template<typename MeshType>
using markededges_t = edges_reference_wrapper_t<MeshType>;
template<typename MeshType>
using boundaryedges_t = edges_reference_wrapper_t<MeshType>;
template<typename MeshType>
using internaledges_t = edges_reference_wrapper_t<MeshType>;

template<typename MeshType>
using points_reference_wrapper_t = boost::tuple<mpl::size_t<MESH_POINTS>,
                                    typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
                                    typename MeshTraits<MeshType>::point_reference_wrapper_const_iterator,
                                    typename MeshTraits<MeshType>::points_reference_wrapper_ptrtype>;
template<typename MeshType>
using allpoints_t = points_reference_wrapper_t<MeshType>;
template<typename MeshType>
using points_pid_t = points_reference_wrapper_t<MeshType>;
template<typename MeshType>
using markedpoints_t = points_reference_wrapper_t<MeshType>;
template<typename MeshType>
using boundarypoints_t = points_reference_wrapper_t<MeshType>;
template<typename MeshType>
using internalpoints_t = points_reference_wrapper_t<MeshType>;

namespace detail
{
template <typename RangeType>
struct submeshrangetype
{
    typedef typename mpl::if_< boost::is_std_list<RangeType>,
                               mpl::identity<RangeType>,
                               mpl::identity<std::list<RangeType> > >::type::type::value_type type;
};
}
template<typename IteratorRangeT>
using submeshrange_t = typename Feel::detail::submeshrangetype<IteratorRangeT>::type;


template<typename IteratorType>
using filter_enum_t = typename boost::tuples::element<0,typename meta::remove_all<IteratorType>::type >::type;
template<typename IteratorType>
using filter_iterator_t = typename boost::tuples::element<1,typename meta::remove_all<IteratorType>::type >::type;
template<typename IteratorType>
using filter_entity_t = typename boost::unwrap_reference<typename filter_iterator_t<IteratorType>::value_type>::type;
template<typename IteratorType>
using ext_entities_from_iterator_t = boost::tuple<filter_enum_t<IteratorType>,
                                                  typename std::vector<boost::reference_wrapper<filter_entity_t<IteratorType> const> >::const_iterator,
                                                  typename std::vector<boost::reference_wrapper<filter_entity_t<IteratorType> const> >::const_iterator,
                                                  std::shared_ptr<std::vector<boost::reference_wrapper<filter_entity_t<IteratorType> const> > >
                                                  >;

template<typename MeshType>
using elements_wrapper_t = typename MeshTraits<MeshType>::elements_reference_wrapper_type;
template<typename MeshType>
using elements_wrapper_ptr_t = typename MeshTraits<MeshType>::elements_reference_wrapper_ptr_type;



} // Feel
#endif /* __Traits_H */
