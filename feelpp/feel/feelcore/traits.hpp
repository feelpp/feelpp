/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Thibaut Metivet <thibaut.metivet@inria.fr>
       Date: 2005-07-28

  Copyright (C) 2007,2009 Université de Grenoble 1
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-28
 */
#ifndef FEELPP_TRAITS_HPP
#define FEELPP_TRAITS_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/typetraits.hpp>

#include <Eigen/Core>

namespace Feel
{
template <typename... T>
using strongest_numeric_type = std::common_type_t<T...>;

template <class T>
struct is_shared_ptr : std::bool_constant<false> {};

template <class T>
struct is_shared_ptr<std::shared_ptr<T> > : std::bool_constant<true> {};

template <class T>
struct remove_shared_ptr
{
    typedef T type;
};

template <class T>
struct remove_shared_ptr<std::shared_ptr<T> >
{
    typedef T type;
};

template <class T>
struct remove_shared_ptr<std::unique_ptr<T>>
{
    typedef T type;
};

template<class T>
constexpr bool is_shared_ptr_v = is_shared_ptr<T>::value;

template<typename T>
struct is_ptr_or_shared_ptr : std::bool_constant<is_shared_ptr_v<T>||std::is_pointer_v<T>>  {};

template<typename T>
using remove_shared_ptr_type = typename remove_shared_ptr<T>::type;

template <typename T>
using remove_shared_ptr_t = typename remove_shared_ptr<T>::type;

template<typename T>
using decay_type = std::decay_t<remove_shared_ptr_type<std::decay_t<T>>>;

template<typename T>
decltype(auto) remove_shared_ptr_f( T&& e )
{
    return hana::if_( hana::bool_<is_shared_ptr_v<T>>{},
                     []( auto&& x ) { return *x; },
                     []( auto&& x ) { return x; } )( std::forward<T>(e) );

}

template< typename T >
struct type_identity
{
    using type = T;
};
template< typename T>
using type_identity_t = typename type_identity<T>::type;


template <typename T, typename = void>
struct is_iterable : std::false_type {};
template <typename T>
struct is_iterable<T, std::void_t<decltype(std::declval<T>().begin()),decltype(std::declval<T>().end())>>
    : std::true_type {};
template <typename T>
constexpr bool is_iterable_v = is_iterable<T>::value;


template <typename T, typename V, typename = void>
struct is_iterable_of : std::false_type {};
template <typename T, typename V>
struct is_iterable_of< T, V, std::void_t<
        decltype(std::declval<T>().begin()), decltype(std::declval<T>().end()), 
        std::is_same<decltype(*std::declval<T>().begin()), V>
    > > : std::true_type {};
template <typename T, typename V>
constexpr bool is_iterable_of_v = is_iterable_of<T, V>::value;

template <class T>
struct is_std_vector : std::bool_constant<false> {};
template <class T>
struct is_std_vector<std::vector<T> > : std::bool_constant<true> {};
template <class T>
using is_std_vector_t = is_std_vector<T>;
template <class T>
inline constexpr  bool is_std_vector_v = is_std_vector_t<T>::value;

template <class T>
struct remove_std_vector
{
    typedef T type;
};
template <class T>
struct remove_std_vector<std::vector<T> >
{
    typedef T type;
};
template<typename T>
using remove_std_vector_t = typename remove_std_vector<T>::type;


template <typename T>
inline constexpr bool is_eigen_matrix_v = std::is_base_of_v<Eigen::MatrixBase<std::decay_t<T>>,std::decay_t<T> >;

/**
 * @brief Utility class to get the element type contained in ObjectType
 * 
 * @tparam ObjectType class
 */
template <typename ObjectType, typename = void>
struct element_type_helper;  // Primary template left undefined on purpose.

template <typename ObjectType>
using element_t = typename element_type_helper<ObjectType>::type;

template <typename ObjectType>
using element_ptr_t = typename element_type_helper<ObjectType>::ptrtype;

// Trait to get the value type for a mesh
template <typename Type, typename Enable = void>
struct value_type_trait
{
}; // The default case is empty.

// Helper type alias to extract the value type
template <typename Type>
using value_t = typename value_type_trait<Type>::type;

} // namespace Feel
#endif
