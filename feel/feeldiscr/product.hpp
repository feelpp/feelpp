/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-03-23

  Copyright (C) 2014 Feel++ Consortium

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
/**
   \file product.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-03-23
 */
#if !defined(FEELPP_PRODUCT_HPP)
#define FEELPP_PRODUCT_HPP 1

#include <boost/fusion/view/joint_view.hpp>

namespace Feel {

namespace detail {

namespace product {

template<typename T>
struct clean
{
    typedef typename boost::remove_reference<typename boost::remove_const<typename boost::remove_pointer<T>::type>::type >::type type;
};
template<typename Space>
struct GetMesh
{
    typedef typename clean<typename Space::element_type>::type::mesh_type type;
};
template<typename Space>
struct GetBasis
{
    typedef typename clean<typename Space::element_type>::type::basis_type type;
};
template<typename Space>
struct GetPeriodicity
{
    typedef typename clean<typename Space::element_type>::type::periodicity_type type;
};
template<typename Space>
struct GetMortar
{
    typedef typename clean<typename Space::element_type>::type::mortar_0_type type;
};

template<typename... SpaceList>
struct Product
{
    typedef typename mpl::transform<fusion::vector<SpaceList...>, GetMesh<mpl::_1>, mpl::back_inserter<meshes<> > >::type mesh_type;
    typedef typename mpl::transform<fusion::vector<SpaceList...>, GetBasis<mpl::_1>, mpl::back_inserter<bases<> > >::type basis_type;
    typedef typename mpl::transform<fusion::vector<SpaceList...>, GetPeriodicity<mpl::_1>, mpl::back_inserter<Periodicity<> > >::type periodicity_type;
    typedef typename mpl::transform<fusion::vector<SpaceList...>, GetMortar<mpl::_1>, mpl::back_inserter<mortars<> > >::type mortar_type;

    typedef FunctionSpace<typename mpl::front<mesh_type>::type,basis_type,double,periodicity_type,mortar_type> type;
    typedef boost::shared_ptr<type> ptrtype;
};
} }

template<typename... SpaceList>
typename Feel::detail::product::Product<SpaceList...>::ptrtype
product( SpaceList... spaces )
{
    typedef typename Feel::detail::product::Product<SpaceList...> product_type;
    typedef typename product_type::type space_type;
    typedef typename product_type::ptrtype space_ptrtype;
    std::cout << "sizeof " << sizeof...(SpaceList) << "\n";
    space_ptrtype Rh( space_type::NewFromList( spaces... ) );
    return Rh;
}

} // Feel++
#endif
