/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013 Feel++ Consortium

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
   \file mesh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_DETAIL_MESH_HPP)
#define FEELPP_DETAIL_MESH_HPP 1

namespace Feel {

/// \cond DETAIL
namespace detail
{
template<typename Args, typename Tag=tag::mesh>
struct mesh
{
    typedef typename boost::remove_pointer<
        typename boost::remove_const<
            typename boost::remove_reference<
                typename parameter::binding<Args, Tag>::type
                >::type
            >::type
    >::type _type;
typedef typename mpl::if_<is_shared_ptr<_type>,
                          mpl::identity<typename _type::element_type>,
                          mpl::identity<_type> >::type::type type;
typedef boost::shared_ptr<type> ptrtype;
};

} }

#endif /* FEELPP_DETAIL_MESH_HPP */
