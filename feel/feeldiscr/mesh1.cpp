/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-07-07

  Copyright (C) 2009,2010 Université Joseph Fourier (Grenoble I)

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
   \file mesh1.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-07-07
 */
#define FEELPP_MESH_IMPL_NOEXTERN 1
#include <feel/feeldiscr/meshimpl.hpp>

namespace Feel
{
//
// Explicit instatiations
//
#if defined( FEELPP_INSTANTIATION_MODE )


BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_SIMPLEX_OP, 3, ( DIMS1, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS1 ), RDIMS1 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_HYPERCUBE_OP, 3, ( DIMS1, BOOST_PP_LIST_FIRST_N( FEELPP_MESH_MAX_ORDER, ORDERS1 ), RDIMS1 ) )



#endif // FEELPP_INSTANTIATION_MODE

}
