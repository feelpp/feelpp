/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-07-07

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file mesh3.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-07-07
 */
#include <feel/feeldiscr/meshimpl.hpp>

namespace Feel
{
//
// Explicit instatiations
//
#if defined( FEEL_INSTANTIATION_MODE )


template class Mesh<SimplexProduct<3, 1, 3> >;

# define DIMS BOOST_PP_TUPLE_TO_LIST(1,(3))
//# define RDIMS BOOST_PP_TUPLE_TO_LIST(2,(1,2))
# define RDIMS BOOST_PP_TUPLE_TO_LIST(1,(3))
# define ORDERS BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))

# define FACTORY(LDIM,LORDER,RDIM) template class Mesh<Simplex<LDIM, LORDER, RDIM> >;

# define FACTORY_OP(_, GDO) FACTORY GDO

// only up to 4 for mesh data structure nit supported for higher order in Gmsh
BOOST_PP_LIST_FOR_EACH_PRODUCT(FACTORY_OP, 3, (DIMS, BOOST_PP_LIST_FIRST_N(BOOST_PP_MIN(4,FEEL_MESH_MAX_ORDER), ORDERS), RDIMS))




#endif // FEEL_INSTANTIATION_MODE

}
