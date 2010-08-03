/* -*- mode: c++ -*-

  This file is part of the Life library

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
#include <life/lifediscr/meshimpl.hpp>

namespace Life
{
//
// Explicit instatiations
//
#if defined( LIFE_INSTANTIATION_MODE )

template class Mesh<Simplex<3, 1, 3> >;
template class Mesh<SimplexProduct<3, 1, 3> >;


#if BOOST_PP_GREATER_EQUAL( LIFE_MESH_MAX_ORDER, 2 )
template class Mesh<Simplex<3, 2, 3> >;
template class Mesh<SimplexProduct<3, 2, 3> >;
#if BOOST_PP_GREATER_EQUAL( LIFE_MESH_MAX_ORDER, 3 )
template class Mesh<Simplex<3, 3, 3> >;
template class Mesh<SimplexProduct<3, 3, 3> >;
#if BOOST_PP_GREATER_EQUAL( LIFE_MESH_MAX_ORDER, 4 )
template class Mesh<Simplex<3, 4, 3> >;
template class Mesh<SimplexProduct<3, 4, 3> >;
#if BOOST_PP_GREATER_EQUAL( LIFE_MESH_MAX_ORDER, 5 )
template class Mesh<Simplex<3, 5, 3> >;
template class Mesh<SimplexProduct<3, 5, 3> >;
#endif // 5
#endif // 4
#endif // 3
#endif // 2

#endif // LIFE_INSTANTIATION_MODE

}
