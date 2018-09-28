/*
  This file is part of the Feel library

  Copyright (C) 2007,2008 University of Coimbra
  Copyright (C) 2010 Université de Grenoble 1

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
#include <boost/preprocessor/comparison/greater_equal.hpp>
#include <feel/feeldiscr/meshhighorderimpl.hpp>

namespace Feel
{
#if defined( FEELPP_INSTANTIATION_MODE )

template class MeshHighOrder< Simplex<2,1> >;
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
template class MeshHighOrder< Simplex<2,2> >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 3 )
template class MeshHighOrder< Simplex<2,3> >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 4 )
template class MeshHighOrder< Simplex<2,4> >;
#endif // FEELPP_MESH_MAX_ORDER
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 5 )
template class MeshHighOrder< Simplex<2,5> >;
#endif // FEELPP_MESH_MAX_ORDER

#endif // FEELPP_INSTANTIATION_MODE
} // Feel
