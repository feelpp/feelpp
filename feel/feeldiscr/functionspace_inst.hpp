/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Tue Feb 25 08:45:06 2014

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
#ifndef FEELPP_FUNCTIONSPACE_INST_HPP
#define FEELPP_FUNCTIONSPACE_INST_HPP 1

namespace Feel {

# define FEDIMS1 BOOST_PP_TUPLE_TO_LIST(1,(1))
# define FEORDERS1 BOOST_PP_TUPLE_TO_LIST(2,(1,2))

# define FEDIMS2 BOOST_PP_TUPLE_TO_LIST(1,(2))
# define FEORDERS2 BOOST_PP_TUPLE_TO_LIST(2,(1,2))

# define FEDIMS3 BOOST_PP_TUPLE_TO_LIST(1,(3))
# define FEORDERS3 BOOST_PP_TUPLE_TO_LIST(2,(1,2))

# define FACTORY_PCH(LDIM,LORDER) template class FunctionSpace<Mesh<Simplex<LDIM> >,bases<Lagrange<LORDER,Scalar> > >;
# define FACTORY_QCH(LDIM,LORDER) template class FunctionSpace<Mesh<Hypercube<LDIM> >,bases<Lagrange<LORDER,Scalar> > >;

# define FACTORY_PCH_OP(_, GDO) FACTORY_PCH GDO
# define FACTORY_QCH_OP(_, GDO) FACTORY_QCH GDO

# define FACTORY_PCH_OP_E(_, GDO) extern FACTORY_PCH GDO
# define FACTORY_QCH_OP_E(_, GDO) extern FACTORY_QCH_E GDO

#if !defined( FEELPP_MESH_IMPL_NOEXTERN )

#if 0
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_SIMPLEX_OP_E, 3, ( DIMS1, ORDERS1, RDIMS1 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_HYPERCUBE_OP_E, 3, ( DIMS1, ORDERS1, RDIMS1 ) )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_SIMPLEX_OP_E, 3, ( DIMS2, ORDERS2, RDIMS2 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_HYPERCUBE_OP_E, 3, ( DIMS2, ORDERS2, RDIMS2 ) )
#endif

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_PCH_OP_E, 2, ( FEDIMS3, FEORDERS3 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_QCH_OP_E, 2, ( FEDIMS3, FEORDERS3 ) )

#endif // FEELPP_MESH_IMPL_NOEXTERN

typedef bases<Lagrange<1, Scalar>> scalar_p1_basis_type;
typedef bases<Lagrange<1, Vectorial>> vectorial_p1_basis_type;

//2d
typedef Mesh< Simplex<2> > triangle_p1_mesh_type;
typedef FunctionSpace<triangle_p1_mesh_type, scalar_p1_basis_type, Periodicity <NoPeriodicity>> pch_2d_type;

#if !defined(FEELPP_FUNCTIONSPACE_NOEXTERN)
extern template class pch_2d_type;

#endif /* FEELPP_FUNCTIONSPACE_NOEXTERN */

}

#endif // FEELPP_FUNCTIONSPACE_INST_HPP
