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

#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

# define FEDIMS1 BOOST_PP_TUPLE_TO_LIST(1,(1))
# define FEDIMS2 BOOST_PP_TUPLE_TO_LIST(1,(2))
# define FEDIMS3 BOOST_PP_TUPLE_TO_LIST(1,(3))

# define FEORDERS BOOST_PP_TUPLE_TO_LIST(4,(0,1,2,3))
# define FETENSOR BOOST_PP_TUPLE_TO_LIST(3,(Scalar,Vectorial,Tensor2))
# define FECONTINUITY   BOOST_PP_TUPLE_TO_LIST(2,(Continuous,Discontinuous))
# define FECONTINUOUS   BOOST_PP_TUPLE_TO_LIST(1,(Continuous))
# define FEDISCONTINOUS BOOST_PP_TUPLE_TO_LIST(1,(Discontinuous))

# define FACTORY_PCH(LDIM,LORDER,TENSOR,CONT) template class FunctionSpace<Mesh<Simplex<LDIM> >,bases<Lagrange<LORDER,TENSOR,CONT,PointSetEquiSpaced> >, Periodicity<NoPeriodicity> >;
# define FACTORY_QCH(LDIM,LORDER,TENSOR,CONT) template class FunctionSpace<Mesh<Hypercube<LDIM> >,bases<Lagrange<LORDER,TENSOR,CONT,PointSetEquiSpaced> >, Periodicity<NoPeriodicity> >;

# define FACTORY_PCH_OP(_, GDO) FACTORY_PCH GDO
# define FACTORY_QCH_OP(_, GDO) FACTORY_QCH GDO

# define FACTORY_PCH_OP_E(_, GDO) extern FACTORY_PCH GDO
# define FACTORY_QCH_OP_E(_, GDO) extern FACTORY_QCH GDO

#if !defined( FEELPP_FUNCTIONSPACE_NOEXTERN )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_PCH_OP_E, 4, ( FEDIMS1, FEORDERS, FETENSOR, FECONTINUITY ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_QCH_OP_E, 4, ( FEDIMS1, FEORDERS, FETENSOR, FECONTINUITY ) )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_PCH_OP_E, 4, ( FEDIMS2, FEORDERS, FETENSOR, FECONTINUITY ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_QCH_OP_E, 4, ( FEDIMS2, FEORDERS, FETENSOR, FECONTINUITY ) )

BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_PCH_OP_E, 4, ( FEDIMS3, FEORDERS, FETENSOR, FECONTINUITY ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_QCH_OP_E, 4, ( FEDIMS3, FEORDERS, FETENSOR, FECONTINUITY ) )

#endif // FEELPP_FUNCTIONSPACE_NOEXTERN

} // Feel

#endif // FEELPP_FUNCTIONSPACE_INST_HPP
