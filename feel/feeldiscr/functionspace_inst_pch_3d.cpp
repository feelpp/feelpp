/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Tue Feb 25 09:00:42 2014

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
#define FEELPP_FUNCTIONSPACE_NOEXTERN 1
#include <feel/feeldiscr/functionspace_inst.hpp>

namespace Feel {
//
// Explicit instatiations
//
#if defined( FEELPP_INSTANTIATION_MODE )


BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY_PCH_OP, 4, ( FEDIMS3, FEORDERS, FETENSOR, FECONTINUOUS ) )




#endif // FEELPP_INSTANTIATION_MODE

}
