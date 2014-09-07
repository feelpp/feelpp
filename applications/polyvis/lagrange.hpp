/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-23

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium

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
   \file lagrange.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-23
 */
#ifndef __plagrange_H
#define __plagrange_H 1

#include "basisgen.hpp"

namespace Feel
{
namespace detail
{
template<typename Convex, uint16_type Order>
PolyvisBase*
createLagrange()
{
    std::cout << "Creating Lagrange(" << Convex::nDim << "," << Order << ")\n";
    return new Polyvis<basis_type<Lagrange<Order> >,convex_type<Convex> >;
}
# /* Generates code for all dim and order. */
# define LAGRANGE_FACTORY_OP(_, GDO)           \
    LAGRANGE_FACTORY GDO                        \
    /**/
#
# define GEN_LAGRANGE(LDIM,LORDER,CONVEX)              \
    "lagrange(" #LDIM "," #LORDER "," #CONVEX ")"       \
    /**/
#
#
#define LAGRANGE_FACTORY(LDIM,LORDER,CONVEX)                     \
    const bool BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT(lag, CONVEX), LDIM), LORDER ) = \
        PolyvisBase::factory_type::instance()                           \
        .registerProduct( GEN_LAGRANGE(LDIM,LORDER,CONVEX), &createLagrange<CONVEX< LDIM >,LORDER> );

} // detail
} // Feel
#endif /* __lagrange_H */
