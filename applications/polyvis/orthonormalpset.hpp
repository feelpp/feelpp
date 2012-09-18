/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-30

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file dubiner.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-30
 */
#ifndef __pdubiner_H
#define __pdubiner_H 1

#include "basisgen.hpp"

namespace Feel
{
namespace detail
{
template<typename Convex,int Order>
PolyvisBase*
createOrthonormalpset()
{
    return new Polyvis<convex_type<Convex>, basis_type<Feel::OrthonormalPolynomialSet<Order> > >;
}
# /* Generates code for all dim and order. */
# define ORTHONORMALPSET_FACTORY_OP(_, GDO)            \
    ORTHONORMALPSET_FACTORY GDO                        \
    /**/
#
# define GEN_ORTHONORMALPSET(LDIM,LORDER,CONVEX)            \
    "orthonormalpset(" #LDIM "," #LORDER "," #CONVEX ")"     \
    /**/
#
#
#define ORTHONORMALPSET_FACTORY(LDIM,LORDER,CONVEX)              \
    const bool BOOST_PP_CAT( BOOST_PP_CAT(BOOST_PP_CAT(orthnorm, CONVEX), LDIM), LORDER ) =        \
        PolyvisBase::factory_type::instance()                           \
        .registerProduct( GEN_ORTHONORMALPSET(LDIM,LORDER,CONVEX), &createOrthonormalpset<CONVEX<LDIM>,LORDER> );

} // detail
} // Feel
#endif /* __orthonormalpset_H */

