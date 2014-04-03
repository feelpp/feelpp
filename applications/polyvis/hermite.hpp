/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-07-02

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium

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
   \file hermite.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-07-02
 */
#ifndef __phermite_H
#define __phermite_H 1

#include "basisgen.hpp"
#include <feel/feelpoly/hermite.hpp>

namespace Feel
{
namespace detail
{
template<typename Convex, uint16_type Order>
PolyvisBase*
createHermite()
{
    std::cout << "Creating Hermite(" << Convex::nDim << "," << Order << ")\n";
    return new Polyvis<basis_type<Hermite<Order> >,convex_type<Convex> >;
}
# /* Generates code for all dim and order. */
# define HERMITE_FACTORY_OP(_, GDO)           \
    HERMITE_FACTORY GDO                        \
    /**/
#
# define GEN_HERMITE(LDIM,LORDER,CONVEX)              \
    "hermite(" #LDIM "," #LORDER "," #CONVEX ")"       \
    /**/
#
#
#define HERMITE_FACTORY(LDIM,LORDER,CONVEX)                     \
    const bool BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT(lag, CONVEX), LDIM), LORDER ) = \
        PolyvisBase::factory_type::instance()                           \
        .registerProduct( GEN_HERMITE(LDIM,LORDER,CONVEX), &createHermite<CONVEX< LDIM >,LORDER> );

} // detail
} // Feel
#endif /* __hermite_H */
