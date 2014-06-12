/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-05-09

  Copyright (C) 2010-2011 Universite Joseph Fourier (Grenoble I)
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
#ifndef __pnedelec_H
#define __pnedelec_H 1

#include "basisgen.hpp"
#include <feel/feelpoly/nedelec.hpp>

namespace Feel
{
namespace detail
{
template<typename Convex, uint16_type Order>
PolyvisBase*
createNedelec()
{
    std::cout << "Creating Nedelec(" << Convex::nDim << "," << Order << ")\n";
    return new Polyvis<basis_type<Nedelec<Order> >,convex_type<Convex> >;
}
# /* Generates code for all dim and order. */
# define NEDELEC_FACTORY_OP(_, GDO)           \
    NEDELEC_FACTORY GDO                        \
    /**/
#
# define GEN_NEDELEC(LDIM,LORDER,CONVEX)              \
    "nedelec(" #LDIM "," #LORDER "," #CONVEX ")"       \
    /**/
#
#
#define NEDELEC_FACTORY(LDIM,LORDER,CONVEX)                     \
    const bool BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT(lag, CONVEX), LDIM), LORDER ) = \
        PolyvisBase::factory_type::instance()                           \
        .registerProduct( GEN_NEDELEC(LDIM,LORDER,CONVEX), &createNedelec<CONVEX< LDIM >,LORDER> );

} // detail
} // Feel
#endif /* __nedelec_H */
