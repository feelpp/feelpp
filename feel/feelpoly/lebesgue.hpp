/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

This file is part of the Feel library

Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2006-03-06

Copyright (C) 2006 EPFL

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
   \file lebesgue.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2006-09-26
*/
#ifndef __lebesgueConstant_H
#define __lebesgueConstant_H 1

#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/lagrange.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

template< class Convex,
          uint16_type Order,
          template<class, uint16_type, class> class PointSetType,
          typename value_type >
value_type lebesgueConstant()
{
    static const uint16_type Dim = Convex::nDim;

    typedef typename mpl::if_< mpl::bool_< Convex::is_simplex >, Simplex<Dim, 1>, Hypercube<Dim, 1> >::type convex_type;

    typedef Lagrange<Order,Scalar,PointSetType>::template apply<Dim, value_type, Convex>::type basis_type;
    basis_type M_basis;

    PointSetEquiSpaced<convex_type, 100 - 80*( Dim-2 ), value_type> test;

    return ublas::norm_1( M_basis.evaluate( test.points() ) );
};

} // Feel
#endif /* __LebesgueConstant_H */
