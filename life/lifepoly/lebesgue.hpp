/* -*- mode: c++ -*-

This file is part of the Life library

Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
Date: 2006-03-06

Copyright (C) 2006 EPFL

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
   \file lebesgue.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2006-09-26
*/
#ifndef __lebesgueConstant_H
#define __lebesgueConstant_H 1

#include <life/lifepoly/equispaced.hpp>
#include <life/lifepoly/lagrange.hpp>

namespace Life
{
namespace ublas = boost::numeric::ublas;

template< class Convex,
          uint16_type Order,
          template<class, uint16_type, class> class PointSetType,
          typename value_type >
value_type lebesgueConstant()
{
    static const uint16_type Dim = Convex::nDim;

    typedef typename mpl::if_< mpl::bool_< Convex::is_simplex >, Simplex<Dim, 1>, SimplexProduct<Dim, 1> >::type convex_type;

    typedef typename mpl::if_< mpl::bool_< Convex::is_simplex >,
        fem::Lagrange<Dim, Order, Scalar, Continuous, value_type, Simplex, PointSetType >,
        fem::Lagrange<Dim, Order, Scalar, Continuous, value_type, SimplexProduct, PointSetType > >::type basis_type;

    basis_type _M_basis;

    PointSetEquiSpaced<convex_type, 100 - 80*(Dim-2), value_type> test;

    return ublas::norm_1( _M_basis.evaluate( test.points() ) );
};

} // Life
#endif /* __LebesgueConstant_H */
