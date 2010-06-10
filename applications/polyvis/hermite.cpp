/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-07-02

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file hermite.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-07-02
 */
#include "hermite.hpp"

namespace Life
{
namespace detail
{
# define ORDERS_HERMITE BOOST_PP_TUPLE_TO_LIST(1,(3))
# define DIMS_HERMITE BOOST_PP_TUPLE_TO_LIST(2,(1, 2))

BOOST_PP_LIST_FOR_EACH_PRODUCT(HERMITE_FACTORY_OP, 3, (DIMS_HERMITE, ORDERS_HERMITE, SIMPLEX))
}
}



