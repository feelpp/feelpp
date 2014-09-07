/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-30

  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)
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
   \file basisgen.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-30
 */
#ifndef __basisgen_H
#define __basisgen_H 1

#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

#include "polyvis.hpp"

namespace Feel
{
namespace detail
{

# define DIMS                                               \
    BOOST_PP_TUPLE_TO_LIST(                                 \
                           3,                               \
                           (                                \
                            1,2,3                           \
                                                            )   \
                                                            )
# define CONVEX                                 \
    BOOST_PP_TUPLE_TO_LIST(                                 \
                           2,                               \
                           (                                \
                            Simplex,Hypercube              \
                                                            )   \
                                                            )
#
# define SIMPLEX BOOST_PP_TUPLE_TO_LIST(1, (Simplex))
#
# define SIMPLEX_PRODUCT BOOST_PP_TUPLE_TO_LIST(1, (Hypercube))
#
# define DIM1 BOOST_PP_TUPLE_TO_LIST(1, (1))
# define DIM2 BOOST_PP_TUPLE_TO_LIST(1, (2))
# define DIM3 BOOST_PP_TUPLE_TO_LIST(1, (3))
#
#if 1
# define ORDER0 BOOST_PP_TUPLE_TO_LIST(1,(0))
# define ORDER1 BOOST_PP_TUPLE_TO_LIST(1,(1))
# define ORDER01 BOOST_PP_TUPLE_TO_LIST(2,(0,1))
# define ORDER0123 BOOST_PP_TUPLE_TO_LIST(4,(0,1,2,3))

# define ORDERS                                             \
    BOOST_PP_TUPLE_TO_LIST(                                 \
                           10,                              \
                           (                                \
                            1,2,3,4,5,6,7,8,9,10            \
                                                            )   \
                                                            )
#else
# define ORDERS BOOST_PP_TUPLE_TO_LIST(3,(1,2,3))
#endif
} // detail
} // Feel
#endif /* __basisgen_H */
