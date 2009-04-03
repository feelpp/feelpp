/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-04-03

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
   \file parameter.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-04-03
 */
#ifndef __lifecore_parameter_H
#define __lifecore_parameter_H 1


#if !defined(BOOST_PARAMETER_MAX_ARITY)
#define BOOST_PARAMETER_MAX_ARITY 10
#else
#undef BOOST_PARAMETER_MAX_ARITY
#define BOOST_PARAMETER_MAX_ARITY 10
#endif

#include <boost/parameter.hpp>

#if 0
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/aux_/maybe.hpp>
#include <boost/parameter/name.hpp>
#include <boost/parameter/preprocessor.hpp>
#endif

namespace Life
{
BOOST_PARAMETER_NAME(matrix)    // Note: no semicolon
BOOST_PARAMETER_NAME(rhs)
BOOST_PARAMETER_NAME(solution)
BOOST_PARAMETER_NAME(prec)
BOOST_PARAMETER_NAME(maxit)
BOOST_PARAMETER_NAME(transpose)
BOOST_PARAMETER_NAME(reuse_prec)
BOOST_PARAMETER_NAME(tolerance)

BOOST_PARAMETER_NAME(test)    // Note: no semicolon
BOOST_PARAMETER_NAME(trial)
BOOST_PARAMETER_NAME(vector)
BOOST_PARAMETER_NAME(pattern)
BOOST_PARAMETER_NAME(do_threshold)
BOOST_PARAMETER_NAME(threshold)
BOOST_PARAMETER_NAME(init)

}
#endif /* __lifecore_parameter_H */
