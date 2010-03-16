/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-04-03

  Copyright (C) 2009 Université de Grenoble 1

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
   \file parameter.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-04-03
 */
#ifndef __lifecore_parameter_H
#define __lifecore_parameter_H 1


#if !defined(BOOST_PARAMETER_MAX_ARITY)
#define BOOST_PARAMETER_MAX_ARITY 10
#endif

#include <boost/parameter.hpp>
#include <boost/type_traits.hpp>
#if 0
#include <boost/parameter/keyword.hpp>
#include <boost/parameter/aux_/maybe.hpp>
#include <boost/parameter/name.hpp>
#include <boost/parameter/preprocessor.hpp>
#endif

namespace Life
{
namespace parameter = boost::parameter;

BOOST_PARAMETER_NAME(matrix)    // Note: no semicolon
BOOST_PARAMETER_NAME(matrixA)
BOOST_PARAMETER_NAME(matrixB)
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
BOOST_PARAMETER_NAME(name)
BOOST_PARAMETER_NAME(nev)
BOOST_PARAMETER_NAME(ncv)
BOOST_PARAMETER_NAME(backend)
BOOST_PARAMETER_NAME(problem)
BOOST_PARAMETER_NAME(solver)
BOOST_PARAMETER_NAME(spectrum)
BOOST_PARAMETER_NAME(transform)
// parameter for description of geometries
BOOST_PARAMETER_NAME(h)
BOOST_PARAMETER_NAME(dim)
BOOST_PARAMETER_NAME(order)
BOOST_PARAMETER_NAME(addmidpoint)
BOOST_PARAMETER_NAME(xmin)
BOOST_PARAMETER_NAME(xmax)
BOOST_PARAMETER_NAME(ymin)
BOOST_PARAMETER_NAME(ymax)
BOOST_PARAMETER_NAME(zmin)
BOOST_PARAMETER_NAME(zmax)
// parameter for xmlParse
BOOST_PARAMETER_NAME(type)
BOOST_PARAMETER_NAME(latex)
BOOST_PARAMETER_NAME(cmdName)
BOOST_PARAMETER_NAME(values)
BOOST_PARAMETER_NAME(dependencies)
BOOST_PARAMETER_NAME(funcs)
BOOST_PARAMETER_NAME(mesh)
BOOST_PARAMETER_NAME(desc)
BOOST_PARAMETER_NAME(shape)
BOOST_PARAMETER_NAME(convex)
// orders
BOOST_PARAMETER_NAME(order_u)
BOOST_PARAMETER_NAME(order_p)

BOOST_PARAMETER_NAME(initial_time)
BOOST_PARAMETER_NAME(final_time)
BOOST_PARAMETER_NAME(time_step)
BOOST_PARAMETER_NAME(strategy)

}
#endif /* __lifecore_parameter_H */
