/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2011 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2012-2014 Feel++ Consortium

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
   \file vf.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-17
 */
#ifndef FEELPP_VF_HPP
#define FEELPP_VF_HPP 1

//#include <boost/numeric/bindings/traits/traits.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/blas/blas.hpp>


#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>

/**
 * \brief allow automatic type naming of complex expression
 */
#define AUTO( a, b ) auto a = (b);

#include <feel/feelcore/feel.hpp>

#include <feel/feelvf/detail/gmc.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/trans.hpp>
#include <feel/feelvf/unary.hpp>
#include <feel/feelvf/one.hpp>
#include <feel/feelvf/pow.hpp>
#include <feel/feelvf/minmax.hpp>


#include <feel/feelvf/operations.hpp>

#include <feel/feelvf/operators.hpp>
//#include <feel/feelvf/operators2.hpp>
//#include <feel/feelvf/operators3.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/stdmathfunctors.hpp>
#include <feel/feelvf/trace.hpp>
#include <feel/feelvf/det.hpp>
#include <feel/feelvf/symm.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/cross.hpp>
#include <feel/feelvf/norm.hpp>
#include <feel/feelvf/norm2.hpp>
#include <feel/feelvf/ones.hpp>
#include <feel/feelvf/inv.hpp>
#include <feel/feelvf/twovalued.hpp>
//#include <feel/feelvf/eye.hpp>
#include <feel/feelvf/val.hpp>
#include <feel/feelvf/function.hpp>
#include <feel/feelvf/matvec.hpp>
//#include <feel/feelvf/integral.hpp>

#include <feel/feelvf/form.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/mean.hpp>
#include <feel/feelvf/measure.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/norml2squared.hpp>
#include <feel/feelvf/normh1.hpp>
#include <feel/feelvf/normsemih1.hpp>
#include <feel/feelvf/on.hpp>
//#include <feel/feelvf/integratordirac.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelvf/evaluator.hpp>
#include <feel/feelvf/evaluatorcontext.hpp>



#include <feel/feelvf/ginac.hpp>

#include <boost/preprocessor/comparison/equal.hpp>

namespace Feel
{
using namespace vf;
}
/// \endcond

#endif /* __VF_H */
