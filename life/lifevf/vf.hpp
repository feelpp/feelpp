/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-01-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-01-17
 */
#ifndef __VF_H
#define __VF_H 1

//#include <boost/numeric/bindings/traits/traits.hpp>
//#include <boost/numeric/bindings/traits/ublas_vector.hpp>
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
//#include <boost/numeric/bindings/blas/blas.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>

/**
 * \brief allow automatic type naming of complex expression
 */
#define AUTO( a, b ) __typeof__( b ) a = (b);

namespace Life
{
//namespace blas = boost::numeric::bindings::blas;
//namespace traits = boost::numeric::bindings::traits;
namespace fusion = boost::fusion;


namespace vf
{
namespace detail
{

/// \cond detail
template<int Index> struct gmc
{
    static const int value = Index;
    typedef mpl::void_ reference_element_type;
} ;

template<typename Geo_t>
struct ExtractGm
{
    typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, detail::gmc<0> >,mpl::identity<detail::gmc<0> >,mpl::identity<detail::gmc<1> > >::type::type key_type;
    typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::pointer gmc_ptrtype;
    typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;

    static gmc_ptrtype get( Geo_t const& geom )
    {
        return fusion::at_key<key_type>( geom ).get();
    }
    static Geo_t clone( Geo_t const& geom )
    {
        Geo_t geom2( geom );
        fusion::at_key<key_type>( geom2 )  = fusion::at_key<key_type>( geom )->clone();
        return geom2;
    }
};
/// \endcond
}
}
}
#include <life/lifecore/life.hpp>
#include <life/lifevf/expr.hpp>

#include <life/lifevf/ppoperators.hpp>

#include <life/lifevf/operators.hpp>
//#include <life/lifevf/operators2.hpp>
//#include <life/lifevf/operators3.hpp>
#include <life/lifevf/geometricdata.hpp>
#include <life/lifevf/stdmathfunctors.hpp>
#include <life/lifevf/trace.hpp>
//#include <life/lifevf/symm.hpp>
#include <life/lifevf/norm.hpp>
#include <life/lifevf/ones.hpp>
#include <life/lifevf/twovalued.hpp>
//#include <life/lifevf/eye.hpp>
#include <life/lifevf/val.hpp>
#include <life/lifevf/function.hpp>
#include <life/lifevf/matvec.hpp>
//#include <life/lifevf/integral.hpp>


#include <life/lifevf/integrator.hpp>
//#include <life/lifevf/integratordirac.hpp>
#include <life/lifevf/projectors.hpp>


#include <life/lifevf/form.hpp>

#include <boost/preprocessor/comparison/equal.hpp>

/// \endcond

#endif /* __VF_H */
