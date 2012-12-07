/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-28

  Copyright (C) 2007,2009 Université de Grenoble 1
  Copyright (C) 2005,2006,2009 EPFL

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
   \file ublas_traits.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-28
 */
#ifndef FEELPP_UBLAS_TRAITS_HPP
#define FEELPP_UBLAS_TRAITS_HPP 1

#if !defined( FEELPP_TRAITS_HPP)
#error feel/feelcore/ublas_traits.hpp must not be used directly, use feel/feelcore/traits.hpp instead
#endif

#include <boost/numeric/ublas/traits.hpp>

namespace boost
{
namespace numeric
{
namespace ublas
{

#if defined( FEELPP_HAS_MPFR )
template<>
struct type_traits<Feel::mp_type>
{
    typedef type_traits<Feel::mp_type> self_type;
    typedef Feel::mp_type value_type;
    typedef const value_type &const_reference;
    typedef value_type &reference;
    typedef value_type real_type;
    typedef Feel::mp_type precision_type;

    static const unsigned plus_complexity = 1;
    static const unsigned multiplies_complexity = 1;

    static
    BOOST_UBLAS_INLINE
    real_type real ( const_reference t )
    {
        return t;
    }
    static
    BOOST_UBLAS_INLINE
    real_type imag ( const_reference /*t*/ )
    {
        return 0;
    }
    static
    BOOST_UBLAS_INLINE
    value_type conj ( const_reference t )
    {
        return t;
    }

    static
    BOOST_UBLAS_INLINE
    real_type abs ( const_reference t )
    {
        return ::abs ( t );
    }

    static
    BOOST_UBLAS_INLINE
    value_type sqrt ( const_reference t )
    {
        return ::sqrt ( t );
    }

    static
    BOOST_UBLAS_INLINE
    value_type type_sqrt ( const_reference t )
    {
        return ::sqrt ( t );
    }

    static
    BOOST_UBLAS_INLINE
    real_type norm_1 ( const_reference t )
    {
        return self_type::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_2 ( const_reference t )
    {
        return self_type::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_inf ( const_reference t )
    {
        return self_type::abs ( t );
    }

    static
    BOOST_UBLAS_INLINE
    bool equals ( const_reference t1, const_reference t2 )
    {
        return self_type::norm_inf ( t1 - t2 ) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
               ( std::max ) ( ( std::max ) ( self_type::norm_inf ( t1 ),
                                             self_type::norm_inf ( t2 ) ),
                              BOOST_UBLAS_TYPE_CHECK_MIN );
    }
};
#endif // FEELPP_HAS_MPFR

#if defined ( FEELPP_HAS_QD_REAL )
template<>
struct type_traits<qd_real>
{
    typedef type_traits<qd_real> self_type;
    typedef qd_real value_type;
    typedef const value_type &const_reference;
    typedef value_type &reference;
    typedef value_type real_type;
    typedef qd_real precision_type;

    static const unsigned plus_complexity = 1;
    static const unsigned multiplies_complexity = 1;

    static
    BOOST_UBLAS_INLINE
    real_type real ( const_reference t )
    {
        return t;
    }
    static
    BOOST_UBLAS_INLINE
    real_type imag ( const_reference /*t*/ )
    {
        return 0;
    }
    static
    BOOST_UBLAS_INLINE
    value_type conj ( const_reference t )
    {
        return t;
    }

    static
    BOOST_UBLAS_INLINE
    real_type abs ( const_reference t )
    {
        return ::abs ( t );
    }

    static
    BOOST_UBLAS_INLINE
    value_type sqrt ( const_reference t )
    {
        return ::sqrt ( t );
    }

    static
    BOOST_UBLAS_INLINE
    value_type type_sqrt ( const_reference t )
    {
        return ::sqrt ( t );
    }

    static
    BOOST_UBLAS_INLINE
    real_type norm_1 ( const_reference t )
    {
        return self_type::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_2 ( const_reference t )
    {
        return self_type::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_inf ( const_reference t )
    {
        return self_type::abs ( t );
    }

    static
    BOOST_UBLAS_INLINE
    bool equals ( const_reference t1, const_reference t2 )
    {
        return self_type::norm_inf ( t1 - t2 ) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
               ( std::max ) ( ( std::max ) ( self_type::norm_inf ( t1 ),
                                             self_type::norm_inf ( t2 ) ),
                              BOOST_UBLAS_TYPE_CHECK_MIN );
    }
};
#endif /* FEELPP_HAS_QD_REAL */
#if defined ( FEELPP_HAS_DD_REAL )
template<>
struct type_traits<dd_real>
{
    typedef type_traits<dd_real> self_type;
    typedef dd_real value_type;
    typedef const value_type &const_reference;
    typedef value_type &reference;
    typedef value_type real_type;
    typedef dd_real precision_type;

    static const unsigned plus_complexity = 1;
    static const unsigned multiplies_complexity = 1;

    static
    BOOST_UBLAS_INLINE
    real_type real ( const_reference t )
    {
        return t;
    }
    static
    BOOST_UBLAS_INLINE
    real_type imag ( const_reference /*t*/ )
    {
        return 0;
    }
    static
    BOOST_UBLAS_INLINE
    value_type conj ( const_reference t )
    {
        return t;
    }

    static
    BOOST_UBLAS_INLINE
    real_type abs ( const_reference t )
    {
        return ::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    value_type sqrt ( const_reference t )
    {
        return ::sqrt ( t );
    }

    static
    BOOST_UBLAS_INLINE
    value_type type_sqrt ( const_reference t )
    {
        return ::sqrt ( t );
    }

    static
    BOOST_UBLAS_INLINE
    real_type norm_1 ( const_reference t )
    {
        return self_type::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_2 ( const_reference t )
    {
        return self_type::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_inf ( const_reference t )
    {
        return self_type::abs ( t );
    }

    static
    BOOST_UBLAS_INLINE
    bool equals ( const_reference t1, const_reference t2 )
    {
        return self_type::norm_inf ( t1 - t2 ) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
               ( std::max ) ( ( std::max ) ( self_type::norm_inf ( t1 ),
                                             self_type::norm_inf ( t2 ) ),
                              BOOST_UBLAS_TYPE_CHECK_MIN );
    }
};
#endif /* FEELPP_HAS_DD_REAL */

#if defined ( FEELPP_HAS_MP_REAL )
template<>
struct type_traits<mp_real>
{
    typedef type_traits<mp_real> self_type;
    typedef mp_real value_type;
    typedef const value_type &const_reference;
    typedef value_type &reference;
    typedef value_type real_type;
    typedef double precision_type;

    static const unsigned plus_complexity = 1;
    static const unsigned multiplies_complexity = 1;

    static
    BOOST_UBLAS_INLINE
    real_type real ( const_reference t )
    {
        return t;
    }
    static
    BOOST_UBLAS_INLINE
    real_type imag ( const_reference /*t*/ )
    {
        return 0;
    }
    static
    BOOST_UBLAS_INLINE
    value_type conj ( const_reference t )
    {
        return t;
    }

    static
    BOOST_UBLAS_INLINE
    real_type abs ( const_reference t )
    {
        return ::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    value_type sqrt ( const_reference t )
    {
        return ::sqrt ( t );
    }

    static
    BOOST_UBLAS_INLINE
    value_type type_sqrt ( const_reference t )
    {
        return ::sqrt ( t );
    }

    static
    BOOST_UBLAS_INLINE
    real_type norm_1 ( const_reference t )
    {
        return self_type::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_2 ( const_reference t )
    {
        return self_type::abs ( t );
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_inf ( const_reference t )
    {
        return self_type::abs ( t );
    }

    static
    BOOST_UBLAS_INLINE
    bool equals ( const_reference t1, const_reference t2 )
    {
        return self_type::norm_inf ( t1 - t2 ) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
               ( std::max ) ( ( std::max ) ( self_type::norm_inf ( t1 ),
                                             self_type::norm_inf ( t2 ) ),
                              BOOST_UBLAS_TYPE_CHECK_MIN );
    }
};
#endif /* FEELPP_HAS_MP_REAL */
}
}
}
#endif /* FEELPP_UBLAS_TRAITS_HPP */
