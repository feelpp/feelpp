/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-07-28

  Copyright (C) 2005,2006 EPFL

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
   \file ublas_traits.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-07-28
 */
#ifndef LIFE_UBLAS_TRAITS_HPP
#define LIFE_UBLAS_TRAITS_HPP 1

#if !defined( LIFE_TRAITS_HPP)
#error life/lifecore/ublas_traits.hpp must not be used directly, use life/lifecore/traits.hpp instead
#endif

#include <boost/numeric/ublas/traits.hpp>

namespace boost { namespace numeric { namespace ublas {

#if defined( HAVE_MPFR )
template<>
struct type_traits<Life::mp_type> {
    typedef type_traits<Life::mp_type> self_type;
    typedef Life::mp_type value_type;
    typedef const value_type &const_reference;
    typedef value_type &reference;
    typedef value_type real_type;
    typedef Life::mp_type precision_type;

    static const unsigned plus_complexity = 1;
    static const unsigned multiplies_complexity = 1;

    static
    BOOST_UBLAS_INLINE
    real_type real (const_reference t) {
        return t;
    }
    static
    BOOST_UBLAS_INLINE
    real_type imag (const_reference /*t*/) {
        return 0;
    }
    static
    BOOST_UBLAS_INLINE
    value_type conj (const_reference t) {
        return t;
    }

    static
    BOOST_UBLAS_INLINE
    real_type abs (const_reference t) {
        return ::abs (t);
    }

    static
    BOOST_UBLAS_INLINE
    value_type sqrt (const_reference t) {
        return ::sqrt (t);
    }

    static
    BOOST_UBLAS_INLINE
    value_type type_sqrt (const_reference t) {
        return ::sqrt (t);
    }

    static
    BOOST_UBLAS_INLINE
    real_type norm_1 (const_reference t) {
        return self_type::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_2 (const_reference t) {
        return self_type::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_inf (const_reference t) {
        return self_type::abs (t);
    }

    static
    BOOST_UBLAS_INLINE
    bool equals (const_reference t1, const_reference t2) {
        return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
            (std::max) ((std::max) (self_type::norm_inf (t1),
                                    self_type::norm_inf (t2)),
                        BOOST_UBLAS_TYPE_CHECK_MIN);
    }
};
#endif // HAVE_MPFR

#if defined ( LIFE_HAVE_QD_REAL )
template<>
struct type_traits<qd_real> {
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
    real_type real (const_reference t) {
        return t;
    }
    static
    BOOST_UBLAS_INLINE
    real_type imag (const_reference /*t*/) {
        return 0;
    }
    static
    BOOST_UBLAS_INLINE
    value_type conj (const_reference t) {
        return t;
    }

    static
    BOOST_UBLAS_INLINE
    real_type abs (const_reference t) {
        return ::abs (t);
    }

    static
    BOOST_UBLAS_INLINE
    value_type sqrt (const_reference t) {
        return ::sqrt (t);
    }

    static
    BOOST_UBLAS_INLINE
    value_type type_sqrt (const_reference t) {
        return ::sqrt (t);
    }

    static
    BOOST_UBLAS_INLINE
    real_type norm_1 (const_reference t) {
        return self_type::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_2 (const_reference t) {
        return self_type::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_inf (const_reference t) {
        return self_type::abs (t);
    }

    static
    BOOST_UBLAS_INLINE
    bool equals (const_reference t1, const_reference t2) {
        return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
            (std::max) ((std::max) (self_type::norm_inf (t1),
                                    self_type::norm_inf (t2)),
                        BOOST_UBLAS_TYPE_CHECK_MIN);
    }
};
#endif /* LIFE_HAVE_QD_REAL */
#if defined ( LIFE_HAVE_DD_REAL )
template<>
struct type_traits<dd_real> {
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
    real_type real (const_reference t) {
        return t;
    }
    static
    BOOST_UBLAS_INLINE
    real_type imag (const_reference /*t*/) {
        return 0;
    }
    static
    BOOST_UBLAS_INLINE
    value_type conj (const_reference t) {
        return t;
    }

    static
    BOOST_UBLAS_INLINE
    real_type abs (const_reference t) {
        return ::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    value_type sqrt (const_reference t) {
        return ::sqrt (t);
    }

    static
    BOOST_UBLAS_INLINE
    value_type type_sqrt (const_reference t) {
        return ::sqrt (t);
    }

    static
    BOOST_UBLAS_INLINE
    real_type norm_1 (const_reference t) {
        return self_type::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_2 (const_reference t) {
        return self_type::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_inf (const_reference t) {
        return self_type::abs (t);
    }

    static
    BOOST_UBLAS_INLINE
    bool equals (const_reference t1, const_reference t2) {
        return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
            (std::max) ((std::max) (self_type::norm_inf (t1),
                                    self_type::norm_inf (t2)),
                        BOOST_UBLAS_TYPE_CHECK_MIN);
    }
};
#endif /* LIFE_HAVE_DD_REAL */

#if defined ( LIFE_HAVE_MP_REAL )
template<>
struct type_traits<mp_real> {
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
    real_type real (const_reference t) {
        return t;
    }
    static
    BOOST_UBLAS_INLINE
    real_type imag (const_reference /*t*/) {
        return 0;
    }
    static
    BOOST_UBLAS_INLINE
    value_type conj (const_reference t) {
        return t;
    }

    static
    BOOST_UBLAS_INLINE
    real_type abs (const_reference t) {
        return ::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    value_type sqrt (const_reference t) {
        return ::sqrt (t);
    }

    static
    BOOST_UBLAS_INLINE
    value_type type_sqrt (const_reference t) {
        return ::sqrt (t);
    }

    static
    BOOST_UBLAS_INLINE
    real_type norm_1 (const_reference t) {
        return self_type::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_2 (const_reference t) {
        return self_type::abs (t);
    }
    static
    BOOST_UBLAS_INLINE
    real_type norm_inf (const_reference t) {
        return self_type::abs (t);
    }

    static
    BOOST_UBLAS_INLINE
    bool equals (const_reference t1, const_reference t2) {
        return self_type::norm_inf (t1 - t2) < BOOST_UBLAS_TYPE_CHECK_EPSILON *
            (std::max) ((std::max) (self_type::norm_inf (t1),
                                    self_type::norm_inf (t2)),
                        BOOST_UBLAS_TYPE_CHECK_MIN);
    }
};
#endif /* LIFE_HAVE_MP_REAL */
}}}
#endif /* LIFE_UBLAS_TRAITS_HPP */
