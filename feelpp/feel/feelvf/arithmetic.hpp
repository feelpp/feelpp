//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 31 Mar 2019
//! @copyright 2019 Feel++ Consortium
//!
#ifndef FEELPP_VF_ARITHMETIC_HPP
#define FEELPP_VF_ARITHMETIC_HPP 1

#if defined( FEELPP_HAS_QD_H ) && defined( FEELPP_HAS_MPFR )
#define VF_CHECK_ARITHMETIC_TYPE( VALUE_TYPE )                                          \
    BOOST_STATIC_ASSERT( ( ::boost::is_arithmetic<VALUE_TYPE>::value ||                 \
                           ::boost::is_same<VALUE_TYPE, std::complex<float>>::value ||  \
                           ::boost::is_same<VALUE_TYPE, std::complex<double>>::value || \
                           ::boost::is_same<VALUE_TYPE, mp_type>::value ||              \
                           ::boost::is_same<VALUE_TYPE, dd_real>::value ||              \
                           ::boost::is_same<VALUE_TYPE, qd_real>::value ) );            \
    /**/
#elif defined( FEELPP_HAS_QD_H )
#define VF_CHECK_ARITHMETIC_TYPE( VALUE_TYPE )                                          \
    BOOST_STATIC_ASSERT( ( ::boost::is_arithmetic<VALUE_TYPE>::value ||                 \
                           ::boost::is_same<VALUE_TYPE, std::complex<float>>::value ||  \
                           ::boost::is_same<VALUE_TYPE, std::complex<double>>::value || \
                           ::boost::is_same<VALUE_TYPE, dd_real>::value ||              \
                           ::boost::is_same<VALUE_TYPE, qd_real>::value ) );            \
    /**/
#elif defined( FEELPP_HAS_MPFR )
#define VF_CHECK_ARITHMETIC_TYPE( VALUE_TYPE )                                          \
    BOOST_STATIC_ASSERT( ( ::boost::is_arithmetic<VALUE_TYPE>::value ||                 \
                           ::boost::is_same<VALUE_TYPE, std::complex<float>>::value ||  \
                           ::boost::is_same<VALUE_TYPE, std::complex<double>>::value || \
                           ::boost::is_same<VALUE_TYPE, mp_type>::value ) );            \
    /**/
#else
#define VF_CHECK_ARITHMETIC_TYPE( VALUE_TYPE )                                            \
    BOOST_STATIC_ASSERT( ( ::boost::is_arithmetic<VALUE_TYPE>::value ||                   \
                           ::boost::is_same<VALUE_TYPE, std::complex<float>>::value ||    \
                           ::boost::is_same<VALUE_TYPE, std::complex<double>>::value ) ); \
    /**/
#endif




#endif

