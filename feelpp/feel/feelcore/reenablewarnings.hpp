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
//! @date 29 Apr 2018
//! @copyright 2018 Feel++ Consortium
//!

#ifndef FEELPP_REENABLEWARNINGS_HPP
#define FEELPP_REENABLEWARNINGS_HPP 1

#ifndef FEELPP_PERMANENTLY_DISABLE_WARNINGS
#ifdef _MSC_VER
#pragma warning( pop )
#elif defined __INTEL_COMPILER
#pragma warning pop
#elif defined __clang__
#pragma clang diagnostic pop
#elif defined __GNUC__ && __GNUC__>=6
#pragma GCC diagnostic pop
#endif

#if defined __NVCC__
//    Don't reenable the diagnostic messages, as it turns out these messages need
//    to be disabled at the point of the template instantiation (i.e the user code)
//    otherwise they'll be triggered by nvcc.
//    #pragma diag_default code_is_unreachable
//    #pragma diag_default initialization_not_reachable
//    #pragma diag_default 2651
//    #pragma diag_default 2653
#endif

#endif

#endif
