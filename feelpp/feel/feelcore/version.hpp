/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 17 Feb 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_VERSION_HPP
#define FEELPP_VERSION_HPP 1

#include <feel/feelinfo.h>

#define FEELPP_MAKE_VERSION( a,b,c ) (((a) << 16) | ((b) << 8) | (c))
#define FEELPP_IS_VERSION(a,b,c) ( FEELPP_VERSION >= FEELPP_MAKE_VERSION(a,b,c) )

#define FEELPP_VERSION_LESS_THAN(major,minor,micro)                   \
    ((FEELPP_VERSION_MAJOR < (major) ||                                  \
      (FEELPP_VERSION_MAJOR == (major) && (FEELPP_VERSION_MINOR < (minor) || \
                                          (FEELPP_VERSION_MINOR == (minor) && \
                                           FEELPP_VERSION_MICRO < (micro))))) ? 1 : 0)

#define FEELPP_VERSION_GREATER_THAN(major,minor,micro)                \
    ((FEELPP_VERSION_MAJOR > (major) ||                                  \
      (FEELPP_VERSION_MAJOR == (major) && (FEELPP_VERSION_MINOR > (minor) || \
                                          (FEELPP_VERSION_MINOR == (minor) && \
                                           FEELPP_VERSION_MICRO > (micro))))) ? 1 : 0)

#define FEELPP_VERSION_GREATER_OR_EQUAL_THAN(major,minor,micro)       \
    ((FEELPP_VERSION_MAJOR > (major) ||                                  \
      (FEELPP_VERSION_MAJOR == (major) && (FEELPP_VERSION_MINOR > (minor) || \
                                          (FEELPP_VERSION_MINOR == (minor) && \
                                           FEELPP_VERSION_MICRO >= (micro))))) ? 1 : 0)



#endif
