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
//! @date 13 Aug 2014
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_FEELGMSH_HPP
#define FEELPP_FEELGMSH_HPP 1

#include <feel/feelcore/feel.hpp>

#if defined( FEELPP_HAS_GMSH_H )

#if defined( FEELPP_HAS_GMSH_API )
#include <gmsh.h>

#define GMSH_VERSION_LESS_THAN(major,minor,subminor)                   \
    ((GMSH_API_VERSION_MAJOR < (major) ||             \
      (GMSH_API_VERSION_MAJOR == (major) && (GMSH_API_VERSION_MINOR < (minor)))) ? 1 : 0)

#define GMSH_VERSION_GREATER_THAN(major,minor,subminor)                \
    ((GMSH_API_VERSION_MAJOR > (major) ||             \
      (GMSH_API_VERSION_MAJOR == (major) && (GMSH_API_VERSION_MINOR > (minor)))) ? 1 : 0)

#define GMSH_VERSION_GREATER_OR_EQUAL_THAN(major,minor,subminor)       \
    ((GMSH_API_VERSION_MAJOR > (major) ||             \
      ((GMSH_API_VERSION_MAJOR == (major)) && ((GMSH_API_VERSION_MINOR >= (minor))))) ? 1 : 0)

#else

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmacro-redefined"
#endif
#if defined(HAVE_DLOPEN)
#define HAVE_DLOPEN_WAS_DEFINED 1
#undef HAVE_DLOPEN
#endif
#include <Gmsh.h>
#include <GmshVersion.h>
#if !defined( HAVE_DLOPEN )
#if defined( HAVE_DLOPEN_WAS_DEFINED )
#define HAVE_DLOPEN 1
#undef HAVE_DLOPEN_WAS_DEFINED
#endif
#endif
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#define GMSH_VERSION_LESS_THAN(major,minor,subminor)                   \
    ((GMSH_MAJOR_VERSION < (major) ||                                  \
      (GMSH_MAJOR_VERSION == (major) && (GMSH_MINOR_VERSION < (minor) || \
                                          (GMSH_MINOR_VERSION == (minor) && \
                                           GMSH_PATCH_VERSION < (subminor))))) ? 1 : 0)

#define GMSH_VERSION_GREATER_THAN(major,minor,subminor)                \
    ((GMSH_MAJOR_VERSION > (major) ||                                  \
      (GMSH_MAJOR_VERSION == (major) && (GMSH_MINOR_VERSION > (minor) || \
                                          (GMSH_MINOR_VERSION == (minor) && \
                                           GMSH_PATCH_VERSION > (subminor))))) ? 1 : 0)

#define GMSH_VERSION_GREATER_OR_EQUAL_THAN(major,minor,subminor)       \
    ((GMSH_MAJOR_VERSION > (major) ||                                  \
      ((GMSH_MAJOR_VERSION == (major)) && ((GMSH_MINOR_VERSION > (minor)) || \
                                         ((GMSH_MINOR_VERSION == (minor)) && \
                                          ( GMSH_PATCH_VERSION >= (subminor))) ))) ? 1 : 0)

#endif


#else

#define GMSH_VERSION_LESS_THAN(major,minor,subminor)

#define GMSH_VERSION_GREATER_THAN(major,minor,subminor)

#define GMSH_VERSION_GREATER_OR_EQUAL_THAN(major,minor,subminor)

#endif /* FEELPP_HAS_GMSH_H */

#endif /* FEELPP_FEELGMSH_HPP */
