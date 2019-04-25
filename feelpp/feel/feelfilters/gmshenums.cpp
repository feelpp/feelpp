/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 28 Dec 2016

 Copyright (C) 2016 Feel++ Consortium

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
#include <feel/feelfilters/gmshenums.hpp>
#include <feel/feelcore/feelgmsh.hpp>

namespace Feel
{
#if GMSH_VERSION_GREATER_OR_EQUAL_THAN(4,0,0)
const char* FEELPP_GMSH_FORMAT_VERSION = "4.1";
#else
const char* FEELPP_GMSH_FORMAT_VERSION = "2.2";
#endif

#if defined(HAVE_METIS)
const GMSH_PARTITIONER GMSH_PARTITIONER_DEFAULT = GMSH_PARTITIONER_METIS;
#else
const GMSH_PARTITIONER GMSH_PARTITIONER_DEFAULT = GMSH_PARTITIONER_CHACO;
#endif

}
