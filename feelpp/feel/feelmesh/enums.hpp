/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-03-11

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

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
   \file enums.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-03-11
 */
#ifndef __MeshEnums_H
#define __MeshEnums_H 1

#include <boost/mp11/integral.hpp>

namespace  Feel
{

/**
 * entity identifier  for mesh iterators
 */
enum ElementsType
{
    MESH_ELEMENTS = 0,           /**< elements */
    MESH_FACES = 1,              /**< faces */
    MESH_INTERNAL_FACES = 2,     /**< internal faces */
    MESH_EDGES = 3,              /**< edges */
    MESH_INTERNAL_EDGES = 4,     /**< internal edges */
    MESH_POINTS = 5              /**< points */
};
using on_elements_t = boost::mp11::mp_int<MESH_ELEMENTS>;
using on_facets_t = boost::mp11::mp_int<MESH_FACES>;
using on_internal_faces_t = boost::mp11::mp_int<MESH_INTERNAL_FACES>;
using on_edges_t = boost::mp11::mp_int<MESH_FACES>;
using on_points_t = boost::mp11::mp_int<MESH_POINTS>;


enum MESH_CHANGES
{
    MESH_CHANGES_POINTS_COORDINATES = 0,
    MESH_CHANGES_CONNECTIVITY       = 1,
    MESH_CHANGES_PARTITION          = 2
};

enum class EntityProcessType {LOCAL_ONLY,GHOST_ONLY,ALL,IGNORE_ENTITY_ON_INTERPROCESS_FACE};
using entity_process_t = EntityProcessType;

}
#endif /* __MeshEnums_H */
