/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2016 Feel++ Consortium

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
   \file savegmshmesh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#ifndef FEELPP_FILTERS_SAVEGMSHMESH_H
#define FEELPP_FILTERS_SAVEGMSHMESH_H

#include <feel/feelfilters/exportergmsh.hpp>
#include <feel/feelfilters/detail/mesh.hpp>
#include <feel/feelfilters/detail/meshfromgeoentity.hpp>

namespace Feel {
/**
 *
 * \brief save a mesh data structure (hold in a shared_ptr<>) in the GMSH format
 *
 * \arg mesh mesh data structure
 * \arg filename filename string (with extension)
 */
template <typename ... Ts>
void saveGMSHMesh( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && mesh = args.get(_mesh);
    std::string const& filename = args.get(_filename);
    bool parametricnodes = args.get_else(_parametricnodes,false);

    using _mesh_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(mesh)>>>;

    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem().string(), 1,  mesh->worldCommPtr() );
    exporter.saveMesh( filename, mesh, parametricnodes );

}

template <typename ... Ts>
void saveGeoEntityAsGMSHMesh( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && geoentity = args.get(_geoentity);
    std::string const& filename = args.get(_filename);
    using meshFromGeoEntity_helper_type = Feel::detail::meshFromGeoEntity<decltype(geoentity)>;
    auto && pointset = args.get_else_invocable(_pointset, [](){ return typename meshFromGeoEntity_helper_type::pointset_type{}; });

    using _mesh_type = typename meshFromGeoEntity_helper_type::type;

    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem().string(), 1,  Environment::worldComm().subWorldCommSeq() );
    exporter.gmshSaveOneElementAsMesh( filename, geoentity, pointset );
}

}
#endif /* FEELPP_SAVEGMSHMESH_HPP */
