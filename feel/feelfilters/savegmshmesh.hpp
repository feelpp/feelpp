/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013 Feel++ Consortium

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
#if !defined(FEELPP_SAVEGMSHMESH_HPP)
#define FEELPP_SAVEGMSHMESH_HPP 1

#include <feel/feelfilters/exportergmsh.hpp>
#include <feel/feelfilters/detail/meshfromgeoentity.hpp>

namespace Feel {
/**
 *
 * \brief save a mesh data structure (hold in a shared_ptr<>) in the GMSH format
 *
 * \arg mesh mesh data structure
 * \arg filename filename string (with extension)
 */
BOOST_PARAMETER_FUNCTION(
    ( void ),          // return type
    saveGMSHMesh,    // 2. function name
    tag,             // 3. namespace of tag types
    ( required
      ( mesh, * )
      ( filename, * ) ) // 4. one required parameter, and
    ( optional
      ( parametricnodes,          *( boost::is_integral<mpl::_> ), 0 ) )
    )
{
    typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

#if BOOST_FILESYSTEM_VERSION == 3
    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem().string(), 1,  mesh->worldComm() );
#elif BOOST_FILESYSTEM_VERSION == 2
    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem(), 1, mesh->worldComm() );
#endif
    exporter.saveMesh( filename, mesh, parametricnodes );

}

BOOST_PARAMETER_FUNCTION(
    ( void ),  // return type
    saveGeoEntityAsGMSHMesh,    // 2. function name
    tag,             // 3. namespace of tag types
    ( required
      ( geoentity, * )
      ( filename, * ) ) // 4. one required parameter, and
    ( optional
      ( pointset, *, typename Feel::detail::meshFromGeoEntity<Args>::pointset_type() ) )
    )
{
    typedef typename Feel::detail::meshFromGeoEntity<Args>::type _mesh_type;

#if BOOST_FILESYSTEM_VERSION == 3
    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem().string(), 1,  Environment::worldComm().subWorldCommSeq() );
#elif BOOST_FILESYSTEM_VERSION == 2
    ExporterGmsh<_mesh_type,1> exporter( fs::path( filename ).stem(), 1, Environment::worldComm().subWorldCommSeq() );
#endif
    exporter.gmshSaveOneElementAsMesh( filename, geoentity, pointset );
}

}
#endif /* FEELPP_SAVEGMSHMESH_HPP */
