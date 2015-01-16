/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s:

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
  \file loadmesh.hpp
  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  \date 2013-12-24
  */
#if !defined(FEELPP_LOADMESH_HPP)
#define FEELPP_LOADMESH_HPP 1

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feelfilters/domain.hpp>

namespace Feel {

/**
 *
 * \brief load a mesh data structure (hold in a shared_ptr<>) using GMSH
 *
 * \arg mesh mesh data structure
 * \arg filename filename string (with extension)
 * \arg refine optionally refine with \p refine levels the mesh (default: 0)
 * \arg update update the mesh data structure (build internal faces and edges) (default : true)
 * \arg physical_are_elementary_regions boolean to load specific meshes formats (default : false)
 */
BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::mesh<Args>::ptrtype ), // return type
    loadMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( mesh, *)

        ) // 4. one required parameter, and

    ( optional
      ( filename, *( boost::is_convertible<mpl::_,std::string> ), soption(_name="gmsh.filename") )
      ( desc, *,boost::shared_ptr<gmsh_type>() )  // geo() can't be used here as default !!

      ( h,              *( boost::is_arithmetic<mpl::_> ), doption(_name="gmsh.hsize") )
      ( straighten,          (bool), boption(_name="gmsh.straighten") )
      ( refine,          *( boost::is_integral<mpl::_> ), ioption(_name="gmsh.refine") )
      ( update,          *( boost::is_integral<mpl::_> ), MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES )
      ( physical_are_elementary_regions,		   (bool), boption(_name="gmsh.physical_are_elementary_regions") )
      ( worldcomm,       (WorldComm), Environment::worldComm() )
      ( force_rebuild,   *( boost::is_integral<mpl::_> ), boption(_name="gmsh.rebuild") )
      ( respect_partition,	(bool), boption(_name="gmsh.respect_partition") )
      ( rebuild_partitions,	(bool), boption(_name="gmsh.partition") )
      ( rebuild_partitions_filename, *( boost::is_convertible<mpl::_,std::string> )	, filename )
      ( partitions,      *( boost::is_integral<mpl::_> ), worldcomm.globalSize() )
      ( partitioner,     *( boost::is_integral<mpl::_> ), ioption(_name="gmsh.partitioner") )
      ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
      ( depends, *( boost::is_convertible<mpl::_,std::string> ), soption(_name="gmsh.depends") )
        )
    )
{
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunsequenced"
#endif
    typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

    // look for mesh_name in various directories (executable directory, current directory. ...)
    // return an empty string if the file is not found

    std::string filenameExpand = Environment::expand(filename);
    fs::path mesh_name=fs::path(Environment::findFile(filenameExpand));
    LOG_IF( WARNING, mesh_name.extension() != ".geo" && mesh_name.extension() != ".msh" )
        << "Invalid filename " << filenameExpand << " it should have either the .geo or .msh extension\n";


    if ( mesh_name.extension() == ".geo" )
    {

        return createGMSHMesh(
            _mesh=mesh,
            _desc= (!desc) ? geo( _filename=mesh_name.string(),
                                  _h=h,
                                  _depends=depends,
                                  _worldcomm=worldcomm  ) : desc ,
            _h=h,
            _straighten=straighten,
            _refine=refine,
            _update=update,
            _physical_are_elementary_regions=physical_are_elementary_regions,
            _force_rebuild=force_rebuild,
            _worldcomm=worldcomm,
            _respect_partition=respect_partition,
            _rebuild_partitions=rebuild_partitions,
            _rebuild_partitions_filename=rebuild_partitions_filename,
            _partitions=partitions,
            _partitioner=partitioner,
            _partition_file=partition_file
            );
    }

    if ( mesh_name.extension() == ".msh"  )
    {
        return loadGMSHMesh( _mesh=mesh,
                             _filename=mesh_name.string(),
                             _straighten=straighten,
                             _refine=refine,
                             _update=update,
                             _physical_are_elementary_regions=physical_are_elementary_regions,
                             _worldcomm=worldcomm,
                             _respect_partition=respect_partition,
                             _rebuild_partitions=rebuild_partitions,
                             _rebuild_partitions_filename=rebuild_partitions_filename,
                             _partitions=partitions,
                             _partitioner=partitioner,
                             _partition_file=partition_file
            );

    }

    LOG(WARNING) << "File " << mesh_name << " not found, generating instead an hypercube in " << _mesh_type::nDim << "D geometry and mesh...";
    return createGMSHMesh(_mesh=mesh,
                          _desc=domain( _name=soption(_name="gmsh.domain.shape"), _h=h, _worldcomm=worldcomm ),
                          _h=h,
                          _refine=refine,
                          _update=update,
                          _physical_are_elementary_regions=physical_are_elementary_regions,
                          _force_rebuild=force_rebuild,
                          _worldcomm=worldcomm,
                          _respect_partition=respect_partition,
                          _rebuild_partitions=rebuild_partitions,
                          _rebuild_partitions_filename=rebuild_partitions_filename,
                          _partitions=partitions,
                          _partitioner=partitioner,
                          _partition_file=partition_file );
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
} // loadMesh

} // Feel namespace

#endif /* FEELPP_LOADMESH_HPP */
