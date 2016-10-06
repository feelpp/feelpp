/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s:

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
  \file loadmesh.hpp
  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  \date 2013-12-24
  */
#if !defined(FEELPP_LOADMESH_HPP)
#define FEELPP_LOADMESH_HPP 1

#ifdef FEELPP_HAS_GMSH

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/importeracusimrawmesh.hpp>
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
      ( prefix,(std::string), "" )
      ( filename, *( boost::is_convertible<mpl::_,std::string> ), soption(_prefix=prefix,_name="gmsh.filename") )
      ( desc, *,boost::shared_ptr<gmsh_type>() )  // geo() can't be used here as default !!

      ( h,              *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.hsize") )
      ( scale,          *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.scale") )
      ( straighten,          (bool), boption(_prefix=prefix,_name="gmsh.straighten") )
      ( refine,          *( boost::is_integral<mpl::_> ), ioption(_prefix=prefix,_name="gmsh.refine") )
      ( update,          *( boost::is_integral<mpl::_> ), MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES )
      ( physical_are_elementary_regions,		   (bool), boption(_prefix=prefix,_name="gmsh.physical_are_elementary_regions") )
      ( worldcomm,       (WorldComm), mesh->worldComm() )
      ( force_rebuild,   *( boost::is_integral<mpl::_> ), boption(_prefix=prefix,_name="gmsh.rebuild") )
      ( respect_partition,	(bool), boption(_prefix=prefix,_name="gmsh.respect_partition") )
      ( rebuild_partitions,	(bool), boption(_prefix=prefix,_name="gmsh.partition") )
      ( rebuild_partitions_filename, *( boost::is_convertible<mpl::_,std::string> )	, filename )
      ( partitions,      *( boost::is_integral<mpl::_> ), worldcomm.globalSize() )
      ( partitioner,     *( boost::is_integral<mpl::_> ), ioption(_prefix=prefix,_name="gmsh.partitioner") )
      ( savehdf5,        *( boost::is_integral<mpl::_> ), boption(_prefix=prefix,_name="gmsh.savehdf5") )
      ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
      ( depends, *( boost::is_convertible<mpl::_,std::string> ), soption(_prefix=prefix,_name="gmsh.depends") )
      ( verbose,   (int), ioption(_prefix=prefix,_name="gmsh.verbosity") )
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
    int proc_rank = worldcomm.globalRank();
    //Environment::isMasterRank()

    LOG_IF( WARNING,
            mesh_name.extension() != ".geo" &&
            mesh_name.extension() != ".json" &&
            mesh_name.extension() != ".msh" &&
            mesh_name.extension() != ".arm" )
        << "Invalid filename " << filenameExpand << " it should have either the .geo. .json or .msh extension\n";


    if ( mesh_name.extension() == ".geo" )
    {
#if defined(FEELPP_HAS_HDF5)
        if ( boption(_name="mesh.load.enable") && soption(_name="mesh.load.format") == "json+h5" )
        {
            auto json_fname = mesh_name.stem().string()+".json";
            if( fs::exists(json_fname) )
            {
                if ( worldcomm.isMasterRank() )
                    std::cout << "[loadMesh] Loading mesh in format json+h5: " << fs::system_complete(json_fname) << "\n";
                LOG(INFO) << " Loading mesh in format json+h5: " << json_fname;
                CHECK( mesh ) << "Invalid mesh pointer to load " << json_fname;
                _mesh_ptrtype m( mesh );
                m->setWorldComm( worldcomm );
                m->loadHDF5( json_fname );
                return m;
            }
        }
#endif
        if ( worldcomm.isMasterRank() )
        {
            std::cout << "[loadMesh] Loading mesh in format geo+msh: " << fs::system_complete(mesh_name) << "\n";
            if ( !desc )
                std::cout << "[loadMesh] Use default geo desc: " << mesh_name.string() << " " << h << " " << depends << "\n";
        }
        auto m = createGMSHMesh(
            _mesh=mesh,
            _desc= (!desc) ? geo( _filename=mesh_name.string(),
                                  _h=h,
                                  _depends=depends,
                                  _worldcomm=worldcomm  ) : desc ,
            _h=h,
            _scale=scale,
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
            _partition_file=partition_file,
            _verbose=verbose
                                );

#if defined(FEELPP_HAS_HDF5)
        if ( savehdf5 )
            m->saveHDF5( mesh_name.stem().string()+".json" );
#endif
        return m;
    }

    if ( mesh_name.extension() == ".msh"  )
    {
        if ( worldcomm.isMasterRank() )
            std::cout << "[loadMesh] Loading mesh in format msh: " << fs::system_complete(mesh_name) << "\n";
        auto m = loadGMSHMesh( _mesh=mesh,
                               _filename=mesh_name.string(),
                               _straighten=straighten,
                               _refine=refine,
                               _scale=scale,
                               _update=update,
                               _physical_are_elementary_regions=physical_are_elementary_regions,
                               _worldcomm=worldcomm,
                               _respect_partition=respect_partition,
                               _rebuild_partitions=rebuild_partitions,
                               _rebuild_partitions_filename=rebuild_partitions_filename,
                               _partitions=partitions,
                               _partitioner=partitioner,
                               _partition_file=partition_file,
                               _verbose=verbose
                               );
#if defined(FEELPP_HAS_HDF5)
        if ( savehdf5 )
            m->saveHDF5( mesh_name.stem().string()+".json" );
#endif
        return m;
    }
#if defined(FEELPP_HAS_HDF5)
    if ( mesh_name.extension() == ".json"  )
    {
        if ( worldcomm.isMasterRank() )
            std::cout << "[loadMesh] Loading mesh in format json+h5: " << fs::system_complete(mesh_name) << "\n";
        LOG(INFO) << " Loading mesh in json+h5 format " << fs::system_complete(mesh_name);
        CHECK( mesh ) << "Invalid mesh pointer to load " << mesh_name;
        _mesh_ptrtype m( mesh );
        m->setWorldComm( worldcomm );
        m->loadHDF5( mesh_name.string(), update );
        return m;
    }
#endif

    // Acusim Raw Mesh
    if ( mesh_name.extension() == ".arm"  )
    {
        if ( worldcomm.isMasterRank() )
            std::cout << "[loadMesh] Loading mesh in format arm(acusolve)h5: " << fs::system_complete(mesh_name) << "\n";
        LOG(INFO) << " Loading mesh in arm(acusolve) format " << fs::system_complete(mesh_name);
        CHECK( mesh ) << "Invalid mesh pointer to load " << mesh_name;
        _mesh_ptrtype m( mesh );
        m->setWorldComm( worldcomm );
        ImporterAcusimRawMesh<_mesh_type> i( mesh_name.string(), worldcomm );
        i.visit( m.get() );
        m->components().reset();
        m->components().set( update );
        m->updateForUse();
#if defined(FEELPP_HAS_HDF5)
        if ( savehdf5 )
            m->saveHDF5( fs::path(filenameExpand).stem().string()+".json" );
#endif
        return m;
    }

    mesh_name = soption(_name="gmsh.domain.shape");
    if ( worldcomm.isMasterRank() )
        std::cout << "[loadMesh] no file name or unrecognized extension provided\n"
                  << "[loadMesh] automatically generating amesh from gmsh.domain.shape in format geo+msh: " 
                  << mesh_name << ".geo\n";
    LOG(WARNING) << "File " << mesh_name << " not found, generating instead an hypercube in " << _mesh_type::nDim << "D geometry and mesh...";
    auto m = createGMSHMesh(_mesh=mesh,
                            _desc=domain( _name=mesh_name.string(), _h=h, _worldcomm=worldcomm ),
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
                            _partition_file=partition_file,
                            _verbose=verbose);

#if defined(FEELPP_HAS_HDF5)
    if ( savehdf5 )
        m->saveHDF5( fs::path(filenameExpand).stem().string()+".json" );
#endif
    return m;
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
} // loadMesh

} // Feel namespace

#endif

#endif /* FEELPP_LOADMESH_HPP */
