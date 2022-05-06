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
   \file creategmshmesh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#ifndef FEELPP_FILTERS_CREATEGMSHMESH_H
#define FEELPP_FILTERS_CREATEGMSHMESH_H

#include <feel/feelcore/parameter.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/detail/mesh.hpp>

#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/importergmsh.hpp>

#include <feel/feelmesh/partitionmesh.hpp>

namespace Feel {

/**
 *
 * \brief create a mesh data structure (hold in a shared_ptr<>) using GMSH
 *
 * \arg mesh mesh data structure
 * \arg descprition
 * \arg h (float, optional, default = 0.1)
 * \arg order (integer, optional, default = 1)
 * \arg parametricnodes (boolean, optional, default = 0)
 * \arg refine (boolean, optional, default = 0)
 * \arg update (boolean, optional, default = 0)
 * \arg force_rebuild boolean (boolean, optional, default = 0)
 * \arg physical_are_elementary_regions change file format (optional, default = false)
 */
template <typename T>
using args_createGMSHMesh_type = NA::arguments<
    typename na::mesh::template required_as_t<std::shared_ptr<T>>,
    typename na::desc::template required_as_t<std::shared_ptr<gmsh_type>>,
    typename na::prefix::template required_as_t<std::string const&>,
    typename na::vm::template required_as_t<po::variables_map const&>,
    typename na::format::template required_as_t<int>,
    typename na::h::template required_as_t<double>,
    typename na::scale::template required_as_t<double>,
    typename na::parametricnodes::template required_as_t<bool>,
    typename na::in_memory::template required_as_t<bool>,
    typename na::straighten::template required_as_t<bool>,
    typename na::refine::template required_as_t<int>,
    typename na::structured::template required_as_t<int>,
    typename na::update::template required_as_t<size_type>,
    typename na::force_rebuild::template required_as_t<bool>,
    typename na::physical_are_elementary_regions::template required_as_t<bool>,
    typename na::periodic::template required_as_t<PeriodicEntities const&>,
    typename na::respect_partition::template required_as_t<bool>,
    typename na::rebuild_partitions::template required_as_t<bool>,
    typename na::rebuild_partitions_filename::template required_as_t<std::string const&>,
    typename na::worldcomm::template required_as_t<worldcomm_ptr_t>,
    typename na::partitions::template required_as_t<rank_type>,
    typename na::partition_file::template required_as_t<int>,
    typename na::partitioner::template required_as_t<int>,
    typename na::verbose::template required_as_t<int>,
    typename na::directory::template required_as_t<std::string const&>
    >;

template <typename MeshType>
std::shared_ptr<MeshType>
createGMSHMesh( args_createGMSHMesh_type<MeshType> && args )
{
    auto && [mesh,desc,prefix,vm,format,h,scale,parametricnodes,in_memory,
             straighten,refine,structured,update,force_rebuild,physical_are_elementary_regions,
             periodic,respect_partition,rebuild_partitions,rebuild_partitions_filename,
             worldcomm,partitions,partition_file,partitioner,verbose,directory] = args.get_all();

    using _mesh_type = unwrap_ptr_t<std::decay_t<decltype(mesh)>>;
    using _mesh_ptrtype = std::shared_ptr<_mesh_type>;

    _mesh_ptrtype _mesh{ mesh };
    _mesh->setWorldComm( worldcomm );

    if ( worldcomm->isActive() )
    {
        desc->setDimension( mesh->nDim );
        desc->setOrder( mesh->nOrder );
        desc->setWorldComm( Environment::worldCommSeqPtr()/*worldcomm*/ );
        desc->setNumberOfPartitions( 1/*partitions*/ );
        desc->setPartitioner( (GMSH_PARTITIONER) partitioner );
        desc->setMshFileByPartition( partition_file );
        desc->setRefinementLevels( refine );
        desc->setFileFormat( (GMSH_FORMAT)format );
        desc->setStructuredMesh( structured );
        desc->setPeriodic( periodic );
        desc->setInMemory( in_memory );
        desc->setGModelName( "feelpp_gmsh_model" );
        desc->setVerbosity( verbose );

        std::string fname;
        if ( worldcomm->isMasterRank() )
        {
            bool generated_or_modified;
            boost::tie( fname, generated_or_modified ) = desc->generate( desc->prefix(), desc->description(), force_rebuild, parametricnodes, true, directory );

            // refinement if option is enabled to a value greater or equal to 1
            // do not refine if the mesh/geo file was previously generated or modified
            if ( refine && !generated_or_modified )
            {
                VLOG(1) << "Refine mesh ( level: " << refine << ")\n";
                fname = desc->refine( fname, refine, parametricnodes );
            }

            ImporterGmsh<_mesh_type> import( fname, FEELPP_GMSH_FORMAT_VERSION, (partitions > 1)? Environment::worldCommSeqPtr() : worldcomm );

            // need to replace physical_regions by elementary_regions for specific meshes
            if ( physical_are_elementary_regions )
            {
                import.setElementRegionAsPhysicalRegion( physical_are_elementary_regions );
            }
            import.setRespectPartition( respect_partition );

            if ( in_memory )
            {
                import.setGModelName( desc->gModelName() );
                import.setDeleteGModelAfterUse( true );
                import.setInMemory( in_memory );
            }

            if ( partitions > 1 )
            {
                _mesh_ptrtype _meshSeq = std::make_shared<_mesh_type>( Environment::worldCommSeqPtr() );
                _meshSeq->accept( import );
                _meshSeq->components().reset();
                _meshSeq->components().set( size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_UPDATE_FACES|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED) );
                _meshSeq->updateForUse();
#if defined(FEELPP_HAS_HDF5)
                using io_t = PartitionIO<_mesh_type>;
                std::string fnamePartitioned = rebuild_partitions_filename;
                if ( fnamePartitioned.empty() )
                    fname = fs::path( fname ).replace_extension( ".json" ).string();
                else
                    fname = fs::path( fnamePartitioned ).replace_extension( ".json" ).string();
                io_t io( fname );
                std::vector<elements_reference_wrapper_t<_mesh_type>> partitionByRange;
                io.write( partitionMesh( _meshSeq, partitions, partitionByRange ) );
#endif
            }
            else
            {
                import.setScaling( scale );
                _mesh->accept( import );
                _mesh->components().reset();
                _mesh->components().set( update );
                _mesh->updateForUse();
            }
        }
#if defined(FEELPP_HAS_HDF5)
        if ( partitions > 1 )
        {
            mpi::broadcast( worldcomm->globalComm(), fname, worldcomm->masterRank() );
            _mesh->loadHDF5( fname, update, scale );
        }
#endif
        if ( straighten && _mesh_type::nOrder > 1 )
            return straightenMesh( _mesh, worldcomm->subWorldCommPtr() );
    }
    return _mesh;
}

template <typename ... Ts>
auto
createGMSHMesh( Ts && ... v )
{
    auto args0 = NA::make_arguments( std::forward<Ts>(v)... )
        .add_default_arguments( NA::make_default_argument( _prefix, "" ),
                                NA::make_default_argument( _vm, Environment::vm() ),
                                NA::make_default_argument( _parametricnodes, 0 ),
                                NA::make_default_argument( _update, MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK ),
                                NA::make_default_argument( _force_rebuild, false ),
                                NA::make_default_argument( _periodic, PeriodicEntities() ),
                                NA::make_default_argument( _rebuild_partitions_filename, "" ),
                                NA::make_default_argument( _partition_file, false ),
                                NA::make_default_argument( _directory, "" )
                                );

    auto && mesh = args0.get(_mesh);
    std::string const& prefix = args0.get(_prefix );
    po::variables_map const& vm = args0.get(_vm );

    auto args1 = std::move( args0 ).add_default_arguments( NA::make_default_argument_invocable( _format, [&prefix,&vm](){ return ioption(_prefix=prefix,_name="gmsh.format",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _h, [&prefix,&vm](){ return doption(_prefix=prefix,_name="gmsh.hsize",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _scale, [&prefix,&vm](){ return doption(_prefix=prefix,_name="mesh.scale",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _in_memory, [&prefix,&vm](){ return boption(_prefix=prefix,_name="gmsh.in-memory",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _straighten, [&prefix,&vm](){ return boption(_prefix=prefix,_name="gmsh.straighten",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _refine, [&prefix,&vm](){ return ioption(_prefix=prefix,_name="gmsh.refine",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _structured, [&prefix,&vm](){ return ioption(_prefix=prefix,_name="gmsh.structured",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _physical_are_elementary_regions, [&prefix,&vm](){ return boption(_prefix=prefix,_name="gmsh.physical_are_elementary_regions",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _respect_partition, [&prefix,&vm](){ return boption(_prefix=prefix,_name="gmsh.respect_partition",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _rebuild_partitions, [&prefix,&vm](){ return boption(_prefix=prefix,_name="gmsh.partition",_vm=vm); } ),
                                                           NA::make_default_argument( _worldcomm, mesh->worldCommPtr() ),
                                                           NA::make_default_argument_invocable( _partitioner, [&prefix,&vm](){ return ioption(_prefix=prefix,_name="gmsh.partitioner",_vm=vm); } ),
                                                           NA::make_default_argument_invocable( _verbose, [&prefix,&vm](){ return ioption(_prefix=prefix,_name="gmsh.verbosity",_vm=vm); } )
                                                           );
    auto && worldcomm = args1.get( _worldcomm );

    auto args = std::move( args1 ).add_default_arguments( NA::make_default_argument( _partitions, (worldcomm)?worldcomm->globalSize():1 ) );

    using mesh_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(mesh)>>>;
    return createGMSHMesh<mesh_type>( std::move( args ) );
}

}

#endif /* FEELPP_CREATEGMSHMESH_HPP */
