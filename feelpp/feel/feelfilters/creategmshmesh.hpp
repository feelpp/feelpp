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
#if !defined(FEELPP_CREATEGMSHMESH_HPP)
#define FEELPP_CREATEGMSHMESH_HPP 1

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
BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::mesh<Args>::ptrtype ), // return type
    createGMSHMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( mesh, * )
      ( desc, * )
        ) // 4. one required parameter, and

    ( optional
      ( prefix,(std::string), "" )
      ( vm, ( po::variables_map const& ), Environment::vm() )
      ( format,         *, ioption(_prefix=prefix,_name="gmsh.format",_vm=vm) )
      ( h,              *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.hsize",_vm=vm) )
      //( geo_parameters,  *( boost::icl::is_map<mpl::_> ), Gmsh::gpstr2map("") )
      ( scale,          *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="mesh.scale",_vm=vm) )
      ( parametricnodes, *( boost::is_integral<mpl::_> ), 0 )
      ( in_memory,       *( boost::is_integral<mpl::_> ), boption(_prefix=prefix,_name="gmsh.in-memory",_vm=vm) )
      ( straighten,      *( boost::is_integral<mpl::_> ), boption(_prefix=prefix,_name="gmsh.straighten",_vm=vm) )
      ( refine,          *( boost::is_integral<mpl::_> ), ioption(_prefix=prefix,_name="gmsh.refine",_vm=vm) )
      ( structured,          *( boost::is_integral<mpl::_> ), ioption(_prefix=prefix,_name="gmsh.structured",_vm=vm) )
      ( update,          *( boost::is_integral<mpl::_> ), MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK )
      ( force_rebuild,   *( boost::is_integral<mpl::_> ), 0 )
      ( physical_are_elementary_regions,           *,boption(_prefix=prefix,_name="gmsh.physical_are_elementary_regions",_vm=vm))
      ( periodic,        *, PeriodicEntities() )
      ( respect_partition,	(bool), boption(_prefix=prefix,_name="gmsh.respect_partition",_vm=vm) )
      ( rebuild_partitions,	(bool), boption(_prefix=prefix,_name="gmsh.partition",_vm=vm) )
      ( rebuild_partitions_filename, *( boost::is_convertible<mpl::_,std::string> )	, "" )
      ( worldcomm,      (worldcomm_ptr_t), mesh->worldCommPtr() )
      ( partitions,   *( boost::is_integral<mpl::_> ), worldcomm->localSize() )
      ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
      ( partitioner,   *( boost::is_integral<mpl::_> ), ioption(_prefix=prefix,_name="gmsh.partitioner",_vm=vm) )
      ( verbose,   (int), ioption(_prefix=prefix,_name="gmsh.verbosity",_vm=vm) )
      ( directory,(std::string), "" )
        )
    )
{
    typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

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

}

#endif /* FEELPP_CREATEGMSHMESH_HPP */
