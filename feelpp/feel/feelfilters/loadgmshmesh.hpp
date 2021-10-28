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
   \file loadgmshmesh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined(FEELPP_LOADGMSHMESH_HPP)
#define FEELPP_LOADGMSHMESH_HPP 1

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/detail/mesh.hpp>
#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelmesh/partitionmesh.hpp>

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
#if 0
BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::mesh<Args>::ptrtype ), // return type
    loadGMSHMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( mesh, * )
      ( filename, * )
      ) // 4. one required parameter, and

    ( optional
      ( prefix,(std::string), "" )
      ( vm, ( po::variables_map const& ), Environment::vm() )
      ( scale,          *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="mesh.scale",_vm=vm) )
      ( straighten,          *( boost::is_integral<mpl::_> ), boption(_prefix=prefix,_name="gmsh.straighten",_vm=vm) )
      ( refine,          *( boost::is_integral<mpl::_> ), ioption(_prefix=prefix,_name="gmsh.refine",_vm=vm) )
      ( update,          *( boost::is_integral<mpl::_> ), 0 )
      ( physical_are_elementary_regions,		   *, boption(_prefix=prefix,_name="gmsh.physical_are_elementary_regions",_vm=vm) )
      ( worldcomm,       (worldcomm_ptr_t), mesh->worldCommPtr() )
      ( respect_partition,	(bool), boption(_prefix=prefix,_name="gmsh.respect_partition",_vm=vm) )
      ( rebuild_partitions,	(bool), boption(_prefix=prefix,_name="gmsh.partition",_vm=vm) )
      ( rebuild_partitions_filename,	*, ""/*filename*/ )
      ( partitions,      *( boost::is_integral<mpl::_> ), worldcomm->globalSize() )
      ( partitioner,     *( boost::is_integral<mpl::_> ), ioption(_prefix=prefix,_name="gmsh.partitioner",_vm=vm) )
      ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
      ( verbose,   (int), ioption(_prefix=prefix,_name="gmsh.verbosity",_vm=vm) )
      )
                         )
#endif

template <typename MeshType>
std::shared_ptr<MeshType>
loadGMSHMeshImpl( std::shared_ptr<MeshType> mesh, std::string const& filename, std::string const& prefix, po::variables_map const& vm,
                  double scale, bool straighten, int refine, size_type update, bool physical_are_elementary_regions, worldcomm_ptr_t const& worldcomm,
                  bool respect_partition, bool rebuild_partitions, std::string const& rebuild_partitions_filename, int partitions,int partitioner, int partition_file, int verbose )
{
    using _mesh_type = MeshType;
    using _mesh_ptrtype = std::shared_ptr<_mesh_type>;
    //typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    //typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

    _mesh_ptrtype _mesh( mesh );
    _mesh->setWorldComm( worldcomm );

    bool allProcessLoadMeshFile = !(rebuild_partitions && partitions > 1) && (refine == 0);

    std::string fnamePartitioned = rebuild_partitions_filename;

    if ( worldcomm->isMasterRank() || allProcessLoadMeshFile )
    {
        std::string filename_with_path = Environment::findFile( filename );
        if ( filename_with_path.empty() )
        {
            std::vector<std::string> plist = Environment::geoPathList();
            std::ostringstream ostr;
            std::for_each( plist.begin(), plist.end(), [&ostr]( std::string s ) { ostr << " - " << s << "\n"; } );
            CHECK( !filename_with_path.empty() ) << "File " << filename << " cannot be found in the following paths list:\n " << ostr.str();
        }

        if ( refine > 0 )
        {
            Gmsh gmsh( _mesh_type::nDim,_mesh_type::nOrder, Environment::worldCommSeqPtr() );
            gmsh.setRefinementLevels( refine );
            gmsh.setNumberOfPartitions( 1/*partitions*/ );
            gmsh.setPartitioner( (GMSH_PARTITIONER)partitioner );
            gmsh.setMshFileByPartition( partition_file );
            gmsh.setVerbosity(verbose);
            filename_with_path = gmsh.refine( filename_with_path, refine );
        }

        ImporterGmsh<_mesh_type> import( filename_with_path, FEELPP_GMSH_FORMAT_VERSION, allProcessLoadMeshFile? worldcomm : Environment::worldCommSeqPtr() );
        fs::path p_fname( filename_with_path );

#if defined (FEELPP_HAS_GMSH_H)
        if ( p_fname.extension() == ".med" ||
             p_fname.extension() == ".bdf" ||
             p_fname.extension() == ".cgns" ||
             p_fname.extension() == ".p3d" ||
             p_fname.extension() == ".mesh"
             )
        {
            tic();
            std::string gmodelName = "feelpp_gmsh_model";
#if defined( FEELPP_HAS_GMSH_API )
            gmsh::model::add( gmodelName );
            // load msh file
            gmsh::open( filename_with_path );
#else
            int status = GmshReaderFactory::instance().at(p_fname.extension().string())( filename_with_path, gmodelName );
            if( status > 1)
                throw std::logic_error( "read  failed: " + filename_with_path );
#endif
            import.setGModelName( gmodelName );
            import.setDeleteGModelAfterUse( true );
            import.setInMemory( true );
            using namespace std::string_literals;
            toc("loadGMSHMesh.reader"s+p_fname.extension().string(), FLAGS_v>0);
        }
#endif // FEELPP_HAS_GMSH_H

        // need to replace physical_region by elementary_region while reading
        if ( physical_are_elementary_regions )
            import.setElementRegionAsPhysicalRegion( physical_are_elementary_regions );
        import.setRespectPartition( respect_partition );
        if ( rebuild_partitions && partitions > 1 )
        {
            _mesh_ptrtype _meshSeq = std::make_shared<_mesh_type>( Environment::worldCommSeqPtr() );
            _meshSeq->accept( import );
            _meshSeq->components().reset();
            _meshSeq->components().set( size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED) );
            _meshSeq->updateForUse();
#if defined(FEELPP_HAS_HDF5)
            using io_t = PartitionIO<_mesh_type>;
            if ( fnamePartitioned.empty() )
                fnamePartitioned = (fs::current_path() / fs::path( filename_with_path ).filename().replace_extension( ".json" )).string();
            else
                fnamePartitioned = fs::path( fnamePartitioned ).replace_extension( ".json" ).string();

            io_t io( fnamePartitioned );
            std::vector<elements_reference_wrapper_t<_mesh_type>> partitionByRange;
            io.write( partitionMesh( _meshSeq, partitions, partitionByRange ) );
#endif
        }
        else
        {
            tic();
            import.setScaling( scale );
            _mesh->accept( import );
            toc("loadGMSHMesh.readmesh", FLAGS_v>0);

            tic();
            _mesh->components().reset();
            _mesh->components().set( update );
            _mesh->updateForUse();
            toc("loadGMSHMesh.update", FLAGS_v>0);
        }
    }
#if defined(FEELPP_HAS_HDF5)
    if ( rebuild_partitions && partitions > 1 )
    {
        mpi::broadcast( worldcomm->globalComm(), fnamePartitioned, worldcomm->masterRank() );
        _mesh->loadHDF5( fnamePartitioned, update, scale );
    }
#endif
    if ( straighten && _mesh_type::nOrder > 1 )
        return straightenMesh( _mesh, worldcomm->subWorldCommPtr() );

    return _mesh;
}

template <typename ... Ts>
auto loadGMSHMesh( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && mesh = args.get(_mesh);
    std::string const& filename = args.get( _filename );
    std::string const& prefix = args.get_else( _prefix, "" );
    po::variables_map const& vm = args.get_else( _vm, Environment::vm() );
    double scale = args.get_else( _scale, doption(_prefix=prefix,_name="mesh.scale",_vm=vm) );
    bool straighten = args.get_else( _straighten, boption(_prefix=prefix,_name="gmsh.straighten",_vm=vm) );
    int refine = args.get_else( _refine, ioption(_prefix=prefix,_name="gmsh.refine",_vm=vm) );
    size_type update = args.get_else( _update, MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    bool physical_are_elementary_regions = args.get_else( _physical_are_elementary_regions, boption(_prefix=prefix,_name="gmsh.physical_are_elementary_regions",_vm=vm) );
    worldcomm_ptr_t worldcomm = args.get_else( _worldcomm, mesh->worldCommPtr() );
    bool respect_partition = args.get_else(_respect_partition, boption(_prefix=prefix,_name="gmsh.respect_partition",_vm=vm) );
    bool rebuild_partitions = args.get_else(_rebuild_partitions, boption(_prefix=prefix,_name="gmsh.partition",_vm=vm) );
    std::string const& rebuild_partitions_filename = args.get_else(_rebuild_partitions_filename, "");
    int partitions = args.get_else(_partitions, worldcomm->globalSize() );
    int partitioner = args.get_else(_partitioner, ioption(_prefix=prefix,_name="gmsh.partitioner",_vm=vm) );
    int partition_file = args.get_else(_partition_file, 0 );
    int verbose = args.get_else(_verbose,ioption(_prefix=prefix,_name="gmsh.verbosity",_vm=vm) );

    using mesh_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(mesh)>>>;
    return loadGMSHMeshImpl( std::shared_ptr<mesh_type>( mesh ), filename, prefix, vm,
                             scale, straighten, refine, update, physical_are_elementary_regions, worldcomm,
                             respect_partition, rebuild_partitions, rebuild_partitions_filename, partitions, partitioner, partition_file, verbose );
}


}
#endif /* FEELPP_LOADGMSHMESH_HPP */
