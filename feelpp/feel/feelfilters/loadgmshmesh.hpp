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
    loadGMSHMesh,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( mesh, * )
      ( filename, * )
      ) // 4. one required parameter, and

    ( optional
      ( prefix,(std::string), "" )
      ( scale,          *( boost::is_arithmetic<mpl::_> ), doption(_prefix=prefix,_name="gmsh.scale") )
      ( straighten,          *( boost::is_integral<mpl::_> ), boption(_prefix=prefix,_name="gmsh.straighten") )
      ( refine,          *( boost::is_integral<mpl::_> ), ioption(_prefix=prefix,_name="gmsh.refine") )
      ( update,          *( boost::is_integral<mpl::_> ), 0 )
      ( physical_are_elementary_regions,		   *, boption(_prefix=prefix,_name="gmsh.physical_are_elementary_regions") )
      ( worldcomm,       (worldcomm_ptr_t), Environment::worldCommPtr() )
      ( respect_partition,	(bool), boption(_prefix=prefix,_name="gmsh.respect_partition") )
      ( rebuild_partitions,	(bool), boption(_prefix=prefix,_name="gmsh.partition") )
      ( rebuild_partitions_filename,	*, filename )
      ( partitions,      *( boost::is_integral<mpl::_> ), worldcomm->globalSize() )
      ( partitioner,     *( boost::is_integral<mpl::_> ), ioption(_prefix=prefix,_name="gmsh.partitioner") )
      ( partition_file,   *( boost::is_integral<mpl::_> ), 0 )
      ( verbose,   (int), ioption(_prefix=prefix,_name="gmsh.verbosity") )
      )
                         )
{
    typedef typename Feel::detail::mesh<Args>::type _mesh_type;
    typedef typename Feel::detail::mesh<Args>::ptrtype _mesh_ptrtype;

    _mesh_ptrtype _mesh( mesh );
    _mesh->setWorldComm( worldcomm );

    std::string filename_with_path = Environment::findFile( filename );
    if ( filename_with_path.empty() )
    {
        std::vector<std::string> plist = Environment::geoPathList();
        std::ostringstream ostr;
        std::for_each( plist.begin(), plist.end(), [&ostr]( std::string s ) { ostr << " - " << s << "\n"; } );
        CHECK( !filename_with_path.empty() ) << "File " << filename << " cannot be found in the following paths list:\n " << ostr.str();
    }

    if ( refine || ( rebuild_partitions && partitions > 1 ) )
    {
        Gmsh gmsh( _mesh_type::nDim,_mesh_type::nOrder, worldcomm );
        gmsh.setRefinementLevels( refine );
        gmsh.setNumberOfPartitions( partitions );
        gmsh.setPartitioner( (GMSH_PARTITIONER)partitioner );
        gmsh.setMshFileByPartition( partition_file );
        gmsh.setVerbosity(verbose);

        // refinement if option is enabled to a value greater or equal to 1
        if ( refine )
        {
            filename_with_path = gmsh.refine( filename_with_path, refine );
        }
        if ( rebuild_partitions && partitions > 1 )
        {
            std::string fnamePartitioned = rebuild_partitions_filename;
            if ( fnamePartitioned.empty() )
                fnamePartitioned = filename_with_path;
            gmsh.rebuildPartitionMsh(filename_with_path,rebuild_partitions_filename);
            filename_with_path = rebuild_partitions_filename;
        }
    }

    ImporterGmsh<_mesh_type> import( filename_with_path, FEELPP_GMSH_FORMAT_VERSION, worldcomm );
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
        CHECK( false ) << "NOT IMPLEMENTED";
#else
        int status = GmshReaderFactory::instance().at(p_fname.extension().string())( filename_with_path, gmodelName );
        if( status > 1)
            throw std::logic_error( "read  failed: " + filename_with_path );
#endif
        import.setGModelName( gmodelName );
        import.setDeleteGModelAfterUse( true );
        import.setInMemory(true);
        using namespace std::string_literals;
        toc("loadGMSHMesh.reader"s+p_fname.extension().string(), FLAGS_v>0);
    }
#endif // FEELPP_HAS_GMSH_H

    // need to replace physical_region by elementary_region while reading
    if ( physical_are_elementary_regions )
        import.setElementRegionAsPhysicalRegion( physical_are_elementary_regions );
    import.setScaling( scale );
    import.setRespectPartition( respect_partition );
    tic();
    _mesh->accept( import );
    toc("loadGMSHMesh.readmesh", FLAGS_v>0);

    tic();
    _mesh->components().reset();
    _mesh->components().set( update );
    _mesh->updateForUse();
    toc("loadGMSHMesh.update", FLAGS_v>0);

    if ( straighten && _mesh_type::nOrder > 1 )
        return straightenMesh( _mesh, worldcomm->subWorldCommPtr() );

    return _mesh;
}

}
#endif /* FEELPP_LOADGMSHMESH_HPP */
