/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2011-08-11

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file fsimesh.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2011-08-11
 */

#include <feel/feelmodels/modelmesh/fsimesh.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelmesh/partitionmesh.hpp>


namespace Feel
{
namespace FeelModels
{


template< class ConvexType >
FSIMesh<ConvexType>::FSIMesh( std::string prefix, worldcomm_ptr_t const& worldcomm )
    :
    M_prefix( prefix ),
    M_worldComm( worldcomm ),
    M_meshSize( 0.1 ),
    M_nPartitions( this->worldComm().localSize() ),
    M_partitioner( ioption(_name="gmsh.partitioner") ),
    M_forceRebuild( false )
{}

//-------------------------------------------------------------------------------------------------//

template< class ConvexType >
void
FSIMesh<ConvexType>::buildFSIMeshFromMsh()
{
#if 0
    if ( this->worldComm().isMasterRank() )
        std::cout << "[FSIMesh] : buildFSIMeshFromMsh start" << std::endl;

    CHECK( !this->mshPathFSI().empty() ) << "mshPathFSI must be specified";

    if ( this->worldComm().isMasterRank() &&
         ( !fs::exists( M_mshfilepathFluidPart1 ) || !fs::exists( M_mshfilepathSolidPart1 ) || this->forceRebuild() ) )
    {
        std::cout << "[FSIMesh] : load fsi mesh ....\n";
        auto meshFSI = loadMesh( _mesh=new mesh_type(this->worldComm().subWorldCommSeqPtr()),
                                 _prefix=this->prefix(),
                                 _filename=this->mshPathFSI().string(),
                                 _straighten=false,
                                 _worldcomm=this->worldComm().subWorldCommSeqPtr(),
                                 _rebuild_partitions=false,
                                 _update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_UPDATE_EDGES)
                                 );
        std::cout << "[FSIMesh] : load fsi mesh finish\n";

        this->buildSubMesh( meshFSI );
    }
    // wait writing meshes done
    this->worldComm().globalComm().barrier();

    this->buildMeshesPartitioning();

    if ( this->worldComm().isMasterRank() )
        std::cout << "[FSIMesh] : buildFSIMeshFromMsh finish" << std::endl;
#endif
}

//-------------------------------------------------------------------------------------------------//

template< class ConvexType >
void
FSIMesh<ConvexType>::buildFSIMeshFromGeo()
{
#if 0
    if ( this->worldComm().isMasterRank() )
        std::cout << "[FSIMesh] : buildFSIMeshFromGeo start" << std::endl;

    CHECK( !this->mshPathFSI().empty() ) << "mshPathFSI must be specified";
    fs::path meshesdirectories = this->mshPathFSI().parent_path();

#if 0
    // go in output path directory for launch loadMesh
    fs::path curPath=fs::current_path();
    bool hasChangedRep=false;
    if ( curPath != meshesdirectories )
    {
        if ( this->worldComm().isMasterRank() )
            std::cout << "[FSIMesh] change repository (temporary) for build mesh from geo : " << meshesdirectories.string() << "\n";
        bool hasChangedRep=true;
        Environment::changeRepository( _directory=boost::format(meshesdirectories.string()), _subdir=false,
                                       _worldcomm=this->worldComm() );
    }
#endif

    if ( this->worldComm().isMasterRank() &&
         (!fs::exists( M_mshfilepathFluidPart1 ) || !fs::exists( M_mshfilepathSolidPart1 ) || this->forceRebuild() ) )
    {
        if ( !fs::exists( meshesdirectories ) ) fs::create_directories( meshesdirectories );



#if defined(FEELPP_HAS_GMSH_H)
        auto/*gmsh_ptrtype*/ geodesc = geo( _filename=this->geoPathFSI().string(),
                                            _prefix=this->prefix(),
                                            _worldcomm=this->worldComm().subWorldCommSeqPtr() );
        // allow to have a geo and msh file with a filename equal to prefix
        geodesc->setPrefix(this->mshPathFSI().stem().string());//this->prefix());
        auto meshFSI = createGMSHMesh(_mesh=new mesh_type(this->worldComm().subWorldCommSeqPtr()),
                                      _desc=geodesc,
                                      _prefix=this->prefix(),
                                      _worldcomm=this->worldComm().subWorldCommSeqPtr(),
                                      _h=this->meshSize(),
                                      _straighten=false,
                                      _partitions=1/*this->worldComm().subWorldCommSeq().localSize()*/,
                                      _rebuild_partitions=false,
                                      _update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_UPDATE_EDGES),
                                      _directory=meshesdirectories.string() );
#else
        
#endif
#if 0
        
        // copy geofile
        std::string nameMeshFile = this->mshPathFSI().stem().string() + ".geo";
        fs::path geoPathFSIcopy = meshesdirectories / fs::path(nameMeshFile);
#if BOOST_FILESYSTEM_VERSION == 3
        boost::system::error_code ec;
        fs::copy_file( this->geoPathFSI(), geoPathFSIcopy, fs::copy_option::overwrite_if_exists, ec );
#elif BOOST_FILESYSTEM_VERSION == 2
        fs::copy_file( this->geoPathFSI(), geoPathFSIcopy, fs::copy_option::overwrite_if_exists );
#endif
        std::cout << "[FSIMesh] : build fsi mesh ....\n";
        auto meshFSI = loadMesh( _mesh=new mesh_type(this->worldComm().subWorldCommSeq()),
                                 _filename=geoPathFSIcopy.string(),
                                 _h=this->meshSize(),
                                 _straighten=false,
                                 _worldcomm=this->worldComm().subWorldCommSeq(),
                                 _rebuild_partitions=false,
                                 _update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_UPDATE_EDGES)
                                 );
#endif
        std::cout << "[FSIMesh] : build fsi mesh finish\n";
        CHECK( fs::exists( this->mshPathFSI() ) ) << "mesh file not exist : " << this->mshPathFSI();
        this->buildSubMesh( meshFSI );
    }

#if 0
        // go back to previous repository
        if ( hasChangedRep )
            Environment::changeRepository( _directory=boost::format(curPath.string()), _subdir=false,
                                           _worldcomm=this->worldComm() );
#endif

    // wait writing meshes done
    this->worldComm().globalComm().barrier();

    this->buildMeshesPartitioning();

    if ( this->worldComm().isMasterRank() )
        std::cout << "[FSIMesh] : buildFSIMeshFromGeo finish" << std::endl;
#endif
} // buildFSImesh


//-------------------------------------------------------------------------------------------------//

template< class ConvexType >
void
FSIMesh<ConvexType>::buildSubMesh( mesh_ptrtype const& fsimesh )
{
    std::cout << "[FSIMesh] : build fluid and structure submesh ......\n";

    // create fluid submesh
    std::list<std::string> myFluidMarkers( this->markersNameFluidVolume().begin(),this->markersNameFluidVolume().end() );
    auto fluidMesh_temp = createSubmesh(_mesh=fsimesh,_range=markedelements(fsimesh,myFluidMarkers),
                                        _context=size_type(EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT),
                                        _update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_UPDATE_EDGES));
    std::cout << "fluid mesh -> numGlobalElements : " << fluidMesh_temp->numGlobalElements() << "\n";
    std::cout << "fluid mesh -> save sequential mesh in : " << this->mshPathFluidPart1() << "\n";
    fs::path meshesdirectoriesFluid = this->mshPathFluidPart1().parent_path();
    if ( !fs::exists( meshesdirectoriesFluid ) ) fs::create_directories( meshesdirectoriesFluid );
#ifdef FEELPP_HAS_HDF5
    PartitionIO<mesh_fluid_type> iofluid( this->mshPathFluidPart1().string() );
    iofluid.write( fluidMesh_temp );
#else
    saveGMSHMesh(_mesh=fluidMesh_temp,_filename=this->mshPathFluidPart1().string());
#endif
    // create struct submesh
    std::list<std::string> mySolidMarkers( this->markersNameSolidVolume().begin(),this->markersNameSolidVolume().end() );
    auto solidMesh_temp = createSubmesh(_mesh=fsimesh,_range=markedelements(fsimesh,mySolidMarkers),
                                        _context=size_type(EXTRACTION_KEEP_MARKERNAMES_ONLY_PRESENT),
                                        _update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_UPDATE_EDGES));
    std::cout << "solid mesh -> numGlobalElements : " << solidMesh_temp->numGlobalElements()<<"\n";
    std::cout << "solid mesh -> save sequential mesh in : " << this->mshPathSolidPart1() << "\n";
    fs::path meshesdirectoriesSolid = this->mshPathSolidPart1().parent_path();
    if ( !fs::exists( meshesdirectoriesSolid ) ) fs::create_directories( meshesdirectoriesSolid );
#ifdef FEELPP_HAS_HDF5
    PartitionIO<mesh_solid_type> iosolid( this->mshPathSolidPart1().string() );
    iosolid.write( solidMesh_temp );
#else
    saveGMSHMesh(_mesh=solidMesh_temp,_filename=this->mshPathSolidPart1().string());
#endif
    std::cout << "[FSIMesh] : build fluid and structure submesh finish\n";
}

//-------------------------------------------------------------------------------------------------//

template< class ConvexType >
void
FSIMesh<ConvexType>::buildMeshesPartitioning()
{
#ifdef FEELPP_HAS_HDF5
    if ( this->worldComm().isMasterRank() )
    {
        if ( !fs::exists( this->mshPathFluidPartN() ) || this->forceRebuild() )
        {
            std::cout << "partitioning fluid submesh in " << this->nPartitions() << " part......\n"
                      << "From : " << this->mshPathFluidPart1() << "\n"
                      << "Write : " << this->mshPathFluidPartN() <<"\n";
            auto fluidmeshSeq = loadMesh(_mesh=new mesh_type(this->worldComm().subWorldCommSeqPtr()), _savehdf5=0,
                                         _filename=this->mshPathFluidPart1().string(),
                                         _update=size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES),
                                         _straighten=false );
            PartitionIO<mesh_fluid_type> iofluid( this->mshPathFluidPartN().string() );
            iofluid.write( partitionMesh( fluidmeshSeq,  this->nPartitions() ) );
        }

        if ( !fs::exists( this->mshPathSolidPartN() ) || this->forceRebuild() )
        {
            std::cout << "partitioning structure submesh in " << this->nPartitions() << " part......\n"
                      << "From : " << this->mshPathSolidPart1() << "\n"
                      << "Write : " << this->mshPathSolidPartN() <<"\n";
            auto solidmeshSeq = loadMesh(_mesh=new mesh_type(this->worldComm().subWorldCommSeqPtr()), _savehdf5=0,
                                         _filename=this->mshPathSolidPart1().string(),
                                         _update=size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES),
                                         _straighten=false );
            PartitionIO<mesh_solid_type> iosolid( this->mshPathSolidPartN().string() );
            iosolid.write( partitionMesh( solidmeshSeq,  this->nPartitions() ) );
        }
    }
    this->worldComm().globalComm().barrier();
#else
    // partitioning if no exist
    if ( !fs::exists( this->mshPathFluidPartN() ) || this->forceRebuild() )
    {
        this->worldComm().globalComm().barrier();

        if ( this->worldComm().isMasterRank() )
            std::cout << "partitioning fluid submesh in " << this->nPartitions() << " part......\n"
                      << "From : " << this->mshPathFluidPart1() << "\n"
                      << "Write : " << this->mshPathFluidPartN() <<"\n";

        // partioning mesh base
        Gmsh gmsh( mesh_type::nDim,
                   mesh_type::nOrder,
                   this->worldComm() );
        gmsh.setNumberOfPartitions( this->nPartitions() );
        gmsh.setPartitioner( (GMSH_PARTITIONER)this->partitioner() );
        gmsh.rebuildPartitionMsh( this->mshPathFluidPart1().string(), this->mshPathFluidPartN().string() );
    }
    if ( !fs::exists( this->mshPathSolidPartN() ) || this->forceRebuild() )
    {
        this->worldComm().globalComm().barrier();
        if ( this->worldComm().isMasterRank() )
            std::cout << "partitioning structure submesh in " << this->nPartitions() << " part......\n"
                      << "From : " << this->mshPathSolidPart1() << "\n"
                      << "Write : " << this->mshPathSolidPartN() <<"\n";
        // partioning mesh base
        Gmsh gmsh( mesh_type::nDim,
                   mesh_type::nOrder,
                   this->worldComm() );
        gmsh.setNumberOfPartitions( this->nPartitions() );
        gmsh.setPartitioner( (GMSH_PARTITIONER)this->partitioner() );
        gmsh.rebuildPartitionMsh( this->mshPathSolidPart1().string(), this->mshPathSolidPartN().string() );
    }
#endif
}

//-------------------------------------------------------------------------------------------------//

template class FSIMesh< Simplex<2,1> >;
template class FSIMesh< Simplex<3,1> >;
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
template class FSIMesh< Simplex<2,2> >;
template class FSIMesh< Simplex<3,2> >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 3 )
template class FSIMesh< Simplex<2,3> >;
template class FSIMesh< Simplex<3,3> >;
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 4 )
template class FSIMesh< Simplex<2,4> >;
template class FSIMesh< Simplex<3,4> >;
#endif
#if 0 // miss straigten instantiation
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 5 )
template class FSIMesh< Simplex<2,5> >;
#endif
#endif


} // namespace FeelModels
} // namespace Feel

