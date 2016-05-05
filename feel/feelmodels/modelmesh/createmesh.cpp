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
 \file createmesh.cpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-08-11
 */

#include <feel/feelmodels/modelmesh/createmesh.hpp>

#include <feel/feelfilters/loadmesh.hpp>
//#include <feel/feelfilters/geo.hpp>
//#include <feel/feelfilters/creategmshmesh.hpp>
//#include <feel/feelfilters/loadgmshmesh.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename MeshType>
boost::shared_ptr<MeshType>
reloadMesh( std::string const& nameFile, WorldComm const& worldComm, int straighten )
{
    typedef MeshType mesh_type;

    // reload mesh path stored in file
    std::ifstream file( nameFile.c_str() );
    if ( !file )
    {
        CHECK( false ) << "Fail to open the txt file containing path of msh file : " << nameFile << "\n";
    }
    std::string mshfile;
    if ( !( file >> mshfile ) )
    {
        CHECK( false ) << "Fail to read the msh path in file : " << nameFile << "\n";
    }
    file.close();

#if 0
    return loadGMSHMesh(_mesh=new mesh_type,
                        _filename=mshfile,
                        _worldcomm=worldComm,
                        _straighten=straighten,
                        _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
#else
    return loadMesh( _mesh = new mesh_type( worldComm ),
                     _filename = mshfile,
                     _worldcomm = worldComm,
                     //_prefix=model.prefix(),
                     _rebuild_partitions = false,
                     _savehdf5 = 0,
                     _straighten = straighten,
                     _update = MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
#endif

} // reloadFluidMesh()

template <typename MeshType>
void createMeshModel( ModelNumerical& model, boost::shared_ptr<MeshType>& mesh, std::string const& modelMeshRestartFile )
{
    typedef MeshType mesh_type;
    std::string fmpath = ( fs::path( model.rootRepository() ) / fs::path( modelMeshRestartFile /*model.fileNameMeshPath()*/ ) ).string();
    if ( model.doRestart() )
    {
        model.log( "createMeshModel", "", "restart with : " + fmpath );

        if ( !model.restartPath().empty() )
        {
            fmpath = ( fs::path( model.restartPath() ) / fs::path( modelMeshRestartFile /*model.fileNameMeshPath()*/ ) ).string();
        }
        mesh = reloadMesh<mesh_type>( fmpath, model.worldComm() );
    }
    else
    {
        if ( model.hasMshfileStr() )
        {
            std::string rootpath = model.rootRepository();
            std::string mshfileRebuildPartitions = rootpath + "/" + model.prefix() + ".msh";

            model.log( "createMeshModel", "", "load mesh file : " + model.mshfileStr() );
            std::string meshFileExt = fs::path( model.mshfileStr() ).extension().string();
            bool rebuildPartition = boption( _prefix = model.prefix(), _name = "gmsh.partition" );
            if ( rebuildPartition && meshFileExt != ".msh" )
                CHECK( false ) << "Can not rebuild at this time the mesh partitionining with other format than .msh : TODO";
#if 0
            mesh = loadGMSHMesh(_mesh=new mesh_type,
                                _filename=model.mshfileStr(),
                                _worldcomm=model.worldComm(),
                                _prefix=model.prefix(),
                                _rebuild_partitions=rebuildPartition,//model.rebuildMeshPartitions(),
                                _rebuild_partitions_filename=mshfileRebuildPartitions,
                                _partitions=model.worldComm().localSize(),
                                _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK);
#else
            mesh = loadMesh( _mesh = new mesh_type( model.worldComm() ),
                             _filename = model.mshfileStr(),
                             _worldcomm = model.worldComm(),
                             _prefix = model.prefix(),
                             _rebuild_partitions = rebuildPartition,
                             _rebuild_partitions_filename = mshfileRebuildPartitions,
                             _partitions = model.worldComm().localSize(),
                             _savehdf5 = 0,
                             _update = MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
#endif

            if ( rebuildPartition /*model.rebuildMeshPartitions()*/ ) model.setMshfileStr( mshfileRebuildPartitions );
        }
        else if ( model.hasGeofileStr() )
        {
            std::string path = model.rootRepository();
            std::string mshfile = path + "/" + model.prefix() + ".msh";
            model.setMshfileStr( mshfile );

            fs::path curPath = fs::current_path();
            bool hasChangedRep = false;
            if ( curPath != fs::path( model.rootRepository() ) )
            {
                model.log( "createMeshModel", "", "change repository (temporary) for build mesh from geo : " + model.rootRepository() );
                bool hasChangedRep = true;
                Environment::changeRepository( _directory = boost::format( model.rootRepository() ), _subdir = false );
            }

            gmsh_ptrtype geodesc = geo( _filename = model.geofileStr(),
                                        _prefix = model.prefix(),
                                        _worldcomm = model.worldComm() );
            // allow to have a geo and msh file with a filename equal to prefix
            geodesc->setPrefix( model.prefix() );
            mesh = createGMSHMesh( _mesh = new mesh_type, _desc = geodesc,
                                   _prefix = model.prefix(), _worldcomm = model.worldComm(),
                                   _partitions = model.worldComm().localSize() );

            // go back to previous repository
            if ( hasChangedRep )
                Environment::changeRepository( _directory = boost::format( curPath.string() ), _subdir = false );
        }
#if 0
        else
        {
            std::string geotoolSavePath;
            if ( model.geotoolSaveDirectory()!=model.appliShortRepository() )
            {
                model.log("createMeshModel","", "change rep -> "+ model.geotoolSaveDirectory() );
                Environment::changeRepository( _directory=boost::format(model.geotoolSaveDirectory()), _subdir=false );
                geotoolSavePath = Environment::rootRepository()+"/"+ model.geotoolSaveDirectory();
            }
            else
            {
                geotoolSavePath = model.rootRepository();
            }

            std::string geotoolSaveName = model.geotoolSaveName();
            std::string mshfile = geotoolSavePath + "/" + geotoolSaveName + ".msh";
            std::string geofilename = geotoolSavePath + "/" + geotoolSaveName;// without .geo
            model.setMshfileStr(mshfile);

            model.loadConfigMeshFile(geofilename);

            if ( model.geotoolSaveDirectory()!=model.appliShortRepository() )
            {
                model.log("createMeshModel","", "change rep -> " + model.rootRepository() );
                Environment::changeRepository( _directory=boost::format(model.appliShortRepository()), _subdir=true );
            }
        }
#endif
        model.saveMSHfilePath( fmpath );
    } //not restart
}

template boost::shared_ptr<Mesh<Simplex<2, 1>>> reloadMesh<Mesh<Simplex<2, 1>>>( std::string const&, WorldComm const&, int );
template boost::shared_ptr<Mesh<Simplex<3, 1>>> reloadMesh<Mesh<Simplex<3, 1>>>( std::string const&, WorldComm const&, int );
template boost::shared_ptr<Mesh<Simplex<1, 1, 2>>> reloadMesh<Mesh<Simplex<1, 1, 2>>>( std::string const&, WorldComm const&, int );
template boost::shared_ptr<Mesh<Simplex<1, 1, 3>>> reloadMesh<Mesh<Simplex<1, 1, 3>>>( std::string const&, WorldComm const&, int );
template void createMeshModel<Mesh<Simplex<2, 1>>>( ModelNumerical&, boost::shared_ptr<Mesh<Simplex<2, 1>>>&, std::string const& );
template void createMeshModel<Mesh<Simplex<3, 1>>>( ModelNumerical&, boost::shared_ptr<Mesh<Simplex<3, 1>>>&, std::string const& );
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
template boost::shared_ptr<Mesh<Simplex<2, 2>>> reloadMesh<Mesh<Simplex<2, 2>>>( std::string const&, WorldComm const&, int );
template boost::shared_ptr<Mesh<Simplex<3, 2>>> reloadMesh<Mesh<Simplex<3, 2>>>( std::string const&, WorldComm const&, int );
template void createMeshModel<Mesh<Simplex<2, 2>>>( ModelNumerical&, boost::shared_ptr<Mesh<Simplex<2, 2>>>&, std::string const& );
template void createMeshModel<Mesh<Simplex<3, 2>>>( ModelNumerical&, boost::shared_ptr<Mesh<Simplex<3, 2>>>&, std::string const& );
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 3 )
template boost::shared_ptr<Mesh<Simplex<2, 3>>> reloadMesh<Mesh<Simplex<2, 3>>>( std::string const&, WorldComm const&, int );
template boost::shared_ptr<Mesh<Simplex<3, 3>>> reloadMesh<Mesh<Simplex<3, 3>>>( std::string const&, WorldComm const&, int );
template void createMeshModel<Mesh<Simplex<2, 3>>>( ModelNumerical&, boost::shared_ptr<Mesh<Simplex<2, 3>>>&, std::string const& );
template void createMeshModel<Mesh<Simplex<3, 3>>>( ModelNumerical&, boost::shared_ptr<Mesh<Simplex<3, 3>>>&, std::string const& );
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 4 )
template boost::shared_ptr<Mesh<Simplex<2, 4>>> reloadMesh<Mesh<Simplex<2, 4>>>( std::string const&, WorldComm const&, int );
template boost::shared_ptr<Mesh<Simplex<3, 4>>> reloadMesh<Mesh<Simplex<3, 4>>>( std::string const&, WorldComm const&, int );
template void createMeshModel<Mesh<Simplex<2, 4>>>( ModelNumerical&, boost::shared_ptr<Mesh<Simplex<2, 4>>>&, std::string const& );
template void createMeshModel<Mesh<Simplex<3, 4>>>( ModelNumerical&, boost::shared_ptr<Mesh<Simplex<3, 4>>>&, std::string const& );
#endif

} // namespace FeelModels
} // namespace Feel
