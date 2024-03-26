/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2016-07-18

  Copyright (C) 2016 Feel++ Consortium

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
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
// #include <feel/feelmesh/partitionmesh.hpp>
// #include <feel/feelfilters/partitionio.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
//#include <feel/feelfilters/exporter.hpp>

using namespace Feel;


std::string
getSubMeshOutputPath( fs::path const& inputPathMesh, std::string const& submeshtype )
{
    fs::path inputDir = inputPathMesh.parent_path();
    std::string inputFilenameWithoutExt = inputPathMesh.stem().string();

    std::string meshFormat = soption(_name="format");
    std::string meshExtension = ( meshFormat == "gmsh" )? ".msh" : ".json";
    std::string outputFilenameWithoutExt = (boost::format("%1%_%2%")%inputFilenameWithoutExt % submeshtype ).str();
    std::string outputFilenameWithExt = outputFilenameWithoutExt + meshExtension;

    std::string outputPathMesh = ( fs::current_path() / fs::path(outputFilenameWithExt) ).string();
    if ( Environment::vm().count("ofile") )
    {
        outputPathMesh = fs::canonical( soption("ofile") ).string();
        if ( true )
        {
            fs::path givenPath = fs::path(outputPathMesh);
            std::string newfilename = (boost::format("%1%_%2%%3%")%givenPath.stem().string() % submeshtype %meshExtension ).str();
            outputPathMesh = (givenPath.parent_path() / fs::path(newfilename)).string();
        }
        CHECK( fs::path(outputPathMesh).extension() == ".json" || fs::path(outputPathMesh).extension() == ".msh" ) << "invalid mesh extension";
    }
    else if ( Environment::vm().count("odir") )
    {
        fs::path odir = fs::canonical( soption("odir") ).string();
        outputPathMesh = (odir / fs::path(outputFilenameWithExt) ).string();
    }
    return outputPathMesh;
}

template <typename MeshType>
void
runCreateSubmeshAndSaveElement( MeshType const& mesh, std::list<std::string> const& markers,
                                std::string const& outputPathMesh, mpl::false_ ) {}
template <typename MeshType>
void
runCreateSubmeshAndSaveElement( std::shared_ptr<MeshType> const& mesh, std::list<std::string> const& markers,
                                std::string const& outputPathMesh, mpl::true_ )
{
    Feel::cout << "-----------------------------\n"
               << "extract submesh of elements\n";
    Feel::cout << "markers :";
    for (std::string const& marker : markers )
        Feel::cout << " " << marker;
    Feel::cout << "\n";

    tic();
    auto submesh = createSubmesh( _mesh=mesh, _range=markedelements(mesh,markers) );
    toc("extract submesh",true);

    Feel::cout << "output mesh path : " << outputPathMesh << "\n";
    tic();
    if ( fs::path(outputPathMesh).extension() == ".json" )
    {
        using io_t = PartitionIO<mesh_t<decltype(submesh)>>;
        io_t io( outputPathMesh );
        io.write( submesh );
    }
    else if ( fs::path(outputPathMesh).extension() == ".msh" )
    {
        saveGMSHMesh(_mesh=submesh,_filename=outputPathMesh );
    }
    toc("save markedelements submesh on disk done",true);
}
template <typename MeshType>
void
runCreateSubmeshAndSaveFace( MeshType const& mesh, std::list<std::string> const& markers,
                             bool extractBoundaryFaces, std::string const& outputPathMesh, mpl::false_ ) {}
template <typename MeshType>
void
runCreateSubmeshAndSaveFace( std::shared_ptr<MeshType> const& mesh, std::list<std::string> const& markers,
                             bool extractBoundaryFaces, std::string const& outputPathMesh, mpl::true_ )
{
    Feel::cout << "-----------------------------\n"
               << "extract submesh of faces\n";
    if ( extractBoundaryFaces )
        Feel::cout << "range : boundary faces\n";
    else
    {
        Feel::cout << "markers :";
        for (std::string const& marker : markers )
            Feel::cout << " " << marker;
        Feel::cout << "\n";
    }

    tic();
    trace_mesh_ptr_t<MeshType> submesh;
    if ( extractBoundaryFaces )
    {
        submesh = createSubmesh( _mesh=mesh, _range=boundaryfaces(mesh), _only_on_boundary_faces=true );
    }
    else
    {
        submesh = createSubmesh( _mesh=mesh, _range=markedfaces(mesh,markers), _update=0 );
    }
    toc("extract submesh",true);

    Feel::cout << "output mesh path : " << outputPathMesh << "\n";

    tic();
    if ( fs::path(outputPathMesh).extension() == ".json" )
    {
        using io_t = PartitionIO<mesh_t<decltype(submesh)>>;
        io_t io( outputPathMesh );
        io.write( submesh );
    }
    else if ( fs::path(outputPathMesh).extension() == ".msh" )
    {
        saveGMSHMesh(_mesh=submesh,_filename=outputPathMesh );
    }
    toc("save face submesh on disk done",true);
}
template <typename MeshType>
void
runCreateSubmeshAndSaveEdge( MeshType const& mesh, std::list<std::string> const& markers,
                             std::string const& outputPathMesh, mpl::false_ ) {}
template <typename MeshType>
void
runCreateSubmeshAndSaveEdge( std::shared_ptr<MeshType> const& mesh, std::list<std::string> const& markers,
                             std::string const& outputPathMesh, mpl::true_ )
{
    Feel::cout << "-----------------------------\n"
               << "extract submesh of edges\n";
    Feel::cout << "markers :";
    for (std::string const& marker : markers )
        Feel::cout << " " << marker;
    Feel::cout << "\n";

    tic();
    auto submesh = createSubmesh( _mesh=mesh, _range=markededges(mesh,markers), _update=0 );
    toc("extract submesh",true);

    Feel::cout << "output mesh path : " << outputPathMesh << "\n";
    tic();
    if ( fs::path(outputPathMesh).extension() == ".json" )
    {
        using io_t = PartitionIO<mesh_t<decltype(submesh)>>;
        io_t io( outputPathMesh );
        io.write( submesh );
    }
    else if ( fs::path(outputPathMesh).extension() == ".msh" )
    {
        saveGMSHMesh(_mesh=submesh,_filename=outputPathMesh );
    }
    toc("save face submesh on disk done",true);
}

template <typename ShapeType>
void
run( std::vector<std::string> const& markers, bool extractBoundaryFaces )
{
    typedef Mesh<ShapeType> mesh_type;

    fs::path inputPathMesh = fs::canonical( soption("ifile") );
    tic();
    auto mesh = loadMesh(_mesh=new mesh_type(Environment::worldCommPtr()), _savehdf5=0,
                         _filename=inputPathMesh.string(),
                         //_update=size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES));
                         _update=size_type(MESH_UPDATE_FACES|MESH_UPDATE_EDGES));
    toc("load parent mesh",true);

    std::list<std::string> markerListElt, markerListFace, markerListEdge, markerListPoint;
    std::set<std::string> markerDone;
    for ( std::string const& marker : markers )
    {
        if ( markerDone.find( marker ) != markerDone.end() )
            continue;
        if ( mesh->hasPointMarker( marker ) )
            markerListPoint.push_back( marker );
        else if ( mesh->hasEdgeMarker( marker ) )
            markerListEdge.push_back( marker );
        else if ( mesh->hasFaceMarker( marker ) )
            markerListFace.push_back( marker );
        else if ( mesh->hasMarker( marker ) )
            markerListElt.push_back( marker );
        else
        {
            if ( Environment::isMasterRank() )
                std::cout << "marker " << marker << " is not in mesh\n";
        }
        markerDone.insert( marker );
    }


    std::string outputPathMeshElts = getSubMeshOutputPath( inputPathMesh, "elements" );
    std::string outputPathMeshFaces = (extractBoundaryFaces)? getSubMeshOutputPath( inputPathMesh, "boundaryfaces" ) : getSubMeshOutputPath( inputPathMesh, "faces" );
    std::string outputPathMeshEdges = getSubMeshOutputPath( inputPathMesh, "edges" );
    std::string outputPathMeshPoints = getSubMeshOutputPath( inputPathMesh, "points" );

    fs::path outputDir = fs::path(outputPathMeshElts).parent_path();
    if ( Environment::isMasterRank() && !fs::exists( outputDir ) )
        fs::create_directories( outputDir );
    Environment::worldComm().barrier();

    if ( !markerListElt.empty() )
    {
        static const bool meshok = mesh_type::nDim > 0;
        runCreateSubmeshAndSaveElement( mesh, markerListElt, outputPathMeshElts, mpl::bool_<meshok>() );
    }
    if ( !markerListFace.empty() || extractBoundaryFaces )
    {
        static const bool meshok = mesh_type::nDim > 1;
        runCreateSubmeshAndSaveFace( mesh, markerListFace, extractBoundaryFaces, outputPathMeshFaces , mpl::bool_<meshok>() );
    }
    if ( !markerListEdge.empty() )
    {
        static const bool meshok = mesh_type::nDim == 3;
        runCreateSubmeshAndSaveEdge( mesh, markerListEdge, outputPathMeshEdges, mpl::bool_<meshok>() );
    }
    if ( !markerListPoint.empty() )
    {
        CHECK( false ) << "implement Mesh 0d";
    }

}

int main( int argc, char** argv )
{
	po::options_description submeshoptions( "Extract SubMesh options" );
	submeshoptions.add_options()
        ( "dim", po::value<int>()->default_value( 3 ), "mesh dimension" )
        ( "realdim", po::value<int>(), "mesh real dimension [optional (default is value of dim)]" )
        ( "shape", po::value<std::string>()->default_value( "simplex" ), "mesh entity" )
        ( "ifile", po::value<std::string>(), "input mesh filename" )
        ( "odir", po::value<std::string>(), "output directory [optional (default is current dir)]" )
        ( "ofile", po::value<std::string>(), "output mesh filename [optional (default name deduced from input mesh filename)]" )
        ( "format", po::value<std::string>()->default_value( "json+h5" ), "output mesh format : json+h5,gmsh" )
        ( "markers", po::value<std::vector<std::string> >()->multitoken(), "markers to extract" )
        ( "boundaryfaces", po::value<bool>()->default_value( false ), "extract boundaryfaces (ignore markers if given)" )
		;

    Environment env( _argc=argc, _argv=argv,
                     _desc=submeshoptions,
                     _about=about( _name="mesh_submesh" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ),
                     _directory=".",_subdir=false );

    // init/check config from options
    int dim = ioption(_name="dim");
    std::string shape = soption(_name="shape");

    int nRealDim = dim;
    if( Environment::vm().count("realdim") )
        nRealDim = ioption(_name="realdim");
    if ( nRealDim < dim )
    {
        std::cout << "do nothing because invalid realdim"<< nRealDim << "\n";
        return 0;
    }
    std::vector<std::string> markers;
    if ( Environment::vm().count("markers"))
        markers = Environment::vm()["markers"].as<std::vector<std::string> >();

    bool extractBoundaryFaces = boption(_name="boundaryfaces");
    if ( markers.empty() && !extractBoundaryFaces )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because no marker (use --markers mark1 mark2)\n";
        return 0;
    }

    if ( !Environment::vm().count("ifile") )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because no input mesh (use --ifile mymesh.msh)\n";
        return 0;
    }

    fs::path pathInputMesh = fs::canonical( soption("ifile") );
    if ( !fs::exists( pathInputMesh ) )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --ifile : " << pathInputMesh.string() << " does not exist\n";
        return 0;
    }

    // run application by template type selection from the mesh shape
    if ( dim == 1 )
    {
        switch ( nRealDim )
        {
        case 1 : run<Simplex<1,1,1>>( markers,extractBoundaryFaces );break;
        case 2 : run<Simplex<1,1,2>>( markers,extractBoundaryFaces );break;
        case 3 : run<Simplex<1,1,3>>( markers,extractBoundaryFaces );break;
        }
    }
    else
    {
        if ( shape == "simplex" )
        {
            switch ( dim )
            {
            case 2 : ( nRealDim==2 )? run<Simplex<2,1,2>>( markers,extractBoundaryFaces ) : run<Simplex<2,1,3>>( markers,extractBoundaryFaces ); break;
            case 3 : run<Simplex<3,1>>( markers,extractBoundaryFaces ); break;
            }
        }
        else if ( shape == "hypercube" )
        {
            switch ( dim )
            {
            case 2 : ( nRealDim==2 )? run<Hypercube<2,1,2>>( markers,extractBoundaryFaces ) : run<Hypercube<2,1,3>>( markers,extractBoundaryFaces ); break;
            case 3 : run<Hypercube<3>>( markers,extractBoundaryFaces ); break;
            }
        }
    }

    return 0;
}
