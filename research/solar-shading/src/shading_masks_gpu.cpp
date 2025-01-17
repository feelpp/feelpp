#include <iostream>




#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>



//#include "shading_bvh_hip_beta.hpp"

#include "shading_mask.hpp"


using namespace Feel;



int main(int argc, char **argv)
{
    using mesh_t = Mesh<Simplex<FEELPP_TOP_DIM,1,3>>;//<Dim,Order,RDim>
    //using mesh_t = Mesh<Simplex<3,1,3>>;

    po::options_description shadingoptions( "Shading mask options" );
    shadingoptions.add_options()
        ( "json-filename", po::value<std::string>()->default_value(""), "json file containing the buildings names" )
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _desc=shadingoptions
                           );

    // Load the json and recover the building names
    auto jsonfile = removeComments( readFromFile( Environment::expand( soption("json-filename" ) ) ) );
    std::istringstream astr( jsonfile );
    json json_buildings = json::parse( astr );

    // Load the mesh
    tic();

#if FEELPP_TOP_DIM==2
    auto mesh = loadMesh( _mesh = new mesh_t, _update=MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES|MESH_GEOMAP_NOT_CACHED );
#else
    auto mesh = loadMesh( _mesh = new mesh_t );
#endif

    toc("Mesh loaded");

     // Compute the shading masks
    ShadingMask<mesh_t> sm(mesh,json_buildings);

    sm.numTypeThread=2; //1: Mode std::async 2: Mode Specx
    sm.QSaveTypeThreadDotON=false;
    sm.numModeTaskUsed=2;

    sm.computeMasksMaster(); 

    std::cout << fmt::format("End Shading mask example\n");


    return 0; 
}
