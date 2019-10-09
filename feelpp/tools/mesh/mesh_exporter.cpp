#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    using Feel::cout;
    	po::options_description meshpartoptions( "Mesh Partitioner options" );
	meshpartoptions.add_options()
        ( "dim", po::value<int>()->default_value( 3 ), "mesh dimension" )
        ( "shape", po::value<std::string>()->default_value( "simplex" ), "mesh dimension" )
        ( "order", po::value<int>()->default_value( 1 ), "mesh geometric order" );
    Environment env( _argc=argc, _argv=argv,
                     _desc=meshpartoptions,
                     _about=about( _name="mesh_exporter" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" )
                     );
    
    using mesh_t = Mesh<Simplex<3, 1>>;
    auto mesh = loadMesh( _mesh = new mesh_t );
    auto ex = exporter( _mesh=mesh );
    ex->addRegions();

    ex->save();
    return 0;

}

