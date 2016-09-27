#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    using namespace Feel::vf;
    // initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="mymesh" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );
    // create a mesh with GMSH using Feel++ geometry tool
    auto mesh = loadMesh(_mesh=new  Mesh<Simplex<2>>);

    LOG(INFO) << "mesh " << soption(_name="gmsh.filename") << " loaded";

    LOG(INFO) << "volume =" << integrate( _range=elements( mesh ), _expr=cst( 1. ) ).evaluate();
    LOG(INFO) << "surface = " << integrate( _range=boundaryfaces( mesh ), _expr=cst( 1. ) ).evaluate();
}
