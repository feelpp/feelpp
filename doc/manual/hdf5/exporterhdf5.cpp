#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/partitionio.hpp>
#include <feel/feelcore/hdf5.hpp>

using namespace Feel;
int main (int argc, char ** argv) { 

    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="exporterhdf5",
                                  _author="Benjamin Vanthong",
                                  _email="benjamin.vanthong@icloud.com") );

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);

    h5partition (_mesh = mesh) ;

    return 0 ;
}


