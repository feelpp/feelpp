#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporterhdf5.hpp>
#include <feel/feelcore/hdf5.hpp>
#include <feel/feeldiscr/mesh.hpp>

using namespace Feel;
int main (int argc, char ** argv) { 

    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="exporterhdf5",
                                  _author="Benjamin Vanthong",
                                  _email="benjamin.vanthong@icloud.com") );

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
/*
   Exporterhdf5 <Mesh<Simplex<2>>, 1> expo ("thomas", Environment::worldComm ()) ;
   expo.write (mesh) ;
*/

    return 0 ;
}


