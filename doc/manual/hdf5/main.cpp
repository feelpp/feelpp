#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelcore/hdf5.hpp>
#include <iostream>

using namespace Feel;
int main (int argc, char ** argv) { 
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="main",
                                  _author="Benjamin Vanthong",
                                  _email="benjamin.vanthong@icloud.com") );
/*
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    std::cout << mesh->numGlobalElements () << std::endl ;
    std::cout << mesh->numFaces () << std::endl ;
    std::cout << mesh->point(0).id() << std::endl ;
    //std::cout << mesh->faceList[0].point(0).localId() << std::endl ;
    std::cout << mesh->edgeList[0].id() << std::endl ;
*/ 
    HDF5 hdf5 ("boom.h5", Feel::detail::Environment::worldComm(), false) ; 
    hdf5.createGroup ("./N/") ;

   hsize_t currentSpacesDims[2] = {4,5} ;
   hdf5.createTableInGroup ("./N/", "elements", H5T_IEEE_F64BE, currentSpacesDims) ;

   double tab[20] ;
   for (int i = 0 ; i < 4 ; i++)
       for (int j = 0 ; j < 5 ; j ++)
           tab[i*5 + j] = i + j  ; 

   hsize_t currentOffset[2] = {0,0} ;
   hdf5.write ("./N/elements", H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, tab) ; 

   hdf5.closeTable ("./N/elements") ;
    
    hdf5.closeFile () ;
    return 0 ;
}


