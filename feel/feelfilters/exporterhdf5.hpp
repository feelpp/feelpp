#ifndef __Exporterhdf5_H
#define __Exporterhdf5_H 1

#include <iostream>

#if defined(FEELPP_HAS_HDF5)
#include <feel/feelcore/hdf5.hpp>

namespace Feel 
{

template <typename MeshType>
class Exporterhdf5
{
    public: 
        typedef MeshType mesh_type ;
        typedef typename mesh_type::value_type value_type ;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype ;

        Exporterhdf5 () {}

        Exporterhdf5 (const std::string& fileName,
                      const WorldComm& comm) ;
        
        virtual ~Exporterhdf5 () {}   

        void write (const mesh_ptrtype& mesh) ;

        void read (mesh_ptrtype& mesh) ;

    private :

        void writePoints () ;
        void writeFaces () ;
        void writeElements () ;
        void writeStats () ;

        WorldComm M_comm ;
        std::string M_fileName ;
        HDF5 M_HDF5 ;
        mesh_ptrtype M_meshOut ;
        mesh_ptrtype M_meshIn ;

        // Mesh geometry
        size_type M_elementNodes ;
        size_type M_maxNumElements ;
        size_type M_maxNumPoints ;

        std::vector<size_type> M_uintBuffer ;
        std::vector<value_type> M_realBuffer ;
};

template<typename MeshType>
inline Exporterhdf5<MeshType>::Exporterhdf5 (const std::string& fileName,
                                             const WorldComm& comm) :
    M_comm (comm),
    M_fileName (fileName)
{

}

template <typename MeshType>
void Exporterhdf5<MeshType>::write (const mesh_ptrtype& mesh)
{
    std::cout << mesh_type::nRealDim << std::endl ;
    M_meshOut = mesh ;
    M_HDF5.openFile (M_fileName, M_comm, false) ;
    writePoints () ;
    writeFaces () ;
    writeElements () ;
    writeStats () ;
    M_HDF5.closeFile () ;

    M_meshOut.reset () ;
}

template <typename MeshType>
void Exporterhdf5<MeshType>::read (mesh_ptrtype& mesh)
{
   mesh.reset () ;
   M_meshIn.reset (new mesh_type) ;
   
   mesh = M_meshIn ;
}

template <typename MeshType>
void Exporterhdf5<MeshType>::writePoints () 
{
    auto pt_it = M_meshOut->beginPointWithProcessId () ;
    auto const pt_en = M_meshOut->endPointWithProcessId () ;
    M_maxNumPoints= std::distance (pt_it, pt_en) ;
    std::cout << "nombre de Points : " << M_maxNumPoints << std::endl ;

    hsize_t currentSpaceDims [2] ;
    hsize_t currentCount [2] ;

    currentSpaceDims[0] = 1;
    currentSpaceDims[1] = M_maxNumPoints ;

    currentCount[0] = 3 ;
    currentCount[1] = M_maxNumPoints ;

    M_HDF5.createTable ("point_coords", H5T_IEEE_F64BE, currentCount) ;
    M_HDF5.createTable ("point_ids", H5T_STD_U32BE, currentSpaceDims) ;

    M_uintBuffer.resize (currentSpaceDims[0]*currentSpaceDims[1], 0) ;
    M_realBuffer.resize (currentCount[0]*currentCount[1], 0) ;

    for (size_type i = 0 ; i < M_maxNumPoints ; i++ , pt_it++) 
    {
        M_uintBuffer[i] = pt_it->id () ;

        M_realBuffer[i] = pt_it->node()[0] ;
        if (mesh_type::nRealDim >= 2)
            M_realBuffer[M_maxNumPoints + i] = pt_it->node()[1] ;
        if (mesh_type::nRealDim >= 3)
            M_realBuffer[2*M_maxNumPoints + i] = pt_it->node()[2] ;
    }

    hsize_t currentOffset[2] = {0, 0} ;

    M_HDF5.write ("point_coords", H5T_NATIVE_DOUBLE, currentCount, currentOffset, &M_realBuffer[0]) ;
    M_HDF5.write ("point_ids", H5T_NATIVE_LLONG, currentSpaceDims, currentOffset , &M_uintBuffer[0]) ;

    
    M_HDF5.closeTable("point_coords") ;
    M_HDF5.closeTable("point_ids") ;
}

template <typename MeshType>
void Exporterhdf5<MeshType>::writeFaces () 
{
/*
    //auto face_it = M_meshOut->beginFaceWithProcessId () ;
    auto face_it = M_meshOut->beginFace () ;
    //auto const face_en = M_meshOut->endFaceWithProcessId () ;
    auto const face_en = M_meshOut->endFace () ;
    size_type maxNumFaces= std::distance (face_it, face_en) ;
    std::cout << "nombre de Face : " << maxNumFaces << std::endl ;

    std::cout << face_it->id() << std::endl ;
    face_it++ ;
    std::cout << face_it->id() << std::endl ;
*/  
}

template <typename MeshType>
void Exporterhdf5<MeshType>::writeElements () 
{
    auto elt_it = M_meshOut->beginElementWithMarker (M_meshOut->beginParts()->first) ;
    auto elt_en = M_meshOut->endElementWithMarker (M_meshOut->beginParts()->first) ;
    M_maxNumElements = std::distance (elt_it, elt_en) ;
    std::cout << "nombre d'Elements : " << M_maxNumElements << std::endl ;
    
    M_elementNodes = M_meshOut-> numLocalVertices () ;
    std::cout << " nombre de Points par element : " << M_elementNodes << std::endl ;

    hsize_t currentSpacesDims [2] ;

    currentSpacesDims [0] = M_elementNodes + 1 ;
    currentSpacesDims [1] = M_maxNumElements ;

    M_HDF5.createTable ("elements", H5T_STD_U32BE, currentSpacesDims) ;

    M_uintBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], -1) ;
    for ( size_type i = 0 ; elt_it != elt_en ; ++elt_it, i++ ) 
    {
        M_uintBuffer[i] = elt_it->id() ;
        for ( size_type j = 1 ; j < M_elementNodes + 1; j++ ) 
        {
            M_uintBuffer[j*M_maxNumElements + i] = elt_it->point(j-1).id() ;
        }
    }

    hsize_t currentOffset[2] = {0, 0} ;
    M_HDF5.write ( "elements", H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &M_uintBuffer[0] ) ;

    M_HDF5.closeTable ("elements") ;
}

template <typename MeshType>
void Exporterhdf5<MeshType>::writeStats () 
{
    size_type numStats = 3 ;
    hsize_t currentSpacesDims[2] = {1, 3} ;
    
    M_HDF5.createTable ( "stats", H5T_STD_U32BE, currentSpacesDims ) ;
    M_uintBuffer.resize (numStats) ;

    M_uintBuffer[0] = M_maxNumPoints ;
    M_uintBuffer[1] = M_maxNumElements ;
    M_uintBuffer[2] = M_elementNodes ;
    
    hsize_t currentOffset [2] = {0, 0} ;
    M_HDF5.write ("stats", H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &M_uintBuffer[0]) ;
    M_HDF5.closeTable ("stats") ;
}

}
#endif /* FEELL_HAS_HDF5 */
#endif /* __Exporterhdf5_H */

