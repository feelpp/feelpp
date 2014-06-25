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

        WorldComm M_comm ;
        std::string M_fileName ;
        HDF5 M_HDF5 ;
        mesh_ptrtype M_meshOut ;
        mesh_ptrtype M_meshIn ;

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
    std::cout << "write mesh" << std::endl ;
    std::cout << mesh_type::nRealDim << std::endl ;
    M_meshOut = mesh ;
    M_HDF5.openFile (M_fileName, M_comm, false) ;
    writePoints () ;
    M_HDF5.closeFile () ;

    M_meshOut.reset () ;
}

template <typename MeshType>
void Exporterhdf5<MeshType>::read (mesh_ptrtype& mesh)
{
    
}

template <typename MeshType>
void Exporterhdf5<MeshType>::writePoints () 
{
    auto pt_it = M_meshOut->beginPointWithProcessId () ;
    auto const pt_en = M_meshOut->endPointWithProcessId () ;
    size_type maxNumPoints= std::distance (pt_it, pt_en) ;
    std::cout << "nombre de Points : " << maxNumPoints << std::endl ;

    hsize_t currentSpaceDims [2] ;
    hsize_t currentCount [2] ;

    currentSpaceDims[0] = 1 ;
    currentSpaceDims[1] = maxNumPoints ;

    currentCount[0] = 3 ;
    currentCount[1] = maxNumPoints ;

    M_HDF5.createTable ("point_ids", H5T_STD_U32BE, currentSpaceDims) ;
    M_HDF5.createTable ("point_coords", H5T_IEEE_F64BE, currentCount) ;

    M_uintBuffer.resize (maxNumPoints, 0) ;
    M_realBuffer.resize (currentCount[0]*currentCount[1], 0) ;

    for (size_type i = 0 ; i < maxNumPoints ; i ++ , pt_it++) 
    {
        M_uintBuffer[i] = pt_it->id () + 1 ;
        std::cout << M_uintBuffer[i] << std::endl ;

        M_realBuffer[i] = pt_it->node()[0] ;
        if (mesh_type::nRealDim >= 2)
            M_realBuffer[maxNumPoints + i] = pt_it->node()[1] ;
        if (mesh_type::nRealDim >= 3)
            M_realBuffer[2*maxNumPoints + i] = pt_it->node()[2] ;
    }

    hsize_t currentOffset[2] = {0, 0} ;

    M_HDF5.write ("point_ids", H5T_NATIVE_INT, currentSpaceDims, currentOffset , &M_uintBuffer[0]) ;
    M_HDF5.write ("point_coords", H5T_NATIVE_DOUBLE, currentCount, currentOffset, &M_realBuffer[0]) ;

    
    M_HDF5.closeTable("point_ids") ;
    M_HDF5.closeTable("point_coords") ;
}

}
#endif /* FEELL_HAS_HDF5 */
#endif /* __Exporterhdf5_H */

