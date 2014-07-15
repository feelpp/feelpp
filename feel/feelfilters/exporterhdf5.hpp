#ifndef __Exporterhdf5_H
#define __Exporterhdf5_H 1

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/timeset.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>

#if defined(FEELPP_HAS_HDF5)
#include <feel/feelcore/hdf5.hpp>

namespace Feel 
{

template <typename MeshType, int N>
class Exporterhdf5
    : 
public Exporter <MeshType, N>
{
    typedef Exporter<MeshType, N> super ;
    public: 
        typedef MeshType mesh_type ;
        typedef typename mesh_type::value_type value_type ;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype ;    
        typedef typename super::timeset_type timeset_type;
        typedef typename super::timeset_ptrtype timeset_ptrtype;
        typedef typename super::timeset_iterator timeset_iterator;
        typedef typename super::timeset_const_iterator timeset_const_iterator;

        Exporterhdf5 () {}

        Exporterhdf5 (const std::string& fileName,
                      const WorldComm& comm) ;
        Exporterhdf5 ( WorldComm const& worldComm = Environment::worldComm() ) ;
        Exporterhdf5 ( std::string const& __p = "default", int freq = 1, WorldComm const& worldComm = Environment::worldComm() ) ;
        Exporterhdf5( po::variables_map const& vm=Environment::vm(), std::string const& exp_prefix = "", WorldComm const& worldComm = Environment::worldComm() );

        Exporterhdf5 ( Exporterhdf5 const & __ex ) ;

        virtual ~Exporterhdf5 () {}   

        void write (const mesh_ptrtype& mesh)  const ;

        void read (mesh_ptrtype& mesh)  const ;

    private :

        void init() ;
        void save () const ;
        void visit ( mesh_type* mesh) ;
        
        void writePoints () const ;
        void writeFaces () const ;
        void writeElements () const ;
        void writeElements1 () const ;
        void writeStats () const ;
        void writeDataNodes () const ;
        template<typename Iterator>
        void saveNodal ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const ;

        void write_xdmf_xml () const ;

        void bubbleSort (size_type * ids, value_type * coords, size_type n) const ;

        mutable WorldComm M_comm ;
        mutable std::string M_fileName ;
        mutable HDF5 M_HDF5 ;
        mutable mesh_ptrtype M_meshOut ;
        mutable mesh_ptrtype M_meshIn ;

        // Mesh geometry
        mutable size_type M_elementNodes ;
        mutable size_type M_maxNumElements ;
        mutable size_type M_maxNumPoints ;
        mutable size_type M_numParts ;
        mutable std::string M_element_type ;

        mutable std::vector<size_type> M_uintBuffer ;
        mutable std::vector<value_type> M_realBuffer ;
};

template<typename MeshType, int N>
inline Exporterhdf5<MeshType, N>::Exporterhdf5 (const std::string& fileName,
                                             const WorldComm& worldComm) :
    M_comm (worldComm),
    M_fileName (fileName),
    super ( worldComm )
{    
    init() ;
}

template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::Exporterhdf5( WorldComm const& worldComm )
:
super( worldComm ),
M_element_type()

{
    init();
}
template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::Exporterhdf5( std::string const& __p, int freq, WorldComm const& worldComm )
    :
    super( "hdf5", __p, freq, worldComm ),
    M_element_type()
{
    init();
}
template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::Exporterhdf5( po::variables_map const& vm, std::string const& exp_prefix, WorldComm const& worldComm )
    :
    super( vm, exp_prefix, worldComm )
{
    init();
}

template<typename MeshType, int N>
Exporterhdf5<MeshType,N>::Exporterhdf5( Exporterhdf5 const & __ex )
    :
    super( __ex ),
    M_element_type( __ex.M_element_type )
{
}

template<typename MeshType, int N>
void
Exporterhdf5<MeshType,N>::init()
{
    if ( mesh_type::nDim == 1 )
        if ( mesh_type::Shape == SHAPE_LINE )
            M_element_type = ( mesh_type::nOrder == 1 )?"Polyline":"Edge_3";

    if ( mesh_type::nDim == 2 )
    {
        if ( mesh_type::Shape == SHAPE_TRIANGLE )
            M_element_type = ( mesh_type::nOrder == 1 )?"Triangle":"Tri_6";

        else if ( mesh_type::Shape == SHAPE_QUAD )
            M_element_type = ( mesh_type::nOrder == 1 )?"Quadrilateral":"Quad_8";
    }

    if ( mesh_type::nDim == 3 )
    {
        if ( mesh_type::Shape == SHAPE_TETRA )
            M_element_type = ( mesh_type::nOrder == 1 )?"Tetrahedron":"Tet_10";

        else if ( mesh_type::Shape == SHAPE_HEXA )
            M_element_type = ( mesh_type::nOrder == 1 )?"Hexahedron":"Hex_20";
    }
}


template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::write (const mesh_ptrtype& mesh) const 
{
    std::cout << "mesh_type::nRealDim : " << mesh_type::nRealDim << std::endl ;
    M_meshOut = mesh ;
    M_fileName = "yo" ;
    M_HDF5.openFile (M_fileName+".h5", M_comm, false) ;
    writePoints () ;
    writeFaces () ;

    writeElements1 () ;
    writeStats () ;

    write_xdmf_xml () ;
    M_HDF5.closeFile () ;

    M_meshOut.reset () ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::read (mesh_ptrtype& mesh) const 
{
   mesh.reset () ;
   M_meshIn.reset (new mesh_type) ;
   
   mesh = M_meshIn ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writePoints () const 
{
    auto pt_it = M_meshOut->beginPointWithProcessId () ;
    auto const pt_en = M_meshOut->endPointWithProcessId () ;
    M_maxNumPoints= std::distance (pt_it, pt_en) ;
    std::cout << "nombre de Points : " << M_maxNumPoints << std::endl ;

    hsize_t currentSpaceDims [2] ;
    hsize_t currentCount [2] ;

    currentSpaceDims[0] = 1;
    currentSpaceDims[1] = M_maxNumPoints ;

    currentCount[0] = M_maxNumPoints ;
    currentCount[1] = 3 ;

    M_HDF5.createTable ("point_coords", H5T_IEEE_F64BE, currentCount) ;
    M_HDF5.createTable ("point_ids", H5T_STD_U32BE, currentSpaceDims) ;

    M_uintBuffer.resize (currentSpaceDims[0]*currentSpaceDims[1], 0) ;
    M_realBuffer.resize (currentCount[0]*currentCount[1], 0) ;

    for (size_type i = 0 ; i < M_maxNumPoints ; i++ , pt_it++) 
    {
        M_uintBuffer[i] = pt_it->id () ;

        M_realBuffer[3*i] = pt_it->node()[0] ;
        if (mesh_type::nRealDim >= 2)
            M_realBuffer[3*i + 1] = pt_it->node()[1] ;
        if (mesh_type::nRealDim >= 3)
            M_realBuffer[3*i + 2] = pt_it->node()[2] ;
    }

    bubbleSort (&M_uintBuffer[0], &M_realBuffer[0], M_maxNumPoints) ;

    hsize_t currentOffset[2] = {0, 0} ;

    M_HDF5.write ("point_coords", H5T_NATIVE_DOUBLE, currentCount, currentOffset, &M_realBuffer[0]) ;
    M_HDF5.write ("point_ids", H5T_NATIVE_LLONG, currentSpaceDims, currentOffset , &M_uintBuffer[0]) ;

    
    M_HDF5.closeTable("point_coords") ;
    M_HDF5.closeTable("point_ids") ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeFaces () const 
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

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeElements1 () const 
{
    typename mesh_type::parts_const_iterator_type p_it = M_meshOut->beginParts();
    typename mesh_type::parts_const_iterator_type p_en = M_meshOut->endParts();
    M_numParts = std::distance (p_it, p_en) ;
    std::cout << "M_numParts : " << M_numParts << std::endl ;
    M_maxNumElements = 0 ;
    for (int i = 0 ; i < M_numParts ; i++, p_it ++) 
    {
        auto elt_it = M_meshOut->beginElementWithMarker (p_it->first) ;
        auto elt_en = M_meshOut->endElementWithMarker (p_it->first) ;
        M_maxNumElements += std::distance (elt_it, elt_en) ;
        std::cout << "std::distance (elt_it, elt_en) : " << std::distance(elt_it, elt_en) << std::endl ;
    }
    std::cout << "M_numMaxElements : " << M_maxNumElements << std::endl ;
    M_elementNodes = M_meshOut-> numLocalVertices () ;
    std::cout << "nombre de Points par element : " << M_elementNodes << std::endl ;

    hsize_t currentSpacesDims [2] ;
    hsize_t currentSpacesDims2 [2] ;

    currentSpacesDims [0] = M_maxNumElements ;
    currentSpacesDims [1] = M_elementNodes ;

    currentSpacesDims2 [0] = 1 ;
    currentSpacesDims2 [1] = M_maxNumElements ;

    M_HDF5.createTable ("element_ids", H5T_STD_U32BE, currentSpacesDims2) ;
    M_HDF5.createTable ("element_nodes", H5T_STD_U32BE, currentSpacesDims) ;

    M_uintBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0) ;
    std::vector<size_type> idsBuffer ;
    idsBuffer.resize (currentSpacesDims2[1], 0) ;

    M_uintBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0) ;
    
    size_type numElementPerPart = 0 ;
    size_type k = 0 ;
    size_type i = 0 ;
    for (p_it = M_meshOut->beginParts (); k < M_numParts ;  p_it++ , k++)
    {
        auto elt_it = M_meshOut->beginElementWithMarker (p_it->first) ;
        auto elt_en = M_meshOut->endElementWithMarker (p_it->first) ;
        numElementPerPart = std::distance (elt_it, elt_en) ;
        std::cout << "numElementPerPart : " << numElementPerPart << std::endl ;
        for ( ; elt_it != elt_en ; ++elt_it , i ++)
        {
            idsBuffer[i] = elt_it->id () ;
            for ( size_type j = 0 ; j < M_elementNodes ; j ++ )
               M_uintBuffer[j + M_elementNodes*i] = elt_it->point(j).id() -1 ; 
        }
    }

    hsize_t currentOffset[2] = {0, 0} ;
    M_HDF5.write ( "element_ids", H5T_NATIVE_LLONG, currentSpacesDims2, currentOffset, &idsBuffer[0] ) ;
    M_HDF5.write ( "element_nodes", H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &M_uintBuffer[0] ) ;

    M_HDF5.closeTable ("element_ids") ;
    M_HDF5.closeTable ("element_nodes") ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeElements () const  
{   
    typename mesh_type::parts_const_iterator_type p_it = M_meshOut->beginParts();
    typename mesh_type::parts_const_iterator_type p_en = M_meshOut->endParts();
    M_numParts = std::distance (p_it, p_en) ;
    std::cout << "M_numParts : " << M_numParts << std::endl ;

    auto elt_it = M_meshOut->beginElement () ;
    auto elt_en = M_meshOut->endElement () ;
    M_maxNumElements = std::distance (elt_it, elt_en) ;
    std::cout<<"nombre d'Elements : " << std::distance (elt_it, elt_en) << std::endl ;
    
    M_elementNodes = M_meshOut-> numLocalVertices () ;
    std::cout << "nombre de Points par element : " << M_elementNodes << std::endl ;

    hsize_t currentSpacesDims [2] ;
    hsize_t currentSpacesDims2 [2] ;

    currentSpacesDims [0] = M_maxNumElements ;
    currentSpacesDims [1] = M_elementNodes ;

    currentSpacesDims2 [0] = 1 ;
    currentSpacesDims2 [1] = M_maxNumElements ;

    std::string str_ids = std::string("element_ids") ;
    std::string str_nodes = std::string ("element_nodes") ;
    std::cout << "element_ids " << str_ids << "  element_nodes " << str_nodes <<std::endl ; 

    //M_HDF5.createTable ("element_ids", H5T_STD_U32BE, currentSpacesDims2) ;
    //M_HDF5.createTable ("element_nodes", H5T_STD_U32BE, currentSpacesDims) ;

    M_HDF5.createTable (str_ids.c_str(), H5T_STD_U32BE, currentSpacesDims2) ;
    M_HDF5.createTable (str_nodes.c_str(), H5T_STD_U32BE, currentSpacesDims) ;

    M_uintBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0) ;
    std::vector<size_type> idsBuffer ;
    idsBuffer.resize (currentSpacesDims2[1], 0) ;

    M_uintBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0) ;
    for (size_type i= 0  ; elt_it != elt_en ; ++elt_it , i ++) 
    {
            idsBuffer[i] = elt_it->id() ;
        for ( size_type j = 0 ; j < M_elementNodes ; j++ ) 
        {
                M_uintBuffer[j + M_elementNodes*i] = elt_it->point(j).id() - 1 ;
        }
    }

    hsize_t currentOffset[2] = {0, 0} ;
    M_HDF5.write ( str_ids.c_str(), H5T_NATIVE_LLONG, currentSpacesDims2, currentOffset, &idsBuffer[0] ) ;
    M_HDF5.write ( str_nodes.c_str(), H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &M_uintBuffer[0] ) ;

    M_HDF5.closeTable (str_ids.c_str()) ;
    M_HDF5.closeTable (str_nodes.c_str()) ;
}

template <typename MeshType, int N>
template <typename Iterator>
void Exporterhdf5<MeshType, N>::saveNodal ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const 
{   
    typename mesh_type::parts_const_iterator_type p_it = M_meshOut->beginParts();
    typename mesh_type::parts_const_iterator_type p_en = M_meshOut->endParts();
    M_numParts = std::distance (p_it, p_en) ;


    for ( size_type i = 0 ; i < M_numParts ; i ++ , p_it++ ) 
    {
        uint16_type nComponents = __var -> second.nComponents ;
        if ( __var->second.is_vectorial )
            nComponents = 3 ;

        auto r = markedelements (M_meshOut, *p_it, EntityProcessType::ALL) ;
        auto elt_it = r.template get<1>() ;
        auto elt_en = r.template get<2>() ;

    }





   /* 
    auto d_it = M_meshOut -> beginNodalScalar () ;
    auto d_en = M_meshOut -> beginNodalScalar () ;

    auto __u = d_it -> second ;
    hsize_t currentSpacesDims [2] ;

    currentSpacesDims [0] = 1 ;
    currentSpacesDims [1] = M_maxNumPoints ;

    M_HDF5.createTable ("dataNodes", H5T_IEEE_F64BE, currentSpacesDims) ;
    M_realBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0) ;
 
    auto pt_it = M_meshOut->beginPoint () ;
    auto pt_en = M_meshOut->endPoint () ;

    for ( size_type j = 0 ; j < M_maxNumPoints ; j++, ++pt_it ) 
    {
        M_realBuffer[j] = __u( pt_it->id() ) ;
    }
    hsize_t currentOffset[2] = {0, 0} ;
    M_HDF5.write ( "dataNodes", H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] ) ;

    M_HDF5.closeTable ("dataNodes") ;
    */
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeDataNodes () const 
{
    hsize_t currentSpacesDims [2] ;

    currentSpacesDims [0] = 1 ;
    currentSpacesDims [1] = M_maxNumPoints ;

    M_HDF5.createTable ("dataNodes", H5T_IEEE_F64BE, currentSpacesDims) ;

    M_realBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0) ;
    for ( size_type j = 0 ; j < M_maxNumPoints ; j++ ) 
    {
        //M_realBuffer[j] =  j/M_maxNumPoints ;
        //M_realBuffer[j] =  1.*j/M_maxNumPoints ;
        //M_realBuffer[j] = j < 610 ? 1 : 0.5*j/M_maxNumPoints ;
        M_realBuffer[j] = j < 610 ? 1 : 0.5 ;
    }
    M_realBuffer[M_maxNumPoints-1] = 0 ;
    hsize_t currentOffset[2] = {0, 0} ;
    M_HDF5.write ( "dataNodes", H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] ) ;

    M_HDF5.closeTable ("dataNodes") ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeStats () const 
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

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::bubbleSort (size_type * ids, value_type * coords, size_type n) const 
{
    size_type int_tmp ;
    value_type float_tmp ;
    bool swapped = false ;
    do 
    {
        swapped = false ;
        for (size_type j = 0 ; j < n-1 ; j ++) 
        {
            if (ids[j] > ids[j+1]) 
            {
                int_tmp = ids[j] ;
                ids[j] = ids[j+1] ;
                ids[j+1] = int_tmp ;

                float_tmp = coords[3*j] ;
                coords[3*j] = coords [3*(j+1)] ;
                coords[3*(j+1)] = float_tmp ;

                float_tmp = coords[3*j+1] ;
                coords[3*j+1] = coords [3*(j+1)+1] ;
                coords[3*(j+1)+1] = float_tmp ;

                float_tmp = coords[3*j+2] ;
                coords[3*j+2] = coords [3*(j+1)+2] ;
                coords[3*(j+1)+2] = float_tmp ;

                swapped = true ;
            }
        }
        n = n -1 ;
    }
    while (swapped) ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::save () const 
{
    std::cout << "hdf5 exporter" << std::endl ;

    timeset_const_iterator __ts_it = this->beginTimeSet () ;
    timeset_const_iterator __ts_en = this->endTimeSet () ;
    
    timeset_ptrtype __ts = *__ts_it ;
    typename timeset_type::step_const_iterator __it = __ts->beginStep () ;
    typename timeset_type::step_const_iterator __end = __ts->endStep () ;
    __it = boost::prior ( __end ) ;

    typename timeset_type::step_ptrtype __step = * __it ;
    

    write (__step->mesh()) ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::visit ( mesh_type* mesh) 
{
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::write_xdmf_xml () const
{
    std::cout << "M_element_type : " << M_element_type << std::endl ;
    FILE * xmf = 0 ;
    xmf = fopen ((M_fileName+".xmf").c_str(), "w") ;
    fprintf (xmf, "<?xml version=\"1.0\" ?>\n") ;
    fprintf (xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n") ;
    fprintf (xmf, "<Xdmf Version=\"2.0\">\n") ;
    fprintf (xmf, " <Domain>\n") ;
    fprintf (xmf, "     <Grid Name=\"%s\" GridType=\"Uniform\">\n", M_fileName.c_str()) ;
    fprintf (xmf, "         <Topology TopologyType=\"%s\" NumberOfElements=\"%zu\" NodesPerElement=\"%zu\">\n", M_element_type.c_str(), M_maxNumElements, M_elementNodes) ;
    fprintf (xmf, "             <DataItem Dimensions=\"%zu %zu\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n", M_maxNumElements, M_elementNodes) ;
    fprintf (xmf, "             %s.h5:/element_nodes\n", M_fileName.c_str()) ;
    fprintf (xmf, "             </DataItem>\n") ;
    fprintf (xmf, "         </Topology>\n") ;
    fprintf (xmf, "         <Geometry GeometryType=\"XYZ\">\n") ;
    fprintf (xmf, "             <DataItem Dimensions=\"%zu 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n", M_maxNumPoints) ;
    fprintf (xmf, "             %s.h5:/point_coords\n", M_fileName.c_str()) ;
    fprintf (xmf, "             </DataItem>\n") ;
    fprintf (xmf, "         </Geometry>\n") ;
    fprintf (xmf, "     </Grid>\n") ;
    fprintf (xmf, " </Domain>\n") ;
    fprintf (xmf, "</Xdmf>\n") ;
    fclose(xmf) ;
}

}
#endif /* FEELL_HAS_HDF5 */
#endif /* __Exporterhdf5_H */

