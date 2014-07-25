#ifndef __Exporterhdf5_H
#define __Exporterhdf5_H 1

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/timeset.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <map>

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

        void save () const ;
        void visit ( mesh_type* mesh) ;

    private :

        void init() ;
        void write ()  const ;
        
        void writePoints () const ;
        void writeElements () const ;
        void writeStats () const ;
        template<typename Iterator>
        void saveNodal ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const ;
        template<typename Iterator>
        void saveElement ( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const ;

        void write_xdmf () const ;
        void open_xdmf_xml () const ;
        void close_xdmf_xml () const ;

        void bubbleSort (size_type * ids, value_type * coords, size_type n) const ;

        mutable WorldComm M_comm ;
        mutable std::string M_fileName ;
        mutable std::string M_fileNameStep ;
        mutable HDF5 M_HDF5 ;
        mutable mesh_ptrtype M_meshOut ;

        // Mesh geometry
        mutable size_type M_elementNodes ;
        mutable size_type M_maxNumElements ;
        mutable size_type M_maxNumPoints ;
        mutable size_type M_numParts ;
        mutable std::string M_element_type ;

        mutable size_type M_step = 0 ;
        mutable std::ofstream M_xmf  ;

        mutable std::vector<size_type> M_uintBuffer ;
        mutable std::vector<value_type> M_realBuffer ;
        mutable std::map<size_type, size_type> M_newPointId ;
        mutable std::map<size_type, size_type> M_newElementId ;
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
void Exporterhdf5<MeshType, N>::write () const 
{

    std::ostringstream str ;
    str <<  this->prefix () << "-" << Environment::worldComm().globalSize()<<"_"<<Environment::worldComm().globalRank()  ;
    M_fileName = str.str () ; 

    open_xdmf_xml () ;

    timeset_const_iterator __ts_it = this->beginTimeSet () ;
    timeset_const_iterator __ts_en = this->endTimeSet () ;
    
    timeset_ptrtype __ts = *__ts_it ;

    M_xmf << "       <Grid Name=\"Simulation over time\" GridType=\"Collection\" CollectionType=\"Temporal\">" << "\n" ;
    M_xmf << "           <Time TimeType=\"HyperSlab\">\n" ;
    M_xmf << "               <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\">\n" ;
    M_xmf << "               " << (*(__ts->beginStep()))->time() << " " << __ts->timeIncrement() << " " << __ts->numberOfTotalSteps() << "\n" ;
    M_xmf << "               </DataItem>\n" ;
    M_xmf << "           </Time>\n" ;

    while ( __ts_it != __ts_en )
    {
        __ts = *__ts_it ;
        typename timeset_type::step_const_iterator __it = __ts->beginStep () ;
        typename timeset_type::step_const_iterator __end = __ts->endStep () ;
        __it = boost::prior ( __end ) ;

        while ( __it != __end )
        {
            typename timeset_type::step_ptrtype __step = * __it ;

            if ( __step->isInMemory() )
            {
                M_meshOut = __step->mesh () ;
                std::ostringstream filestr ;
                filestr << "-" << M_step++ ;
                M_fileNameStep = M_fileName+filestr.str() ;
                M_HDF5.openFile (M_fileNameStep+".h5", Environment::worldCommSeq(), false) ;
                M_xmf << "           <Grid Name=\"" << M_fileNameStep << "\" GridType=\"Uniform\">\n" ;

                writePoints () ;
                writeElements () ;

                std::cout << "time                          : " << __step->time () << std::endl ; 
                std::cout << "time increment                : " << __ts->timeIncrement () << std::endl ; 
                std::cout << "numberOfSteps                 : " << __ts->numberOfSteps () << std::endl ; 
                std::cout << "numberOfTotalSteps            : " << __ts->numberOfTotalSteps () << std::endl ;

                saveNodal (__step, __step->beginNodalScalar(), __step->endNodalScalar() ) ;
                saveNodal (__step, __step->beginNodalVector(), __step->endNodalVector() ) ;
                
                saveElement (__step, __step->beginElementScalar(), __step->endElementScalar() ) ;
                saveElement (__step, __step->beginElementVector(), __step->endElementVector() ) ;

                M_xmf << "           </Grid>\n" ;
                writeStats() ; 
                M_HDF5.closeFile () ;
            }
            ++__it ;
        }
        ++__ts_it ;
    }
    close_xdmf_xml () ;
    M_meshOut.reset () ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writePoints () const 
{
    auto pt_it = M_meshOut->beginPointWithProcessId () ;
    auto const pt_en = M_meshOut->endPointWithProcessId () ;
    M_maxNumPoints= std::distance (pt_it, pt_en) ;

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

    for (size_type i = 0 ; i < M_maxNumPoints ; i ++) 
        M_newPointId[M_uintBuffer[i]] = i ;

    hsize_t currentOffset[2] = {0, 0} ;

    M_HDF5.write ("point_coords", H5T_NATIVE_DOUBLE, currentCount, currentOffset, &M_realBuffer[0]) ;
    M_HDF5.write ("point_ids", H5T_NATIVE_LLONG, currentSpaceDims, currentOffset , &M_uintBuffer[0]) ;

    M_HDF5.closeTable("point_coords") ;
    M_HDF5.closeTable("point_ids") ;    

    M_xmf << "           <Geometry GeometryType=\"XYZ\">\n" ;
    M_xmf << "               <DataItem Dimensions=\"" << M_maxNumPoints << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n" ;
    M_xmf << "               " << M_fileNameStep << ".h5:/point_coords\n" ;
    M_xmf << "               </DataItem>\n" ;
    M_xmf << "           </Geometry>\n" ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeElements () const 
{
    typename mesh_type::parts_const_iterator_type p_it = M_meshOut->beginParts();
    typename mesh_type::parts_const_iterator_type p_en = M_meshOut->endParts();
    M_numParts = std::distance (p_it, p_en) ;
    M_maxNumElements = 0 ;
    for (int i = 0 ; i < M_numParts ; i++, p_it ++) 
    {
        auto elt_it = M_meshOut->beginElementWithMarker (p_it->first) ;
        auto elt_en = M_meshOut->endElementWithMarker (p_it->first) ;
        M_maxNumElements += std::distance (elt_it, elt_en) ;
    }
    M_elementNodes = M_meshOut-> numLocalVertices () ;

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
    
    size_type k = 0 ;
    size_type i = 0 ;
    M_newElementId.clear () ;
    for (p_it = M_meshOut->beginParts (); k < M_numParts ;  p_it++ , k++)
    {
        auto elt_it = M_meshOut->beginElementWithMarker (p_it->first) ;
        auto elt_en = M_meshOut->endElementWithMarker (p_it->first) ;
        for ( ; elt_it != elt_en ; ++elt_it , i ++)
        {
            idsBuffer[i] = elt_it->id () ;
            M_newElementId[idsBuffer[i]] = i ;
            for ( size_type j = 0 ; j < M_elementNodes ; j ++ )
               M_uintBuffer[j + M_elementNodes*i] = M_newPointId[elt_it->point(j).id()]  ; 
        }
    }


    hsize_t currentOffset[2] = {0, 0} ;
    M_HDF5.write ( "element_ids", H5T_NATIVE_LLONG, currentSpacesDims2, currentOffset, &idsBuffer[0] ) ;
    M_HDF5.write ( "element_nodes", H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &M_uintBuffer[0] ) ;

    M_HDF5.closeTable ("element_ids") ;
    M_HDF5.closeTable ("element_nodes") ;

    M_xmf << "           <Topology TopologyType=\"" << M_element_type << "\" NumberOfElements=\"" << M_maxNumElements << "\" NodesPerElement=\"" << M_elementNodes << "\">\n" ;
    M_xmf << "               <DataItem Dimensions=\"" <<M_maxNumElements << " " << M_elementNodes << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n" ;
    M_xmf << "               " << M_fileNameStep << ".h5:/element_nodes\n" ;
    M_xmf << "               </DataItem>\n" ;
    M_xmf << "           </Topology>\n" ;
}

template <typename MeshType, int N>
template <typename Iterator>
void Exporterhdf5<MeshType, N>::saveNodal ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const 
{   
    while ( __var != en )
    {    
        std::string attributeType ("Scalar") ;

        std::string solutionName = __var->first ;

        uint16_type nComponents = __var -> second.nComponents ;

        if ( __var->second.is_scalar )
        {
             solutionName += ".scl" ; 
             attributeType = "Scalar" ;
             nComponents = 1 ;
        }
        else if ( __var->second.is_vectorial )
        {
            solutionName += ".vec" ;
            attributeType = "Vector" ;
            nComponents = 3 ;
        }
        else if ( __var->second.is_tensor2 )
        {
            solutionName += ".tsr" ;
            attributeType = "Tensor" ;
            nComponents = 9 ;
        }

        solutionName += ".node" ;

        std::cout << "solution name                 : " << solutionName << std::endl ;

        hsize_t currentSpacesDims [2] ;

        currentSpacesDims [0] = nComponents ;
        currentSpacesDims [1] = M_maxNumPoints ;

        M_HDF5.createTable (solutionName.c_str(), H5T_IEEE_F64BE, currentSpacesDims) ;

        M_realBuffer.resize (M_maxNumPoints*nComponents, 0) ;

        typename mesh_type::parts_const_iterator_type p_it = M_meshOut->beginParts();
        typename mesh_type::parts_const_iterator_type p_en = M_meshOut->endParts();
        M_numParts = std::distance (p_it, p_en) ;

        for ( ; p_it != p_en ; p_it++) 
        {

            auto r = markedelements (M_meshOut, p_it->first, EntityProcessType::LOCAL_ONLY) ;
            auto elt_it = r.template get<1>() ;
            auto elt_en = r.template get<2>() ;

            Feel::detail::MeshPoints<float> mp ( __step->mesh().get(), elt_it, elt_en, true, true, true ) ;

            size_type e = 0 ; 
            
            for ( ; elt_it != elt_en ; ++elt_it )
            {
                for ( uint16_type c = 0 ; c < nComponents ; ++c )
                {
                    for ( uint16_type p = 0 ; p < __step->mesh()->numLocalVertices() ; ++p, ++e )
                    {
                        size_type  ptid = M_newPointId[elt_it->get().point(p).id()]  ;
                        size_type global_node_id = mp.ids.size()*c + ptid ;
                        if ( c < __var->second.nComponents ) 
                        {
                            size_type dof_id = boost::get<0>( __var->second.functionSpace()->dof()->localToGlobal ( elt_it->get().id(), p, c ) ) ;
                            M_realBuffer[global_node_id] = __var->second.globalValue ( dof_id ) ;
                        }
                        else
                        {
                            M_realBuffer[global_node_id] = 0.0 ;
                        }
                    }
                }
            }
        }    
        hsize_t currentOffset[2] = {0, 0} ;

        M_HDF5.write ( solutionName.c_str(), H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] ) ;

        M_HDF5.closeTable (solutionName.c_str()) ;    

        M_xmf << "           <Attribute AttributeType=\""<< attributeType << "\" Name=\"" << solutionName << "\" Center=\"Node\">\n" ;
        M_xmf << "               <DataItem Dimensions=\""<< nComponents << " " << M_maxNumPoints << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n" ;
        M_xmf << "               "<< M_fileNameStep <<".h5:/"<< solutionName <<"\n" ;
        M_xmf << "               </DataItem>\n" ;    
        M_xmf << "           </Attribute>\n" ;
        ++__var ;
    }
}

template<typename MeshType, int N>
template<typename Iterator>
void Exporterhdf5<MeshType, N>::saveElement ( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const 
{
    while ( __evar != __evaren ) 
    {
        std::string attributeType ("Scalar") ;

        std::string solutionName = __evar->first ;
        std::cout << "solution Name                 : " << solutionName << std::endl ;
        uint16_type nComponents = __evar->second.nComponents ;

        if ( __evar->second.is_scalar )
        {
             solutionName += ".scl" ; 
             attributeType = "Scalar" ;
             nComponents = 1 ;
        }
        else if ( __evar->second.is_vectorial )
        {
            solutionName += ".vec" ;
            attributeType = "Vector" ;
            nComponents = 3 ;
        }
        else if ( __evar->second.is_tensor2 )
        {
            solutionName += ".tsr" ;
            attributeType = "Tensor" ;
            nComponents = 9 ;
        }
        solutionName += ".element" ;
        hsize_t currentSpacesDims [2] ;

        currentSpacesDims [0] = nComponents ;
        currentSpacesDims [1] = M_maxNumElements ;

        M_HDF5.createTable (solutionName.c_str(), H5T_IEEE_F64BE, currentSpacesDims) ;

        M_realBuffer.resize (M_maxNumElements*nComponents, 0) ;

        typename mesh_type::parts_const_iterator_type p_it = M_meshOut->beginParts() ;
        typename mesh_type::parts_const_iterator_type p_en = M_meshOut->endParts() ;
        for ( ; p_it != p_en ; p_it++ ) 
        {
            typename mesh_type::marker_element_const_iterator elt_st ;
            typename mesh_type::marker_element_const_iterator elt_en ;
            boost::tie( elt_st, elt_en ) = __step->mesh()->elementsWithMarker( p_it->first, __evar->second.worldComm().localRank() ) ;

            if ( !__evar->second.areGlobalValuesUpdated() )
                __evar->second.updateGlobalValues() ;

            size_type ncells = std::distance ( elt_st, elt_en ) ;
            for ( int c = 0 ; c < nComponents ; ++c )
            {
                size_type e = 0 ;
                for ( auto elt_it = elt_st ; elt_it != elt_en ; ++elt_it, ++e )
                {
                    size_type global_node_id = c*ncells+e ;
                    if ( c < __evar->second.nComponents )
                    {
                        size_type dof_id = boost::get<0>( __evar->second.functionSpace()->dof()->localToGlobal( elt_it->id(), 0, c ) ) ;
                        M_realBuffer[global_node_id] = __evar->second.globalValue ( dof_id ) ;
                    }
                    else 
                        M_realBuffer[global_node_id] = 0 ;
                }
            }
        }   
        hsize_t currentOffset[2] = {0, 0} ;

        M_HDF5.write ( solutionName.c_str(), H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] ) ;

        M_HDF5.closeTable (solutionName.c_str()) ;    

        M_xmf << "           <Attribute AttributeType=\"" << attributeType << "\" Name=\"" << solutionName << "\" Center=\"Cell\">\n" ;
        M_xmf << "               <DataItem Dimensions=\""<< nComponents <<" "<< M_maxNumElements << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n" ;
        M_xmf << "               "<< M_fileNameStep << ".h5:/" << solutionName << "\n" ;
        M_xmf << "               </DataItem>\n" ;    
        M_xmf << "           </Attribute>\n" ;
        __evar++ ;
    }
}


template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeStats () const 
{
    size_type numStats = 3 ;
    hsize_t currentSpacesDims[2] = {1, numStats} ;
    
    M_HDF5.createTable ( "stats", H5T_STD_U32BE, currentSpacesDims ) ;
    M_uintBuffer.resize (numStats) ;

    M_uintBuffer[0] = M_maxNumPoints ;
    M_uintBuffer[1] = M_maxNumElements ;
    M_uintBuffer[2] = M_elementNodes ;

    std::cout << "nombre de Points              : " << M_maxNumPoints << std::endl ;
    std::cout << "M_numMaxElements              : " << M_maxNumElements << std::endl ;
    std::cout << "nombre de Points par element  : " << M_elementNodes << std::endl ;
    std::cout << "M_numParts                    : " << M_numParts << std::endl ;
    std::cout << "mesh_type::nRealDim           : " << mesh_type::nRealDim << std::endl ;
    std::cout << "fileNameStep                  : " << M_fileNameStep << ".h5" << std::endl ;
    std::cout << "fileName                      : " << M_fileName << ".xmf" << std::endl ;
    
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
    std::cout << "+---------------+ hdf5 exporter +---------------+" << std::endl ;
    write () ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::visit ( mesh_type* mesh) 
{
    mesh->beginParts () ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::write_xdmf () const 
{

}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::open_xdmf_xml () const
{
    M_xmf.open ((M_fileName+".xmf").c_str(), std::ofstream::out) ;
    M_xmf << "<?xml version=\"1.0\" ?>\n" ;
    M_xmf << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" ;
    M_xmf << "<Xdmf Version=\"2.0\">\n" ;
    M_xmf << "   <Domain>\n" ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::close_xdmf_xml () const
{
    M_xmf << "       </Grid>\n" ;
    M_xmf << "   </Domain>\n" ;
    M_xmf << "</Xdmf>\n"  ;
    M_xmf.close () ;
}
}
#endif /* FEELL_HAS_HDF5 */
#endif /* __Exporterhdf5_H */

