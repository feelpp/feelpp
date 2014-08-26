#ifndef __Exporterhdf5_H
#define __Exporterhdf5_H 1

#include <feel/feelcore/feel.hpp>

#include <iostream>
#include <math.h>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/timeset.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/detail/fileindex.hpp>
#include <map>
#include <fstream>
#include <boost/lambda/lambda.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#if defined(FEELPP_HAS_HDF5)
#include <feel/feelcore/hdf5.hpp>

/*!
 * \file exporterhdf5.hpp
 * \brief HDF5 and XDMF exporter
 * \author VANTHONG Benjamin
 */

namespace Feel 
{
namespace fs = boost::filesystem;

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
        /*!
         * \brief Fonction used in almost all constructor to initialize the element's type
         */
        void init() ;
        
        /*!
         * \brief write .xmf and .h5 files for each process
         */
        void write ()  const ;
        
        /*!
         * \brief only one write .xmf file and one .h5 file for each time step  
         */
        void writeMerge ()  const ;
        
        /*!
         * \brief write points' coordonates   
         */
        void writePoints () const ;

        /*!
         * \brief write points' coordonates (merge version)
         */
        void writePointsMerge () const ;
        
        /*!
         * \brief write elements (an element is formed by several nodes) 
         */
        void writeElements () const ;
        
        /*!
         * \brief write elements (merge version) 
         */
        void writeElementsMerge () const ;

        /*!
         * \brief write informations of the mesh in .h5 file (unused for now)
         */
        void writeStats () const ;

        /*!
         * \brief save solutions on nodes 
         * \param __step a time step
         * \param __var  iterator on solutions (begin)
         * \param en     iterator on solutions (end)
         */
        template<typename Iterator>
        void saveNodal ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const ;
        template<typename Iterator>

        /*!
         * \brief save solutions on nodes (merge version)
         * \param __step a time step
         * \param __var  iterator on solutions (begin)
         * \param en     iterator on solutions (end)
         */
        void saveNodalMerge ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const ;
        template<typename Iterator>
        
        /*!
         * \brief save solutions on elements  )
         * \param __step   a time step
         * \param __evar   iterator on solutions (begin)
         * \param __evaren iterator on solutions (end)
         */
        void saveElement ( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const ;
        template<typename Iterator>

        /*!
         * \brief save solutions on elements (merge version) 
         * \param __step   a time step
         * \param __evar   iterator on solutions (begin)
         * \param __evaren iterator on solutions (end)
         */
        void saveElementMerge ( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const ;

        /*!
         * \brief open XDMF file (used only in version without merge)
         */
        void open_xdmf_xml () const ;
        
        /*!
         * \brief close properly XDMF file (used only in version without merge)
         */
        void close_xdmf_xml () const ;

        
        /*!
         * \brief a bubble sort of 2 arrays in the same time
         * \param ids    array of points' identifier
         * \param coords coordinates of points
         * \param n      size of those array
         */
        void bubbleSort (size_type * ids, value_type * coords, size_type n) const ;

        mutable WorldComm M_comm ;              /*!< MPI worldComm */
        mutable std::string M_fileName ;        /*!< file name */
        mutable std::string M_fileNameStep ;    /*!< file name + time step */
        mutable HDF5 M_HDF5 ;                   /*!< HDF5 IO */
        mutable mesh_ptrtype M_meshOut ;        /*!< pointer on current mesh in a time step */

        // Mesh geometry
        mutable size_type M_elementNodes ;      /*!< number of nodes for one element */ 
        mutable size_type M_maxNumElements ;    /*!< number of elements for the current process */
        mutable size_type M_maxNumPoints ;      /*!< number of points for the current process */
        mutable size_type M_numParts ;          /*!< number of partitions */
        mutable std::string M_element_type ;    /*!< element's type */

        mutable size_type M_step = 0 ;          /*!< number of the current step */
        mutable std::ofstream M_xmf  ;          /*!< Out stream to write the .xmf file */

        mutable std::vector<size_type> M_uintBuffer ;           /*!< buffer of integer */
        mutable std::vector<value_type> M_realBuffer ;          /*!< buffer of double */
        mutable std::map<size_type, size_type> M_newPointId ;   /*!< new point identifier after sort */
        mutable MPI_File fh ;                                   /*!< file descriptor of the .xmf file (merge version only) */
        mutable std::ostringstream M_str ;                      /*!< buffer of string */
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
void Exporterhdf5<MeshType, N>::writeMerge () const 
{

    MPI_Status status;
    int size ;
    M_str <<  this->prefix () << ".xmf" ;
    M_fileName = M_str.str () ; 

    char * strTmp = strdup (M_fileName.c_str()) ;
    if(this->worldComm().isMasterRank() && fs::exists(strTmp))
    {
        MPI_File_delete(strTmp, MPI_INFO_NULL);
    }
    MPI_Barrier( Environment::worldComm().comm() );
    MPI_File_open( this->worldComm().comm(), strTmp, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
    free (strTmp) ;

    size_type rank = Environment::worldComm().globalRank () ;

    timeset_const_iterator __ts_it = this->beginTimeSet () ;
    timeset_const_iterator __ts_en = this->endTimeSet () ;
    
    timeset_ptrtype __ts = *__ts_it ;

    M_str.str("") ;
    M_str << "<?xml version=\"1.0\" ?>\n" ;
    M_str << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" ;
    M_str << "<Xdmf Version=\"2.0\">\n" ;
    M_str << "   <Domain>\n" ;
    M_str << "       <Grid Name=\"Simulation\" GridType=\"Collection\" CollectionType=\"Spatial\">\n" ;

    char * buffer  = (char *) malloc (M_str.str().length()*sizeof(char) + 1) ;
    if ( this->worldComm().isMasterRank() ) 
    { size = M_str.str().length() ; }
    else
    { size = 0 ; }
    M_str.str().copy (buffer, size, 0) ;
    MPI_File_write_ordered (fh, buffer, size, MPI_CHAR, &status) ;

    if ( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
        std::cout << "file generated                : " << M_fileName << std::endl ;

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
                M_fileNameStep = this->prefix() + filestr.str() ;
                M_HDF5.openFile (M_fileNameStep+".h5", Environment::worldComm(), false) ;

                M_str.str("") ;
                M_str << "           <Grid Name=\"" << M_fileNameStep << "\" GridType=\"Uniform\">\n" ;
                M_str << "               <Time Value=\"" << __step->time() << "\"/>\n" ;  

                writePointsMerge () ;
                writeElementsMerge () ;
                
                if ( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
                {
                    //std::cout << "time                          : " << __step->time () << std::endl ; 
                    //std::cout << "time increment                : " << __ts->timeIncrement () << std::endl ; 
                    //std::cout << "numberOfSteps                 : " << __ts->numberOfSteps () << std::endl ; 
                    //std::cout << "numberOfTotalSteps            : " << __ts->numberOfTotalSteps () << std::endl ;
                    std::cout << "file generated                : " << M_fileNameStep  << ".h5" << std::endl ;
                }

                saveNodalMerge (__step, __step->beginNodalScalar(), __step->endNodalScalar() ) ;
                saveNodalMerge (__step, __step->beginNodalVector(), __step->endNodalVector() ) ;
                saveNodalMerge( __step, __step->beginNodalTensor2(), __step->endNodalTensor2() );

                saveElementMerge (__step, __step->beginElementScalar(), __step->endElementScalar() ) ;
                saveElementMerge (__step, __step->beginElementVector(), __step->endElementVector() ) ;
                saveElementMerge( __step, __step->beginElementTensor2(), __step->endElementTensor2() );

                M_str << "           </Grid>\n" ;
                free (buffer) ;
                buffer  = (char *) malloc (M_str.str().length()*sizeof(char) + 1) ;
                size = M_str.str().length() ; 
                M_str.str().copy (buffer, size, 0) ;
                MPI_File_write_ordered (fh, buffer, size, MPI_CHAR, &status) ;

                M_HDF5.closeFile () ;
            }
            ++__it ;
        }
        ++__ts_it ;
    }
    M_str.str("") ;
    M_str << "       </Grid>\n" ;
    M_str << "   </Domain>\n" ;
    M_str << "</Xdmf>\n"  ;

    free (buffer) ;
    buffer  = (char *) malloc (M_str.str().length()*sizeof(char) + 1) ;
    if ( this->worldComm().isMasterRank() ) 
    { size = M_str.str().length() ; }
    else
    { size = 0 ; }
    M_str.str().copy (buffer, size, 0) ;
    MPI_File_write_ordered (fh, buffer, size, MPI_CHAR, &status) ;

    MPI_File_close(&fh) ;
    free (buffer) ;
}


template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::write () const 
{

    std::ostringstream str ;
    str <<  this->prefix () << "-" << Environment::worldComm().globalSize()<<"_"<<Environment::worldComm().globalRank()  ;
    M_fileName = str.str () ; 

    std::cout << "file generated                : " << M_fileName << ".xmf" << std::endl ;
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

                if ( Environment::worldCommSeq().globalRank() == Environment::worldCommSeq().masterRank() )
                {
                    //std::cout << "time                          : " << __step->time () << std::endl ; 
                    //std::cout << "time increment                : " << __ts->timeIncrement () << std::endl ; 
                    //std::cout << "numberOfSteps                 : " << __ts->numberOfSteps () << std::endl ; 
                    //std::cout << "numberOfTotalSteps            : " << __ts->numberOfTotalSteps () << std::endl ;
                    std::cout << "file generated                : " << M_fileNameStep  << ".h5" << std::endl ;
                }

                saveNodal (__step, __step->beginNodalScalar(), __step->endNodalScalar() ) ;
                saveNodal (__step, __step->beginNodalVector(), __step->endNodalVector() ) ;
                saveNodal( __step, __step->beginNodalTensor2(), __step->endNodalTensor2() );

                saveElement (__step, __step->beginElementScalar(), __step->endElementScalar() ) ;
                saveElement (__step, __step->beginElementVector(), __step->endElementVector() ) ;
                saveElement( __step, __step->beginElementTensor2(), __step->endElementTensor2() );

                M_xmf << "           </Grid>\n" ;
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
void Exporterhdf5<MeshType, N>::writePointsMerge () const 
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
    std::ostringstream filestr0 ;
    filestr0 << Environment::worldComm().globalRank () ;
    std::vector <size_type> globalNumPoint ;
    globalNumPoint.resize (Environment::worldComm().globalSize(), 0) ;
    boost::mpi::all_gather(Environment::worldComm(), M_maxNumPoints, globalNumPoint) ;
    for (size_type i = 0 ; i < Environment::worldComm().globalSize () ; i++)
    {
        std::ostringstream filestr ;
        filestr << i  ;
        currentCount[0] = globalNumPoint[i] ;
        M_HDF5.createTable (filestr.str(), "point_coords", H5T_IEEE_F64BE, currentCount, false) ;
    }
    currentCount[0] = M_maxNumPoints ;
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

    Environment::worldComm().barrier() ;
    M_HDF5.write (filestr0.str()+"point_coords", H5T_NATIVE_DOUBLE, currentCount, currentOffset, &M_realBuffer[0]) ;
    Environment::worldComm().barrier() ;

    for (size_type i = 0 ; i < Environment::worldComm().globalSize() ; i++)
    {
        std::ostringstream filestr1 ;
        filestr1 << i ;
        M_HDF5.closeTable(filestr1.str()+"point_coords") ;
    }

    M_str << "               <Geometry GeometryType=\"XYZ\">\n" ;
    M_str << "                   <DataItem Dimensions=\"" << M_maxNumPoints << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n" ;
    M_str << "                   " << M_fileNameStep << ".h5:/" << filestr0.str() << "/point_coords\n" ;
    M_str << "                   </DataItem>\n" ;
    M_str << "               </Geometry>\n" ;
}

template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::writeElementsMerge () const 
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

    std::ostringstream filestr0 ;
    filestr0 << Environment::worldComm().globalRank () ;
    std::vector <size_type> globalNumElements ;
    globalNumElements.resize (Environment::worldComm().globalSize(), 0) ;
    boost::mpi::all_gather(Environment::worldComm(), M_maxNumElements, globalNumElements) ;

    for (size_type i = 0 ; i < Environment::worldComm().globalSize () ; i++)
    {
        std::ostringstream filestr ;
        filestr << i  ;
        currentSpacesDims [0] = globalNumElements[i] ;

        M_HDF5.createTable (filestr.str(), "element_nodes", H5T_STD_U32BE, currentSpacesDims, true) ;
    }
        currentSpacesDims [0] = M_maxNumElements ;

    M_uintBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0) ;
    std::vector<size_type> idsBuffer ;
    idsBuffer.resize (currentSpacesDims2[1], 0) ;

    M_uintBuffer.resize (currentSpacesDims[0]*currentSpacesDims[1], 0) ;

    size_type k = 0 ;
    size_type i = 0 ;
    for (p_it = M_meshOut->beginParts (); k < M_numParts ;  p_it++ , k++)
    {
        auto elt_it = M_meshOut->beginElementWithMarker (p_it->first) ;
        auto elt_en = M_meshOut->endElementWithMarker (p_it->first) ;
        for ( ; elt_it != elt_en ; ++elt_it , i ++)
        {
            idsBuffer[i] = elt_it->id () ;
            for ( size_type j = 0 ; j < M_elementNodes ; j ++ )
                M_uintBuffer[j + M_elementNodes*i] = M_newPointId[elt_it->point(j).id()]  ; 
        }
    }


    hsize_t currentOffset[2] = {0, 0} ;
    M_HDF5.write ( filestr0.str()+"element_nodes", H5T_NATIVE_LLONG, currentSpacesDims, currentOffset, &M_uintBuffer[0] ) ;

    for (size_type i = 0 ; i < Environment::worldComm().globalSize() ; i++)
    {
        std::ostringstream filestr1 ;
        filestr1 << i ;
        M_HDF5.closeTable(filestr1.str()+"element_nodes") ;
    }
    M_str << "               <Topology TopologyType=\"" << M_element_type << "\" NumberOfElements=\"" << M_maxNumElements << "\" NodesPerElement=\"" << M_elementNodes << "\">\n" ;
    M_str << "                   <DataItem Dimensions=\"" <<M_maxNumElements << " " << M_elementNodes << "\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n" ;
    M_str << "                   " << M_fileNameStep << ".h5:/" << filestr0.str() << "/element_nodes\n" ;
    M_str << "                   </DataItem>\n" ;
    M_str << "               </Topology>\n" ;
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
    for (p_it = M_meshOut->beginParts (); k < M_numParts ;  p_it++ , k++)
    {
        auto elt_it = M_meshOut->beginElementWithMarker (p_it->first) ;
        auto elt_en = M_meshOut->endElementWithMarker (p_it->first) ;
        for ( ; elt_it != elt_en ; ++elt_it , i ++)
        {
            idsBuffer[i] = elt_it->id () ;
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
    M_xmf << "               <DataItem Dimensions=\"" <<M_maxNumElements << " " << M_elementNodes << "\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n" ;
    M_xmf << "               " << M_fileNameStep << ".h5:/element_nodes\n" ;
    M_xmf << "               </DataItem>\n" ;
    M_xmf << "           </Topology>\n" ;
}

template <typename MeshType, int N>
template <typename Iterator>
void Exporterhdf5<MeshType, N>::saveNodalMerge ( typename timeset_type::step_ptrtype __step, Iterator __var, Iterator en ) const 
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

//        if ( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
//            std::cout << "solution name                 : " << solutionName << std::endl ;

        hsize_t currentSpacesDims [2] ;

        currentSpacesDims [0] = nComponents ;
        currentSpacesDims [1] = M_maxNumPoints ;

        std::ostringstream filestr0 ;
        filestr0 << Environment::worldComm().globalRank () ;

        std::vector <size_type> globalNumPoint ;
        globalNumPoint.resize (Environment::worldComm().globalSize(), 0) ;
        boost::mpi::all_gather(Environment::worldComm(), M_maxNumPoints, globalNumPoint) ;

        for (size_type i = 0 ; i < Environment::worldComm().globalSize () ; i++)
        {
            std::ostringstream filestr ;
            filestr << i  ;
            currentSpacesDims [1] = globalNumPoint[i] ;

            M_HDF5.createTable (filestr.str(), solutionName.c_str() , H5T_IEEE_F64BE, currentSpacesDims, true) ;
        }
        
        currentSpacesDims [1] = M_maxNumPoints ;

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

        M_HDF5.write ( filestr0.str() + solutionName, H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] ) ;

        for (size_type i = 0 ; i < Environment::worldComm().globalSize() ; i++)
        {
            std::ostringstream filestr1 ;
            filestr1 << i ;
            M_HDF5.closeTable(filestr1.str()+solutionName) ;
        }

        M_str << "               <Attribute AttributeType=\""<< attributeType << "\" Name=\"" << solutionName << "\" Center=\"Node\">\n" ;
        M_str << "                   <DataItem Dimensions=\""<< nComponents << " " << M_maxNumPoints << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n" ;
        M_str << "                   "<< M_fileNameStep <<".h5:/"<< filestr0.str() << "/" <<solutionName <<"\n" ;
        M_str << "                   </DataItem>\n" ;    
        M_str << "               </Attribute>\n" ;
        ++__var ;
    }
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
void Exporterhdf5<MeshType, N>::saveElementMerge ( typename timeset_type::step_ptrtype __step, Iterator __evar, Iterator __evaren ) const 
{
    while ( __evar != __evaren ) 
    {
        std::string attributeType ("Scalar") ;

        std::string solutionName = __evar->first ;
//        if ( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
//            std::cout << "solution Name                 : " << solutionName << std::endl ;
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

        std::ostringstream filestr0 ;
        filestr0 << Environment::worldComm().globalRank () ;
        std::vector <size_type> globalNumElements ;
        globalNumElements.resize (Environment::worldComm().globalSize(), 0) ;
        boost::mpi::all_gather(Environment::worldComm(), M_maxNumElements, globalNumElements) ;

        for (size_type i = 0 ; i < Environment::worldComm().globalSize () ; i++)
        {
            std::ostringstream filestr ;
            filestr << i  ;
            currentSpacesDims [1] = globalNumElements[i] ;

            M_HDF5.createTable (filestr.str(), solutionName, H5T_IEEE_F64BE, currentSpacesDims, true) ;
        }
        currentSpacesDims [1] = M_maxNumElements ;

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

        M_HDF5.write ( filestr0.str() + solutionName, H5T_NATIVE_DOUBLE, currentSpacesDims, currentOffset, &M_realBuffer[0] ) ;

        for (size_type i = 0 ; i < Environment::worldComm().globalSize() ; i++)
        {
            std::ostringstream filestr1 ;
            filestr1 << i ;
            M_HDF5.closeTable(filestr1.str()+solutionName) ;
        }

        M_str << "               <Attribute AttributeType=\"" << attributeType << "\" Name=\"" << solutionName << "\" Center=\"Cell\">\n" ;
        M_str << "                   <DataItem Dimensions=\""<< nComponents <<" "<< M_maxNumElements << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Endian=\"Big\">\n" ;
        M_str << "                   "<< M_fileNameStep << ".h5:/" << filestr0.str() << "/" << solutionName << "\n" ;
        M_str << "                   </DataItem>\n" ;    
        M_str << "               </Attribute>\n" ;
        __evar++ ;
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

    if ( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
    {
        std::cout << "nombre de Points              : " << M_maxNumPoints << std::endl ;
        std::cout << "M_numMaxElements              : " << M_maxNumElements << std::endl ;
        std::cout << "nombre de Points par element  : " << M_elementNodes << std::endl ;
        std::cout << "M_numParts                    : " << M_numParts << std::endl ;
        std::cout << "mesh_type::nRealDim           : " << mesh_type::nRealDim << std::endl ;
        std::cout << "fileNameStep                  : " << M_fileNameStep << ".h5" << std::endl ;
        std::cout << "fileName                      : " << M_fileName << ".xmf" << std::endl ;
    }

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
    if ( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
        std::cout << "exporter.merge                : " << (boption (_name = "exporter.merge" ) ? "true" : "false")  << std::endl ;
    MPI_Barrier( Environment::worldComm().comm() );
    if ( boption ( _name = "exporter.merge" ) )
        writeMerge() ;
    else 
        write () ;
}

    template <typename MeshType, int N>
void Exporterhdf5<MeshType, N>::visit ( mesh_type* mesh) 
{
    mesh->beginParts () ;
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

