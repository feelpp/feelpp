/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-10-16

  Copyright (C) 2013-2015 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef FEELPP_PARTITIONIO_HPP
#define FEELPP_PARTITIONIO_HPP 1

#if defined(FEELPP_HAS_HDF5)
#include <boost/algorithm/string/split.hpp>
#include <feel/feelcore/hdf5.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp> 

namespace Feel
{
namespace pt =  boost::property_tree;
/**
   \brief Class that handles I/O of mesh parts (for offline partitioning mode)

  This class is used to write the mesh parts produced by (offline) mesh
  partitioning into a single HDF5 container. This part is done with a single
  MPI process (usually on a workstation).

  The class is later used during an online simulation (multiple MPI processes,
  on a cluster, supercomputer etc.) to read the mesh parts.

  Usage:
      - call the constructor (a default one is also available)
      - if the default constructor was used, call setup method to initialize
        the object
      - call write method to store a set of mesh parts in the HDF5 file
      - call read method to load the mesh part associate with the current
        MPI process from the HDF5 file

  Creating, opening and closing the HDF5 file is done automatically by the
  object.

  Description of the storage format:
  N - number of mesh parts

  The HDF5 container contains 6 tables (num. of rows is the leading dimension):

  1. stats - N x 15 (unsigned integer) - each row belongs to a mesh part and
     contains the following values, in order:
         - N
         - graph part size (unused)
         - number of points
         - number of boundary points
         - number of vertices
         - number of boundary vertices
         - number of global vertices
         - number of edges
         - number of boundary edges
         - number of global edges
         - number of faces
         - number of boundary faces
         - number of global faces
         - number of volumes
         - number of global volumes
  2. point_ids - (3 * N) x (maximum number of points in mesh parts) (unsigned
     integer) - each mesh part uses three consecutive rows, with the following
     contents:
         - point markers
         - point global ids
         - point flags
  3. point_coords - (3 * N) x (maximum number of points in mesh parts) (double)
     each mesh part uses three consecutive rows, with the following contents:
         - point x coordinates
         - point y coordinates
         - point z coordinates
  4. edges - (5 * N) x (maximum number of edges in mesh parts) (unsigned
     integer) - each mesh part uses three consecutive rows, with the following
     contents:
         - local ids of first point of edges
         - local ids of second point of edges
         - edge markers
         - edge global ids
         - edge flags
  5. faces ((7 + num of face nodes) * N) x (maximum number of faces in mesh
     parts) (unsigned integer) - each mesh part uses (7 + num of face nodes)
     consecutive rows, with the following contents:
         - one row for the local ids of each face node of the faces in the part
         - face markers
         - face global ids
         - id of first neighbour element
         - id of second neighbour element
         - position of first neighbour element
         - position of second neighbour element
         - face flags
  6. elements ((3 + num of element nodes) * N) x (maximum number of elements
     in mesh parts) (unsigned integer) - each mesh part uses (2 + num of
     element nodes) consecutive rows, with the following contents:
         - one row for the local ids of each element node of the elements in
           the mesh part
         - element markers
         - element global ids
         - element flag

 * @author Radu Popescu <radu.popescu@epfl.ch>
 * @author Christophe Prud'homme (adaptation from LifeV to Feel++)
 * @see
 */
template<typename MeshType>
class PartitionIO
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{
    typedef MeshType mesh_type;
    typedef typename mesh_type::value_type value_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef mesh_type meshparts_type;
    typedef boost::shared_ptr<meshparts_type> meshparts_ptrtype;

    using node_type = typename mesh_type::node_type;
    using point_type = typename mesh_type::point_type;
    using edge_type = typename mesh_type::edge_type;
    using face_type = typename mesh_type::face_type;
    using element_type = typename mesh_type::element_type;


    //@}

    /** @name Constructors & Destructor
     */
    //@{

    //! Default empty constructor
    PartitionIO() {}

    //! Constructor
    /*!
     * Non-empty constructor
     * \param fileName the name of the HDF5 file to be used
     * \param comm WorldComm
     * \param transposeInFile (default false) write/read tables transposed
     *        in the HDF5. This option is available because the transposed
     *        format may be faster on certain machines. Leave as default
     *        unless certain that it works for a given machine.
     */
    PartitionIO (const std::string& fileName,
                 const bool transposeInFile = false);

    //! Empty destructor
    virtual ~PartitionIO() {}
    //@}

    /** @name  Public Methods
     */
    //@{

    //! Initialization method
    /*!
     * Initialization method that is used with the default constructor
     * \param fileName the name of the HDF5 file to be used
     * \param comm WorldComm
     * \param transposeInFile (default false) write/read tables transposed
     *        in the HDF5. This option is available because the transposed
     *        format may be faster on certain machines. Leave as default
     *        unless certain that it works for a given machine.
     */
    void setup (const std::string& fileName,
                const bool transposeInFile = false);
    //! Write method
    /*!
     * Call this method to write the mesh parts to disk
     * \param mesh pointer to a vector containing pointers to the mesh
     *        parts (RegionMesh objects). This pointer is released after
     *        writing, so the mesh parts can be cleared from memory
     */
    void write ( mesh_ptrtype mesh);
    //! Read method
    /*!
     * Call this method to read from the HDF5 file the mesh part associated
     * with the rank of the current MPI process
     * \param mesh pointer to a mesh . If the RegionMesh has been initialized and contains any
     *        state, it will be destroyed before reading
     */
    void read (mesh_ptrtype mesh);


    //@}



protected:

private:
    // Copy constructor and assignment operator are disabled
    PartitionIO (const PartitionIO&);
    PartitionIO& operator= (const PartitionIO&);

    //! Private Methods
    //@{
    
    /**
     * write meta data 
     */
    void writeMetaData( mesh_ptrtype mesh );
    
    /**
     * read meta data
     */
    void readMetaData( mesh_ptrtype mesh );

    // Methods for writing
    void writeStats();
    void writePoints();
    template<typename T> void writeEdges(typename std::enable_if<is_3d<T>::value>::type* = nullptr);
    template<typename T> void writeEdges(typename std::enable_if<mpl::or_<is_2d<T>,is_1d<T>>::value>::type* = nullptr) {}
    void writeFaces();
    void writeElements();
    // Methods for reading
    void readStats();
    void readPoints();
    template<typename T> void readEdges(typename std::enable_if<is_3d<T>::value>::type*  = nullptr);
    template<typename T> void readEdges(typename std::enable_if<mpl::or_<is_2d<T>,is_1d<T>>::value>::type* = nullptr) {}
    void readFaces();
    void readElements();
    
    template<typename IteratorType>
    void writeGhosts( std::string const& tablename, 
                      IteratorType begin, IteratorType end,
                      size_type maxNumEntities );
    template<typename IteratorType, typename ContainerType>
    void readGhosts( std::string const& tablename, 
                     IteratorType begin, IteratorType end,
                     ContainerType& container );
    //@}

    //! Private Data Members
    //@{
    WorldComm M_comm;
    size_type M_myRank;
    std::string M_filename;
    std::string M_h5_filename;
    bool M_transposeInFile;
    meshparts_ptrtype M_meshPartsOut;
    mesh_ptrtype M_meshPartIn;
    // Mesh geometry
    size_type M_elementNodes;
    size_type M_faceNodes;
    size_type M_edgeNodes;
    size_type M_maxNumPoints;
    size_type M_maxNumEdges;
    size_type M_maxNumFaces;
    size_type M_maxNumElements;
    size_type M_numParts;
    // Counters when loading the mesh part
    size_type M_numPoints,M_numGPoints;
    size_type M_numEdges,M_numGEdges;
    size_type M_numFaces,M_numGFaces;
    size_type M_numElements,M_numGElements;
    //HDF5 I/O filter
    HDF5 M_HDF5IO;
    // Buffers for reading/writing
    std::vector<unsigned int> M_uintBuffer;
    std::vector<value_type> M_realBuffer;
//@}
};

template<typename MeshType>
inline PartitionIO<MeshType>::PartitionIO (const std::string& fileName,
                                           const bool transposeInFile) :
    M_filename (fileName),
    M_h5_filename (),
    M_transposeInFile (transposeInFile),
    M_elementNodes(0),
    M_faceNodes(0),
    M_maxNumPoints (0),
    M_maxNumEdges (0),
    M_maxNumFaces (0),
    M_maxNumElements (0)
{

}

template<typename MeshType>
inline void PartitionIO<MeshType>::setup (const std::string& fileName,
                                          const bool transposeInFile)
{
    M_filename = fileName;
    M_h5_filename;
    M_transposeInFile = transposeInFile;
    M_maxNumPoints = 0;
    M_maxNumEdges = 0;
    M_maxNumFaces = 0;
    M_maxNumElements = 0;

    M_elementNodes = 0;
    M_faceNodes = 0;
}

template<typename MeshType>
void PartitionIO<MeshType>::write (mesh_ptrtype meshParts)
{
    if ( Environment::isMasterRank() )
        writeMetaData( meshParts );
    M_h5_filename = fs::path(M_filename).stem().string() + ".h5";
    LOG(INFO) << "writing mesh in HDF5 format in " << M_h5_filename;
    M_meshPartsOut = meshParts;
    // 1 partition per process
    // in the future we have to handle many partitions in one process
    M_numParts = 1;//M_meshParts->size();

    tic();
    M_HDF5IO.openFile (M_h5_filename, meshParts->worldComm(), false);
    writeStats();
    tic();
    writePoints();
    toc("writing points",FLAGS_v>0);
    writeEdges<MeshType>();
    tic();
    writeFaces();
    toc("writing faces",FLAGS_v>0);
    tic();
    writeElements();
    toc("writing elements",FLAGS_v>0);
    
    M_HDF5IO.closeFile();
    toc("writing hdf5 file",FLAGS_v>0);

}
template<typename MeshType>
void PartitionIO<MeshType>::writeMetaData (mesh_ptrtype meshParts)
{
    pt::ptree pt;
    fs::path p(M_filename);
    pt.put("mesh.h5",p.stem().string()+".h5");
    for( auto m : meshParts->markerNames() )
    {
        pt.put("mesh.physicals."+m.first, m.second );
    }
    
    pt::write_json(M_filename, pt);
}
template<typename MeshType>
void PartitionIO<MeshType>::read (mesh_ptrtype meshParts)
{
    readMetaData( meshParts );
    M_meshPartIn = meshParts;
    M_numParts = meshParts->worldComm().localSize();
    tic();
    M_HDF5IO.openFile (M_h5_filename, meshParts->worldComm(), true);
    readStats();
    tic();
    readPoints();
    toc("reading points",FLAGS_v>0);

    readEdges<MeshType>();
    
    tic();
    readFaces();
    toc("reading faces",FLAGS_v>0);
    
    tic();
    readElements();
    toc("reading elements",FLAGS_v>0);

    M_HDF5IO.closeFile();
    toc("reading hdf5 file",FLAGS_v>0);
    tic();
    M_meshPartIn->components().set( MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_CHECK );
    M_meshPartIn->updateForUse();
    toc("mesh update for use",FLAGS_v>0);
}
template<typename MeshType>
void PartitionIO<MeshType>::readMetaData (mesh_ptrtype meshParts)
{
    pt::ptree pt;
    pt::read_json(M_filename, pt);
    M_h5_filename = pt.get<std::string>("mesh.h5");
    auto physicals = pt.get_child_optional("mesh.physicals");
    if ( physicals )
    {
        for (auto& item : pt.get_child("mesh.physicals"))
        {
            std::string v = item.second.data();//item.get<std::string>();
            LOG(INFO) << "name: " << item.first << " " << v;
            typedef std::vector< std::string > split_vector_type;
    
            split_vector_type SplitVec; 
            boost::split( SplitVec, v, boost::is_any_of(" "), boost::token_compress_on ); 
            std::vector<size_type> m;
            std::for_each( SplitVec.begin(), SplitVec.end(), 
                           [&m]( std::string const& s )
                           {
                              m.push_back( std::stoi( s ) );  
                           } );
            meshParts->addMarkerName( std::make_pair( item.first, m ) );
            
        }
    }
}
template<typename MeshType>
void PartitionIO<MeshType>::writeStats()
{
    M_numParts = M_meshPartsOut->worldComm().localSize();
    // Write mesh partition stats (N = number of parts)
    // This is an N x 15 table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = M_numParts;
        currentSpaceDims[1] = 15;
        currentCount[0] = 1;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = 15;
        currentSpaceDims[1] = M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 1;
    }

    // Create new table
    M_HDF5IO.createTable ("stats", H5T_STD_U32BE, currentSpaceDims);

    // Fill buffer
    M_uintBuffer.resize (15);
    //for (size_type i = 0; i < M_numParts; ++i)
    {
        //mesh_type& currentPart = (* (*M_meshPartsOut) [i]);
        auto & currentPart = *M_meshPartsOut;
        M_uintBuffer[0] = M_meshPartsOut->worldComm().localRank();
        M_uintBuffer[1] = 0;
        M_uintBuffer[2] = currentPart.numPoints();
        M_uintBuffer[3] = currentPart.numGlobalPoints();
        M_uintBuffer[4] = currentPart.numVertices();
        M_uintBuffer[5] = 0;
        M_uintBuffer[6] = currentPart.numGlobalVertices();
        M_uintBuffer[7] = currentPart.numEdges();
        M_uintBuffer[8] = 0;
        M_uintBuffer[9] = currentPart.numGlobalEdges();
        M_uintBuffer[10] = currentPart.numFaces();
        M_uintBuffer[11] = 0;
        M_uintBuffer[12] = currentPart.numGlobalFaces();
        M_uintBuffer[13] = currentPart.numElements();
        M_uintBuffer[14] = currentPart.numGlobalElements();

        hsize_t currentOffset[2];
        if (! M_transposeInFile)
        {
            currentOffset[0] = currentPart.worldComm().localRank();
            currentOffset[1] = 0;
        }
        else
        {
            currentOffset[0] = 0;
            currentOffset[1] = currentPart.worldComm().localRank();
        }

        M_HDF5IO.write ("stats", H5T_NATIVE_UINT, currentCount, currentOffset,
                        &M_uintBuffer[0]);
    }
    M_HDF5IO.closeTable ("stats");
}

template<typename MeshType>
void PartitionIO<MeshType>::writePoints()
{

    // Write points (N = number of parts)
    // Two tables: one is (3 * N) x max_num_points of int
    //             second is (d * N) x max_num_points of real
    hsize_t currentSpaceDims[2],currentSpaceDims1[2],currentSpaceDimsGhost[2];
    hsize_t currentCount[2],currentCount1[2],currentCountGhost[2];
    int d = M_meshPartsOut->nRealDim;
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = d * M_numParts;
        currentSpaceDims[1] = M_meshPartsOut->maxNumPoints();
        currentSpaceDims1[0] = 5 * M_numParts;
        currentSpaceDims1[1] = M_meshPartsOut->maxNumPoints();
        currentCount[0] = d;

        currentCount[1] = currentSpaceDims[1];
        currentCount1[0] = 5;
        currentCount1[1] = currentSpaceDims[1];
        
        
    }
    else
    {
        currentSpaceDims[0] = M_meshPartsOut->maxNumPoints();
        currentSpaceDims[1] = d * M_numParts;
        currentSpaceDims1[0] = M_meshPartsOut->maxNumPoints();
        currentSpaceDims1[1] = 5 * M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = d;
        currentCount1[0] = currentSpaceDims[0];
        currentCount1[1] = 5;
    }

    M_HDF5IO.createTable ("point_ids", H5T_STD_U32BE, currentSpaceDims1);
    M_HDF5IO.createTable ("point_coords", H5T_IEEE_F64BE, currentSpaceDims);


    // Fill buffer
    M_uintBuffer.resize (currentCount1[0] * currentCount1[1], 0);
    M_realBuffer.resize (currentCount[0] * currentCount[1], 0);

    size_type stride = currentCount[1], strideIds = currentCount1[1];
    
    
    if (! M_transposeInFile)
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            //mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            auto & currentPart = *M_meshPartsOut;
            auto pt_it = currentPart.beginPoint();
            auto pt_en = currentPart.endPoint();
            for( int j = 0 ; pt_it != pt_en; ++pt_it, ++j )
            {
                M_uintBuffer[0 * strideIds + j] = pt_it->id();
                M_uintBuffer[1 * strideIds + j] = pt_it->marker().value();
                M_uintBuffer[2 * strideIds + j] = pt_it->processId();
                M_uintBuffer[3 * strideIds + j] = pt_it->isOnBoundary();
                M_uintBuffer[4 * strideIds + j] = pt_it->numberOfPartitions();

                M_realBuffer[j] = pt_it->operator[](0);
                if ( pt_it->nRealDim >= 2 )
                    M_realBuffer[stride + j] = pt_it->operator[](1);
                if ( pt_it->nRealDim >= 3 )
                    M_realBuffer[2 * stride + j] = pt_it->operator[](2);
                
            }
            
            hsize_t currentOffset1[2] = {currentPart.worldComm().localRank()* currentCount1[0], 0};
            hsize_t currentOffset2[2] = {currentPart.worldComm().localRank()* currentCount[0], 0};

            M_HDF5IO.write ("point_ids", H5T_NATIVE_UINT, currentCount1,
                            currentOffset1, &M_uintBuffer[0]);
            M_HDF5IO.write ("point_coords", H5T_NATIVE_DOUBLE, currentCount,
                            currentOffset2, &M_realBuffer[0]);

            writeGhosts( "points_ghosts", currentPart.beginPoint(), currentPart.endPoint(),
                         M_meshPartsOut->maxNumPoints());
        }
    }
    else
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            auto & currentPart = *M_meshPartsOut;
            auto pt_it = currentPart.beginPoint();
            auto pt_en = currentPart.endPoint();
            for( int j = 0; pt_it != pt_en; ++pt_it, ++j )
            {
                M_uintBuffer[stride * j] = pt_it->marker().value();
                M_uintBuffer[stride * j + 1] = pt_it->id();
                M_uintBuffer[stride * j + 2] = 0;
                M_realBuffer[stride * j] = pt_it->operator[](0);
                if ( pt_it->nRealDim >= 2 )
                    M_realBuffer[stride * j + 1] = pt_it->operator[](1);
                if ( pt_it->nRealDim >= 3 )
                    M_realBuffer[stride * j + 2] = pt_it->operator[](2);
            }

            hsize_t currentOffset[2] = {0, currentPart.worldComm().localRank()* currentCount[1]};
            hsize_t currentOffset1[2] = {0, currentPart.worldComm().localRank()* currentCount1[1]};
            M_HDF5IO.write ("point_ids", H5T_NATIVE_UINT, currentCount1,
                            currentOffset1, &M_uintBuffer[0]);
            M_HDF5IO.write ("point_coords", H5T_NATIVE_DOUBLE, currentCount,
                            currentOffset, &M_realBuffer[0]);
        }
    }

    M_realBuffer.resize (0);
    M_uintBuffer.resize (0);



    M_HDF5IO.closeTable ("point_coords");
    M_HDF5IO.closeTable ("point_ids");

}

template<typename MeshType>
template<typename T>
void PartitionIO<MeshType>::writeEdges( typename std::enable_if<is_3d<T>::value>::type* )
{
    // Write edges (N = number of parts)
    // Table is (5 * N) x max_num_edges of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    int I = 7;
    M_edgeNodes = edge_type::numPoints;
    M_maxNumEdges = M_meshPartsOut->maxNumEdges();
    LOG(INFO)<< "writeEdges " << M_edgeNodes << "  " << M_maxNumEdges;
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = (I+M_edgeNodes) * M_numParts;
        currentSpaceDims[1] = M_maxNumEdges;
        currentCount[0] = I+M_edgeNodes;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = M_maxNumEdges;
        currentSpaceDims[1] = (I+M_edgeNodes) * M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = (I+M_edgeNodes);
    }

    M_HDF5IO.createTable ("edges", H5T_STD_U32BE, currentSpaceDims);

    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    size_type stride = currentCount[1];
    if (! M_transposeInFile)
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            auto& currentPart = *M_meshPartsOut;
            auto elt_it = currentPart.beginEdge();
            auto elt_en = currentPart.endEdge();
            CHECK(currentPart.maxNumEdges() >= std::distance( elt_it, elt_en ) )
                << "invalid num edges " << M_maxNumEdges
                << " dist : " << std::distance( elt_it, elt_en ) << std::endl;
            for( int j = 0 ; elt_it != elt_en; ++ elt_it, ++j )
            {
                M_uintBuffer[ 0 * stride + j] = elt_it->id();
                M_uintBuffer[ 1 * stride + j] = elt_it->marker().value();
                M_uintBuffer[ 2 * stride + j] = elt_it->marker2().value();
                M_uintBuffer[ 3 * stride + j] = elt_it->marker3().value();
                M_uintBuffer[ 4 * stride + j] = elt_it->isOnBoundary();
                M_uintBuffer[ 5 * stride + j] = elt_it->processId();
                M_uintBuffer[ 6 * stride + j] = elt_it->numberOfPartitions();

                // TODO : add ghost data

                for (size_type k = 0; k < M_edgeNodes; ++k)
                {
                    M_uintBuffer[ (I+k) * stride + j] = elt_it->point (k).id();
                }
            }
            hsize_t currentOffset[2] = {M_meshPartsOut->worldComm().localRank()* currentCount[0], 0};
            M_HDF5IO.write ("edges", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
    else
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {

            hsize_t currentOffset[2] = {0, M_meshPartsOut->worldComm().localRank()* currentCount[1]};
            M_HDF5IO.write ("edges", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }

    M_HDF5IO.closeTable ("edges");
    //writeGhosts( "edges_ghosts", M_meshPartsOut->beginEdge(), M_meshPartsOut->endEdge(),
    //M_meshPartsOut->maxNumEdges());
        
}

template<typename MeshType>
void PartitionIO<MeshType>::writeFaces()
{
    // Write faces (N = number of parts)
    // Table is ((I + num_face_nodes * N) x max_num_faces of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    int I = 7;
    M_faceNodes = face_type::numPoints;
    M_maxNumFaces = M_meshPartsOut->maxNumFaces();
    LOG(INFO)<< "writeFaces " << M_faceNodes << "  " << M_maxNumFaces;
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = (I + M_faceNodes) * M_numParts;
        currentSpaceDims[1] = M_maxNumFaces;
        currentCount[0] = I + M_faceNodes;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = M_maxNumFaces;
        currentSpaceDims[1] = (I + M_faceNodes) * M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = I + M_faceNodes;
    }

    M_HDF5IO.createTable ("faces", H5T_STD_U32BE, currentSpaceDims);
    LOG(INFO) << "creating face table of size " << currentSpaceDims[0] << "x" << currentSpaceDims[1];
    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    size_type stride = currentCount[1];
    if (! M_transposeInFile)
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            auto& currentPart = *M_meshPartsOut;
            auto elt_it = currentPart.beginFace();
            auto elt_en = currentPart.endFace();
            CHECK(currentPart.maxNumFaces() >= std::distance( elt_it, elt_en ) )
                << "invalid num element " << M_maxNumFaces
                << " dist : " << std::distance( elt_it, elt_en ) << std::endl;
            for( int j = 0 ; elt_it != elt_en; ++ elt_it, ++j )
            {
                M_uintBuffer[ 0 * stride + j] = elt_it->id();
                M_uintBuffer[ 1 * stride + j] = elt_it->marker().value();
                M_uintBuffer[ 2 * stride + j] = elt_it->marker2().value();
                M_uintBuffer[ 3 * stride + j] = elt_it->marker3().value();
                M_uintBuffer[ 4 * stride + j] = elt_it->isOnBoundary();
                M_uintBuffer[ 5 * stride + j] = elt_it->processId();
                M_uintBuffer[ 6 * stride + j] = elt_it->numberOfPartitions();

                // TODO : add ghost data

                for (size_type k = 0; k < M_faceNodes; ++k)
                {
                    M_uintBuffer[ (I+k) * stride + j] = elt_it->point (k).id();
                }

            }

            
            hsize_t currentOffset[2] = {M_meshPartsOut->worldComm().localRank()* currentCount[0], 0};
            M_HDF5IO.write ("faces", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
    else
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            hsize_t currentOffset[2] = {0, M_meshPartsOut->worldComm().localRank()* currentCount[1]};
            M_HDF5IO.write ("faces", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }

    M_HDF5IO.closeTable ("faces");
    //writeGhosts( "faces_ghosts", M_meshPartsOut->beginFace(), M_meshPartsOut->endFace(),
    //M_meshPartsOut->maxNumFaces());
}

template<typename MeshType>
void PartitionIO<MeshType>::writeElements()
{
    // Write elements (N = number of parts)
    // Table is ((2 + num_element_nodes * N) x max_num_elements of int
    M_elementNodes = mesh_type::element_type::numLocalVertices;
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    const int nptinfo = 7;
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = (nptinfo + M_elementNodes) * M_numParts;
        currentSpaceDims[1] = M_meshPartsOut->maxNumElements();
        currentCount[0] = nptinfo + M_elementNodes;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = M_meshPartsOut->maxNumElements();
        currentSpaceDims[1] = (nptinfo + M_elementNodes) * M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = nptinfo + M_elementNodes;
    }

    M_HDF5IO.createTable ("elements", H5T_STD_U32BE, currentSpaceDims);

    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    size_type stride = currentCount[1];

    if (! M_transposeInFile)
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            auto& currentPart = *M_meshPartsOut;
            auto elt_it = currentPart.beginElement();
            auto elt_en = currentPart.endElement();
            CHECK(currentPart.maxNumElements() >= std::distance( elt_it, elt_en ) )
                << "invalid num element " << M_maxNumElements
                << " dist : " << std::distance( elt_it, elt_en ) << std::endl;
            for( int j = 0 ; elt_it != elt_en; ++ elt_it, ++j )
            {
                M_uintBuffer[ 0 * stride + j] = elt_it->id();
                M_uintBuffer[ 1 * stride + j] = elt_it->marker().value();
                M_uintBuffer[ 2 * stride + j] = elt_it->marker2().value();
                M_uintBuffer[ 3 * stride + j] = elt_it->marker3().value();
                M_uintBuffer[ 4 * stride + j] = elt_it->isOnBoundary();
                M_uintBuffer[ 5 * stride + j] = elt_it->processId();
                M_uintBuffer[ 6 * stride + j] = elt_it->numberOfPartitions();

                // TODO : add ghost data

                for (size_type k = 0; k < M_elementNodes; ++k)
                {
                    M_uintBuffer[ (nptinfo+k) * stride + j] = elt_it->point (k).id();
                }

            }

            hsize_t currentOffset[2] = { currentPart.worldComm().localRank() * currentCount[0], 0};
            M_HDF5IO.write ("elements", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
    else
    {
        for (size_type i = 0; i < M_numParts; ++i)
        {

#if 0
            hsize_t currentOffset[2] = {0, i* currentCount[1]};
            M_HDF5IO.write ("elements", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
#endif
        }

    }

    M_HDF5IO.closeTable ("elements");
    writeGhosts( "elements_ghosts", M_meshPartsOut->beginElement(), M_meshPartsOut->endElement(),
                 M_meshPartsOut->maxNumElements());
}

template<typename MeshType>
void PartitionIO<MeshType>::readStats()
{
    // Write mesh partition stats (N = number of parts)
    // This is an N x 15 table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];

    M_HDF5IO.openTable ("stats", currentSpaceDims);

    if (! M_transposeInFile)
    {
        currentCount[0] = 1;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = M_meshPartIn->worldComm().localRank();
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 1;
        currentOffset[0] = 0;
        currentOffset[1] = M_meshPartIn->worldComm().localRank();
    }

    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("stats", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("stats");

    M_numPoints = M_uintBuffer[2];
    M_numEdges = M_uintBuffer[7];
    M_numFaces = M_uintBuffer[10];
    M_numElements = M_uintBuffer[13];
    
    M_numGPoints = M_uintBuffer[3];
    M_numGEdges = M_uintBuffer[9];
    M_numGFaces = M_uintBuffer[12];
    M_numGElements = M_uintBuffer[14];

    M_meshPartIn->setNumGlobalElements({M_numGPoints,M_numGEdges,M_numGFaces,M_numGElements});
}

template<typename MeshType>
void PartitionIO<MeshType>::readPoints()
{
    // Read mesh points (N = number of parts)
    // There are two tables: a (3 * N) x max_num_points table of int and
    // a (3 * N) x max_num_points table of real
    hsize_t currentSpaceDims[2],currentSpaceDimsIds[2];
    hsize_t currentCount[2],currentCountIds[2];
    hsize_t currentOffset[2],currentOffsetIds[2];

    M_HDF5IO.openTable ("point_ids", currentSpaceDimsIds);
    LOG(INFO) << "loaded points_ids:" << currentSpaceDimsIds[0] << "x" << currentSpaceDimsIds[1];
    M_HDF5IO.openTable ("point_coords", currentSpaceDims);
    LOG(INFO) << "loaded points_coords:" << currentSpaceDims[0] << "x" << currentSpaceDims[1];
    int d = M_meshPartIn->nRealDim;
    if (! M_transposeInFile)
    {
        // Ids
        currentCountIds[0] = 5;
        currentCountIds[1] = currentSpaceDimsIds[1];
        currentOffsetIds[0] = currentCountIds[0] * M_meshPartIn->worldComm().localRank();
        currentOffsetIds[1] = 0;
        // coords
        currentCount[0] = d;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_meshPartIn->worldComm().localRank();
        currentOffset[1] = 0;
    }
    else
    {
        // Ids
        currentCountIds[0] = currentSpaceDimsIds[0];
        currentCountIds[1] = 5;
        currentOffsetIds[0] = 0;
        currentOffsetIds[1] = currentCountIds[1] * M_meshPartIn->worldComm().localRank();
        // coords
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = d;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_meshPartIn->worldComm().localRank();
    }

    M_uintBuffer.resize (currentCountIds[0] * currentCountIds[1], 0);
    M_realBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("point_ids", H5T_NATIVE_UINT, currentCountIds, currentOffsetIds,
                   &M_uintBuffer[0]);
    M_HDF5IO.read ("point_coords", H5T_NATIVE_DOUBLE, currentCount,
                   currentOffset, &M_realBuffer[0]);

    M_HDF5IO.closeTable ("point_ids");
    M_HDF5IO.closeTable ("point_coords");

    size_type stride = currentCount[1], strideIds=currentCountIds[1];
    node_type coords( d );
    if (! M_transposeInFile)
    {
        for (size_type j = 0; j < M_numPoints; ++j)
        {
            //
            int id        = M_uintBuffer[0*strideIds + j];
            int marker    = M_uintBuffer[1*strideIds + j];
            int pid       = M_uintBuffer[2*strideIds + j];
            bool onbdy    = M_uintBuffer[3*strideIds + j];
            int nparts    = M_uintBuffer[4*strideIds + j];

            // ghost table
            //std           = M_uintBuffer[4*stride + j];

            for( int c = 0; c < d; ++c )
            {
                coords[c] = M_realBuffer[c*stride+j];
            }

            point_type pt( id, coords, onbdy );
            pt.setMarker( marker );
            pt.setProcessId( pid );
            pt.setProcessIdInPartition( M_meshPartIn->worldComm().localRank() );
            pt.setNumberOfPartitions( nparts );
            M_meshPartIn->addPoint( pt );
            
        }
    }
    else
    {
        for (size_type j = 0; j < M_numPoints; ++j)
        {

        }
    }
    M_realBuffer.resize (0);
    readGhosts("points_ghosts", M_meshPartIn->beginPoint(), M_meshPartIn->endPoint(), M_meshPartIn->points());
}

template<typename MeshType>
template<typename T>
void PartitionIO<MeshType>::readEdges(typename std::enable_if<is_3d<T>::value>::type*  )
{
    // Read mesh edges (N = number of parts)
    // Read a (5 * N) x max_num_edges table of int and
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];
    int I = 7;
    M_edgeNodes = edge_type::numPoints;

    M_HDF5IO.openTable ("edges", currentSpaceDims);
    LOG(INFO) << "Loading hdf5 mesh edges table " << currentSpaceDims[0] << " x " << currentSpaceDims[1];
    if (! M_transposeInFile)
    {
        currentCount[0] = I+M_edgeNodes;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_meshPartIn->worldComm().localRank();
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = I+M_edgeNodes;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_meshPartIn->worldComm().localRank();
    }

    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("edges", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("edges");


    size_type stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (size_type j = 0; j < M_numEdges; ++j)
        {
            size_type id = M_uintBuffer[0*stride+j];
            int marker   = M_uintBuffer[1*stride+j];
            int marker2  = M_uintBuffer[2*stride+j];
            int marker3  = M_uintBuffer[3*stride+j];
            bool onbdy   = M_uintBuffer[4*stride+j];
            int pid      = M_uintBuffer[5*stride+j];
            int npart    = M_uintBuffer[6*stride+j];


            edge_type e;
            e.setId( id );
            e.setProcessIdInPartition( M_meshPartIn->worldComm().localRank() );
            e.setMarker( marker );
            e.setMarker2( marker2 );
            e.setMarker3( marker3 );
            e.setOnBoundary( onbdy );
            e.setNumberOfPartitions( npart );
            e.setProcessId( pid );

            // TODO: Ghost data
            //e.setNeighborPartitionIds( __e.ghosts );

            for (size_type k = 0; k < M_edgeNodes; ++k)
            {
                //M_uintBuffer[ (4+k) * stride + j] = elt_it->point (k).id();
                DCHECK( M_uintBuffer[ (I+k) * stride + j] <= M_numGPoints ) 
                    << "element " << e.id() 
                    <<  " Invalid Point Id " << M_uintBuffer[ (I+k) * stride + j] 
                    << " with local pt id " << k << " from buffer id " 
                    << (I+k) * stride + j;
                
                DCHECK( M_meshPartIn->point( M_uintBuffer[ (I+k) * stride + j]).node().size() ==  M_meshPartIn->nRealDim ) 
                    << " Invalid Pt " << M_uintBuffer[ (I+k) * stride + j] << " : " 
                    << M_meshPartIn->point( M_uintBuffer[ (I+k) * stride + j]);
                
                e.setPoint( k, M_meshPartIn->point( M_uintBuffer[ (I+k) * stride + j]) );

            }
            M_meshPartIn->addEdge( e );

        }
    }
    else
    {
        for (size_type j = 0; j < M_numEdges; ++j)
        {

        }
    }
    //readGhosts("edges_ghosts", M_meshPartIn->beginEdge(), M_meshPartIn->endEdge(), M_meshPartIn->edges());
}


template<typename MeshType>
void PartitionIO<MeshType>::readFaces()
{
    // read mesh faces (N = number of parts)
    // Read a ((7 + num_face_points) * N) x max_num_faces table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];
    int I = 7;
    M_faceNodes = face_type::numPoints;
    M_maxNumFaces = M_meshPartIn->maxNumFaces();
    LOG(INFO)<< "read Faces " << M_faceNodes << "  " << M_maxNumFaces;
    M_HDF5IO.openTable ("faces", currentSpaceDims);
    LOG(INFO) << "Loading hdf5 mesh faces table " << currentSpaceDims[0] << " x " << currentSpaceDims[1];
    if (! M_transposeInFile)
    {
        currentCount[0] = I + M_faceNodes;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_meshPartIn->worldComm().localRank();
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = I + M_faceNodes;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_meshPartIn->worldComm().localRank();
    }

    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("faces", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("faces");

    size_type stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (size_type j = 0; j < M_numFaces; ++j)
        {
            size_type id = M_uintBuffer[0*stride+j];
            int marker   = M_uintBuffer[1*stride+j];
            int marker2  = M_uintBuffer[2*stride+j];
            int marker3  = M_uintBuffer[3*stride+j];
            bool onbdy   = M_uintBuffer[4*stride+j];
            int pid      = M_uintBuffer[5*stride+j];
            int npart    = M_uintBuffer[6*stride+j];


            face_type e;
            e.setId( id );
            e.setProcessIdInPartition( M_meshPartIn->worldComm().localRank() );
            e.setMarker( marker );
            e.setMarker2( marker2 );
            e.setMarker3( marker3 );
            e.setOnBoundary( onbdy );
            e.setNumberOfPartitions( npart );
            e.setProcessId( pid );

            // TODO: Ghost data
            //e.setNeighborPartitionIds( __e.ghosts );

            for (size_type k = 0; k < M_faceNodes; ++k)
            {
                //M_uintBuffer[ (4+k) * stride + j] = elt_it->point (k).id();
                DCHECK( M_uintBuffer[ (I+k) * stride + j] <= M_numGPoints ) 
                    << "element " << e.id() 
                    <<  " Invalid Point Id " << M_uintBuffer[ (I+k) * stride + j] 
                    << " with local pt id " << k << " from buffer id " 
                    << (I+k) * stride + j;
                
                DCHECK( M_meshPartIn->point( M_uintBuffer[ (I+k) * stride + j]).node().size() ==  M_meshPartIn->nRealDim ) 
                    << " Invalid Pt " << M_uintBuffer[ (I+k) * stride + j] << " : " 
                    << M_meshPartIn->point( M_uintBuffer[ (I+k) * stride + j]);
                
                e.setPoint( k, M_meshPartIn->point( M_uintBuffer[ (I+k) * stride + j]) );

            }
            M_meshPartIn->addFace( e );

        }
    }
    else
    {
        for (size_type j = 0; j < M_numFaces; ++j)
        {
        }
    }
    //readGhosts("faces_ghosts", M_meshPartIn->beginFace(), M_meshPartIn->endFace(), M_meshPartIn->faces());

}

template<typename MeshType>
void PartitionIO<MeshType>::readElements()
{
    // Read mesh elements (N = number of parts)
    // Read a ((3 + num_element_points) * N) x max_num_elements table of int
    M_elementNodes = mesh_type::element_type::numLocalVertices;
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];
    const int nptinfo = 7;

    M_HDF5IO.openTable ("elements", currentSpaceDims);
    LOG(INFO) << "loaded elements:" << currentSpaceDims[0] << "x" << currentSpaceDims[1];
    if (! M_transposeInFile)
    {
        currentCount[0] = nptinfo + M_elementNodes;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_meshPartIn->worldComm().localRank();
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = nptinfo + M_elementNodes;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_meshPartIn->worldComm().localRank();
    }

    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("elements", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("elements");

    size_type stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (size_type j = 0; j < M_numElements; ++j)
        {
            size_type id = M_uintBuffer[0*stride+j];
            int marker   = M_uintBuffer[1*stride+j];
            int marker2  = M_uintBuffer[2*stride+j];
            int marker3  = M_uintBuffer[3*stride+j];
            bool onbdy   = M_uintBuffer[4*stride+j];
            int pid      = M_uintBuffer[5*stride+j];
            int npart    = M_uintBuffer[6*stride+j];

            element_type e;
            e.setId( id );
            e.setProcessIdInPartition( M_meshPartIn->worldComm().localRank() );
            e.setMarker( marker );
            e.setMarker2( marker2 );
            e.setMarker3( marker3 );
            e.setOnBoundary( onbdy );
            e.setNumberOfPartitions( npart );
            e.setProcessId( pid );
            // TODO: Ghost data
            //e.setNeighborPartitionIds( __e.ghosts );

            for (size_type k = 0; k < M_elementNodes; ++k)
            {
                //M_uintBuffer[ (4+k) * stride + j] = elt_it->point (k).id();
                DCHECK( M_uintBuffer[ (nptinfo+k) * stride + j] <= M_numGPoints ) << "element " << e.id() <<  " Invalid Point Id " << M_uintBuffer[ (nptinfo+k) * stride + j] << " with local pt id " << k << " from buffer id " << (nptinfo+k) * stride + j;
                
                DCHECK( M_meshPartIn->point( M_uintBuffer[ (nptinfo+k) * stride + j]).node().size() ==  M_meshPartIn->nRealDim ) 
                    << " Invalid Pt " << M_uintBuffer[ (nptinfo+k) * stride + j] << " : " 
                    << M_meshPartIn->point( M_uintBuffer[ (nptinfo+k) * stride + j]);
                
                e.setPoint( k, M_meshPartIn->point( M_uintBuffer[ (nptinfo+k) * stride + j]) );

            }
            M_meshPartIn->addElement( e, false );
        }
    }
    else
    {
        for (size_type j = 0; j < M_numElements; ++j)
        {
        }
    }
    readGhosts("elements_ghosts", M_meshPartIn->beginElement(), M_meshPartIn->endElement(), M_meshPartIn->elements());

}
template<typename MeshType>
template<typename IteratorType>
void PartitionIO<MeshType>::writeGhosts( std::string const& tablename, 
                                         IteratorType begin, IteratorType end,
                                         size_type maxNumEntities )
{
    using entity_type = typename IteratorType::value_type;
    int topodim = entity_type::nDim;
    auto it = begin;
    size_type max_neighbors = 0;
    size_type max_num_entity = 0;
    DLOG(INFO) << "number of entities for table " << tablename << " :" << std::distance(begin,end) << "  max num entities " << maxNumEntities;
    for( ; it != end; ++ it )
    {
        CHECK( it->numberOfNeighborPartitions() == it->numberOfPartitions()-1 ) << "number of partitions to which the point belong is invalid " << it->numberOfNeighborPartitions() << " vs " << it->numberOfPartitions()-1;
        DLOG_IF(INFO,(topodim==1 || topodim==2) && (it->numberOfNeighborPartitions()  || it->idInOthersPartitions().size() ) ) << " id " << it->id() << " "<< it->numberOfNeighborPartitions() << "  " << it->idInOthersPartitions().size();


        if ( it->numberOfNeighborPartitions() >= 1 )
        {
            max_neighbors = (max_neighbors > it->numberOfNeighborPartitions() )?max_neighbors:it->numberOfNeighborPartitions();
            max_num_entity++;
        }

        
    }
    LOG(INFO) << "max number of neighbor partitions:" << max_neighbors << " max num pt " << max_num_entity;
    std::vector<size_type> globals{max_neighbors,max_num_entity};
    mpi::all_reduce( M_meshPartsOut->worldComm(), mpi::inplace(globals.data()), 2, mpi::maximum<size_type>() );
    hsize_t currentSpaceDimsGhost[2];
    hsize_t currentCountGhost[2];
    //hsize_t currentOffsetGhost[2];
    currentSpaceDimsGhost[0] = (3+3*globals[0])*M_numParts;
    currentSpaceDimsGhost[1] = globals[1];
    currentCountGhost[0] = (3+3*globals[0]);
    currentCountGhost[1] = currentSpaceDimsGhost[1];
    LOG(INFO) << "creating ghost table " << tablename << " dims:" << currentSpaceDimsGhost[0] << " x " << currentSpaceDimsGhost[1] << " globals[0]:" << globals[0];
    M_HDF5IO.createTable (tablename, H5T_STD_U32BE, currentSpaceDimsGhost);
    
    hsize_t currentOffsetGhost[2]={M_meshPartsOut->worldComm().localRank()* currentCountGhost[0], 0};
    std::vector<unsigned int> ghost_info( currentCountGhost[0]*currentCountGhost[1], 0 );
    
    
    size_type strideGhost = currentCountGhost[1];
    it = begin;
    for( int j = 0 ; it != end; ++it )
    {
        DCHECK(it->numberOfNeighborPartitions() == it->numberOfPartitions()-1) 
            << "Invalid partition data for entity id " << it->id() << " neighbord partitions " 
            << it->numberOfNeighborPartitions() << "  number of partitions " << it->numberOfPartitions() 
            << " should be equal, -1 for the number of partitions";
        
        if ( it->numberOfNeighborPartitions() >= 1 )
        {
            DCHECK(it->numberOfNeighborPartitions() >= it->idInOthersPartitions().size() )
                << " Invalid partitioning data size npids "  
                << " ID:" << it->id()
                << " current part:" << it->pidInPartition() << " pid:" << it->processId()
                << " size:"<< it->numberOfNeighborPartitions() << " != " << it->idInOthersPartitions().size()
                << " neigh:" << it->neighborPartitionIds()
                << " iop:" << it->idInOthersPartitions();
            ghost_info[0*strideGhost+j] = it->id();
            ghost_info[1*strideGhost+j] = it->numberOfNeighborPartitions();
            ghost_info[2*strideGhost+j] = it->idInOthersPartitions().size();
            auto const& npids = it->neighborPartitionIds();
            
            for(int n = 0; n < it->numberOfNeighborPartitions();++n)
            {
                ghost_info[(0*globals[0]+n+3)*strideGhost+j] = npids[n];
                
            }
            int n = 0;
            for( auto const& i: it->idInOthersPartitions())
            {
                ghost_info[(1*globals[0]+n+3)*strideGhost+j] = i.first;
                ghost_info[(2*globals[0]+n+3)*strideGhost+j] = i.second;
                n++;
            }
            ++j;
        }
    }
    
    M_HDF5IO.write (tablename, H5T_NATIVE_UINT, currentCountGhost,
                    currentOffsetGhost, &ghost_info[0]);

    M_HDF5IO.closeTable (tablename);

}    

template<typename MeshType>
template<typename IteratorType, typename ContainerType>
void PartitionIO<MeshType>::readGhosts( std::string const& tablename, 
                                        IteratorType begin, IteratorType end,
                                        ContainerType& container )
{
    auto it = begin;
    hsize_t currentSpaceDimsGhost[2]{0,0};
    M_HDF5IO.openTable (tablename, currentSpaceDimsGhost);
    LOG(INFO) << "loading ghost table " << tablename << " dims:" << currentSpaceDimsGhost[0] << " x " << currentSpaceDimsGhost[1] << " for " << M_numParts << "  partitions";
    
    hsize_t currentCountGhost[2]{0,0};
    
    currentCountGhost[0] = currentSpaceDimsGhost[0]/M_numParts;
    currentCountGhost[1] = currentSpaceDimsGhost[1];
    LOG(INFO) << "current  ghost : " << currentCountGhost[0] << " x " << currentCountGhost[1];

    hsize_t maxnumghost = (currentCountGhost[0]-3)/3;
    LOG(INFO) << "readGhosts: max num ghosts : " << maxnumghost;
    hsize_t currentOffsetGhost[2]={M_meshPartIn->worldComm().localRank()* currentCountGhost[0], 0};
    
    std::vector<unsigned int> ghost_info( currentCountGhost[0]*currentCountGhost[1], 0 );
    if ( !ghost_info.size() ) 
    {
        M_HDF5IO.closeTable (tablename);
        return;
    }
    M_HDF5IO.read (tablename, H5T_NATIVE_UINT, currentCountGhost,
                    currentOffsetGhost, &ghost_info[0]);
    M_HDF5IO.closeTable (tablename);
    
    size_type strideGhost = currentCountGhost[1];
    it = begin;
    
    for( int j = 0 ; it != end; ++it )
    {
        if ( it->numberOfPartitions() > 1 )
        {
            int id = ghost_info[0*strideGhost+j];
            int nnpids = ghost_info[1*strideGhost+j];
            int nids = ghost_info[2*strideGhost+j];
            std::vector<rank_type> npids;
            std::map<rank_type,size_type> iop;
            for(int n = 0; n < nnpids;++n)
            {
                
                npids.push_back(ghost_info[(0*maxnumghost+n+3)*strideGhost+j]);
            }
            for(int n = 0; n < nids;++n)
            {
                rank_type p = ghost_info[(1*maxnumghost+n+3)*strideGhost+j];
                size_type io = ghost_info[(2*maxnumghost+n+3)*strideGhost+j];
                DLOG(INFO) << "id " << id << " part "  << n << "/" << it->numberOfPartitions()-1 
                           << " npid: " << ghost_info[(0*maxnumghost+n+3)*strideGhost+j] 
                           << " p: " << p << " io: " << io;
                iop[p]=io;
                DCHECK(npids.size()>=iop.size()) << "Invalid size " << npids.size() << " != " << iop.size()
                                                 << " maxnumghost: " << maxnumghost
                                                 << " id: " << id
                                                 << " n:" << n
                                                 << " npids:" << npids << " iop:" << iop;
            }
            using entity_type = typename IteratorType::value_type;
            container.modify(it,
                             [&npids,&iop]( entity_type& e ) 
                             { 
                                 e.setNeighborPartitionIds( npids ); 
                                 e.setIdInOtherPartitions( iop );
                             } );
            DCHECK(it->numberOfNeighborPartitions() == it->numberOfPartitions()-1) << "[read] Invalid partition data for entity id " << it->id() << " neighbord partitions " << it->numberOfNeighborPartitions() << "  number of partitions " << it->numberOfPartitions() << " should be equal, -1 for the number of partitions";
            DCHECK(it->numberOfNeighborPartitions() >= it->idInOthersPartitions().size() )
                << "[read] Invalid partitioning data size npids "  << it->numberOfNeighborPartitions() << " != " << it->idInOthersPartitions().size();
            ++j;
        }
        DCHECK(j <= currentCountGhost[1] ) << " Invalid element index " << j << " > " << currentCountGhost[1];
    }
    

}    


} // Feel



#endif /* FEELPP_HAS_HDF5 */

#if 0

BOOST_PARAMETER_FUNCTION( ( void ),
                          h5partition,                                       // 2. name of the function template
                          tag,                                        // 3. namespace of tag types
                          ( required                                  // 4. one required parameter, and
                            ( mesh, * )
                              ) // required
                          ( optional                                  // 4. one required parameter, and
                            ( name,  *, Environment::about().appName() )
                            ( worldcomm, *, Environment::worldComm() )
                            ( path, *( boost::is_convertible<mpl::_,std::string> ), std::string(".") )
                          ) )
{
    typedef typename Feel::detail::partition<Args>::type mesh_type;
    std::string filename = prefix + "_mesh.h5";
    PartitionIO<mesh_type> p( filename, worldcomm );
    p.write(mesh);
}

BOOST_PARAMETER_FUNCTION( ( void ),
                          h5read,                                       // 2. name of the function template
                          tag,                                        // 3. namespace of tag types
                          ( required                                  // 4. one required parameter, and
                            ( in_out(mesh), * )
                              ) // required
                          ( optional                                  // 4. one required parameter, and
                            ( prefix,  *, Environment::about().appName() )
                            ( worldcomm, *, Environment::worldComm() )
                            ( path, *( boost::is_convertible<mpl::_,std::string> ), std::string(".") )
                          ) )
{
    typedef typename Feel::detail::partition<Args>::type mesh_type;
    std::string filename = prefix + "_mesh.h5";
    PartitionIO<mesh_type> p( filename, worldcomm );
    p.read(mesh);
}
#endif

#endif /* __PartitionIO_H */
