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
#include <feel/feelcore/hdf5.hpp>

namespace Feel
{
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
    // Methods for writing
    void writeStats();
    void writePoints();
    void writeEdges();
    void writeFaces();
    void writeElements();
    // Methods for reading
    void readStats();
    void readPoints();
    void readEdges();
    void readFaces();
    void readElements();
    //@}

    //! Private Data Members
    //@{
    WorldComm M_comm;
    size_type M_myRank;
    std::string M_fileName;
    bool M_transposeInFile;
    meshparts_ptrtype M_meshPartsOut;
    mesh_ptrtype M_meshPartIn;
    // Mesh geometry
    size_type M_elementNodes;
    size_type M_faceNodes;
    size_type M_maxNumPoints;
    size_type M_maxNumEdges;
    size_type M_maxNumFaces;
    size_type M_maxNumElements;
    size_type M_numParts;
    // Counters when loading the mesh part
    size_type M_numPoints;
    size_type M_numEdges;
    size_type M_numFaces;
    size_type M_numElements;
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
    M_fileName (fileName),
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
    M_fileName = fileName;
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
    M_meshPartsOut = meshParts;
    // 1 partition per process
    // in the future we have to handle many partitions in one process
    M_numParts = 1;//M_meshParts->size();

    M_HDF5IO.openFile (M_fileName, meshParts->worldComm(), false);
    writeStats();
    writePoints();
    //writeEdges();
    //writeFaces();
    writeElements();
    M_HDF5IO.closeFile();


}

template<typename MeshType>
void PartitionIO<MeshType>::read (mesh_ptrtype meshParts)
{
    M_meshPartIn = meshParts;

    M_HDF5IO.openFile (M_fileName, meshParts->worldComm(), true);
    readStats();
    readPoints();
    readEdges();
    readFaces();
    readElements();
    M_HDF5IO.closeFile();


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
        M_uintBuffer[3] = 0;//currentPart.numBPoints();
        M_uintBuffer[4] = currentPart.numVertices();
        M_uintBuffer[5] = 0;//currentPart.numBVertices();
        M_uintBuffer[6] = currentPart.numGlobalVertices();
        M_uintBuffer[7] = currentPart.numEdges();
        M_uintBuffer[8] = 0;//currentPart.numBEdges();
        M_uintBuffer[9] = 0;//currentPart.numGlobalEdges();
        M_uintBuffer[10] = currentPart.numFaces();
        M_uintBuffer[11] = 0;//currentPart.numBFaces();
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
    hsize_t currentSpaceDims[2],currentSpaceDims1[2];
    hsize_t currentCount[2],currentCount1[2];
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
void PartitionIO<MeshType>::writeEdges()
{
    // Write edges (N = number of parts)
    // Table is (5 * N) x max_num_edges of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = 5 * M_numParts;
        currentSpaceDims[1] = M_maxNumEdges;
        currentCount[0] = 5;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = M_maxNumEdges;
        currentSpaceDims[1] = 5 * M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 5;
    }

    M_HDF5IO.createTable ("edges", H5T_STD_U32BE, currentSpaceDims);
#if 0
    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    size_type stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (size_type i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (size_type j = 0; j < currentPart.numEdges(); ++j)
            {
                M_uintBuffer[j] = currentPart.edgeList[j].point (0).localId();
                M_uintBuffer[stride + j] =
                    currentPart.edgeList[j].point (1).localId();
                M_uintBuffer[2 * stride + j] =
                    currentPart.edgeList[j].markerID();
                M_uintBuffer[3 * stride + j] = currentPart.edgeList[j].id();
                M_uintBuffer[4 * stride + j] =
                    static_cast<int> (currentPart.edgeList[j].flag() );
            }

            hsize_t currentOffset[2] = {i* currentCount[0], 0};
            M_HDF5IO.write ("edges", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
    else
    {
        for (size_type i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (size_type j = 0; j < currentPart.numEdges(); ++j)
            {
                M_uintBuffer[stride * j] =
                    currentPart.edgeList[j].point (0).localId();
                M_uintBuffer[stride * j + 1] =
                    currentPart.edgeList[j].point (1).localId();
                M_uintBuffer[stride * j + 2] =
                    currentPart.edgeList[j].markerID();
                M_uintBuffer[stride * j + 3] = currentPart.edgeList[j].id();
                M_uintBuffer[stride * j + 4] =
                    static_cast<int> (currentPart.edgeList[j].flag() );
            }

            hsize_t currentOffset[2] = {0, i* currentCount[1]};
            M_HDF5IO.write ("edges", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
#endif
    M_HDF5IO.closeTable ("edges");
}

template<typename MeshType>
void PartitionIO<MeshType>::writeFaces()
{
    // Write faces (N = number of parts)
    // Table is ((7 + num_face_nodes * N) x max_num_faces of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = (7 + M_faceNodes) * M_numParts;
        currentSpaceDims[1] = M_maxNumFaces;
        currentCount[0] = 7 + M_faceNodes;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = M_maxNumFaces;
        currentSpaceDims[1] = (7 + M_faceNodes) * M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 7 + M_faceNodes;
    }

    M_HDF5IO.createTable ("faces", H5T_STD_U32BE, currentSpaceDims);
#if 0
    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    size_type stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (size_type i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (size_type j = 0; j < currentPart.numFaces(); ++j)
            {
                for (size_type k = 0; k < M_faceNodes; ++k)
                {
                    M_uintBuffer[k * stride + j] =
                        currentPart.faceList[j].point (k).localId();
                }
                M_uintBuffer[M_faceNodes * stride + j] =
                    currentPart.faceList[j].markerID();
                M_uintBuffer[ (M_faceNodes + 1) * stride + j] =
                    currentPart.faceList[j].id();
                M_uintBuffer[ (M_faceNodes + 2) * stride + j] =
                    currentPart.faceList[j].firstAdjacentElementIdentity();
                M_uintBuffer[ (M_faceNodes + 3) * stride + j] =
                    currentPart.faceList[j].secondAdjacentElementIdentity();
                M_uintBuffer[ (M_faceNodes + 4) * stride + j] =
                    currentPart.faceList[j].firstAdjacentElementPosition();
                M_uintBuffer[ (M_faceNodes + 5) * stride + j] =
                    currentPart.faceList[j].secondAdjacentElementPosition();
                M_uintBuffer[ (M_faceNodes + 6) * stride + j] =
                    static_cast<int> (currentPart.faceList[j].flag() );
            }

            hsize_t currentOffset[2] = {i* currentCount[0], 0};
            M_HDF5IO.write ("faces", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
    else
    {
        for (size_type i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (size_type j = 0; j < currentPart.numFaces(); ++j)
            {
                for (size_type k = 0; k < M_faceNodes; ++k)
                {
                    M_uintBuffer[stride * j + k] =
                        currentPart.faceList[j].point (k).localId();
                }
                M_uintBuffer[stride * j + M_faceNodes] =
                    currentPart.faceList[j].markerID();
                M_uintBuffer[stride * j + M_faceNodes + 1] =
                    currentPart.faceList[j].id();
                M_uintBuffer[stride * j + M_faceNodes + 2] =
                    currentPart.faceList[j].firstAdjacentElementIdentity();
                M_uintBuffer[stride * j + M_faceNodes + 3] =
                    currentPart.faceList[j].secondAdjacentElementIdentity();
                M_uintBuffer[stride * j + M_faceNodes + 4] =
                    currentPart.faceList[j].firstAdjacentElementPosition();
                M_uintBuffer[stride * j + M_faceNodes + 5] =
                    currentPart.faceList[j].secondAdjacentElementPosition();
                M_uintBuffer[stride * j + M_faceNodes + 6] =
                    static_cast<int> (currentPart.faceList[j].flag() );
            }

            hsize_t currentOffset[2] = {0, i* currentCount[1]};
            M_HDF5IO.write ("faces", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
#endif
    M_HDF5IO.closeTable ("faces");
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
                    M_uintBuffer[ (4+k) * stride + j] = elt_it->point (k).id();
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

    // Insert stats into mesh partition object
    M_numPoints = M_uintBuffer[2];
    //M_meshPartIn->setMaxNumPoints (M_uintBuffer[2], true);
    //M_meshPartIn->setNumBPoints (M_uintBuffer[3]);

    // Vertices
    //M_meshPartIn->setNumVertices (M_uintBuffer[4]);
    //M_meshPartIn->setNumBVertices (M_uintBuffer[5]);
    //M_meshPartIn->setNumGlobalVertices (M_uintBuffer[6]);

    // Edges
    M_numEdges = M_uintBuffer[7];
    //M_meshPartIn->setNumEdges (M_uintBuffer[7]);
    //M_meshPartIn->setMaxNumEdges (M_uintBuffer[7], true);
    //M_meshPartIn->setNumBEdges (M_uintBuffer[8]);
    //M_meshPartIn->setMaxNumGlobalEdges (M_uintBuffer[9]);

    // Faces
    M_numFaces = M_uintBuffer[10];
    //M_meshPartIn->setNumFaces (M_uintBuffer[10]);
    //M_meshPartIn->setMaxNumFaces (M_uintBuffer[10], true);
    //M_meshPartIn->setNumBFaces (M_uintBuffer[11]);
    //M_meshPartIn->setMaxNumGlobalFaces (M_uintBuffer[12]);

    // Volumes
    M_numElements = M_uintBuffer[13];
    //M_meshPartIn->setMaxNumVolumes (M_uintBuffer[13], true);
    //M_meshPartIn->setMaxNumGlobalVolumes (M_uintBuffer[14]);

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

    //hid_t realDataset = H5Dopen(M_fileId, "point_coords", H5P_DEFAULT);

    M_HDF5IO.openTable ("point_ids", currentSpaceDimsIds);
    M_HDF5IO.openTable ("point_coords", currentSpaceDims);
    int d = M_meshPartsOut->nRealDim;
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
    M_realBuffer.resize (currentCountIds[0] * currentCountIds[1], 0);
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

}

template<typename MeshType>
void PartitionIO<MeshType>::readEdges()
{
    // Read mesh edges (N = number of parts)
    // Read a (5 * N) x max_num_edges table of int and
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];

    M_HDF5IO.openTable ("edges", currentSpaceDims);

    if (! M_transposeInFile)
    {
        currentCount[0] = 5;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_myRank;
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 5;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_myRank;
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

        }
    }
    else
    {
        for (size_type j = 0; j < M_numEdges; ++j)
        {

        }
    }

}


template<typename MeshType>
void PartitionIO<MeshType>::readFaces()
{
    // read mesh faces (N = number of parts)
    // Read a ((7 + num_face_points) * N) x max_num_faces table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];

    M_HDF5IO.openTable ("faces", currentSpaceDims);

    if (! M_transposeInFile)
    {
        currentCount[0] = 7 + M_faceNodes;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_myRank;
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 7 + M_faceNodes;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_myRank;
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

        }
    }
    else
    {
        for (size_type j = 0; j < M_numFaces; ++j)
        {
        }
    }


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
                e.setPoint( k, M_meshPartIn->point( M_uintBuffer[ (nptinfo+k) * stride + j]) );

            }
            M_meshPartIn->addElement( e );
        }
    }
    else
    {
        for (size_type j = 0; j < M_numElements; ++j)
        {
        }
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
