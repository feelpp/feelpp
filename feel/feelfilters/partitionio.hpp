/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-10-16

  Copyright (C) 2013-2016 Feel++ Consortium

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

  1. stats - N x 19 (unsigned integer) - each row belongs to a mesh part and
     contains the following values, in order:
         - partition id
         - number of global points (in mesh file)
         - number of local points (in current partition)
         - points offset (in current partition)
         - number of global active elements (in mesh file)
         - number of local active elements (in current partition)
         - active elements offset (in current partition)
         - number of global ghost elements (in mesh file)
         - number of local ghost elements (in current partition)
         - ghost elements offset (in current partition)
         - number of global marked faces (in mesh file)
         - number of local marked faces (in current partition)
         - marked faces offset (in current partition)
         - number of global marked edges (in mesh file)
         - number of local marked edges (in current partition)
         - marked edges offset (in current partition)
         - number of global marked points (in mesh file)
         - number of local marked points (in current partition)
         - marked points offset (in current partition)
  2. point_ids - 1 x (number of global points) (unsigned integer):
         - point id
  3. point_coords - 1 x (3*number of global points) (double):
         - point x coordinates
         - point y coordinates
         - point z coordinates
  4. elements - 1 x ( (2 + num of element nodes)* number of global active elements ):
         - element id
         - element markers
         - points ids
  5. elements_ghosts - 1 x (3* number of global ghost elements ):
         - element id in current partition
         - process id
         - element id in active partition
  5. elements_ghosts - 1 x ( ((1+num of face node)*number of global ghost elements ):
  6. marked_subentities - 1 x ( (1+num of face node)* number of global marked faces
     + (1+num of face node)* number of global marked edges + (1+1)*number of global marked points ):
         - marker
         - points ids

 * @author Christophe Prud'homme (adaptation from LifeV to Feel++)
 * @author Vincent Chabannes
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
    void writeElements();
    void writeGhostElements();
    void writeMarkedSubEntities();
    // Methods for reading
    void readStats();
    void readPoints();
    void readElements();
    void readGhostElements();
    void readMarkedSubEntities();
    void prepareUpdateForUse();

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
    size_type M_numParts;
    size_type M_numGlobalPoints, M_numLocalPoints, M_pointsOffSet;
    size_type M_numGlobalElements, M_numLocalElements, M_elementsOffSet;
    size_type M_numGlobalGhostElements, M_numLocalGhostElements, M_ghostElementsOffSet;
    size_type M_numGlobalMarkedFaces, M_numLocalMarkedFaces, M_markedFacesOffSet;
    size_type M_numGlobalMarkedEdges, M_numLocalMarkedEdges, M_markedEdgesOffSet;
    size_type M_numGlobalMarkedPoints, M_numLocalMarkedPoints, M_markedPointsOffSet;

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
    M_transposeInFile (transposeInFile)
{
    std::string h5filename = fs::path(M_filename).stem().string() + ".h5";
    M_h5_filename = (fs::path(M_filename).parent_path() / fs::path( h5filename )).string();
}

template<typename MeshType>
inline void PartitionIO<MeshType>::setup (const std::string& fileName,
                                          const bool transposeInFile)
{
    M_filename = fileName;
    M_h5_filename;
    M_transposeInFile = transposeInFile;
}

template<typename MeshType>
void PartitionIO<MeshType>::write (mesh_ptrtype meshParts)
{
    if ( Environment::isMasterRank() )
        writeMetaData( meshParts );

    //std::string h5filename = fs::path(M_filename).stem().string() + ".h5";
    //M_h5_filename = (fs::path(M_filename).parent_path() / fs::path( h5filename )).string();

    LOG(INFO) << "writing mesh in HDF5 format in " << M_h5_filename;
    M_meshPartsOut = meshParts;
    // 1 partition per process
    // in the future we have to handle many partitions in one process
    M_numParts = 1;//M_meshParts->size();

    tic();
    M_HDF5IO.openFile (M_h5_filename, meshParts->worldComm(), false);
    tic();
    writeStats();
    toc("PartitionIO writing stats",FLAGS_v>0);
    tic();
    writePoints();
    toc("PartitionIO writing points",FLAGS_v>0);
    tic();
    writeElements();
    toc("PartitionIO writing elements",FLAGS_v>0);
    tic();
    writeGhostElements();
    toc("PartitionIO writing ghost_elements",FLAGS_v>0);
    tic();
    writeMarkedSubEntities();
    toc("PartitionIO writing marked_subentities",FLAGS_v>0);

    M_HDF5IO.closeFile();
    toc("PartitionIO writing hdf5 file",FLAGS_v>0);

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

    tic();
    readStats();
    toc("PartitionIO reading stats",FLAGS_v>0);
    tic();
    readPoints();
    toc("PartitionIO reading points",FLAGS_v>0);
    tic();
    readElements();
    toc("PartitionIO reading elements",FLAGS_v>0);
    tic();
    readGhostElements();
    toc("PartitionIO reading ghost_elements",FLAGS_v>0);
    tic();
    readMarkedSubEntities();
    toc("PartitionIO reading marked_subentities",FLAGS_v>0);

    M_HDF5IO.closeFile();
    toc("PartitionIO reading hdf5 file",FLAGS_v>0);

    prepareUpdateForUse();

    tic();
    M_meshPartIn->components().set( MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_CHECK );
    M_meshPartIn->updateForUse();
    toc("PartitionIO mesh update for use",FLAGS_v>0);
}
template<typename MeshType>
void PartitionIO<MeshType>::readMetaData (mesh_ptrtype meshParts)
{
    pt::ptree pt;
    pt::read_json(M_filename, pt);
    //M_h5_filename = pt.get<std::string>("mesh.h5");
    //if ( fs::exist( fs::path(M_filename).parent_path()/ fs::path(M_h5_filename) ) )

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
    // This is an N x 19 table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = M_numParts;
        currentSpaceDims[1] = 19;
        currentCount[0] = 1;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = 19;
        currentSpaceDims[1] = M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 1;
    }

    M_numGlobalPoints = 0;
    M_numLocalPoints = M_meshPartsOut->statNumPointsAll();
    M_pointsOffSet = 0;
    M_numGlobalElements = 0;
    M_numLocalElements = M_meshPartsOut->statNumElementsAll();
    M_elementsOffSet = 0;
    M_numGlobalGhostElements = 0;
    M_numLocalGhostElements = M_meshPartsOut->statNumElementsAll() - M_meshPartsOut->statNumElementsActive();
    M_ghostElementsOffSet = 0;
    M_numGlobalMarkedFaces = 0;
    M_numLocalMarkedFaces = M_meshPartsOut->statNumFacesMarkedAll();
    M_markedFacesOffSet = 0;
    M_numGlobalMarkedEdges = 0;
    M_numLocalMarkedEdges = M_meshPartsOut->statNumEdgesMarkedAll();
    M_markedEdgesOffSet = 0;
    M_numGlobalMarkedPoints = 0;
    M_numLocalMarkedPoints = M_meshPartsOut->statNumPointsMarkedAll();
    M_markedPointsOffSet = 0;
    for ( rank_type p=0 ; p< M_numParts ; ++p )
    {
        size_type nPtAll = M_meshPartsOut->statNumPointsAll(p);
        M_numGlobalPoints += nPtAll;
        size_type nEltAll = M_meshPartsOut->statNumElementsAll(p);
        M_numGlobalElements += nEltAll;
        size_type nGhostEltAll = M_meshPartsOut->statNumElementsAll(p) - M_meshPartsOut->statNumElementsActive(p);
        M_numGlobalGhostElements += nGhostEltAll;
        size_type nMarkedFacesAll = M_meshPartsOut->statNumFacesMarkedAll(p);
        M_numGlobalMarkedFaces += nMarkedFacesAll;
        size_type nMarkedEdgesAll = M_meshPartsOut->statNumEdgesMarkedAll(p);
        M_numGlobalMarkedEdges += nMarkedEdgesAll;
        size_type nMarkedPointsAll = M_meshPartsOut->statNumPointsMarkedAll(p);
        M_numGlobalMarkedPoints += nMarkedPointsAll;
        if ( p < M_meshPartsOut->worldComm().localRank() )
        {
            M_pointsOffSet += nPtAll;
            M_elementsOffSet += nEltAll;
            M_ghostElementsOffSet += nGhostEltAll;
            M_markedFacesOffSet += nMarkedFacesAll;
            M_markedEdgesOffSet += nMarkedEdgesAll;
            M_markedPointsOffSet += nMarkedPointsAll;
        }
    }
    // Create new table
    M_HDF5IO.createTable ("stats", H5T_STD_U32BE, currentSpaceDims);

    // Fill buffer
    M_uintBuffer.resize( 19 );
    //for (size_type i = 0; i < M_numParts; ++i)
    {
        //mesh_type& currentPart = (* (*M_meshPartsOut) [i]);
        auto & currentPart = *M_meshPartsOut;
        M_uintBuffer[0] = M_meshPartsOut->worldComm().localRank();
        M_uintBuffer[1] = M_numGlobalPoints;
        M_uintBuffer[2] = M_numLocalPoints;
        M_uintBuffer[3] = M_pointsOffSet;
        M_uintBuffer[4] = M_numGlobalElements;
        M_uintBuffer[5] = M_numLocalElements;
        M_uintBuffer[6] = M_elementsOffSet;
        M_uintBuffer[7] = M_numGlobalGhostElements;
        M_uintBuffer[8] = M_numLocalGhostElements;
        M_uintBuffer[9] = M_ghostElementsOffSet;
        M_uintBuffer[10] = M_numGlobalMarkedFaces;
        M_uintBuffer[11] = M_numLocalMarkedFaces;
        M_uintBuffer[12] = M_markedFacesOffSet;
        M_uintBuffer[13] = M_numGlobalMarkedEdges;
        M_uintBuffer[14] = M_numLocalMarkedEdges;
        M_uintBuffer[15] = M_markedEdgesOffSet;
        M_uintBuffer[16] = M_numGlobalMarkedPoints;
        M_uintBuffer[17] = M_numLocalMarkedPoints;
        M_uintBuffer[18] = M_markedPointsOffSet;

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
    hsize_t globalDimsIds[2], globalDimsCoords[2];
    hsize_t localDimsIds[2], localDimsCoords[2];
    hsize_t offsetIds[2], offsetCoords[2];

    int d = M_meshPartsOut->nRealDim;
    int nValIds = 1;
    if ( !M_transposeInFile )
    {
        globalDimsIds[0] = 1;
        globalDimsIds[1] = nValIds*M_numGlobalPoints;
        globalDimsCoords[0] = 1;
        globalDimsCoords[1] = d*M_numGlobalPoints;

        localDimsIds[0] = 1;
        localDimsIds[1] = nValIds*M_numLocalPoints;
        localDimsCoords[0] = 1;
        localDimsCoords[1] = d*M_numLocalPoints;
        offsetIds[0] = 0;
        offsetIds[1] = nValIds*M_pointsOffSet;
        offsetCoords[0] = 0;
        offsetCoords[1] = d*M_pointsOffSet;
    }
    else
    {
        CHECK( false ) << "TODO";
    }

    M_HDF5IO.createTable ("point_ids", H5T_STD_U32BE, globalDimsIds);
    M_HDF5IO.createTable ("point_coords", H5T_IEEE_F64BE, globalDimsCoords);


    // Fill buffer
    M_uintBuffer.resize( localDimsIds[0]*localDimsIds[1], 0 );
    M_realBuffer.resize( localDimsCoords[0]*localDimsCoords[1], 0 );

    if ( !M_transposeInFile )
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            //mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            auto & currentPart = *M_meshPartsOut;
            size_type currentBufferIndexIds = 0,currentBufferIndexCoords = 0;
            auto pt_it = currentPart.beginPoint();
            auto pt_en = currentPart.endPoint();
            for( int j = 0 ; pt_it != pt_en; ++pt_it, ++j )
            {
                M_uintBuffer[currentBufferIndexIds++] = pt_it->id();

                M_realBuffer[currentBufferIndexCoords++] = pt_it->operator[](0);
                if ( mesh_type::nRealDim >= 2 )
                    M_realBuffer[currentBufferIndexCoords++] = pt_it->operator[](1);
                if ( mesh_type::nRealDim >= 3 )
                    M_realBuffer[currentBufferIndexCoords++] = pt_it->operator[](2);
            }

            M_HDF5IO.write ("point_ids", H5T_NATIVE_UINT, localDimsIds, offsetIds, &M_uintBuffer[0]);
            M_HDF5IO.write ("point_coords", H5T_NATIVE_DOUBLE, localDimsCoords, offsetCoords, &M_realBuffer[0]);
        }
    }
    else
    {
    }

    M_realBuffer.resize (0);
    M_uintBuffer.resize (0);

    M_HDF5IO.closeTable ("point_coords");
    M_HDF5IO.closeTable ("point_ids");
}

template<typename MeshType>
void PartitionIO<MeshType>::writeElements()
{
    hsize_t globalDims[2];
    hsize_t localDims[2];
    hsize_t offset[2];

    int nVal = 2 + element_type::numPoints;// id+marker+ points
    if (! M_transposeInFile)
    {
        globalDims[0] = 1;
        globalDims[1] = nVal*M_numGlobalElements;
        localDims[0] = 1;
        localDims[1] = nVal*M_numLocalElements;
        offset[0] = 0;
        offset[1] = nVal*M_elementsOffSet;
    }
    else
    {
        CHECK( false ) << "TODO";
    }

    M_HDF5IO.createTable ("elements", H5T_STD_U32BE, globalDims);

    // Fill buffer
    M_uintBuffer.resize( localDims[0]*localDims[1], 0 );

    if (! M_transposeInFile)
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            auto& currentPart = *M_meshPartsOut;
            size_type currentBufferIndex = 0;
            auto elt_it = currentPart.beginElement();
            auto elt_en = currentPart.endElement();
            for( ; elt_it != elt_en; ++ elt_it )
            {
                M_uintBuffer[currentBufferIndex++] = elt_it->id();
                M_uintBuffer[currentBufferIndex++] = elt_it->marker().value();
                for (size_type k = 0; k < element_type::numPoints; ++k)
                {
                    M_uintBuffer[currentBufferIndex++] = elt_it->point(k).id();
                }

            }
            M_HDF5IO.write ("elements", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);
        }
    }
    else
    {
    }

    M_HDF5IO.closeTable ("elements");
}

template<typename MeshType>
void PartitionIO<MeshType>::writeGhostElements()
{
    hsize_t globalDims[2];
    hsize_t localDims[2];
    hsize_t offset[2];

    int nVal = 3;//(idGhost,pid,idActive)
    if (! M_transposeInFile)
    {
        globalDims[0] = 1;
        globalDims[1] = nVal*M_numGlobalGhostElements;
        localDims[0] = 1;
        localDims[1] = nVal*M_numLocalGhostElements;
        offset[0] = 0;
        offset[1] = nVal*M_ghostElementsOffSet;
    }
    else
    {
        CHECK( false ) << "TODO";
    }

    M_HDF5IO.createTable ("elements_ghosts", H5T_STD_U32BE, globalDims);

    M_uintBuffer.resize( localDims[0]*localDims[1], 0 );

    if ( !M_transposeInFile )
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            size_type currentBufferIndex = 0;
            auto& currentPart = *M_meshPartsOut;
            auto elt_it = currentPart.beginGhostElement();
            auto elt_en = currentPart.endGhostElement();
            for( ; elt_it != elt_en; ++ elt_it )
            {
                M_uintBuffer[currentBufferIndex++] = elt_it->id();
                M_uintBuffer[currentBufferIndex++] = elt_it->processId();
                M_uintBuffer[currentBufferIndex++] = elt_it->idInOthersPartitions( elt_it->processId() );
            }
        }
    }
    M_HDF5IO.write ("elements_ghosts", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("elements_ghosts");

}

namespace detail
{
template<typename MeshType>
void updateMarkedSubEntitiesBuffer( MeshType const& mesh,std::vector<unsigned int> & buffer, mpl::int_<1> /**/ )
{
    size_type currentBufferIndex = 0;
    auto point_it = mesh.beginPoint();
    auto point_en = mesh.endPoint();
    for( ; point_it != point_en; ++point_it )
    {
        if ( point_it->marker().isOff() ) continue;
        buffer[currentBufferIndex++] = point_it->marker().value();
        buffer[currentBufferIndex++] = point_it->id();
    }
}
template<typename MeshType>
void updateMarkedSubEntitiesBuffer( MeshType const& mesh,std::vector<unsigned int> & buffer, mpl::int_<2> /**/ )
{
    size_type currentBufferIndex = 0;
    auto face_it = mesh.beginFace();
    auto face_en = mesh.endFace();
    for( ; face_it != face_en; ++face_it )
    {
        if ( face_it->marker().isOff() ) continue;
        buffer[currentBufferIndex++] = face_it->marker().value();
        for ( uint16_type vLocId = 0; vLocId < MeshType::face_type::numPoints ; ++vLocId )
        {
            buffer[currentBufferIndex++] = face_it->point( vLocId ).id();
        }
    }
    auto point_it = mesh.beginPoint();
    auto point_en = mesh.endPoint();
    for( ; point_it != point_en; ++point_it )
    {
        if ( point_it->marker().isOff() ) continue;
        buffer[currentBufferIndex++] = point_it->marker().value();
        buffer[currentBufferIndex++] = point_it->id();
    }
}
template<typename MeshType>
void updateMarkedSubEntitiesBuffer( MeshType const& mesh,std::vector<unsigned int> & buffer, mpl::int_<3> /**/ )
{
    size_type currentBufferIndex = 0;
    auto face_it = mesh.beginFace();
    auto face_en = mesh.endFace();
    for( ; face_it != face_en; ++face_it )
    {
        if ( face_it->marker().isOff() ) continue;
        buffer[currentBufferIndex++] = face_it->marker().value();
        for ( uint16_type vLocId = 0; vLocId < MeshType::face_type::numPoints ; ++vLocId )
        {
            buffer[currentBufferIndex++] = face_it->point( vLocId ).id();
        }
    }
    auto edge_it = mesh.beginEdge();
    auto edge_en = mesh.endEdge();
    for( ; edge_it != edge_en; ++edge_it )
    {
        if ( edge_it->marker().isOff() ) continue;
        buffer[currentBufferIndex++] = edge_it->marker().value();
        for ( uint16_type vLocId = 0; vLocId < MeshType::edge_type::numPoints ; ++vLocId )
        {
            buffer[currentBufferIndex++] = edge_it->point( vLocId ).id();
        }
    }
    auto point_it = mesh.beginPoint();
    auto point_en = mesh.endPoint();
    for( ; point_it != point_en; ++point_it )
    {
        if ( point_it->marker().isOff() ) continue;
        buffer[currentBufferIndex++] = point_it->marker().value();
        buffer[currentBufferIndex++] = point_it->id();
    }
}

template<typename MeshType>
void updateMarkedSubEntitiesBuffer( MeshType const& mesh,std::vector<unsigned int> & buffer )
{
    updateMarkedSubEntitiesBuffer( mesh, buffer, mpl::int_<MeshType::nDim>() );
}

} // namespace detail

template<typename MeshType>
void PartitionIO<MeshType>::writeMarkedSubEntities()
{
    hsize_t globalDims[2];
    hsize_t localDims[2];
    hsize_t offset[2];

    int d = M_meshPartsOut->nDim;
    int nValFace = 1 + face_type::numPoints;// (marker,ptId1,ptId2,...)
    int nValEdge = 1 + edge_type::numPoints;// (marker,ptId1,ptId2,...)
    int nValPoint = 1 + 1;// (marker,ptId1)
    if (! M_transposeInFile)
    {
        globalDims[0] = 1;
        globalDims[1] = nValFace*M_numGlobalMarkedFaces + nValEdge*M_numGlobalMarkedEdges + nValPoint*M_numGlobalMarkedPoints ;
        localDims[0] = 1;
        localDims[1] = nValFace*M_numLocalMarkedFaces + nValEdge*M_numLocalMarkedEdges + nValPoint*M_numLocalMarkedPoints;
        offset[0] = 0;
        offset[1] = nValFace*M_markedFacesOffSet + nValEdge*M_markedEdgesOffSet + nValPoint*M_markedPointsOffSet;
    }
    else
    {
        CHECK( false ) << "TODO";
    }

    M_HDF5IO.createTable ("marked_subentities", H5T_STD_U32BE, globalDims);

    M_uintBuffer.resize( localDims[0]*localDims[1], 0 );

    if ( !M_transposeInFile )
    {
        //for (size_type i = 0; i < M_numParts; ++i)
        {
            auto& currentPart = *M_meshPartsOut;
            detail::updateMarkedSubEntitiesBuffer( currentPart, M_uintBuffer );
        }
    }
    M_HDF5IO.write ("marked_subentities", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("marked_subentities");

}


template<typename MeshType>
void PartitionIO<MeshType>::readStats()
{
    // Write mesh partition stats (N = number of parts)
    // This is an N x 19 table of int
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

    M_numGlobalPoints = M_uintBuffer[1];
    M_numLocalPoints = M_uintBuffer[2];
    M_pointsOffSet = M_uintBuffer[3];
    LOG(INFO) << "read stat points " << M_numGlobalPoints << ","<<M_numLocalPoints<<","<<M_pointsOffSet<<"\n";

    M_numGlobalElements = M_uintBuffer[4];
    M_numLocalElements = M_uintBuffer[5];
    M_elementsOffSet = M_uintBuffer[6];
    LOG(INFO) << "read stat elements " << M_numGlobalElements << ","<<M_numLocalElements<<","<<M_elementsOffSet<<"\n";

    M_numGlobalGhostElements = M_uintBuffer[7];
    M_numLocalGhostElements = M_uintBuffer[8];
    M_ghostElementsOffSet = M_uintBuffer[9];
    LOG(INFO) << "read stat ghost_elements " << M_numGlobalGhostElements << ","<<M_numLocalGhostElements<<","<<M_ghostElementsOffSet<<"\n";

    M_numGlobalMarkedFaces = M_uintBuffer[10];
    M_numLocalMarkedFaces = M_uintBuffer[11];
    M_markedFacesOffSet = M_uintBuffer[12];
    LOG(INFO) << "read stat marked_faces " << M_numGlobalMarkedFaces << ","<<M_numLocalMarkedFaces<<","<<M_markedFacesOffSet<<"\n";

    M_numGlobalMarkedEdges = M_uintBuffer[13];
    M_numLocalMarkedEdges = M_uintBuffer[14];
    M_markedEdgesOffSet = M_uintBuffer[15];
    LOG(INFO) << "read stat marked_edges " << M_numGlobalMarkedEdges << ","<<M_numLocalMarkedEdges<<","<<M_markedEdgesOffSet<<"\n";

    M_numGlobalMarkedPoints = M_uintBuffer[16];
    M_numLocalMarkedPoints = M_uintBuffer[17];
    M_markedPointsOffSet = M_uintBuffer[18];
    LOG(INFO) << "read stat marked_points " << M_numGlobalMarkedPoints << ","<<M_numLocalMarkedPoints<<","<<M_markedPointsOffSet<<"\n";
}

template<typename MeshType>
void PartitionIO<MeshType>::readPoints()
{
    hsize_t globalDimsIds[2], globalDimsCoords[2];
    hsize_t localDimsIds[2], localDimsCoords[2];
    hsize_t offsetIds[2], offsetCoords[2];


    M_HDF5IO.openTable ("point_ids", globalDimsIds);
    LOG(INFO) << "loaded points_ids:" << globalDimsIds[0] << "x" << globalDimsIds[1];
    M_HDF5IO.openTable ("point_coords", globalDimsCoords);
    LOG(INFO) << "loaded points_coords:" << globalDimsCoords[0] << "x" << globalDimsCoords[1];
    int d = mesh_type::nRealDim;//M_meshPartIn->nRealDim;
    int nValIds = 1;
    if (! M_transposeInFile)
    {
        CHECK( globalDimsIds[0] == 1 && globalDimsIds[1] == nValIds*M_numGlobalPoints ) << "BUG";
        CHECK( globalDimsCoords[0] == 1 && globalDimsCoords[1] == d*M_numGlobalPoints ) << "BUG";

        localDimsIds[0] = 1;
        localDimsIds[1] = nValIds*M_numLocalPoints;
        localDimsCoords[0] = 1;
        localDimsCoords[1] = d*M_numLocalPoints;
        offsetIds[0] = 0;
        offsetIds[1] = nValIds*M_pointsOffSet;
        offsetCoords[0] = 0;
        offsetCoords[1] = d*M_pointsOffSet;
    }
    else
    {
        CHECK( false ) << "TODO";
    }

    M_uintBuffer.resize( localDimsIds[0]*localDimsIds[1], 0 );
    M_realBuffer.resize( localDimsCoords[0]*localDimsCoords[1], 0 );
    M_HDF5IO.read("point_ids", H5T_NATIVE_UINT, localDimsIds, offsetIds, &M_uintBuffer[0]);
    M_HDF5IO.read("point_coords", H5T_NATIVE_DOUBLE, localDimsCoords, offsetCoords, &M_realBuffer[0]);

    M_HDF5IO.closeTable("point_ids");
    M_HDF5IO.closeTable("point_coords");

    node_type coords( d );
    if (! M_transposeInFile)
    {
        rank_type partId = M_meshPartIn->worldComm().localRank();
        size_type currentBufferIndexIds = 0, currentBufferIndexCoords = 0;
        for (size_type j = 0; j < M_numLocalPoints; ++j)
        {
            int id = M_uintBuffer[currentBufferIndexIds++];
            for( int c = 0; c < mesh_type::nRealDim; ++c )
            {
                coords[c] = M_realBuffer[currentBufferIndexCoords++];
            }

            point_type pt( id, coords, false/*onbdy*/ );
            pt.setProcessIdInPartition( partId );
            M_meshPartIn->addPoint( pt );
        }
    }
    else
    {
    }
    M_realBuffer.resize(0);

}

template<typename MeshType>
void PartitionIO<MeshType>::readElements()
{
    int nVal = 2 + element_type::numPoints;// id+marker+ points

    hsize_t globalDims[2];
    hsize_t localDims[2];
    hsize_t offset[2];

    M_HDF5IO.openTable ("elements", globalDims);
    LOG(INFO) << "loaded elements:" << globalDims[0] << "x" << globalDims[1];
    if (! M_transposeInFile)
    {
        CHECK( globalDims[0] == 1 && globalDims[1] == nVal*M_numGlobalElements ) << "BUG";
        localDims[0] = 1;
        localDims[1] = nVal*M_numLocalElements;
        offset[0] = 0;
        offset[1] = nVal*M_elementsOffSet;
    }
    else
    {
    }

    M_uintBuffer.resize( localDims[0]*localDims[1], 0 );
    M_HDF5IO.read("elements", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("elements");

    element_type e;
    rank_type local_rank = M_meshPartIn->worldComm().localRank();
    if (! M_transposeInFile)
    {
        size_type currentBufferIndex = 0;
        for (size_type j = 0; j < M_numLocalElements; ++j)
        {
            size_type id = M_uintBuffer[currentBufferIndex++];
            int marker   = M_uintBuffer[currentBufferIndex++];
            e.setId( id );
            e.setProcessIdInPartition( local_rank );
            e.setMarker( marker );
            e.setProcessId( local_rank );// update correctlty for ghost in preparUpdateForUse()

            for (size_type k = 0; k < element_type::numPoints; ++k)
            {
                e.setPoint( k, M_meshPartIn->point( M_uintBuffer[ currentBufferIndex++]) );
            }
            M_meshPartIn->addElement( e, false );
        }
    }
    else
    {
    }
}

template<typename MeshType>
void PartitionIO<MeshType>::readGhostElements()
{
    int nVal = 3; //(idGhost,pid,idActive)

    hsize_t globalDims[2];
    hsize_t localDims[2];
    hsize_t offset[2];

    M_HDF5IO.openTable ("elements_ghosts", globalDims);
    LOG(INFO) << "loaded ghost_elements:" << globalDims[0] << "x" << globalDims[1];
    if ( !M_transposeInFile )
    {
        CHECK( globalDims[0] == 1 && globalDims[1] == nVal*M_numGlobalGhostElements ) << "BUG";
        localDims[0] = 1;
        localDims[1] = nVal*M_numLocalGhostElements;
        offset[0] = 0;
        offset[1] = nVal*M_ghostElementsOffSet;
    }
    else
    {
    }

    M_uintBuffer.resize( localDims[0]*localDims[1], 0 );
    M_HDF5IO.read("elements_ghosts", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);

    M_HDF5IO.closeTable("elements_ghosts");

    if ( !M_transposeInFile )
    {
        size_type currentBufferIndex = 0;
        for (size_type j = 0; j < M_numLocalGhostElements; ++j)
        {
            size_type id = M_uintBuffer[currentBufferIndex++];
            int pid = M_uintBuffer[currentBufferIndex++];
            size_type idInActivePart = M_uintBuffer[currentBufferIndex++];

            auto it = M_meshPartIn->elementIterator( id, M_meshPartIn->worldComm().localRank() );
            M_meshPartIn->elements().modify( it, Feel::detail::UpdateProcessId(pid) );
            M_meshPartIn->elements().modify( it, Feel::detail::updateIdInOthersPartitions( pid, idInActivePart ) );
            M_meshPartIn->elements().modify( it, Feel::detail::UpdateNeighborPartition( M_meshPartIn->worldComm().localRank() ) );
        }
    }
}


namespace detail
{
template<typename MeshType>
void updateMarkedSubEntitiesMesh( std::vector<unsigned int> const& buffer, std::tuple<size_type,size_type,size_type> const& numLocalSubEntities, MeshType & mesh, mpl::int_<1> /**/ )
{
    size_type numLocalMarkedPoints = std::get<2>( numLocalSubEntities );
    size_type currentBufferIndex = 0;
    for (size_type j = 0; j < numLocalMarkedPoints; ++j)
    {
        size_type marker = buffer[currentBufferIndex++];
        size_type id = buffer[currentBufferIndex++];
        mesh.points().modify( mesh.pointIterator( id ), Feel::detail::UpdateMarker( marker ) );
    }
}
template<typename MeshType>
void updateMarkedSubEntitiesMesh( std::vector<unsigned int> const& buffer, std::tuple<size_type,size_type,size_type> const& numLocalSubEntities, MeshType & mesh, mpl::int_<2> /**/ )
{
    rank_type rank = mesh.worldComm().localRank();
    size_type numLocalMarkedFaces = std::get<0>( numLocalSubEntities );
    size_type numLocalMarkedPoints = std::get<2>( numLocalSubEntities );
    size_type currentBufferIndex = 0;

    typename MeshType::face_type newFace;
    for (size_type j = 0; j < numLocalMarkedFaces; ++j)
    {
        size_type marker = buffer[currentBufferIndex++];
        newFace.setId( mesh.numFaces() );
        newFace.setProcessIdInPartition( rank );
        newFace.setMarker( marker );
        newFace.setOnBoundary( true );
        for ( uint16_type vLocId = 0; vLocId < MeshType::face_type::numPoints ; ++vLocId )
            newFace.setPoint( vLocId, mesh.point( buffer[ currentBufferIndex++ ]) );
        mesh.addFace( newFace );
    }
    for (size_type j = 0; j < numLocalMarkedPoints; ++j)
    {
        size_type marker = buffer[currentBufferIndex++];
        size_type id = buffer[currentBufferIndex++];
        mesh.points().modify( mesh.pointIterator( id ), Feel::detail::UpdateMarker( marker ) );
    }
}
template<typename MeshType>
void updateMarkedSubEntitiesMesh( std::vector<unsigned int> const& buffer, std::tuple<size_type,size_type,size_type> const& numLocalSubEntities, MeshType & mesh, mpl::int_<3> /**/ )
{
    rank_type rank = mesh.worldComm().localRank();
    size_type numLocalMarkedFaces = std::get<0>( numLocalSubEntities );
    size_type numLocalMarkedEdges = std::get<1>( numLocalSubEntities );
    size_type numLocalMarkedPoints = std::get<2>( numLocalSubEntities );
    size_type currentBufferIndex = 0;

    typename MeshType::face_type newFace;
    for (size_type j = 0; j < numLocalMarkedFaces; ++j)
    {
        size_type marker = buffer[currentBufferIndex++];
        newFace.setId( mesh.numFaces() );
        newFace.setProcessIdInPartition( rank );
        newFace.setMarker( marker );
        newFace.setOnBoundary( true );
        for ( uint16_type vLocId = 0; vLocId < MeshType::face_type::numPoints ; ++vLocId )
            newFace.setPoint( vLocId, mesh.point( buffer[ currentBufferIndex++ ]) );
        mesh.addFace( newFace );
    }
    typename MeshType::edge_type newEdge;
    for (size_type j = 0; j < numLocalMarkedEdges; ++j)
    {
        size_type marker = buffer[currentBufferIndex++];
        newEdge.setId( mesh.numEdges() );
        newEdge.setProcessIdInPartition( rank );
        newEdge.setMarker( marker );
        newEdge.setOnBoundary( true );
        for ( uint16_type vLocId = 0; vLocId < MeshType::edge_type::numPoints ; ++vLocId )
            newEdge.setPoint( vLocId, mesh.point( buffer[ currentBufferIndex++ ]) );
        mesh.addEdge( newEdge );
    }
    for (size_type j = 0; j < numLocalMarkedPoints; ++j)
    {
        size_type marker = buffer[currentBufferIndex++];
        size_type id = buffer[currentBufferIndex++];
        mesh.points().modify( mesh.pointIterator( id ), Feel::detail::UpdateMarker( marker ) );
    }
}
template<typename MeshType>
void updateMarkedSubEntitiesMesh( std::vector<unsigned int> const& buffer, std::tuple<size_type,size_type,size_type> const& numLocalSubEntities, MeshType & mesh )
{
    updateMarkedSubEntitiesMesh( buffer, numLocalSubEntities, mesh, mpl::int_<MeshType::nDim>() );
}

} // namespace detail

template<typename MeshType>
void PartitionIO<MeshType>::readMarkedSubEntities()
{
    int nValFace = 1 + face_type::numPoints;//(marker,ptId1,ptId2,...)
    int nValEdge = 1 + edge_type::numPoints;// (marker,ptId1,ptId2,...)
    int nValPoint = 1 + 1;// (marker,ptId)

    hsize_t globalDims[2];
    hsize_t localDims[2];
    hsize_t offset[2];

    M_HDF5IO.openTable ("marked_subentities", globalDims);
    LOG(INFO) << "loaded marked_subentities:" << globalDims[0] << "x" << globalDims[1];
    if ( !M_transposeInFile )
    {
        CHECK( globalDims[0] == 1 && globalDims[1] == (nValFace*M_numGlobalMarkedFaces + nValEdge*M_numGlobalMarkedEdges + nValPoint*M_numGlobalMarkedPoints ) ) << "BUG";
        localDims[0] = 1;
        localDims[1] = nValFace*M_numLocalMarkedFaces + nValEdge*M_numLocalMarkedEdges + nValPoint*M_numLocalMarkedPoints;
        offset[0] = 0;
        offset[1] = nValFace*M_markedFacesOffSet + nValEdge*M_markedEdgesOffSet + nValPoint*M_markedPointsOffSet;
    }
    else
    {
    }

    M_uintBuffer.resize( localDims[0]*localDims[1], 0 );
    M_HDF5IO.read("marked_subentities", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);

    M_HDF5IO.closeTable("marked_subentities");

    detail::updateMarkedSubEntitiesMesh( M_uintBuffer, std::make_tuple( M_numLocalMarkedFaces,M_numLocalMarkedEdges,M_numLocalMarkedPoints ), *M_meshPartIn );
}


template<typename MeshType>
void PartitionIO<MeshType>::prepareUpdateForUse()
{
    rank_type rank = M_meshPartIn->worldComm().localRank();

    // put temporarly the points processId to numPart (which are not a valid part id) in ghost partition
    // and store neigboor part at each point id
    rank_type ghostPointPidDetection = M_numParts;
    std::map<size_type,std::set<rank_type> > pointsToNeihborPart;
    auto ghostelt_it = M_meshPartIn->beginGhostElement();
    auto const ghostelt_en = M_meshPartIn->endGhostElement();
    for ( ; ghostelt_it != ghostelt_en ; ++ghostelt_it )
    {
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::element_type::numPoints; ++vLocId )
        {
            size_type vId = ghostelt_it->point( vLocId ).id();
            pointsToNeihborPart[vId].insert( ghostelt_it->processId() );
            M_meshPartIn->points().modify( M_meshPartIn->pointIterator( vId ), Feel::detail::UpdateProcessId( ghostPointPidDetection ) );
        }
    }

    // update neighbor partion in active elements which touch interprocess boundary
    auto elt_it = M_meshPartIn->beginElementWithProcessId( rank );
    auto const elt_en = M_meshPartIn->endElementWithProcessId( rank );
    for ( ; elt_it != elt_en ; ++elt_it )
    {
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::element_type::numPoints; ++vLocId )
        {
            auto const& thepoint = elt_it->point( vLocId );
            size_type vId = thepoint.id();
            if ( thepoint.processId() == ghostPointPidDetection )
            {
                auto const& neighIds = pointsToNeihborPart.find( vId )->second;
                for ( auto const& neighId : neighIds )
                    M_meshPartIn->elements().modify( elt_it, Feel::detail::UpdateNeighborPartition( neighId ) );
            }
        }
    }

    // update points processId in active elements
    elt_it = M_meshPartIn->beginElementWithProcessId( rank );
    for ( ; elt_it != elt_en ; ++elt_it )
    {
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::element_type::numPoints; ++vLocId )
        {
            auto const& thepoint = elt_it->point( vLocId );
            size_type vId = thepoint.id();
            M_meshPartIn->points().modify( M_meshPartIn->pointIterator( vId ), Feel::detail::UpdateProcessId( rank ) );
        }
    }

    // reset points processId in ghost elements (which do not belong the active part)
    ghostelt_it = M_meshPartIn->beginGhostElement();
    for ( ; ghostelt_it != ghostelt_en ; ++ghostelt_it )
    {
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::element_type::numPoints; ++vLocId )
        {
            auto const& thepoint = ghostelt_it->point( vLocId );
            if ( thepoint.processId() == ghostPointPidDetection )
            {
                size_type vId = thepoint.id();
                M_meshPartIn->points().modify( M_meshPartIn->pointIterator( vId ), Feel::detail::UpdateProcessId( invalid_rank_type_value ) );
            }
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
