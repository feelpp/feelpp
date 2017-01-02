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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <feel/feelcore/hdf5.hpp>
#include <feel/feelmesh/meshpartitionset.hpp>

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
  4. elements - 1 x ( (1 + num of element nodes)* number of global active elements ):
         - element markers
         - points ids
  5. elements_ghosts - 1 x ( number of global ghost elements ):
         - process id
  5. elements_ghosts - 1 x ( ((1+num of face node)*number of global ghost elements ):
  6. marked_subentities - 1 x ( (1+num of face node)* number of global marked faces
     + (1+num of face node)* number of global marked edges + (1+1)*number of global marked points ):
         - marker
         - points ids

 * @author Christophe Prud'homme
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
    typedef MeshPartitionSet<mesh_type> mesh_partitionset_type;
    typedef std::unique_ptr<mesh_partitionset_type> mesh_partitionset_ptrtype;

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
    PartitionIO() = default;
    PartitionIO (const PartitionIO&) = delete;
    PartitionIO& operator= (const PartitionIO&) = delete;
    
    
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
    void write ( mesh_partitionset_ptrtype&& mesh);
    //! Read method
    /*!
     * Call this method to read from the HDF5 file the mesh part associated
     * with the rank of the current MPI process
     * \param mesh pointer to a mesh . If the RegionMesh has been initialized and contains any
     *        state, it will be destroyed before reading
     */
    void read (mesh_ptrtype mesh, size_type ctxMeshUpdate = MESH_UPDATE_EDGES|MESH_UPDATE_FACES );


    //@}



protected:

private:
    

    //! Private Methods
    //@{

    /**
     * write meta data
     */
    FEELPP_NO_EXPORT void writeMetaData( mesh_ptrtype mesh );

    /**
     * read meta data
     */
    FEELPP_NO_EXPORT void readMetaData( mesh_ptrtype mesh );

    // Methods for writing
    FEELPP_NO_EXPORT void writeStats();
    FEELPP_NO_EXPORT void writePoints();
    FEELPP_NO_EXPORT void writeElements();
    FEELPP_NO_EXPORT void writeGhostElements();
    FEELPP_NO_EXPORT void writeMarkedSubEntities();
    // Methods for reading
    FEELPP_NO_EXPORT void readStats();
    FEELPP_NO_EXPORT void readPoints();
    FEELPP_NO_EXPORT void readElements( std::map<rank_type,std::vector<size_type>> & mapGhostHdf5IdToFeelId );
    FEELPP_NO_EXPORT void readGhostElements( std::map<rank_type,std::vector<size_type>> const& mapGhostHdf5IdToFeelId );
    FEELPP_NO_EXPORT void readMarkedSubEntities();
    FEELPP_NO_EXPORT void prepareUpdateForUseStep1();
    FEELPP_NO_EXPORT void prepareUpdateForUseStep2();

    //@}

    //! Private Data Members
    //@{
    std::string M_filename;
    std::string M_h5_filename;
    bool M_transposeInFile;
    meshparts_ptrtype M_meshPartsOut;
    mesh_ptrtype M_meshPartIn;
    mesh_partitionset_ptrtype M_meshPartitionSet;

    // Mesh geometry
    size_type M_numGlobalPoints,M_numGlobalElements,M_numGlobalGhostElements;
    size_type M_numGlobalMarkedFaces,M_numGlobalMarkedEdges,M_numGlobalMarkedPoints;
    std::map<rank_type,size_type> M_numLocalPoints, M_pointsOffSet;
    std::map<rank_type,size_type> M_numLocalElements, M_elementsOffSet;
    std::map<rank_type,size_type> M_numLocalGhostElements, M_ghostElementsOffSet;
    std::map<rank_type,size_type> M_numLocalMarkedFaces, M_markedFacesOffSet;
    std::map<rank_type,size_type> M_numLocalMarkedEdges, M_markedEdgesOffSet;
    std::map<rank_type,size_type> M_numLocalMarkedPoints, M_markedPointsOffSet;

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
    fs::path filenamefs = fs::path( fileName );
    std::string h5filename = filenamefs.stem().string() + ".h5";
    if ( filenamefs.is_relative() )
        M_h5_filename = (fs::current_path()/fs::path(h5filename)).string();
    else
        M_h5_filename = (filenamefs.parent_path() / fs::path( h5filename )).string();
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
    //M_meshPartitionSet = std::make_shared<mesh_partitionset_type>( meshParts );
    this->write( std::make_unique<mesh_partitionset_type>( meshParts ) );
}
template<typename MeshType>
void PartitionIO<MeshType>::write ( mesh_partitionset_ptrtype&& meshpartset )
{
    M_meshPartitionSet = std::move(meshpartset);
    M_meshPartsOut = M_meshPartitionSet->mesh();

    if ( M_meshPartsOut->worldComm().isMasterRank() )
        writeMetaData( M_meshPartsOut );

    LOG(INFO) << "writing mesh in HDF5 format in " << M_h5_filename;
    tic();
    M_HDF5IO.openFile (M_h5_filename, M_meshPartsOut->worldComm(), false);

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
    pt.put("mesh.partition.n",meshParts->numberOfPartitions());
    for( auto m : meshParts->markerNames() )
    {
        pt.put("mesh.physicals."+m.first, m.second );
    }

    pt::write_json(M_filename, pt);
}
template<typename MeshType>
void PartitionIO<MeshType>::read (mesh_ptrtype meshParts, size_type ctxMeshUpdate)
{
    readMetaData( meshParts );
    M_meshPartIn = meshParts;

    tic();
    M_HDF5IO.openFile (M_h5_filename, meshParts->worldComm(), true);
    tic();
    readStats();
    toc("PartitionIO reading stats",FLAGS_v>0);
    tic();
    readPoints();
    toc("PartitionIO reading points",FLAGS_v>0);
    tic();
    std::map<rank_type,std::vector<size_type>> mapGhostHdf5IdToFeelId;
    readElements( mapGhostHdf5IdToFeelId );
    toc("PartitionIO reading elements",FLAGS_v>0);
    tic();
    readGhostElements( mapGhostHdf5IdToFeelId );
    toc("PartitionIO reading ghost_elements",FLAGS_v>0);
    tic();
    readMarkedSubEntities();
    toc("PartitionIO reading marked_subentities",FLAGS_v>0);

    M_HDF5IO.closeFile();
    toc("PartitionIO reading hdf5 file",FLAGS_v>0);

    prepareUpdateForUseStep1();
    prepareUpdateForUseStep2();

    tic();
    M_meshPartIn->components().reset();
    M_meshPartIn->components().set( ctxMeshUpdate );
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
    meshParts->setNumberOfPartitions(pt.get("mesh.partition.n",meshParts->worldComm().globalSize()));
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
    rank_type numGlobalPartition = M_meshPartitionSet->numGlobalPartition();
    // Write mesh partition stats (N = number of parts)
    // This is an N x 19 table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = numGlobalPartition;
        currentSpaceDims[1] = 19;
        currentCount[0] = 1;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = 19;
        currentSpaceDims[1] = numGlobalPartition;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 1;
    }

    M_numGlobalPoints = 0;
    M_numGlobalElements = 0;
    M_numGlobalGhostElements = 0;
    M_numGlobalMarkedFaces = 0;
    M_numGlobalMarkedEdges = 0;
    M_numGlobalMarkedPoints = 0;

    for ( rank_type partId : M_meshPartitionSet->localPartitionIds() )
    {
        M_numLocalPoints[partId] = M_meshPartitionSet->statNumPointsAll( partId );
        M_pointsOffSet[partId] = 0;
        M_numLocalElements[partId] = M_meshPartitionSet->statNumElementsAll( partId );
        M_elementsOffSet[partId] = 0;
        M_numLocalGhostElements[partId] = M_meshPartitionSet->statNumElementsAll( partId ) - M_meshPartitionSet->statNumElementsActive( partId );
        M_ghostElementsOffSet[partId] = 0;
        M_numLocalMarkedFaces[partId] = M_meshPartitionSet->statNumFacesMarkedAll( partId );
        M_markedFacesOffSet[partId] = 0;
        M_numLocalMarkedEdges[partId] = M_meshPartitionSet->statNumEdgesMarkedAll( partId );
        M_markedEdgesOffSet[partId] = 0;
        M_numLocalMarkedPoints[partId] = M_meshPartitionSet->statNumPointsMarkedAll( partId );
        M_markedPointsOffSet[partId] = 0;
    }

    for ( rank_type p=0 ; p< numGlobalPartition ; ++p )
    {
        size_type nPtAll = M_meshPartitionSet->statNumPointsAll(p);
        M_numGlobalPoints += nPtAll;
        size_type nEltAll = M_meshPartitionSet->statNumElementsAll(p);
        M_numGlobalElements += nEltAll;
        size_type nGhostEltAll = M_meshPartitionSet->statNumElementsAll(p) - M_meshPartitionSet->statNumElementsActive(p);
        M_numGlobalGhostElements += nGhostEltAll;
        size_type nMarkedFacesAll = M_meshPartitionSet->statNumFacesMarkedAll(p);
        M_numGlobalMarkedFaces += nMarkedFacesAll;
        size_type nMarkedEdgesAll = M_meshPartitionSet->statNumEdgesMarkedAll(p);
        M_numGlobalMarkedEdges += nMarkedEdgesAll;
        size_type nMarkedPointsAll = M_meshPartitionSet->statNumPointsMarkedAll(p);
        M_numGlobalMarkedPoints += nMarkedPointsAll;

        for ( rank_type partId : M_meshPartitionSet->localPartitionIds() )
        {
            if ( p < partId )
            {
                M_pointsOffSet[partId] += nPtAll;
                M_elementsOffSet[partId] += nEltAll;
                M_ghostElementsOffSet[partId] += nGhostEltAll;
                M_markedFacesOffSet[partId] += nMarkedFacesAll;
                M_markedEdgesOffSet[partId] += nMarkedEdgesAll;
                M_markedPointsOffSet[partId] += nMarkedPointsAll;
            }
        }
    }

    // Create new table
    M_HDF5IO.createTable ("stats", H5T_STD_U32BE, currentSpaceDims);

    for ( rank_type partId : M_meshPartitionSet->localPartitionIds() )
    {
        // Fill buffer
        M_uintBuffer.resize( 19 );
        M_uintBuffer[0] = partId;
        M_uintBuffer[1] = M_numGlobalPoints;
        M_uintBuffer[2] = M_numLocalPoints[partId];
        M_uintBuffer[3] = M_pointsOffSet[partId];
        M_uintBuffer[4] = M_numGlobalElements;
        M_uintBuffer[5] = M_numLocalElements[partId];
        M_uintBuffer[6] = M_elementsOffSet[partId];
        M_uintBuffer[7] = M_numGlobalGhostElements;
        M_uintBuffer[8] = M_numLocalGhostElements[partId];
        M_uintBuffer[9] = M_ghostElementsOffSet[partId];
        M_uintBuffer[10] = M_numGlobalMarkedFaces;
        M_uintBuffer[11] = M_numLocalMarkedFaces[partId];
        M_uintBuffer[12] = M_markedFacesOffSet[partId];
        M_uintBuffer[13] = M_numGlobalMarkedEdges;
        M_uintBuffer[14] = M_numLocalMarkedEdges[partId];
        M_uintBuffer[15] = M_markedEdgesOffSet[partId];
        M_uintBuffer[16] = M_numGlobalMarkedPoints;
        M_uintBuffer[17] = M_numLocalMarkedPoints[partId];
        M_uintBuffer[18] = M_markedPointsOffSet[partId];

        std::ostringstream ostr;ostr << "write stat :";
        for (int k=0;k<M_uintBuffer.size();++k)
            ostr << " " << M_uintBuffer[k];
        VLOG(1) << ostr.str();

        hsize_t currentOffset[2];
        if ( !M_transposeInFile )
        {
            currentOffset[0] = partId;
            currentOffset[1] = 0;
        }
        else
        {
            currentOffset[0] = 0;
            currentOffset[1] = partId;
        }
        M_HDF5IO.write("stats", H5T_NATIVE_UINT, currentCount, currentOffset, &M_uintBuffer[0]);
    }
    M_HDF5IO.closeTable("stats");
}

template<typename MeshType>
void PartitionIO<MeshType>::writePoints()
{
    int d = M_meshPartsOut->nRealDim;
    int nValIds = 1;
    hsize_t globalDimsIds[2], globalDimsCoords[2];
    globalDimsIds[0] = 1;
    globalDimsIds[1] = nValIds*M_numGlobalPoints;
    globalDimsCoords[0] = 1;
    globalDimsCoords[1] = d*M_numGlobalPoints;
    if ( M_transposeInFile )
    {
        globalDimsIds[0] = globalDimsIds[1];
        globalDimsIds[1] = 1;
        globalDimsCoords[0] = globalDimsCoords[1];
        globalDimsCoords[1] = 1;
    }

    M_HDF5IO.createTable ("point_ids", H5T_STD_U32BE, globalDimsIds);
    M_HDF5IO.createTable ("point_coords", H5T_IEEE_F64BE, globalDimsCoords);

    for ( rank_type partId : M_meshPartitionSet->localPartitionIds() )
    {
        hsize_t localDimsIds[2], localDimsCoords[2];
        hsize_t offsetIds[2], offsetCoords[2];

        localDimsIds[0] = 1;
        localDimsIds[1] = nValIds*M_numLocalPoints[partId];
        localDimsCoords[0] = 1;
        localDimsCoords[1] = d*M_numLocalPoints[partId];
        offsetIds[0] = 0;
        offsetIds[1] = nValIds*M_pointsOffSet[partId];
        offsetCoords[0] = 0;
        offsetCoords[1] = d*M_pointsOffSet[partId];
        if ( M_transposeInFile )
        {
            localDimsIds[0] = localDimsIds[1];
            localDimsIds[1] = 1;
            localDimsCoords[0] = localDimsCoords[1];
            localDimsCoords[1] = 1;
            offsetIds[0] = offsetIds[1];
            offsetIds[1] = 0;
            offsetCoords[0] = offsetCoords[1];
            offsetCoords[1] = 0;
        }
        // Fill buffer
        M_uintBuffer.resize( localDimsIds[0]*localDimsIds[1], 0 );
        M_realBuffer.resize( localDimsCoords[0]*localDimsCoords[1], 0 );

        if ( !M_uintBuffer.empty() && !M_realBuffer.empty() )
        {
            size_type currentBufferIndexIds = 0,currentBufferIndexCoords = 0;
            auto pt_it = M_meshPartitionSet->beginPoint( partId );
            auto pt_en = M_meshPartitionSet->endPoint( partId );
            for(  ; pt_it != pt_en; ++pt_it )
            {
                auto const& thepoint = boost::unwrap_ref(*pt_it);
                M_uintBuffer[currentBufferIndexIds++] = thepoint.id();

                M_realBuffer[currentBufferIndexCoords++] = thepoint.operator[](0);
                if ( mesh_type::nRealDim >= 2 )
                    M_realBuffer[currentBufferIndexCoords++] = thepoint.operator[](1);
                if ( mesh_type::nRealDim >= 3 )
                    M_realBuffer[currentBufferIndexCoords++] = thepoint.operator[](2);
            }

            M_HDF5IO.write ("point_ids", H5T_NATIVE_UINT, localDimsIds, offsetIds, &M_uintBuffer[0]);
            M_HDF5IO.write ("point_coords", H5T_NATIVE_DOUBLE, localDimsCoords, offsetCoords, &M_realBuffer[0]);
        }
    }

    M_realBuffer.resize (0);
    M_uintBuffer.resize (0);

    M_HDF5IO.closeTable ("point_coords");
    M_HDF5IO.closeTable ("point_ids");
}

template<typename MeshType>
void PartitionIO<MeshType>::writeElements()
{
    //int nVal = 2 + element_type::numPoints;// id+marker+ points
    int nVal = 1 + element_type::numPoints;// marker+ points

    hsize_t globalDims[2];
    globalDims[0] = 1;
    globalDims[1] = nVal*M_numGlobalElements;
    if ( M_transposeInFile )
    {
        globalDims[0] = globalDims[1];
        globalDims[1] = 1;
    }

    M_HDF5IO.createTable ("elements", H5T_STD_U32BE, globalDims);

    for ( rank_type partId : M_meshPartitionSet->localPartitionIds() )
    {
        hsize_t localDims[2];
        hsize_t offset[2];
        localDims[0] = 1;
        localDims[1] = nVal*M_numLocalElements[partId];
        offset[0] = 0;
        offset[1] = nVal*M_elementsOffSet[partId];
        if ( M_transposeInFile )
        {
            localDims[0] = localDims[1];
            localDims[1] = 1;
            offset[0] = offset[1];
            offset[1] = 0;
        }
        // Fill buffer
        M_uintBuffer.resize( localDims[0]*localDims[1], 0 );
        if ( !M_uintBuffer.empty() )
        {
            size_type currentBufferIndex = 0;
            auto elt_it = M_meshPartitionSet->beginActiveElement( partId );
            auto elt_en = M_meshPartitionSet->endActiveElement( partId );
            for( ; elt_it != elt_en; ++ elt_it )
            {
                auto const& theelt = boost::unwrap_ref(*elt_it);
                //M_uintBuffer[currentBufferIndex++] = theelt.id();
                M_uintBuffer[currentBufferIndex++] = theelt.marker().value();
                for (size_type k = 0; k < element_type::numPoints; ++k)
                {
                    M_uintBuffer[currentBufferIndex++] = theelt.point(k).id();
                }
            }
            auto ghostelt_it = M_meshPartitionSet->beginGhostElement( partId );
            auto ghostelt_en = M_meshPartitionSet->endGhostElement( partId );
            for( ; ghostelt_it != ghostelt_en; ++ghostelt_it )
            {
                auto const& theghostelt = boost::unwrap_ref(*ghostelt_it);
                //M_uintBuffer[currentBufferIndex++] = theghostelt.id();
                M_uintBuffer[currentBufferIndex++] = theghostelt.marker().value();
                for (size_type k = 0; k < element_type::numPoints; ++k)
                {
                    M_uintBuffer[currentBufferIndex++] = theghostelt.point(k).id();
                }
            }

            M_HDF5IO.write ("elements", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);
        }
    }

    M_HDF5IO.closeTable ("elements");
}

template<typename MeshType>
void PartitionIO<MeshType>::writeGhostElements()
{
    hsize_t globalDims[2];

    //int nVal = 3;//(idGhost,pid,idActive)
    int nVal = 1; // pid
    globalDims[0] = 1;
    globalDims[1] = nVal*M_numGlobalGhostElements;
    if ( M_transposeInFile )
    {
        globalDims[0] = globalDims[1];
        globalDims[1] = 1;
    }

    M_HDF5IO.createTable ("elements_ghosts", H5T_STD_U32BE, globalDims);

    for ( rank_type partId : M_meshPartitionSet->localPartitionIds() )
    {
        hsize_t localDims[2];
        hsize_t offset[2];
        localDims[0] = 1;
        localDims[1] = nVal*M_numLocalGhostElements[partId];
        offset[0] = 0;
        offset[1] = nVal*M_ghostElementsOffSet[partId];
        if ( M_transposeInFile )
        {
            localDims[0] = localDims[1];
            localDims[1] = 1;
            offset[0] = offset[1];
            offset[1] = 0;
        }
        M_uintBuffer.resize( localDims[0]*localDims[1], 0 );

        if ( !M_uintBuffer.empty() )
        {
            size_type currentBufferIndex = 0;
            auto ghostelt_it = M_meshPartitionSet->beginGhostElement( partId );
            auto ghostelt_en = M_meshPartitionSet->endGhostElement( partId );
            for( ; ghostelt_it != ghostelt_en; ++ghostelt_it )
            {
                auto const& theghostelt = boost::unwrap_ref(*ghostelt_it);
                //M_uintBuffer[currentBufferIndex++] = theghostelt.id();
                M_uintBuffer[currentBufferIndex++] = theghostelt.processId();
            }
            M_HDF5IO.write ("elements_ghosts", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);
        }
    }

    M_HDF5IO.closeTable ("elements_ghosts");

}

namespace detail
{
template<typename MeshPartType>
void updateMarkedSubEntitiesBuffer( MeshPartType const& meshPartSet, rank_type partId, std::vector<unsigned int> & buffer, mpl::int_<1> /**/ )
{
    size_type currentBufferIndex = 0;
    auto point_it = meshPartSet.beginMarkedPoint( partId );
    auto point_en = meshPartSet.endMarkedPoint( partId );
    for( ; point_it != point_en; ++point_it )
    {
        auto const& thepoint = boost::unwrap_ref(*point_it);
        //if ( point_it->marker().isOff() ) continue;
        CHECK( thepoint.marker().isOn() ) << "not a marked point";
        buffer[currentBufferIndex++] = thepoint.marker().value();
        buffer[currentBufferIndex++] = thepoint.id();
    }
}
template<typename MeshPartType>
void updateMarkedSubEntitiesBuffer( MeshPartType const& meshPartSet, rank_type partId, std::vector<unsigned int> & buffer, mpl::int_<2> /**/ )
{
    size_type currentBufferIndex = 0;
    auto face_it = meshPartSet.beginMarkedFace( partId );
    auto face_en = meshPartSet.endMarkedFace( partId );
    for( ; face_it != face_en; ++face_it )
    {
        auto const& theface = boost::unwrap_ref(*face_it);
        //if ( face_it->marker().isOff() ) continue;
        CHECK( theface.marker().isOn() ) << "not a marked face";
        buffer[currentBufferIndex++] = theface.marker().value();
        for ( uint16_type vLocId = 0; vLocId < MeshPartType::mesh_type::face_type::numPoints ; ++vLocId )
        {
            buffer[currentBufferIndex++] = theface.point( vLocId ).id();
        }
    }
    auto point_it = meshPartSet.beginMarkedPoint( partId );
    auto point_en = meshPartSet.endMarkedPoint( partId );
    for( ; point_it != point_en; ++point_it )
    {
        auto const& thepoint = boost::unwrap_ref(*point_it);
        //if ( point_it->marker().isOff() ) continue;
        CHECK( thepoint.marker().isOn() ) << "not a marked point";
        buffer[currentBufferIndex++] = thepoint.marker().value();
        buffer[currentBufferIndex++] = thepoint.id();
    }
}
template<typename MeshPartType>
void updateMarkedSubEntitiesBuffer( MeshPartType const& meshPartSet, rank_type partId, std::vector<unsigned int> & buffer, mpl::int_<3> /**/ )
{
    size_type currentBufferIndex = 0;
    auto face_it = meshPartSet.beginMarkedFace( partId );
    auto face_en = meshPartSet.endMarkedFace( partId );
    for( ; face_it != face_en; ++face_it )
    {
        auto const& theface = boost::unwrap_ref(*face_it);
        //if ( face_it->marker().isOff() ) continue;
        CHECK( theface.marker().isOn() ) << "not a marked face";
        buffer[currentBufferIndex++] = theface.marker().value();
        for ( uint16_type vLocId = 0; vLocId < MeshPartType::mesh_type::face_type::numPoints ; ++vLocId )
        {
            buffer[currentBufferIndex++] = theface.point( vLocId ).id();
        }
    }
    auto edge_it = meshPartSet.beginMarkedEdge( partId );
    auto edge_en = meshPartSet.endMarkedEdge( partId );
    for( ; edge_it != edge_en; ++edge_it )
    {
        auto const& theedge = boost::unwrap_ref(*edge_it);
        //if ( edge_it->marker().isOff() ) continue;
        CHECK( theedge.marker().isOn() ) << "not a marked edge";
        buffer[currentBufferIndex++] = theedge.marker().value();
        for ( uint16_type vLocId = 0; vLocId < MeshPartType::mesh_type::edge_type::numPoints ; ++vLocId )
        {
            buffer[currentBufferIndex++] = theedge.point( vLocId ).id();
        }
    }
    auto point_it = meshPartSet.beginMarkedPoint( partId );
    auto point_en = meshPartSet.endMarkedPoint( partId );
    for( ; point_it != point_en; ++point_it )
    {
        auto const& thepoint = boost::unwrap_ref(*point_it);
        //if ( point_it->marker().isOff() ) continue;
        CHECK( thepoint.marker().isOn() ) << "not a marked point";
        buffer[currentBufferIndex++] = thepoint.marker().value();
        buffer[currentBufferIndex++] = thepoint.id();
    }
}

template<typename MeshPartType>
void updateMarkedSubEntitiesBuffer( MeshPartType const& meshPartSet, rank_type partId, std::vector<unsigned int> & buffer )
{
    updateMarkedSubEntitiesBuffer( meshPartSet, partId,buffer, mpl::int_<MeshPartType::mesh_type::nDim>() );
}

} // namespace detail

template<typename MeshType>
void PartitionIO<MeshType>::writeMarkedSubEntities()
{
    int nValFace = 1 + face_type::numPoints;// (marker,ptId1,ptId2,...)
    int nValEdge = 1 + edge_type::numPoints;// (marker,ptId1,ptId2,...)
    int nValPoint = 1 + 1;// (marker,ptId1)

    hsize_t globalDims[2];
    globalDims[0] = 1;
    globalDims[1] = nValFace*M_numGlobalMarkedFaces + nValEdge*M_numGlobalMarkedEdges + nValPoint*M_numGlobalMarkedPoints ;
    if ( M_transposeInFile )
    {
        globalDims[0] = globalDims[1];
        globalDims[1] = 1;
    }

    M_HDF5IO.createTable ("marked_subentities", H5T_STD_U32BE, globalDims);

    for ( rank_type partId : M_meshPartitionSet->localPartitionIds() )
    {
        hsize_t localDims[2];
        hsize_t offset[2];
        localDims[0] = 1;
        localDims[1] = nValFace*M_numLocalMarkedFaces[partId] + nValEdge*M_numLocalMarkedEdges[partId] + nValPoint*M_numLocalMarkedPoints[partId];
        offset[0] = 0;
        offset[1] = nValFace*M_markedFacesOffSet[partId] + nValEdge*M_markedEdgesOffSet[partId] + nValPoint*M_markedPointsOffSet[partId];
        if ( M_transposeInFile )
        {
            localDims[0] = localDims[1];
            localDims[1] = 1;
            offset[0] = offset[1];
            offset[1] = 1;
        }
        M_uintBuffer.resize( localDims[0]*localDims[1], 0 );
        if ( !M_uintBuffer.empty() )
        {
            Feel::detail::updateMarkedSubEntitiesBuffer( *M_meshPartitionSet, partId, M_uintBuffer );

            M_HDF5IO.write ("marked_subentities", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);
        }
    }

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

    rank_type partId = M_meshPartIn->worldComm().localRank();

    M_HDF5IO.openTable ("stats", currentSpaceDims);

    if (! M_transposeInFile)
    {
        currentCount[0] = 1;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = partId;
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 1;
        currentOffset[0] = 0;
        currentOffset[1] = partId;
    }

    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("stats", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("stats");

    M_numGlobalPoints = M_uintBuffer[1];
    M_numLocalPoints[partId] = M_uintBuffer[2];
    M_pointsOffSet[partId] = M_uintBuffer[3];
    LOG(INFO) << "read stat points " << M_numGlobalPoints << ","<<M_numLocalPoints[partId]<<","<<M_pointsOffSet[partId]<<"\n";

    M_numGlobalElements = M_uintBuffer[4];
    M_numLocalElements[partId] = M_uintBuffer[5];
    M_elementsOffSet[partId] = M_uintBuffer[6];
    LOG(INFO) << "read stat elements " << M_numGlobalElements << ","<<M_numLocalElements[partId]<<","<<M_elementsOffSet[partId]<<"\n";

    M_numGlobalGhostElements = M_uintBuffer[7];
    M_numLocalGhostElements[partId] = M_uintBuffer[8];
    M_ghostElementsOffSet[partId] = M_uintBuffer[9];
    LOG(INFO) << "read stat ghost_elements " << M_numGlobalGhostElements << ","<<M_numLocalGhostElements[partId]<<","<<M_ghostElementsOffSet[partId]<<"\n";

    M_numGlobalMarkedFaces = M_uintBuffer[10];
    M_numLocalMarkedFaces[partId] = M_uintBuffer[11];
    M_markedFacesOffSet[partId] = M_uintBuffer[12];
    LOG(INFO) << "read stat marked_faces " << M_numGlobalMarkedFaces << ","<<M_numLocalMarkedFaces[partId]<<","<<M_markedFacesOffSet[partId]<<"\n";

    M_numGlobalMarkedEdges = M_uintBuffer[13];
    M_numLocalMarkedEdges[partId] = M_uintBuffer[14];
    M_markedEdgesOffSet[partId] = M_uintBuffer[15];
    LOG(INFO) << "read stat marked_edges " << M_numGlobalMarkedEdges << ","<<M_numLocalMarkedEdges[partId]<<","<<M_markedEdgesOffSet[partId]<<"\n";

    M_numGlobalMarkedPoints = M_uintBuffer[16];
    M_numLocalMarkedPoints[partId] = M_uintBuffer[17];
    M_markedPointsOffSet[partId] = M_uintBuffer[18];
    LOG(INFO) << "read stat marked_points " << M_numGlobalMarkedPoints << ","<<M_numLocalMarkedPoints[partId]<<","<<M_markedPointsOffSet[partId]<<"\n";
}

template<typename MeshType>
void PartitionIO<MeshType>::readPoints()
{
    hsize_t globalDimsIds[2], globalDimsCoords[2];
    hsize_t localDimsIds[2], localDimsCoords[2];
    hsize_t offsetIds[2], offsetCoords[2];

    rank_type partId = M_meshPartIn->worldComm().localRank();

    M_HDF5IO.openTable ("point_ids", globalDimsIds);
    LOG(INFO) << "loaded points_ids:" << globalDimsIds[0] << "x" << globalDimsIds[1];
    M_HDF5IO.openTable ("point_coords", globalDimsCoords);
    LOG(INFO) << "loaded points_coords:" << globalDimsCoords[0] << "x" << globalDimsCoords[1];
    int d = mesh_type::nRealDim;//M_meshPartIn->nRealDim;
    int nValIds = 1;

    int indexTable0 = 0, indexTable1 = 1;
    if ( M_transposeInFile )
    {
        indexTable0 = 1;
        indexTable1 = 0;
    }

    CHECK( globalDimsIds[indexTable0] == 1 && globalDimsIds[indexTable1] == nValIds*M_numGlobalPoints ) << "BUG";
    CHECK( globalDimsCoords[indexTable0] == 1 && globalDimsCoords[indexTable1] == d*M_numGlobalPoints ) << "BUG";
    if ( globalDimsIds[indexTable1] > 0 && globalDimsCoords[indexTable0] > 0 )
    {
        localDimsIds[indexTable0] = 1;
        localDimsIds[indexTable1] = nValIds*M_numLocalPoints[partId];
        localDimsCoords[indexTable0] = 1;
        localDimsCoords[indexTable1] = d*M_numLocalPoints[partId];
        offsetIds[indexTable0] = 0;
        offsetIds[indexTable1] = nValIds*M_pointsOffSet[partId];
        offsetCoords[indexTable0] = 0;
        offsetCoords[indexTable1] = d*M_pointsOffSet[partId];

        M_uintBuffer.resize( localDimsIds[indexTable0]*localDimsIds[indexTable1], 0 );
        M_realBuffer.resize( localDimsCoords[indexTable0]*localDimsCoords[indexTable1], 0 );

        if ( !M_uintBuffer.empty() && !M_realBuffer.empty() )
        {
            M_HDF5IO.read("point_ids", H5T_NATIVE_UINT, localDimsIds, offsetIds, &M_uintBuffer[0]);
            M_HDF5IO.read("point_coords", H5T_NATIVE_DOUBLE, localDimsCoords, offsetCoords, &M_realBuffer[0]);
        }
        else
        {
            size_type uselessValueIds=0;
            double uselessValueCoords=0;
            M_HDF5IO.read("point_ids", H5T_NATIVE_UINT, localDimsIds, offsetIds, &uselessValueIds);
            M_HDF5IO.read("point_coords", H5T_NATIVE_DOUBLE, localDimsCoords, offsetCoords, &uselessValueCoords);
        }
    }
    M_HDF5IO.closeTable("point_ids");
    M_HDF5IO.closeTable("point_coords");

    node_type coords( d );
    size_type currentBufferIndexIds = 0, currentBufferIndexCoords = 0;
    for (size_type j = 0; j < M_numLocalPoints[partId]; ++j)
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
    M_realBuffer.resize(0);

}

template<typename MeshType>
void PartitionIO<MeshType>::readElements( std::map<rank_type,std::vector<size_type>> & mapGhostHdf5IdToFeelId )
{
    rank_type partId = M_meshPartIn->worldComm().localRank();

    //int nVal = 2 + element_type::numPoints;// id+marker+ points
    int nVal = 1 + element_type::numPoints;// marker+ points

    hsize_t globalDims[2];
    hsize_t localDims[2];
    hsize_t offset[2];

    M_HDF5IO.openTable ("elements", globalDims);
    LOG(INFO) << "loaded elements:" << globalDims[0] << "x" << globalDims[1];

    int indexTable0 = 0, indexTable1 = 1;
    if ( M_transposeInFile )
    {
        indexTable0 = 1;
        indexTable1 = 0;
    }

    CHECK( globalDims[indexTable0] == 1 && globalDims[indexTable1] == nVal*M_numGlobalElements ) << "BUG";
    if ( globalDims[indexTable1] > 0 )
    {
        localDims[indexTable0] = 1;
        localDims[indexTable1] = nVal*M_numLocalElements[partId];
        offset[indexTable0] = 0;
        offset[indexTable1] = nVal*M_elementsOffSet[partId];

        M_uintBuffer.resize( localDims[indexTable0]*localDims[indexTable1], 0 );
        if ( !M_uintBuffer.empty() )
            M_HDF5IO.read("elements", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);
        else
        {
            size_type uselessValue = 0;
            M_HDF5IO.read("elements", H5T_NATIVE_UINT, localDims, offset, &uselessValue);
        }
    }

    M_HDF5IO.closeTable ("elements");


    size_type nActiveElement = M_numLocalElements[partId]-M_numLocalGhostElements[partId];
    mapGhostHdf5IdToFeelId[partId].resize( M_numLocalGhostElements[partId] );
    element_type e;

    size_type currentBufferIndex = 0;
    for (size_type j = 0; j < M_numLocalElements[partId]; ++j)
    {
        //size_type id = M_uintBuffer[currentBufferIndex++];
        int marker   = M_uintBuffer[currentBufferIndex++];
        //e.setId( id );
        e.setProcessIdInPartition( partId );
        e.setMarker( marker );
        e.setProcessId( partId );// update correctlty for ghost in preparUpdateForUse()

        for (size_type k = 0; k < element_type::numPoints; ++k)
        {
            size_type ptId = M_uintBuffer[ currentBufferIndex++];
            CHECK( M_meshPartIn->hasPoint( ptId ) ) << "point id " << ptId << " not present in mesh";
            e.setPoint( k, M_meshPartIn->point( ptId ) );
        }
        auto const& eltInserted = M_meshPartIn->addElement( e, true/*false*/ );

        if ( j >= nActiveElement )
        {
            mapGhostHdf5IdToFeelId[partId][j-nActiveElement] = eltInserted.id();
        }
    }
}

template<typename MeshType>
void PartitionIO<MeshType>::readGhostElements( std::map<rank_type,std::vector<size_type>> const& mapGhostHdf5IdToFeelId )
{

    rank_type partId = M_meshPartIn->worldComm().localRank();

    //int nVal = 3; //(idGhost,pid,idActive)
    int nVal = 1; // pid

    hsize_t globalDims[2];
    hsize_t localDims[2];
    hsize_t offset[2];

    M_HDF5IO.openTable ("elements_ghosts", globalDims);
    LOG(INFO) << "loaded ghost_elements:" << globalDims[0] << "x" << globalDims[1];

    int indexTable0 = 0, indexTable1 = 1;
    if ( M_transposeInFile )
    {
        indexTable0 = 1;
        indexTable1 = 0;
    }

    CHECK( globalDims[indexTable0] == 1 && globalDims[indexTable1] == nVal*M_numGlobalGhostElements ) << "BUG";
    if ( globalDims[indexTable1] > 0 )
    {
        localDims[indexTable0] = 1;
        localDims[indexTable1] = nVal*M_numLocalGhostElements[partId];
        offset[indexTable0] = 0;
        offset[indexTable1] = nVal*M_ghostElementsOffSet[partId];

        M_uintBuffer.resize( localDims[indexTable0]*localDims[indexTable1], 0 );
        if ( !M_uintBuffer.empty() )
            M_HDF5IO.read("elements_ghosts", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);
        else
        {
            size_type uselessValue = 0;
            M_HDF5IO.read("elements_ghosts", H5T_NATIVE_UINT, localDims, offset, &uselessValue);
        }
    }

    M_HDF5IO.closeTable("elements_ghosts");

    size_type currentBufferIndex = 0;
    auto itFindGhostIdPartId = mapGhostHdf5IdToFeelId.find( partId );

    for (size_type j = 0; j < M_numLocalGhostElements[partId]; ++j)
    {
        size_type id = itFindGhostIdPartId->second[j];
        //size_type idBB = M_uintBuffer[currentBufferIndex++];
        int pid = M_uintBuffer[currentBufferIndex++];
        //size_type idInActivePart = M_uintBuffer[currentBufferIndex++];

        auto it = M_meshPartIn->elementIterator( id, partId );
        M_meshPartIn->elements().modify( it, Feel::detail::UpdateProcessId(pid) );
        //M_meshPartIn->elements().modify( it, Feel::detail::updateIdInOthersPartitions( pid, idInActivePart ) );
        M_meshPartIn->elements().modify( it, Feel::detail::UpdateNeighborPartition( partId ) );
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
    rank_type partId = M_meshPartIn->worldComm().localRank();

    int nValFace = 1 + face_type::numPoints;//(marker,ptId1,ptId2,...)
    int nValEdge = 1 + edge_type::numPoints;// (marker,ptId1,ptId2,...)
    int nValPoint = 1 + 1;// (marker,ptId)

    hsize_t globalDims[2];
    hsize_t localDims[2];
    hsize_t offset[2];

    M_HDF5IO.openTable ("marked_subentities", globalDims);
    LOG(INFO) << "loaded marked_subentities:" << globalDims[0] << "x" << globalDims[1];
    int indexTable0 = 0, indexTable1 = 1;
    if ( M_transposeInFile )
    {
        indexTable0 = 1;
        indexTable1 = 0;
    }

    CHECK( globalDims[indexTable0] == 1 && globalDims[indexTable1] == (nValFace*M_numGlobalMarkedFaces + nValEdge*M_numGlobalMarkedEdges + nValPoint*M_numGlobalMarkedPoints ) ) << "BUG";
    if ( globalDims[indexTable1] > 0 )
    {
        localDims[indexTable0] = 1;
        localDims[indexTable1] = nValFace*M_numLocalMarkedFaces[partId] + nValEdge*M_numLocalMarkedEdges[partId] + nValPoint*M_numLocalMarkedPoints[partId];
        offset[indexTable0] = 0;
        offset[indexTable1] = nValFace*M_markedFacesOffSet[partId] + nValEdge*M_markedEdgesOffSet[partId] + nValPoint*M_markedPointsOffSet[partId];

        M_uintBuffer.resize( localDims[indexTable0]*localDims[indexTable1], 0 );
        if ( !M_uintBuffer.empty() )
            M_HDF5IO.read("marked_subentities", H5T_NATIVE_UINT, localDims, offset, &M_uintBuffer[0]);
        else
        {
            size_type uselessValue=0;
            M_HDF5IO.read("marked_subentities", H5T_NATIVE_UINT, localDims, offset, &uselessValue);
        }
    }

    M_HDF5IO.closeTable("marked_subentities");

    if ( globalDims[indexTable1] > 0 && !M_uintBuffer.empty() )
        Feel::detail::updateMarkedSubEntitiesMesh( M_uintBuffer, std::make_tuple( M_numLocalMarkedFaces[partId],M_numLocalMarkedEdges[partId],M_numLocalMarkedPoints[partId] ),
                                                   *M_meshPartIn );
}


template<typename MeshType>
void PartitionIO<MeshType>::prepareUpdateForUseStep1()
{
    rank_type rank = M_meshPartIn->worldComm().localRank();

    // put temporarly the points processId to numPart(=number of process, which are not a valid part id) in ghost partition
    // and store neigboor part at each point id
    rank_type ghostPointPidDetection = M_meshPartIn->worldComm().localSize();
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

template<typename MeshType>
void PartitionIO<MeshType>::prepareUpdateForUseStep2()
{
    auto const& theWorldComm = M_meshPartIn->worldComm();
    const rank_type nProc = theWorldComm.localSize();
    const rank_type partId = theWorldComm.localRank();

    if (nProc == 1)
        return;

    // prepare container to send  : ( pid -> (ptId1,PtId2,..),(ptId1,PtId2,..) )
    std::map< rank_type, std::vector< std::vector<size_type> > > dataToSend;
    std::map< rank_type, std::vector<size_type> > memoryMsgToSend;
    auto ghostelt_it = M_meshPartIn->beginGhostElement();
    auto const ghostelt_en = M_meshPartIn->endGhostElement();
    for ( ; ghostelt_it!=ghostelt_en ; ++ghostelt_it )
    {
        const rank_type pid = ghostelt_it->processId();
        std::vector<size_type> ptIdsInElt( mesh_type::element_type::numPoints );
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::element_type::numPoints; ++vLocId )
            ptIdsInElt[vLocId] = ghostelt_it->point( vLocId ).id();
        dataToSend[pid].push_back( ptIdsInElt );
        memoryMsgToSend[pid].push_back( ghostelt_it->id() );
    }

    // init mpi request
    int nbRequest=0;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        if ( dataToSend.find(proc) != dataToSend.end() )
            nbRequest+=2;
    }
    if ( nbRequest ==0 ) return;
    mpi::request * reqs = new mpi::request[nbRequest];
    int cptRequest=0;

    // first send
    auto itDataToSend = dataToSend.begin();
    auto const enDataToSend = dataToSend.end();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        reqs[cptRequest] = theWorldComm.localComm().isend( itDataToSend->first , 0, itDataToSend->second );
        ++cptRequest;
    }
    // first recv
    std::map< rank_type, std::vector< std::vector<size_type> > > dataToRecv;
    for ( rank_type proc=0; proc<nProc; ++proc )
    {
        // if we send a msg at proc else we recv a msg from proc
        if ( dataToSend.find(proc) != dataToSend.end() )
        {
            reqs[cptRequest] = theWorldComm.localComm().irecv( proc , 0, dataToRecv[proc] );
            ++cptRequest;
        }
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);

    // build map which allow to identify element from point ids (only for elt which touch the interprocess part)
    std::map<std::set<size_type>,size_type> mapPointIdsToEltId;
    auto elt_it = M_meshPartIn->beginElementWithProcessId( partId );
    auto const elt_en = M_meshPartIn->endElementWithProcessId( partId );
    for ( ; elt_it != elt_en ; ++elt_it )
    {
        if ( elt_it->numberOfNeighborPartitions() < 1 )
            continue;
        std::set<size_type> ptIds;
        for ( uint16_type vLocId = 0 ; vLocId < mesh_type::element_type::numPoints; ++vLocId )
            ptIds.insert( elt_it->point( vLocId ).id() );
        mapPointIdsToEltId[ptIds] = elt_it->id();
    }
    // treatment of recv request and prepare containers to re-send
    std::map<rank_type, std::vector<size_type> > dataToReSend;
    auto itDataRecv = dataToRecv.begin();
    auto const enDataRecv = dataToRecv.end();
    for ( ; itDataRecv!=enDataRecv ; ++itDataRecv )
    {
        const rank_type pid = itDataRecv->first;
        const int nDataRecv = itDataRecv->second.size();
        dataToReSend[pid].resize( nDataRecv );
        for ( int k=0; k<nDataRecv; ++k )
        {
            std::vector<size_type> const& ptIdsInEltVec = itDataRecv->second[k];
            std::set<size_type> ptIdsInEltSet(ptIdsInEltVec.begin(),ptIdsInEltVec.end());
            auto itFindElt = mapPointIdsToEltId.find( ptIdsInEltSet );
            CHECK( itFindElt != mapPointIdsToEltId.end() ) << "element not find from point ids";

            M_meshPartIn->elements().modify( M_meshPartIn->elementIterator( itFindElt->second ), Feel::detail::UpdateNeighborPartition( pid ) );//NEW

            dataToReSend[pid][k] = itFindElt->second;
        }
    }
    // send respond to the request
    cptRequest=0;
    auto itDataToReSend = dataToReSend.begin();
    auto const enDataToReSend = dataToReSend.end();
    for ( ; itDataToReSend!=enDataToReSend ; ++itDataToReSend )
    {
        reqs[cptRequest] = theWorldComm.localComm().isend( itDataToReSend->first , 0, itDataToReSend->second );
        ++cptRequest;
    }
    // recv the initial request
    std::map<rank_type, std::vector<size_type> > finalDataToRecv;
    itDataToSend = dataToSend.begin();
    for ( ; itDataToSend!=enDataToSend ; ++itDataToSend )
    {
        const rank_type pid = itDataToSend->first;
        reqs[cptRequest] = theWorldComm.localComm().irecv( pid, 0, finalDataToRecv[pid] );
        ++cptRequest;
    }
    // wait all requests
    mpi::wait_all(reqs, reqs + nbRequest);
    // delete reqs because finish comm
    delete [] reqs;

    // update ghost element with element id in active partition
    auto itFinalDataToRecv = finalDataToRecv.begin();
    auto const enFinalDataToRecv = finalDataToRecv.end();
    for ( ; itFinalDataToRecv!=enFinalDataToRecv ; ++itFinalDataToRecv)
    {
        const rank_type pid = itFinalDataToRecv->first;
        const int nDataRecv = itFinalDataToRecv->second.size();
        for ( int k=0; k<nDataRecv; ++k )
        {
            auto eltToUpdate = M_meshPartIn->elementIterator( memoryMsgToSend[pid][k],pid );
            M_meshPartIn->elements().modify( eltToUpdate, Feel::detail::updateIdInOthersPartitions( pid, itFinalDataToRecv->second[k] ) );
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
