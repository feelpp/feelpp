/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2005-2006 EPFL
    SPDX-FileCopyrightText: 2011-2020, 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
   @file geoentity.hpp
   @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   @date 2005-08-10
   @brief Base class for all geometric entities (C++20/23 modernized)
 */
#ifndef __GeoEntity_H
#define __GeoEntity_H 1

#include <algorithm>
#include <concepts>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/geoentitystorage.hpp>
#include <feel/feelpoly/order.hpp>

namespace Feel
{

/**
   \class GeoEntity
   \brief base class for all geometric entities

   @author Christophe Prud'homme
   @see
*/
template<typename Entity, typename T = double, typename IndexT = uint32_type>
class GeoEntity
    :
        boost::equality_comparable<GeoEntity<Entity,T,IndexT> >,
    boost::less_than_comparable<GeoEntity<Entity,T,IndexT> >,
    boost::less_than_comparable<GeoEntity<Entity,T,IndexT>, IndexT>,
    public Entity
{
    static inline const uint16_type nBitShiftedGeoEntityContext = 0;
    static inline const uint16_type nBitShiftedReferenceGeometry = 2;
    static inline const uint16_type nBitShiftedReferenceShapes = 8;
public:


    /** @name Typedefs
     */
    //@{

    typedef Entity super;
    typedef GeoEntity<Entity,T,IndexT> GeoShape;
    typedef GeoEntity<Entity,T,IndexT> self_type;
    typedef T value_type;
    using index_type = IndexT;
    using size_type = index_type;
    typedef typename super::topological_face_type face_type;
    typedef face_type GeoBShape;
    typedef typename Entity::edge_permutation_type edge_permutation_type;
    typedef typename Entity::face_permutation_type face_permutation_type;

    static constexpr size_type Shape = super::Shape;
    static constexpr size_type Geometry = super::Geometry;

    static constexpr uint16_type nDim = super::nDim;
    static constexpr uint16_type nRealDim = super::nRealDim;

    //! @brief True if Order is known at compile time (i.e., Order >= 0)
    static constexpr bool is_order_static = super::is_order_static;
    //! @brief True if Order is determined at runtime (i.e., Order == Dynamic)
    static constexpr bool is_order_dynamic = super::is_order_dynamic;
    //! @brief Template order parameter value (may be Dynamic = -1)
    static constexpr int nOrder_v = super::nOrder_v;
    //! @brief Static order (or 1 as placeholder for dynamic case)
    static constexpr uint16_type nOrder = super::nOrder;

    static constexpr uint16_type numVertices = super::numVertices;
    static constexpr uint16_type numFaces = super::numFaces;
    static constexpr uint16_type numGeometricFaces = super::numGeometricFaces;
    static constexpr uint16_type numTopologicalFaces = super::numTopologicalFaces;
    static constexpr uint16_type numEdges = super::numEdges;
    static constexpr uint16_type numNormals = super::numNormals;

    static constexpr uint16_type numPoints = super::numPoints;
    static constexpr uint16_type nbPtsPerVertex = super::nbPtsPerVertex;
    static constexpr uint16_type nbPtsPerEdge = super::nbPtsPerEdge;
    static constexpr uint16_type nbPtsPerFace = super::nbPtsPerFace;
    static constexpr uint16_type nbPtsPerVolume = super::nbPtsPerVolume;

    using convex_type = Entity;

    static constexpr bool is_simplex = super::is_simplex;
    static constexpr bool is_hypercube = super::is_hypercube;

    /**
     * helper class to construct the associated reference convex.
     */
    template<typename TT = double>
    struct reference_convex
    {
        typedef Reference<Entity, nDim, nOrder, nRealDim, TT> type;
    };
    template<typename TT = double>
    using reference_convex_type =  Reference<Entity, nDim, nOrder, nRealDim, TT>;

    using marker_type = Marker<flag_type/*uint16_type*/>;
    /** Compact trivially copyable entity record. */
    using core_type = GeoEntityCore<index_type>;
    /** Detached entity storage used before insertion in a mesh. */
    using detached_storage_type =
        GeoEntityDetachedMetadata<index_type, marker_type, value_type>;
    /** Compatibility name for detached sparse metadata. */
    using sparse_metadata_type = detached_storage_type;
    /** Mesh-owned connectivity and sparse metadata storage type. */
    using entity_storage_type =
        GeoEntityStorage<index_type, marker_type, value_type>;
    /** Compatibility name for mesh-owned sparse metadata storage. */
    using sparse_metadata_storage_type = entity_storage_type;
    /** Remote partition identifier mapping type. */
    using remote_id_map_type =
        typename sparse_metadata_storage_type::remote_id_map_type;
    /** Point/edge-to-element incidence set type. */
    using incidence_set_type =
        typename sparse_metadata_storage_type::incidence_set_type;
    /** Marker mapping type. */
    using marker_map_type =
        typename sparse_metadata_storage_type::marker_map_type;
    /** Cached measure metadata type. */
    using measure_metadata_type =
        typename sparse_metadata_storage_type::measure_metadata_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    GeoEntity()
        :
        super(),
        M_core{
            .id = 0,
            .entityContext =
                static_cast<uint32_type>(
                    ( MESH_ENTITY_INTERNAL <<
                      nBitShiftedGeoEntityContext ) |
                    ( Geometry << nBitShiftedReferenceGeometry ) |
                    ( Shape << nBitShiftedReferenceShapes ) )
        }
    {
    }

    explicit GeoEntity( size_type i,
                        size_type geometry = Geometry,
                        size_type shape = Shape,
                        size_type context = MESH_ENTITY_INTERNAL )
        :
        super(),
        M_core{
            .id = i,
            .entityContext =
                static_cast<uint32_type>(
                    ( context << nBitShiftedGeoEntityContext ) |
                    ( geometry << nBitShiftedReferenceGeometry ) |
                    ( shape << nBitShiftedReferenceShapes ) )
        }
    {
    }

    GeoEntity( GeoEntity const& other )
        :
        super( other ),
        M_core( other.M_core )
    {
        copySparseMetadataFrom( other );
    }

    GeoEntity( GeoEntity&& other ) noexcept
        :
        super( std::move( other ) ),
        M_core( other.M_core )
    {
        moveSparseMetadataFrom( std::move( other ) );
    }

    GeoEntity& operator=( GeoEntity const& other )
    {
        if ( this != &other )
        {
            auto* storage = M_sparseMetadataStorage;
            auto const oldId = M_core.id;
            auto const oldOffset = M_core.connectivityOffset;
            auto const oldCount = M_core.connectivityCount;
            super::operator=( other );
            M_core = other.M_core;
            replaceSparseMetadata(
                other.sparseMetadataSnapshot(), storage, oldId,
                oldOffset, oldCount );
        }
        return *this;
    }

    GeoEntity& operator=( GeoEntity&& other )
    {
        if ( this != &other )
        {
            auto* storage = M_sparseMetadataStorage;
            auto const oldId = M_core.id;
            auto const oldOffset = M_core.connectivityOffset;
            auto const oldCount = M_core.connectivityCount;
            super::operator=( std::move( other ) );
            M_core = other.M_core;
            replaceSparseMetadata(
                other.releaseSparseMetadata(), storage, oldId,
                oldOffset, oldCount );
        }
        return *this;
    }

    ~GeoEntity() override
    {}

    //@}

    /** @name Operator overloads
     */
    //@{
    bool operator==( GeoEntity const& e ) const
    {
        return M_core.id == e.id();
    }
    bool operator<( GeoEntity const& e ) const
    {
        return M_core.id < e.id();
    }

    bool operator<( size_type __i ) const
    {
        return M_core.id < __i;
    }

    //@}

    /** @name Accessors
     */
    //@{

    size_type id() const noexcept
    {
        return M_core.id;
    }

    /**
     * @brief Return the compact hot entity record.
     * @return immutable core record
     */
    [[nodiscard]] core_type const& core() const noexcept
    {
        return M_core;
    }

    /**
     * @brief Attach this entity to mesh-owned storage.
     *
     * Any metadata accumulated while the entity was detached is transferred
     * to the storage. Copying an attached entity intentionally creates a
     * detached metadata snapshot so copies cannot alias another mesh.
     *
     * @param storage mesh-owned storage for this entity dimension
     */
    void attachEntityStorage( entity_storage_type& storage )
    {
        if ( M_sparseMetadataStorage == &storage )
            return;

        auto data = releaseSparseMetadata();
        M_sparseMetadataStorage = &storage;
        auto connectivity = std::move( data.connectivity );
        data.connectivity.clear();
        if ( !data.empty() )
            M_sparseMetadataStorage->assign(
                M_core.id, std::move( data ) );
        setConnectivity(
            std::span<index_type const>(
                connectivity.data(), connectivity.size() ) );
    }

    /**
     * @brief Compatibility wrapper for attaching mesh-owned storage.
     * @param storage mesh-owned entity storage
     */
    void attachSparseMetadataStorage(
        sparse_metadata_storage_type& storage )
    {
        attachEntityStorage( storage );
    }

    /**
     * @brief Return true when metadata is owned by a mesh sidecar.
     */
    [[nodiscard]] bool hasSparseMetadataStorage() const noexcept
    {
        return M_sparseMetadataStorage != nullptr;
    }

    /** @brief Return true when connectivity and sidecars are mesh-owned. */
    [[nodiscard]] bool hasEntityStorage() const noexcept
    {
        return M_sparseMetadataStorage != nullptr;
    }

    /**
     * @brief Return a detached copy of all sparse metadata.
     */
    [[nodiscard]] sparse_metadata_type sparseMetadataSnapshot() const
    {
        if ( M_sparseMetadataStorage )
        {
            auto data = M_sparseMetadataStorage->snapshot( id() );
            auto connectivity = connectivityIds();
            data.connectivity.assign(
                connectivity.begin(), connectivity.end() );
            return data;
        }
        return M_detachedSparseMetadata
                   ? *M_detachedSparseMetadata
                   : sparse_metadata_type{};
    }

    /**
     * @brief Return the offset in mesh-owned connectivity storage.
     */
    [[nodiscard]] uint32_type connectivityOffset() const noexcept
    {
        return M_core.connectivityOffset;
    }

    /**
     * @brief Return the number of entries in the connectivity span.
     */
    [[nodiscard]] uint16_type connectivityCount() const noexcept
    {
        return M_core.connectivityCount;
    }

    /**
     * @brief Set the mesh-owned connectivity span.
     * @param offset first connectivity entry
     * @param count number of connectivity entries
     */
    void setConnectivitySpan(
        uint32_type offset, uint16_type count ) noexcept
    {
        M_core.connectivityOffset = offset;
        M_core.connectivityCount = count;
    }

    /**
     * @brief Return connectivity identifiers in local mathematical order.
     */
    [[nodiscard]] std::span<index_type const> connectivityIds() const
    {
        if ( M_sparseMetadataStorage )
            return M_sparseMetadataStorage->connectivity(
                M_core.connectivityOffset,
                M_core.connectivityCount );
        if ( M_detachedSparseMetadata )
            return M_detachedSparseMetadata->connectivity;
        return {};
    }

    /**
     * @brief Replace connectivity with one contiguous mesh-owned span.
     * @param connectivity identifiers in local mathematical order
     */
    void setConnectivity(
        std::span<index_type const> connectivity )
    {
        if ( connectivity.size() >
             std::numeric_limits<uint16_type>::max() )
            throw std::length_error(
                "entity connectivity exceeds 16-bit cardinality" );

        if ( M_sparseMetadataStorage )
        {
            M_sparseMetadataStorage->releaseConnectivity(
                M_core.connectivityOffset,
                M_core.connectivityCount );
            M_core.connectivityOffset =
                M_sparseMetadataStorage->allocateConnectivity(
                    connectivity );
        }
        else
        {
            if ( connectivity.empty() &&
                 !M_detachedSparseMetadata )
            {
                M_core.connectivityOffset =
                    std::numeric_limits<uint32_type>::max();
                M_core.connectivityCount = 0;
                return;
            }
            auto& detached = mutableDetachedSparseMetadata();
            detached.connectivity.assign(
                connectivity.begin(), connectivity.end() );
            M_core.connectivityOffset =
                std::numeric_limits<uint32_type>::max();
        }
        M_core.connectivityCount =
            static_cast<uint16_type>( connectivity.size() );
    }

    /**
     * @brief Update one local connectivity identifier.
     * @param localIndex local mathematical connectivity index
     * @param id external entity identifier
     * @param cardinality complete entity connectivity cardinality
     */
    void setConnectivityEntry(
        uint16_type localIndex, index_type id,
        uint16_type cardinality )
    {
        if ( localIndex >= cardinality )
            throw std::out_of_range(
                "invalid local entity connectivity index" );
        bool const hasCompleteConnectivity =
            M_core.connectivityCount == cardinality &&
            ( M_sparseMetadataStorage
                  ? ( cardinality == 0 ||
                      M_core.connectivityOffset !=
                          std::numeric_limits<uint32_type>::max() )
                  : ( M_detachedSparseMetadata &&
                      M_detachedSparseMetadata->connectivity.size() ==
                          cardinality ) );
        if ( !hasCompleteConnectivity )
        {
            std::vector<index_type> connectivity(
                cardinality, invalid_v<index_type> );
            if ( M_sparseMetadataStorage &&
                 M_core.connectivityOffset !=
                     std::numeric_limits<uint32_type>::max() )
            {
                auto current = connectivityIds();
                std::copy_n(
                    current.begin(),
                    std::min(
                        current.size(), connectivity.size() ),
                    connectivity.begin() );
            }
            else if ( M_detachedSparseMetadata )
            {
                auto const& current =
                    M_detachedSparseMetadata->connectivity;
                std::copy_n(
                    current.begin(),
                    std::min(
                        current.size(), connectivity.size() ),
                    connectivity.begin() );
            }
            connectivity[localIndex] = id;
            setConnectivity( connectivity );
            return;
        }

        if ( M_sparseMetadataStorage )
            M_sparseMetadataStorage->setConnectivityEntry(
                M_core.connectivityOffset,
                M_core.connectivityCount, localIndex, id );
        else
            mutableDetachedSparseMetadata()
                .connectivity[localIndex] = id;
    }

    /**
     * @brief Remove this entity's connectivity and sparse sidecars.
     *
     * Containers call this before erasing the entity record.
     */
    void clearEntityStorage()
    {
        if ( M_sparseMetadataStorage )
        {
            M_sparseMetadataStorage->releaseConnectivity(
                M_core.connectivityOffset,
                M_core.connectivityCount );
            M_sparseMetadataStorage->erase( id() );
            M_sparseMetadataStorage = nullptr;
        }
        M_detachedSparseMetadata.reset();
        M_core.connectivityOffset =
            std::numeric_limits<uint32_type>::max();
        M_core.connectivityCount = 0;
    }


    /**
     * the dimension of the reference shape
     *
     *
     * @return the dimension of the reference shape
     */
    constexpr uint16_type refDim() const
    {
        return super::nDim;
    }

    /**
     * number of points on the reference shape
     *
     * @return the number of points on the reference shape
     */
    [[nodiscard]] constexpr uint16_type nPoints() const
    {
        return super::numPoints;
    }

    /**
     * @brief Get the polynomial order of the geometric entity.
     * For dynamic order, this returns the runtime order from the underlying shape.
     * @return The order value.
     */
    [[nodiscard]] constexpr uint16_type order() const noexcept
    {
        if constexpr ( is_order_static )
            return nOrder;
        else
            return super::order();
    }

    /**
     * number of vertices on the reference shape
     *
     * @return the number of vertices on the reference shape
     */
    constexpr uint16_type nVertices() const
    {
        return super::numVertices;
    }

    /**
     * number of edges on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    constexpr uint16_type nEdges() const
    {
        return super::numEdges;
    }

    /**
     * number of faces on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    constexpr uint16_type nFaces() const
    {
        return super::numFaces;
    }

    /**
     * number of topological faces on the reference shape
     *
     * @return the number of topological faces on the reference shape
     */
    constexpr uint16_type nTopologicalFaces() const
    {
        return super::numTopologicalFaces;
    }

    /**
     * number of faces on the reference shape
     *
     * @return the number of edges on the reference shape
     */
    constexpr uint16_type nGeometricFaces() const
    {
        return super::numGeometricFaces;
    }

    /**
     * number of normals on the reference shape
     *
     * @return the number of normals on the reference shape
     */
    constexpr uint16_type nNormals() const
    {
        return super::numNormals;
    }


    /**
     *
     *
     *
     * @return true if the entoty has the shape \c __shape, false otherwise
     */
    bool hasShape( size_type __shape ) const
    {
        return hasEntityContext(
            __shape << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a volume
     */
    bool isAVolume() const
    {
        return hasEntityContext(
            GEOMETRY_VOLUME << nBitShiftedReferenceGeometry );
    }

    /**
     * @return true of the entity is a surface
     */
    bool isASurface() const
    {
        return hasEntityContext(
            GEOMETRY_SURFACE << nBitShiftedReferenceGeometry );
    }

    /**
     * @return true of the entity is a line
     */
    bool isALine() const
    {
        return hasEntityContext(
            GEOMETRY_LINE << nBitShiftedReferenceGeometry );
    }

    /**
     * @return true of the entity is a point
     */
    bool isAPoint() const
    {
        return hasEntityContext(
            GEOMETRY_POINT << nBitShiftedReferenceGeometry );
    }

    /**
     * @return true of the entity is a shape point
     */
    bool isAPointShape() const
    {
        return hasEntityContext(
            SHAPE_POINT << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a shape line
     */
    bool isALineShape() const
    {
        return hasEntityContext(
            SHAPE_LINE << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a triangle shape
     */
    bool isATriangleShape() const
    {
        return hasEntityContext(
            SHAPE_TRIANGLE << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a quadrangle
     */
    bool isAQuadrangleShape() const
    {
        return hasEntityContext(
            SHAPE_QUAD << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a tetrahedra shape
     */
    bool isATetrahedraShape() const
    {
        return hasEntityContext(
            SHAPE_TETRA << nBitShiftedReferenceShapes );
    }

    /**
     * @return true of the entity is a hexahedra
     */
    bool isAHexahedraShape() const
    {
        return hasEntityContext(
            SHAPE_HEXA << nBitShiftedReferenceShapes );
    }

    /**
     * @return true if the shape is linear, false otherwise
     */
    bool isLinear() const
    {
        return hasEntityContext(
            SHAPE_LINEAR << nBitShiftedReferenceShapes );
    }

    /**
     * @return true if the shape is bilinear, false otherwise
     */
    bool isBilinear() const
    {
        return hasEntityContext(
            SHAPE_BILINEAR << nBitShiftedReferenceShapes );
    }

    /**
     * @return true if the shape is quadratic, false otherwise
     */
    bool isQuadratic() const
    {
        return hasEntityContext(
            SHAPE_QUADRATIC << nBitShiftedReferenceShapes );
    }

    /**
     * @return true if the entity is internal, false otherwise
     */
    bool isInternal() const
    {
        return hasEntityContext(
            MESH_ENTITY_INTERNAL << nBitShiftedGeoEntityContext );
    }


    /**
     * Tells if  item is on the boundary
     * @return true if on boundary, false otherwise
     */
    bool isOnBoundary() const noexcept
    {
        return hasEntityContext(
            MESH_ENTITY_BOUNDARY << nBitShiftedGeoEntityContext );
    }

    /**
     * maximum dimension of the entity of the element touching the boundary
     */
    uint16_type boundaryEntityDimension() const noexcept
    {
        return M_core.boundaryEntityDimension;
    }
    /**
     * \return \c true if ghost cell, \c false otherwise
     */
    bool isGhostCell() const noexcept
    {
        //return (this->worldComm().localRank()!=M_pid);
        //mpi::communicator world;
        //return (world.rank()!=M_pid);
        return M_core.processIdInPartition != M_core.processId;
    }

    /**
     * \return the processor id of the entity
     */
    rank_type processId() const noexcept
    {
        return M_core.processId;
    }

    /**
     * set the processor id of the entity
     & \param pid processor id
     */
    void setProcessId( rank_type pid ) noexcept
    {
        M_core.processId = pid;
    }

    /**
     * \return the processor id of the entity
     */
    rank_type pidInPartition() const noexcept
    {
        return M_core.processIdInPartition;
    }
    /**
     * set the processor id of the entity
     & \param pid processor id
     */
    void setProcessIdInPartition( rank_type pid ) noexcept
    {
        M_core.processIdInPartition = pid;
    }

    /**
     * \return the partition id
     */
    rank_type partitionId() const noexcept
    {
        return M_core.processId;
    }

    /**
     * \return the number of partition the element is linked to including the
     * partition to which it belongs
     */
    rank_type numberOfPartitions() const noexcept
    {
        return static_cast<rank_type>(
            idInOthersPartitions().size() + 1 );
    }

    /**
     * \return the number of partition the element is linked to
     */
    rank_type numberOfNeighborPartitions() const
    {
        return static_cast<rank_type>(
            idInOthersPartitions().size() );
    }

    /**
     * \return true if the element is linked to other partitions through one of
     * more of its faces
     */
    bool isLinkedToOtherPartitions() const
    {
        return !idInOthersPartitions().empty();
    }

    /**
     * \return the id of the entity in a partition pid
     */
    size_type idInOthersPartitions( rank_type pid ) const
    {
        if ( pid == this->pidInPartition() )
            return this->id();
        auto const& remoteIds = idInOthersPartitions();
        auto const it = remoteIds.find( pid );
        DCHECK( it != remoteIds.end() )
            << " local id " << this->id() << " is unknown for this partition " << pid << "\n";
        return it->second;
    }

    /**
     * \return idInOthersPartitions map
     */
    remote_id_map_type const& idInOthersPartitions() const noexcept
    {
        if ( M_sparseMetadataStorage )
            return M_sparseMetadataStorage->remoteIds( id() );
        return M_detachedSparseMetadata
                   ? M_detachedSparseMetadata->remoteIds
                   : emptySparseMetadata().remoteIds;
    }

    remote_id_map_type const& neighborProcessIds() const noexcept
    {
        return idInOthersPartitions();
    }

    /**
     * clear id in others partitions container
     */
    void clearIdInOthersPartitions()
    {
        if ( M_sparseMetadataStorage )
            M_sparseMetadataStorage->clearRemoteIds( id() );
        else if ( M_detachedSparseMetadata )
            M_detachedSparseMetadata->remoteIds.clear();
    }

    /**
     * \return \c true if active, \c false otherwise
     *
     * \note for now it is a dummy function that returns always true,
     * will change when work on AMR starts
     */
    bool active() const
    {
        return true;
    }

    /**
     * \return the measure of the entity
     */
    virtual value_type measure() const { return value_type{0}; }

    //@}

    /** @name  Mutators
     */
    //@{
    void setId( size_type id )
    {
        if ( M_sparseMetadataStorage )
            M_sparseMetadataStorage->rekey( M_core.id, id );
        M_core.id = id;
    }

    /**
     * set the boundary flag
     * @param b true if the item is on the boundary, false otherwise
     */
    void setOnBoundary( bool b, uint16_type ent_d = invalid_uint16_type_value )
    {
        if ( b )
        {
            setEntityContext(
                MESH_ENTITY_BOUNDARY << nBitShiftedGeoEntityContext );
            clearEntityContext(
                MESH_ENTITY_INTERNAL << nBitShiftedGeoEntityContext );
        }
        else
        {
            clearEntityContext(
                MESH_ENTITY_BOUNDARY << nBitShiftedGeoEntityContext );
            setEntityContext(
                MESH_ENTITY_INTERNAL << nBitShiftedGeoEntityContext );
        }
        M_core.boundaryEntityDimension = ent_d;
    }


    /**
     * set id in a partition pid of the entity
     */
    void setIdInOtherPartitions( rank_type pid, size_type id )
    {
        mutableRemoteIds()[pid] = id;
    }

    /**
     * set (partition,id) in other partitions of the entity
     */
    void setIdInOtherPartitions( remote_id_map_type const& iop )
    {
        assignRemoteIds( iop );
    }


    /**
     * set (partition,id) in other partitions of the entity
     */
    void setIdInOtherPartitions( remote_id_map_type&& iop )
    {
        assignRemoteIds( std::move( iop ) );
    }

    template <typename Q>
    void setNeighborProcessIds( Q && pids )
    {
        auto& remoteIds = mutableRemoteIds();
        remoteIds = std::forward<Q>( pids );
        std::for_each(
            remoteIds.begin(), remoteIds.end(),
            []( auto& e )
            {
                std::get<1>( e ) = invalid_v<index_type>;
            } );
    }

    /**
     * \return the number of partition the element is linked to
     */
    void setNeighborPartitionIds( std::vector<rank_type> const& npids )
    {
        clearIdInOthersPartitions();
        auto& remoteIds = mutableRemoteIds();
        for ( rank_type p : npids)
        {
            CHECK( p != pidInPartition() ) << fmt::format("neighor pid {} should not be the current process id",p );
            remoteIds.emplace( p, invalid_v<index_type> );
        }
    }

    void addNeighborPartitionId( rank_type p )
    {
        CHECK( p != pidInPartition() ) << fmt::format("neighor pid {} should not be the current process id",p );
        mutableRemoteIds().emplace( p, invalid_v<index_type> );
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * eToP(i,j) = localId of jth point on ith local edge
     */
    static uint16_type eToP( uint16_type const __localEdge, uint16_type const __point )
    {
        return super::e2p( __localEdge, __point );
    }

    /**
     * fToP(i,j) = localId of jth point on ith local edge
     */
    static uint16_type fToP( uint16_type const __localFace, uint16_type const __point )
    {
        return super::f2p( __localFace, __point );
    }

    /**
     * fToE(i,j) = localId of jth edge on ith local face
     */
    static uint16_type fToE( uint16_type const __localFace, uint16_type const __edge )
    {
        return super::f2e( __localFace, __edge );
    }

    /**
     * add a new element to which the point belongs
     */
    self_type& addElement( size_type e, int id_in_element = 0 )
    {
        mutableElementIncidence().insert(
            std::make_pair( e, id_in_element ) );
        return *this;
    }

    /**
     * \return the number of elements whom the point belongs to
     */
    size_type numberOfElements() const
    {
        return elements().size();
    }

    /**
     * \return the set of ids of elements whom the point belongs to
     */
    incidence_set_type const& elements() const
    {
        if ( M_sparseMetadataStorage )
            return M_sparseMetadataStorage->incidence( id() );
        return M_detachedSparseMetadata
                   ? M_detachedSparseMetadata->elementIncidence
                   : emptySparseMetadata().elementIncidence;
    }
    incidence_set_type& elements()
    {
        return mutableElementIncidence();
    }

    //@}

    /**
     * set the tags associated to the points
     * - tags[0] physical region
     * - tags[1] elementary region
     * - tags[2] particular region
     */
    virtual void setTags( std::vector<int> const& tags )
    {
        auto& markerValues = mutableMarkers();
        markerValues[1].assign( tags[0] );
        if ( tags.size() > 1 )
            markerValues[2].assign( tags[1] );

        if ( tags.size() > 2 )
        {
            this->setProcessId( tags[3] );

            if ( tags[2] > 1 )
            {
                // ghosts
                std::vector<rank_type> p( tags[2] - 1 );

                for ( size_type i = 0; i < p.size(); ++i )
                {
                    p[i] = tags[4 + i];
                }

                this->setNeighborPartitionIds( p );
            }
        }
    }

    //! return all markers
    marker_map_type const&
    markers() const
    {
        if ( M_sparseMetadataStorage )
            return M_sparseMetadataStorage->markers( id() );
        return M_detachedSparseMetadata
                   ? M_detachedSparseMetadata->markers
                   : emptySparseMetadata().markers;
    }
    //! set all markers
    void setMarkers( marker_map_type const& markers )
    {
        assignMarkers( markers );
    }

    //! return true if has marker type id k
    bool hasMarkerType( uint16_type k ) const
    {
        auto const& markerValues = markers();
        auto itFindMarker = markerValues.find( k );
        if ( itFindMarker == markerValues.end() )
            return false;
        if ( itFindMarker->second.isOff() )
            return false;
        return true;
    }
    //! return true if has marker type id k
    bool hasMarker( uint16_type k ) const { return this->hasMarkerType( k ); }

    //! return marker with type id k
    marker_type const& marker( uint16_type k ) const
    {
        DCHECK( this->hasMarkerType( k ) ) << "no marker type " << k;
        return markers().find( k )->second;
    }
    //! set marker type id k with flag v
    void setMarker( uint16_type k, flag_type v )
    {
        mutableMarkers()[k].assign( v );
    }
    //! add flag v to the marker type id k
    void addMarker( uint16_type k, flag_type v )
    {
        mutableMarkers()[k].insert( v );
    }
    //! set marker type id k with flag v
    template <typename TT>
        requires is_iterable_of_v<TT, flag_type>
    void setMarker( uint16_type k,TT const& vs )
    {
        mutableMarkers()[k].assign( vs );
    }
    //! add flags vs to the marker type id k
    template <typename TT>
        requires is_iterable_of_v<TT, flag_type>
    void addMarker( uint16_type k, TT const& vs )
    {
        for ( auto const& v : vs )
            this->addMarker( k, v );
    }

    // methods for marker type id  1
    bool hasMarker() const
    {
        return this->hasMarkerType( 1 );
    }
    marker_type const& marker() const
    {
        return this->marker( 1 );
    }
    void setMarker( flag_type v )
    {
        this->setMarker( 1, v );
    }
    void addMarker( flag_type v )
    {
        this->addMarker( 1, v );
    }
    template <typename TT>
        requires is_iterable_of_v<TT, flag_type>
    void setMarker( TT const& vs )
    {
        this->setMarker( 1, vs );
    }
    template <typename TT>
        requires is_iterable_of_v<TT, flag_type>
    void addMarker( TT const& vs )
    {
        this->addMarker( 1, vs );
    }

    // methods for marker type id  2
    bool hasMarker2() const
    {
        return this->hasMarkerType( 2 );
    }
    marker_type const& marker2() const
    {
        return this->marker( 2 );
    }
    void setMarker2( flag_type v )
    {
        this->setMarker( 2, v );
    }

    // methods for marker type id  3
    bool hasMarker3() const
    {
        return this->hasMarkerType( 3 );
    }
    marker_type const& marker3() const
    {
        return this->marker( 3 );
    }
    void setMarker3( flag_type v )
    {
        this->setMarker( 3, v );
    }

    /**
     * @brief Return cached measure metadata or an immutable empty payload.
     */
    [[nodiscard]] measure_metadata_type const&
    measureMetadata() const noexcept
    {
        if ( M_sparseMetadataStorage )
            return M_sparseMetadataStorage->measures( id() );
        return M_detachedSparseMetadata
                   ? M_detachedSparseMetadata->measures
                   : emptySparseMetadata().measures;
    }

    /**
     * @brief Return mutable cached measure metadata.
     *
     * Only the measure sidecar is created; markers, incidence, and
     * cross-partition metadata remain absent.
     */
    [[nodiscard]] measure_metadata_type& mutableMeasureMetadata() const
    {
        if ( M_sparseMetadataStorage )
            return M_sparseMetadataStorage->mutableMeasures( id() );
        return mutableDetachedSparseMetadata().measures;
    }

protected:

private:

    [[nodiscard]] bool hasEntityContext( size_type bits ) const noexcept
    {
        return ( M_core.entityContext &
                 static_cast<uint32_type>( bits ) ) != 0;
    }

    void setEntityContext( size_type bits ) noexcept
    {
        M_core.entityContext |= static_cast<uint32_type>( bits );
    }

    void clearEntityContext( size_type bits ) noexcept
    {
        M_core.entityContext &=
            ~static_cast<uint32_type>( bits );
    }

    [[nodiscard]] static sparse_metadata_type const&
    emptySparseMetadata() noexcept
    {
        static sparse_metadata_type const value;
        return value;
    }

    [[nodiscard]] sparse_metadata_type&
    mutableDetachedSparseMetadata() const
    {
        if ( !M_detachedSparseMetadata )
            M_detachedSparseMetadata =
                std::make_unique<sparse_metadata_type>();
        return *M_detachedSparseMetadata;
    }

    [[nodiscard]] remote_id_map_type& mutableRemoteIds()
    {
        if ( M_sparseMetadataStorage )
            return M_sparseMetadataStorage->mutableRemoteIds( id() );
        return mutableDetachedSparseMetadata().remoteIds;
    }

    template <typename RemoteIds>
    void assignRemoteIds( RemoteIds&& remoteIds )
    {
        if ( remoteIds.empty() )
        {
            clearIdInOthersPartitions();
            return;
        }
        mutableRemoteIds() = std::forward<RemoteIds>( remoteIds );
    }

    [[nodiscard]] incidence_set_type& mutableElementIncidence()
    {
        if ( M_sparseMetadataStorage )
            return M_sparseMetadataStorage->mutableIncidence( id() );
        return mutableDetachedSparseMetadata().elementIncidence;
    }

    [[nodiscard]] marker_map_type& mutableMarkers()
    {
        if ( M_sparseMetadataStorage )
            return M_sparseMetadataStorage->mutableMarkers( id() );
        return mutableDetachedSparseMetadata().markers;
    }

    template <typename Markers>
    void assignMarkers( Markers&& markers )
    {
        if ( markers.empty() )
        {
            if ( M_sparseMetadataStorage )
                M_sparseMetadataStorage->clearMarkers( id() );
            else if ( M_detachedSparseMetadata )
                M_detachedSparseMetadata->markers.clear();
            return;
        }
        mutableMarkers() = std::forward<Markers>( markers );
    }

    void copySparseMetadataFrom( GeoEntity const& other )
    {
        auto data = other.sparseMetadataSnapshot();
        M_core.connectivityOffset =
            std::numeric_limits<uint32_type>::max();
        M_core.connectivityCount =
            static_cast<uint16_type>( data.connectivity.size() );
        if ( !data.empty() )
            M_detachedSparseMetadata =
                std::make_unique<sparse_metadata_type>(
                    std::move( data ) );
    }

    void moveSparseMetadataFrom( GeoEntity&& other )
    {
        if ( other.M_sparseMetadataStorage )
            M_sparseMetadataStorage =
                std::exchange( other.M_sparseMetadataStorage, nullptr );
        else
            M_detachedSparseMetadata =
                std::move( other.M_detachedSparseMetadata );
    }

    [[nodiscard]] sparse_metadata_type releaseSparseMetadata()
    {
        if ( M_sparseMetadataStorage )
        {
            auto data = sparseMetadataSnapshot();
            M_sparseMetadataStorage->releaseConnectivity(
                M_core.connectivityOffset,
                M_core.connectivityCount );
            M_sparseMetadataStorage->erase( id() );
            M_sparseMetadataStorage = nullptr;
            M_core.connectivityOffset =
                std::numeric_limits<uint32_type>::max();
            M_core.connectivityCount =
                static_cast<uint16_type>(
                    data.connectivity.size() );
            return data;
        }
        if ( !M_detachedSparseMetadata )
            return {};
        auto data = std::move( *M_detachedSparseMetadata );
        M_detachedSparseMetadata.reset();
        return data;
    }

    void replaceSparseMetadata(
        sparse_metadata_type data,
        sparse_metadata_storage_type* storage,
        size_type oldId,
        uint32_type oldConnectivityOffset,
        uint16_type oldConnectivityCount )
    {
        M_detachedSparseMetadata.reset();
        M_sparseMetadataStorage = storage;
        if ( M_sparseMetadataStorage )
        {
            M_sparseMetadataStorage->releaseConnectivity(
                oldConnectivityOffset, oldConnectivityCount );
            M_sparseMetadataStorage->erase( oldId );
        }

        auto connectivity = std::move( data.connectivity );
        data.connectivity.clear();
        M_core.connectivityOffset =
            std::numeric_limits<uint32_type>::max();
        M_core.connectivityCount = 0;
        if ( M_sparseMetadataStorage )
        {
            if ( !data.empty() )
                M_sparseMetadataStorage->assign(
                    id(), std::move( data ) );
            setConnectivity( connectivity );
        }
        else
        {
            data.connectivity = std::move( connectivity );
            M_core.connectivityCount =
                static_cast<uint16_type>(
                    data.connectivity.size() );
            if ( !data.empty() )
            M_detachedSparseMetadata =
                std::make_unique<sparse_metadata_type>(
                    std::move( data ) );
        }
    }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        DVLOG( 2 ) << "Serializing GeoEntity...\n";
        ar& M_core.id;
        ar& M_core.entityContext;
        ar& M_core.processId;
        ar& M_core.processIdInPartition;

        auto remoteIds = idInOthersPartitions();
        auto markerValues = markers();
        ar& remoteIds;
        ar& markerValues;

        if constexpr ( Archive::is_loading::value )
        {
            assignRemoteIds( std::move( remoteIds ) );
            assignMarkers( std::move( markerValues ) );
        }
    }

private:
    /** Compact hot data. */
    core_type M_core;
    /** Mesh-owned sparse metadata storage, null while detached. */
    sparse_metadata_storage_type* M_sparseMetadataStorage = nullptr;
    /** Lazy sparse metadata used only by detached entities. */
    mutable std::unique_ptr<sparse_metadata_type>
        M_detachedSparseMetadata;
};

typedef GeoEntity<Simplex<0, 1> > GeoPoint;

// simplices
typedef GeoEntity<Simplex<1, 1> > LinearLine;
typedef GeoEntity<Simplex<2, 1> > LinearTriangle;
typedef GeoEntity<Simplex<3, 1> > LinearTetra;
typedef GeoEntity<Simplex<1, 2> > QuadraticLine;
typedef GeoEntity<Simplex<2, 2> > QuadraticTriangle;
typedef GeoEntity<Simplex<3, 2> > QuadraticTetra;

// tensor products
typedef GeoEntity<Hypercube<2, 1> > LinearQuad;
typedef GeoEntity<Hypercube<3, 1> > LinearHexa;
typedef GeoEntity<Hypercube<2, 2> > QuadraticQuad;
typedef GeoEntity<Hypercube<3, 2> > QuadraticHexa;


namespace detail
{
//! return true if the entity satisfy the predicate partitions defined by EPT and part
template <entity_process_t EPT, typename EntityType>
bool checkPartitionPredicate( EntityType const& entity, rank_type part )
{
    if constexpr ( EPT == entity_process_t::LOCAL_ONLY )
        return entity.processId() == part;
    else if constexpr ( EPT == entity_process_t::GHOST_ONLY )
        return entity.processId() != part;
    else if constexpr ( EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
    {
        if constexpr ( is_topological_face<EntityType>::value )
            return entity.processId() == part || entity.isInterProcessDomain( part );
        else
        {
            //CHECK( false ) << "TODO";
            // only for element entity, because INTERPROCESS has no impact (TODO check if entity is a mesh element)
            return entity.processId() == part;
        }
    }
    return true;
}

//! return true if the entity satisfy the predicate partitions defined by EPT and part
template <entity_process_t EPT, typename EntityType, typename InterprocessMappingType>
bool checkPartitionPredicate( EntityType const& entity, rank_type part, InterprocessMappingType const& interprocessMapping )
{
    if constexpr ( EPT == entity_process_t::LOCAL_ONLY )
        return entity.processId() == part;
    else if constexpr ( EPT == entity_process_t::GHOST_ONLY )
        return entity.processId() != part;
    else if constexpr ( EPT == entity_process_t::LOCAL_AND_INTERPROCESS_ONLY )
        return entity.processId() == part || interprocessMapping.find( entity.id() ) != interprocessMapping.end();
    else
        return true;
}


}

} // Feel

#endif /* __GeoEntity_H */
