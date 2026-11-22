/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file entityview.hpp
 * @brief Value-semantic non-owning views of mathematical mesh entities.
 */
#if !defined( FEELPP_FEELMESH_ENTITYVIEW_HPP )
#define FEELPP_FEELMESH_ENTITYVIEW_HPP 1

#include <feel/feelcore/feeltypes.hpp>
#include <feel/feelmesh/concepts.hpp>

#include <Eigen/Core>

#include <concepts>
#include <cstdint>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace Feel
{

/**
 * @class GeometricPointView
 * @brief Non-owning point view backed by mesh coordinate storage.
 *
 * The view has value semantics and does not expose pointer identity as a
 * mathematical property. It is invalidated when its point is erased or when
 * the owning mesh is destroyed. Coordinate-map lifetime follows the owning
 * `GeoCoordinateStorage`.
 *
 * @tparam MeshType owning mesh type
 */
template <typename MeshType>
class GeometricPointView
{
  public:
    /** Owning mesh type. */
    using mesh_type = std::remove_cvref_t<MeshType>;
    /** Compatibility point type. */
    using point_type = typename mesh_type::point_type;
    /** Coordinate scalar type. */
    using value_type = typename mesh_type::value_type;
    /** External identifier type. */
    using index_type = typename mesh_type::index_type;
    /** Dense local coordinate-column index. */
    using local_index_type = uint32_type;
    /** Zero-copy point-coordinate map type. */
    using coordinate_map_type =
        Eigen::Map<Eigen::Matrix<
            value_type, mesh_type::nRealDim, 1> const>;

    /** Construct an empty, non-dereferenceable view. */
    GeometricPointView() = default;

    /**
     * @brief Bind a view to a point identifier.
     * @param mesh owning mesh
     * @param id external point identifier
     */
    GeometricPointView( mesh_type const& mesh, index_type id )
        : M_mesh( std::addressof( mesh ) ),
          M_id( id ),
          M_localIndex( mesh.pointLocalIndex( id ) )
    {
    }

    /** @brief Return whether this view is bound to a mesh point. */
    [[nodiscard]] explicit operator bool() const noexcept
    {
        return M_mesh != nullptr;
    }

    /** @brief Return the external point identifier. */
    [[nodiscard]] index_type id() const noexcept
    {
        return M_id;
    }

    /** @brief Return the dense local coordinate-column index. */
    [[nodiscard]] local_index_type localIndex() const noexcept
    {
        return M_localIndex;
    }

    /** @brief Return a zero-copy coordinate-column map. */
    [[nodiscard]] coordinate_map_type coordinates() const
    {
        this->requireBound();
        auto const allCoordinates = M_mesh->pointCoordinates();
        return coordinate_map_type(
            allCoordinates.data() +
            static_cast<Eigen::Index>( M_localIndex ) *
                mesh_type::nRealDim );
    }

    /**
     * @brief Return one physical coordinate component.
     * @param component physical component
     */
    [[nodiscard]] value_type
    coordinate( uint16_type component ) const
    {
        if ( component >= mesh_type::nRealDim )
            throw std::out_of_range(
                "invalid physical coordinate component" );
        return this->coordinates()( component );
    }

    /** @brief Return the compatibility point facade. */
    [[nodiscard]] point_type const& point() const
    {
        this->requireBound();
        return M_mesh->point( M_id );
    }

    /** @brief Return the owning mesh. */
    [[nodiscard]] mesh_type const& mesh() const
    {
        this->requireBound();
        return *M_mesh;
    }

    /** @brief Return whether the point is on the mesh boundary. */
    [[nodiscard]] bool isOnBoundary() const
    {
        return this->point().isOnBoundary();
    }

    /** @brief Return the owning process identifier. */
    [[nodiscard]] rank_type processId() const
    {
        return this->point().processId();
    }

    /** @brief Return all mathematical marker components. */
    [[nodiscard]] decltype(auto) markers() const
    {
        return this->point().markers();
    }

    /**
     * @brief Return one mathematical marker component.
     * @param component marker component
     */
    [[nodiscard]] decltype(auto)
    marker( uint16_type component ) const
    {
        return this->point().marker( component );
    }

    /** @brief Compare view handles by owner and local point. */
    [[nodiscard]] friend bool operator==(
        GeometricPointView const&,
        GeometricPointView const& ) = default;

  private:
    void requireBound() const
    {
        if ( !M_mesh )
            throw std::logic_error(
                "unbound geometric point view" );
    }

    mesh_type const* M_mesh = nullptr;
    index_type M_id{};
    local_index_type M_localIndex{};
};

/**
 * @class GeometricEntityView
 * @brief Mathematical edge, face, or cell view over compact mesh storage.
 *
 * Point traversal uses the entity's compact canonical connectivity and the
 * mesh coordinate owner. The compatibility entity remains available during
 * migration for topology and metadata that have not yet moved to compact
 * sidecars.
 *
 * @tparam MeshType owning mesh type
 * @tparam EntityType compatibility entity type
 */
template <typename MeshType, typename EntityType>
class GeometricEntityView
{
  public:
    /** Owning mesh type. */
    using mesh_type = std::remove_cvref_t<MeshType>;
    /** Compatibility entity type. */
    using entity_type = std::remove_cvref_t<EntityType>;
    /** Point view type. */
    using point_view_type = GeometricPointView<mesh_type>;
    /** Coordinate scalar type. */
    using value_type = typename mesh_type::value_type;
    /** External identifier type. */
    using index_type = typename mesh_type::index_type;
    /** Dense local point index type. */
    using local_index_type = uint32_type;
    /** Allocating convenience coordinate-matrix type. */
    using coordinate_matrix_type =
        Eigen::Matrix<value_type, mesh_type::nRealDim,
                      Eigen::Dynamic, Eigen::ColMajor>;

    /** Topological dimension of the viewed entity. */
    static constexpr uint16_type nDim = entity_type::nDim;
    /** Ambient physical dimension. */
    static constexpr uint16_type nRealDim = entity_type::nRealDim;

    /** Construct an empty, non-dereferenceable view. */
    GeometricEntityView() = default;

    /**
     * @brief Bind one compatibility entity to its owning mesh.
     * @param mesh owning mesh
     * @param entity viewed entity
     */
    GeometricEntityView(
        mesh_type const& mesh, entity_type const& entity )
        : M_mesh( std::addressof( mesh ) ),
          M_entity( std::addressof( entity ) )
    {
    }

    /** @brief Return whether this view is bound to a mesh entity. */
    [[nodiscard]] explicit operator bool() const noexcept
    {
        return M_mesh != nullptr && M_entity != nullptr;
    }

    /** @brief Return the external entity identifier. */
    [[nodiscard]] index_type id() const
    {
        return this->entity().id();
    }

    /** @brief Return the number of geometric points. */
    [[nodiscard]] uint16_type pointCount() const
    {
        return static_cast<uint16_type>(
            this->pointIds().size() );
    }

    /** @brief Return point IDs in canonical local mathematical order. */
    [[nodiscard]] std::span<index_type const> pointIds() const
    {
        return this->entity().pointIds();
    }

    /**
     * @brief Return a lazy random-access range of dense local point indices.
     *
     * No connectivity allocation is performed. Dense external IDs resolve by
     * arithmetic; sparse imported or partition-local IDs use the coordinate
     * owner's fallback lookup.
     */
    [[nodiscard]] auto pointLocalIndices() const
    {
        auto const* mesh = M_mesh;
        return this->pointIds() |
               std::views::transform(
                   [mesh]( index_type id )
                   {
                       return mesh->pointLocalIndex( id );
                   } );
    }

    /**
     * @brief Return one point view in canonical local order.
     * @param localPoint local point index
     */
    [[nodiscard]] point_view_type
    point( uint16_type localPoint ) const
    {
        auto const ids = this->pointIds();
        if ( localPoint >= ids.size() )
            throw std::out_of_range(
                "invalid local geometric point index" );
        return point_view_type( this->mesh(), ids[localPoint] );
    }

    /**
     * @brief Gather all geometric coordinates into caller-owned storage.
     *
     * The destination may be an Eigen or uBLAS matrix supporting
     * `operator()(row,column)`.
     *
     * @tparam Matrix writable matrix type
     * @param destination matrix with `nRealDim` rows and `pointCount` columns
     */
    template <typename Matrix>
        requires requires( Matrix& matrix, uint16_type component,
                           uint16_type point, value_type value )
        {
            matrix( component, point ) = value;
        }
    void updateCoordinates( Matrix& destination ) const
    {
        this->gatherCoordinates(
            destination, this->pointCount() );
    }

    /**
     * @brief Gather vertex coordinates into caller-owned storage.
     * @tparam Matrix writable matrix type
     * @param destination matrix with `nRealDim` rows and `numVertices` columns
     */
    template <typename Matrix>
        requires requires( Matrix& matrix, uint16_type component,
                           uint16_type point, value_type value )
        {
            matrix( component, point ) = value;
        }
    void updateVertexCoordinates( Matrix& destination ) const
    {
        this->gatherCoordinates(
            destination, entity_type::numVertices );
    }

    /** @brief Return an allocating convenience matrix of geometric points. */
    [[nodiscard]] coordinate_matrix_type coordinates() const
    {
        coordinate_matrix_type result(
            nRealDim, this->pointCount() );
        this->updateCoordinates( result );
        return result;
    }

    /** @brief Return the compatibility entity facade. */
    [[nodiscard]] entity_type const& entity() const
    {
        this->requireBound();
        return *M_entity;
    }

    /** @brief Return the owning mesh. */
    [[nodiscard]] mesh_type const& mesh() const
    {
        this->requireBound();
        return *M_mesh;
    }

    /** @brief Return whether the entity is on the mesh boundary. */
    [[nodiscard]] bool isOnBoundary() const
    {
        return this->entity().isOnBoundary();
    }

    /** @brief Return the owning process identifier. */
    [[nodiscard]] rank_type processId() const
    {
        return this->entity().processId();
    }

    /** @brief Return all mathematical marker components. */
    [[nodiscard]] decltype(auto) markers() const
    {
        return this->entity().markers();
    }

    /**
     * @brief Return one mathematical marker component.
     * @param component marker component
     */
    [[nodiscard]] decltype(auto)
    marker( uint16_type component ) const
    {
        return this->entity().marker( component );
    }

    /** @brief Return the cached entity measure when supported. */
    [[nodiscard]] value_type measure() const
        requires requires( entity_type const& entity )
        {
            { entity.measure() } -> std::convertible_to<value_type>;
        }
    {
        return this->entity().measure();
    }

    /** @brief Return the entity diameter when supported. */
    [[nodiscard]] value_type h() const
        requires requires( entity_type const& entity )
        {
            { entity.h() } -> std::convertible_to<value_type>;
        }
    {
        return this->entity().h();
    }

    /**
     * @brief Return one subentity permutation.
     * @tparam Args permutation-selector argument types
     * @param args permutation-selector arguments
     */
    template <typename... Args>
    [[nodiscard]] decltype(auto)
    permutation( Args&&... args ) const
        requires requires( entity_type const& entity )
        {
            entity.permutation(
                std::forward<Args>( args )... );
        }
    {
        return this->entity().permutation(
            std::forward<Args>( args )... );
    }

    /**
     * @brief Return one face as a mathematical view.
     * @param localFace local face number
     */
    [[nodiscard]] auto face( uint16_type localFace ) const
        requires requires( entity_type const& entity,
                           mesh_type const& mesh )
        {
            entity.face( localFace ).id();
            mesh.faceView( entity.face( localFace ).id() );
        }
    {
        auto const id = this->entity().face( localFace ).id();
        return this->mesh().faceView( id );
    }

    /**
     * @brief Return one edge as a mathematical view.
     * @param localEdge local edge number
     */
    [[nodiscard]] auto edge( uint16_type localEdge ) const
        requires requires( entity_type const& entity,
                           mesh_type const& mesh )
        {
            entity.edge( localEdge ).id();
            mesh.edgeView( entity.edge( localEdge ).id() );
        }
    {
        auto const id = this->entity().edge( localEdge ).id();
        return this->mesh().edgeView( id );
    }

    /**
     * @brief Return one adjacent cell when it exists.
     * @param localFace local face number
     */
    [[nodiscard]] auto neighbor( uint16_type localFace ) const
        requires requires( entity_type const& entity,
                           mesh_type const& mesh )
        {
            { entity.neighbor( localFace ) } ->
                std::convertible_to<index_type>;
            mesh.hasElement( entity.neighbor( localFace ) );
            mesh.element( entity.neighbor( localFace ) );
        }
    {
        using cell_view_type =
            typename mesh_type::element_view_type;
        auto const neighborId =
            static_cast<index_type>(
                this->entity().neighbor( localFace ) );
        if ( !this->mesh().hasElement( neighborId ) )
            return std::optional<cell_view_type>{};
        return std::optional<cell_view_type>{
            this->mesh().elementView( neighborId )
        };
    }

    /** @brief Compare view handles by owner and compatibility entity. */
    [[nodiscard]] friend bool operator==(
        GeometricEntityView const&,
        GeometricEntityView const& ) = default;

  private:
    /**
     * @brief Gather a canonical prefix of entity coordinates.
     *
     * Eigen destinations use column expressions so fixed ambient dimensions
     * vectorize. Other matrix APIs, including uBLAS, use the common scalar
     * fallback.
     *
     * @tparam Matrix writable Eigen- or uBLAS-compatible matrix type
     * @param destination caller-owned coordinate matrix
     * @param pointCount number of canonical local points to gather
     */
    template <typename Matrix>
    void gatherCoordinates(
        Matrix& destination, uint16_type pointCount ) const
    {
        auto const rows = [&destination]()
        {
            if constexpr ( requires { destination.rows(); } )
                return static_cast<size_type>(
                    destination.rows() );
            else
                return static_cast<size_type>(
                    destination.size1() );
        }();
        auto const columns = [&destination]()
        {
            if constexpr ( requires { destination.cols(); } )
                return static_cast<size_type>(
                    destination.cols() );
            else
                return static_cast<size_type>(
                    destination.size2() );
        }();
        if ( rows != nRealDim || columns < pointCount )
            throw std::invalid_argument(
                "incompatible geometric coordinate destination" );

        auto const allCoordinates =
            this->mesh().pointCoordinates();
        auto const ids = this->pointIds();
        if constexpr (
            requires( Matrix& matrix,
                      decltype( allCoordinates ) coordinates )
            {
                matrix.col( Eigen::Index{} ) =
                    coordinates.col( Eigen::Index{} );
            } )
        {
            for ( uint16_type localPoint = 0;
                  localPoint < pointCount; ++localPoint )
            {
                auto const pointIndex =
                    this->mesh().pointLocalIndex(
                        ids[localPoint] );
                destination.col( localPoint ) =
                    allCoordinates.col( pointIndex );
            }
        }
        else
        {
            for ( uint16_type localPoint = 0;
                  localPoint < pointCount; ++localPoint )
            {
                auto const pointIndex =
                    this->mesh().pointLocalIndex(
                        ids[localPoint] );
                for ( uint16_type component = 0;
                      component < nRealDim; ++component )
                    destination( component, localPoint ) =
                        allCoordinates(
                            component, pointIndex );
            }
        }
    }

    void requireBound() const
    {
        if ( !M_mesh || !M_entity )
            throw std::logic_error(
                "unbound geometric entity view" );
    }

    mesh_type const* M_mesh = nullptr;
    entity_type const* M_entity = nullptr;
};

/**
 * @brief Construct a point view.
 * @tparam MeshType owning mesh type
 * @param mesh owning mesh
 * @param id external point identifier
 */
template <typename MeshType>
[[nodiscard]] auto
pointView( MeshType const& mesh,
           typename MeshType::index_type id )
{
    return GeometricPointView<MeshType>( mesh, id );
}

/**
 * @brief Construct an edge, face, or cell view.
 * @tparam MeshType owning mesh type
 * @tparam EntityType compatibility entity type
 * @param mesh owning mesh
 * @param entity viewed entity
 */
template <typename MeshType, typename EntityType>
[[nodiscard]] auto
entityView( MeshType const& mesh, EntityType const& entity )
{
    return GeometricEntityView<MeshType, EntityType>(
        mesh, entity );
}

/**
 * @brief Return the lazy mathematical point-view range of a mesh.
 * @tparam MeshLike mesh object, raw pointer, or smart pointer type
 * @param mesh mesh handle
 */
template <typename MeshLike>
[[nodiscard]] auto pointViews( MeshLike const& mesh )
    requires(
        requires { mesh.pointViews(); } ||
        requires { mesh->pointViews(); } )
{
    if constexpr ( requires { mesh.pointViews(); } )
        return mesh.pointViews();
    else
        return mesh->pointViews();
}

/**
 * @brief Return the lazy mathematical edge-view range of a mesh.
 * @tparam MeshLike mesh object, raw pointer, or smart pointer type
 * @param mesh mesh handle
 */
template <typename MeshLike>
[[nodiscard]] auto edgeViews( MeshLike const& mesh )
    requires(
        requires { mesh.edgeViews(); } ||
        requires { mesh->edgeViews(); } )
{
    if constexpr ( requires { mesh.edgeViews(); } )
        return mesh.edgeViews();
    else
        return mesh->edgeViews();
}

/**
 * @brief Return the lazy mathematical face-view range of a mesh.
 * @tparam MeshLike mesh object, raw pointer, or smart pointer type
 * @param mesh mesh handle
 */
template <typename MeshLike>
[[nodiscard]] auto faceViews( MeshLike const& mesh )
    requires(
        requires { mesh.faceViews(); } ||
        requires { mesh->faceViews(); } )
{
    if constexpr ( requires { mesh.faceViews(); } )
        return mesh.faceViews();
    else
        return mesh->faceViews();
}

/**
 * @brief Return the lazy mathematical cell-view range of a mesh.
 * @tparam MeshLike mesh object, raw pointer, or smart pointer type
 * @param mesh mesh handle
 */
template <typename MeshLike>
[[nodiscard]] auto elementViews( MeshLike const& mesh )
    requires(
        requires { mesh.elementViews(); } ||
        requires { mesh->elementViews(); } )
{
    if constexpr ( requires { mesh.elementViews(); } )
        return mesh.elementViews();
    else
        return mesh->elementViews();
}

} // namespace Feel

#endif // FEELPP_FEELMESH_ENTITYVIEW_HPP
