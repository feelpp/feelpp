/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2001-2004 INRIA
 * SPDX-FileCopyrightText: 2008 Université Joseph Fourier
 * SPDX-FileCopyrightText: 2001-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
  @file geond.hpp
  @brief Multi-dimensional geometric entities with C++20/23 and Dynamic order support
*/
#ifndef _GEOND_HH_
#define _GEOND_HH_

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/serialization/version.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <concepts>
#include <array>
#include <memory>
#include <type_traits>
#include <utility>
#include <variant>

#include <feel/feelmesh/geo0d.hpp>
#include <feel/feelmesh/geoentity.hpp>
#include <feel/feelmesh/marker.hpp>
#include <feel/feelmesh/meshbase.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/order.hpp>

namespace Feel
{
template <typename IndexT>
class MeshBase;

template <int Dim, int Order, int RealDim, template <int, int, int> class Entity, typename T>
struct GT_Lagrange;

template <class Convex, uint16_type O, typename T2>
class Gauss;
template <int IMORDER,
          int DIM,
          template <int, int, int> class Entity,
          template <class Convex, uint16_type O, typename T2> class QPS,
          typename T>
struct IMGeneric;

/// \cond detail
namespace detail
{
/**
 * \class ReversePoint
 * An utility to  invert point numbering on a GeoShape
 */
template <typename GeoShape>
class ReversePoint
{
  public:
    uint16_type operate( uint16_type const& point )
    {
        return ( point < GeoShape::numVertices ) ? GeoShape::numVertices - point : GeoShape::numPoints - point + GeoShape::numVertices;
    }
};
} // namespace detail
/// \endcond

/**
 * @class GeoNDCommon
 * @ingroup Mesh
 * @brief Common data shared in a collection of multi-dimensional geometrical entity.
 */
template <typename GeoNDType>
class GeoNDCommon
{
  public:
    typedef GeoNDType geond_type;
    typedef typename geond_type::gm_ptrtype gm_ptrtype;
    typedef typename geond_type::gm1_ptrtype gm1_ptrtype;
    typedef typename geond_type::index_type index_type;
    using point_type = typename geond_type::point_type;

    GeoNDCommon()
        : M_mesh( nullptr )
    {
    }
    explicit GeoNDCommon( MeshBase<index_type> const* mesh )
        : M_mesh( mesh )
    {
    }
    GeoNDCommon( MeshBase<index_type> const* mesh, gm_ptrtype const& gm, gm1_ptrtype const& gm1 )
        : M_mesh( mesh ),
          M_gm( gm ),
          M_gm1( gm1 )
    {
    }
    GeoNDCommon( GeoNDCommon const& e ) = default;
    GeoNDCommon( GeoNDCommon&& e ) = default;

    MeshBase<index_type> const* mesh() const { return M_mesh; }
    gm_ptrtype const& gm() const { return M_gm; }
    gm1_ptrtype const& gm1() const { return M_gm1; }

    void setMesh( MeshBase<index_type> const* mesh ) { M_mesh = mesh; }
    void setGm( gm_ptrtype const& gm, gm1_ptrtype const& gm1 )
    {
        M_gm = gm;
        M_gm1 = gm1;
    }

    /**
     * @brief Bind compact point identifiers to the mesh point facade.
     *
     * One resolver is shared by every entity of the same dimension. Stored
     * entities therefore do not need their own point-pointer vectors.
     *
     * @tparam MeshType concrete mesh type
     * @param mesh owning mesh
     */
    template <typename MeshType>
        requires requires( MeshType& mesh, index_type id )
        {
            static_cast<point_type&>( mesh.point( id ) );
            static_cast<point_type const&>(
                std::as_const( mesh ).point( id ) );
        }
    void bindPointResolver( MeshType& mesh ) noexcept
    {
        M_pointOwner = std::addressof( mesh );
        M_mutablePointResolver =
            []( void* owner, index_type id ) -> point_type&
            {
                return static_cast<point_type&>(
                    static_cast<MeshType*>( owner )->point( id ) );
            };
        M_constPointResolver =
            []( void const* owner, index_type id ) -> point_type const&
            {
                return static_cast<point_type const&>(
                    static_cast<MeshType const*>( owner )->point( id ) );
            };
    }

    /** @brief Return whether compact point identifiers can be resolved. */
    [[nodiscard]] bool hasPointResolver() const noexcept
    {
        return M_pointOwner && M_mutablePointResolver &&
               M_constPointResolver;
    }

    /**
     * @brief Resolve a mutable point facade from its stable identifier.
     * @param id external point identifier
     */
    [[nodiscard]] point_type& point( index_type id )
    {
        CHECK( this->hasPointResolver() )
            << "GeoND common point resolver is not bound";
        return M_mutablePointResolver( M_pointOwner, id );
    }

    /**
     * @brief Resolve an immutable point facade from its stable identifier.
     * @param id external point identifier
     */
    [[nodiscard]] point_type const& point( index_type id ) const
    {
        CHECK( this->hasPointResolver() )
            << "GeoND common point resolver is not bound";
        return M_constPointResolver( M_pointOwner, id );
    }

  private:
    using mutable_point_resolver_type =
        point_type& ( * )( void*, index_type );
    using const_point_resolver_type =
        point_type const& ( * )( void const*, index_type );

    // mesh to which the geond element belongs to
    MeshBase<index_type> const* M_mesh;
    // geometric mapping
    gm_ptrtype M_gm;
    gm1_ptrtype M_gm1;
    /** Mesh-wide point owner used by the compatibility facade. */
    void* M_pointOwner = nullptr;
    /** Mutable point resolver shared by this entity dimension. */
    mutable_point_resolver_type M_mutablePointResolver = nullptr;
    /** Immutable point resolver shared by this entity dimension. */
    const_point_resolver_type M_constPointResolver = nullptr;
};

template< typename T>
class GeoNDMeasuresStorage
{
  public:
    /**
     * @brief Construct the measure-storage policy.
     *
     * Cached values now live in the mesh-owned sparse metadata sidecar; this
     * policy remains only as a compile-time capability tag.
     */
    explicit GeoNDMeasuresStorage( uint16_type ) noexcept {}
};

template< typename T>
class GeoNDNoMeasuresStorage
{
  public:
    /** @brief Construct the no-measure-storage policy. */
    explicit GeoNDNoMeasuresStorage( uint16_type ) noexcept {}
};

/**
 * @class GeoND
 * @ingroup Mesh
 * @brief Base class for Multi-dimensional basis Geometrical Entities.
 *
 * @note C++20/23 modernized with Dynamic order support
 */
template <uint16_type Dim,
          typename GEOSHAPE,
          typename T = double,
          typename IndexT = uint32_type,
          typename POINTTYPE = Geo0D<Dim, T, IndexT>,
          bool UseMeasuresStorage = false >
class GeoND
    : public GeoEntity<GEOSHAPE, T, IndexT>,
      public std::conditional_t<UseMeasuresStorage, GeoNDMeasuresStorage<T>, GeoNDNoMeasuresStorage<T>>
{
    using super = GeoEntity<GEOSHAPE, T, IndexT>;
    using super2 = std::conditional_t<UseMeasuresStorage, GeoNDMeasuresStorage<T>, GeoNDNoMeasuresStorage<T>>;

  public:
    typedef T value_type;
    using index_type = IndexT;
    using size_type = index_type;
    typedef GeoND<Dim, GEOSHAPE, value_type, index_type, POINTTYPE, UseMeasuresStorage> self_type;
    typedef self_type element_type;

    typedef GEOSHAPE GeoShape;
    typedef POINTTYPE PointType;

    //typedef PointType point_type;
    typedef typename super::face_type face_type;
    using point_type = PointType;

    /**
     * @class PointPointerCache
     * @brief Lazy compatibility range for detached point pointers.
     *
     * The cache owns one contiguous pointer array only while an entity is
     * detached or a caller explicitly requests the legacy `points()` range.
     * Stored entities keep the cache empty and resolve compact connectivity
     * through `GeoNDCommon`.
     */
    class PointPointerCache
    {
      public:
        using iterator = point_type**;
        using const_iterator = point_type* const*;

        PointPointerCache() = default;

        /** @brief Construct and materialize a cache of `size` null pointers. */
        explicit PointPointerCache( uint16_type size )
            : M_size( size )
        {
            this->materialize();
        }

        /** @brief Deep-copy a materialized cache, or preserve emptiness. */
        PointPointerCache( PointPointerCache const& other )
            : M_size( other.M_size )
        {
            if ( other.M_data )
            {
                this->materialize();
                std::copy_n(
                    other.M_data.get(), M_size,
                    M_data.get() );
            }
        }

        PointPointerCache( PointPointerCache&& ) noexcept = default;

        PointPointerCache& operator=( PointPointerCache const& other )
        {
            if ( this == std::addressof( other ) )
                return *this;
            M_size = other.M_size;
            M_data.reset();
            if ( other.M_data )
            {
                this->materialize();
                std::copy_n(
                    other.M_data.get(), M_size,
                    M_data.get() );
            }
            return *this;
        }

        PointPointerCache&
        operator=( PointPointerCache&& ) noexcept = default;

        /** @brief Set the mathematical cardinality and clear old storage. */
        void resize( uint16_type size )
        {
            if ( size == M_size )
                return;
            M_size = size;
            M_data.reset();
        }

        /** @brief Release the compatibility pointer payload. */
        void reset() noexcept
        {
            M_data.reset();
        }

        /** @brief Return whether the pointer payload is materialized. */
        [[nodiscard]] explicit operator bool() const noexcept
        {
            return static_cast<bool>( M_data );
        }

        /** @brief Return the mathematical pointer count. */
        [[nodiscard]] uint16_type size() const noexcept
        {
            return M_size;
        }

        /** @brief Return retained pointer capacity. */
        [[nodiscard]] uint16_type capacity() const noexcept
        {
            return M_data ? M_size : 0;
        }

        /** @brief Return one mutable point pointer. */
        [[nodiscard]] point_type*& operator[]( size_type index )
        {
            this->materialize();
            return M_data[index];
        }

        /** @brief Return one immutable point pointer. */
        [[nodiscard]] point_type* const&
        operator[]( size_type index ) const
        {
            this->materialize();
            return M_data[index];
        }

        [[nodiscard]] iterator begin()
        {
            this->materialize();
            return M_data.get();
        }
        [[nodiscard]] const_iterator begin() const
        {
            this->materialize();
            return M_data.get();
        }
        [[nodiscard]] iterator end()
        {
            auto first = this->begin();
            return M_size == 0 ? first : first + M_size;
        }
        [[nodiscard]] const_iterator end() const
        {
            auto first = this->begin();
            return M_size == 0 ? first : first + M_size;
        }

        /** @brief Allocate and null-initialize the pointer payload. */
        void materialize() const
        {
            if ( M_data || M_size == 0 )
                return;
            M_data = std::make_unique<point_type*[]>( M_size );
            std::fill_n( M_data.get(), M_size, nullptr );
        }

      private:
        mutable std::unique_ptr<point_type*[]> M_data;
        uint16_type M_size = 0;
    };

    static inline const size_type Shape = super::Shape;
    static inline const uint16_type numPoints = super::numPoints;
    static inline const uint16_type numVertices = super::numVertices;
    static inline const uint16_type numLocalPoints = super::numPoints;
    static inline const uint16_type numLocalEdges = super::numEdges;
    static inline const uint16_type numLocalVertices = super::numVertices;
    static inline const int numFaces = super::numFaces;
    static inline const int numEdges = super::numEdges;
    static inline const int numTopologicalFaces = super::numTopologicalFaces;
    static inline const uint16_type numNeighbors = super::numTopologicalFaces;

    using point_pointer_container_type = PointPointerCache;
    using point_iterator =
        typename point_pointer_container_type::iterator;
    using point_const_iterator =
        typename point_pointer_container_type::const_iterator;

    typedef typename matrix_node<value_type>::type matrix_node_type;
    typedef typename node<value_type>::type node_type;

    static constexpr uint16_type nDim = super::nDim;
    static constexpr uint16_type nRealDim = super::nRealDim;

    //! @brief True if Order is known at compile time (i.e., Order >= 0)
    static constexpr bool is_order_static = GEOSHAPE::is_order_static;
    //! @brief True if Order is determined at runtime (i.e., Order == Dynamic)
    static constexpr bool is_order_dynamic = GEOSHAPE::is_order_dynamic;
    //! @brief Template order parameter value (may be Dynamic = -1)
    static constexpr int nOrder_v = GEOSHAPE::nOrder_v;
    //! @brief Static order (or 1 as placeholder for dynamic case)
    static constexpr uint16_type nOrder = super::nOrder;

    //-------------------------------------------------------------------------
    // Geometric mapping helper - C++20 modernized
    //-------------------------------------------------------------------------
    template <int GmOrder>
    struct GetGm
    {
        using type = std::conditional_t<GeoShape::is_hypercube,
                                        GT_Lagrange<nDim, GmOrder, nRealDim, Hypercube, T>,
                                        GT_Lagrange<nDim, GmOrder, nRealDim, Simplex, T>>;
        using ptrtype = std::shared_ptr<type>;
    };

    using gm_type = typename GetGm<nOrder>::type;
    using gm_ptrtype = typename GetGm<nOrder>::ptrtype;

    using gm1_type = typename GetGm<1>::type;
    using gm1_ptrtype = typename GetGm<1>::ptrtype;

    using reference_convex_type = typename gm_type::super::reference_convex_type;
    using reference_convex1_type = typename gm1_type::super::reference_convex_type;

    using vertex_permutation_type = typename super::vertex_permutation_type;
    using edge_permutation_type = typename super::edge_permutation_type;
    using face_permutation_type = typename super::face_permutation_type;

    //! @brief Permutation type based on dimension (C++20 style)
    using permutation_type = std::conditional_t<
        nDim == 1,
        vertex_permutation_type,
        std::conditional_t<nDim == 2, edge_permutation_type, face_permutation_type>>;

    //! @brief Quadrature order for measure computation (static value for static order)
    static constexpr uint16_type meas_quad_order = ( nOrder >= 1 ) ? ( nOrder - 1 ) * nDim : 0;

    //-------------------------------------------------------------------------
    // Integration method helper - C++20 modernized
    //-------------------------------------------------------------------------
    template <int GeoOrder>
    struct GetImMeasure
    {
        static constexpr uint16_type quad_order = ( GeoOrder >= 1 ) ? ( GeoOrder - 1 ) * nDim : 0;
        using type = std::conditional_t<GeoShape::is_hypercube,
                                        typename IMGeneric<quad_order, Dim, Hypercube, Gauss, value_type>::type,
                                        typename IMGeneric<quad_order, Dim, Simplex, Gauss, value_type>::type>;
    };

    using quad_meas_type = typename GetImMeasure<nOrder>::type;
    using quad_meas1_type = typename GetImMeasure<1>::type;

    //-------------------------------------------------------------------------
    // Runtime order support
    //-------------------------------------------------------------------------

    /**
     * @brief Get the polynomial order of the geometric entity.
     * @return The order (compile-time for static order, runtime for dynamic).
     */
    [[nodiscard]] constexpr uint16_type order() const noexcept
    {
        if constexpr ( is_order_static )
            return nOrder;
        else
            return M_runtime_order;
    }

    /**
     * @brief Get the quadrature order for measure computation.
     * @return The quadrature order (accounts for runtime order if dynamic).
     */
    [[nodiscard]] constexpr uint16_type measureQuadOrder() const noexcept
    {
        if constexpr ( is_order_static )
            return meas_quad_order;
        else
            return ( M_runtime_order >= 1 ) ? ( M_runtime_order - 1 ) * nDim : 0;
    }

    /**
     * @brief Get the number of points in the element (runtime-aware).
     * @return The number of points.
     */
    [[nodiscard]] uint16_type nPointsRuntime() const noexcept
    {
        if constexpr ( is_order_static )
            return numPoints;
        else
        {
            if constexpr ( GeoShape::is_simplex )
                return static_cast<uint16_type>( Feel::detail::simplexTotal( nDim, M_runtime_order ) );
            else
                return static_cast<uint16_type>( Feel::detail::hypercubeTotal( nDim, M_runtime_order ) );
        }
    }


    /**
     * default constructor (static order only)
     */
    GeoND() requires( is_order_static )
        : super( 0 ),
          super2( numTopologicalFaces ),
          M_pointPointers( numPoints ),
          M_commonData( nullptr )
    {
    }

    /**
     * constructor with runtime order (dynamic order only)
     *
     * @param runtime_order the polynomial order specified at runtime
     */
    explicit GeoND( RuntimeOrder runtime_order ) requires( is_order_dynamic )
        : super( 0 ),
          super2( numTopologicalFaces ),
          M_pointPointers(
              computeNumPoints( runtime_order.value ) ),
          M_commonData( nullptr ),
          M_runtime_order( runtime_order.value )
    {
    }

  private:
    //! @brief Helper to compute number of points for dynamic order (used in initializer list)
    static constexpr size_type computeNumPoints( uint16_type order ) noexcept
    {
        if constexpr ( GeoShape::is_simplex )
            return Feel::detail::simplexTotal( nDim, order );
        else
            return Feel::detail::hypercubeTotal( nDim, order );
    }

  public:

    /**
     * constructor from an id (static order only)
     *
     * @param id identifier for the element to store
     *
     */
    explicit GeoND( size_type id ) requires( is_order_static )
        : super( id ),
          super2( numTopologicalFaces ),
          M_pointPointers( numPoints ),
          M_commonData( nullptr )
    {
    }

    /**
     * constructor from an id with runtime order (dynamic order only)
     *
     * @param id identifier for the element to store
     * @param runtime_order the polynomial order specified at runtime
     */
    GeoND( size_type id, RuntimeOrder runtime_order ) requires( is_order_dynamic )
        : super( id ),
          super2( numTopologicalFaces ),
          M_pointPointers(
              computeNumPoints( runtime_order.value ) ),
          M_commonData( nullptr ),
          M_runtime_order( runtime_order.value )
    {
    }

    GeoND( GeoND const& e )
        : super( e ),
          super2( e ),
          M_pointPointers( e.M_pointPointers ),
          M_neighbors( e.M_neighbors ),
          M_neighborsInitialized( e.M_neighborsInitialized ),
          M_commonData( e.M_commonData ),
          M_runtime_order( e.M_runtime_order )
    {
    }
    GeoND( GeoND&& e ) noexcept
        : super( std::move( e ) ),
          super2( std::move( e ) ),
          M_pointPointers( std::move( e.M_pointPointers ) ),
          M_neighbors( std::move( e.M_neighbors ) ),
          M_neighborsInitialized( e.M_neighborsInitialized ),
          M_commonData( std::move( e.M_commonData ) ),
          M_runtime_order( std::move( e.M_runtime_order ) )
    {
    }

    GeoND& operator=( GeoND const& e )
    {
        if ( this == std::addressof( e ) )
            return *this;
        super::operator=( e );
        super2::operator=( e );
        M_pointPointers = e.M_pointPointers;
        M_neighbors = e.M_neighbors;
        M_neighborsInitialized = e.M_neighborsInitialized;
        M_commonData = e.M_commonData;
        if constexpr ( is_order_dynamic )
            M_runtime_order = e.M_runtime_order;
        return *this;
    }
    GeoND& operator=( GeoND&& e )
    {
        super::operator=( std::move( e ) );
        super2::operator=( std::move( e ) );
        M_pointPointers = std::move( e.M_pointPointers );
        M_neighbors = std::move( e.M_neighbors );
        M_neighborsInitialized = e.M_neighborsInitialized;
        M_commonData = std::move( e.M_commonData );
        if constexpr ( is_order_dynamic )
            M_runtime_order = std::move( e.M_runtime_order );
        return *this;
    }

    /**
     * destructor, make it virtual for derived classes
     */
    ~GeoND() override
    {
    }

    /**
     * set common data
     */
    void setCommonData( GeoNDCommon<self_type> * commonData )
    {
        M_commonData = commonData;
        if ( this->canResolveCompactPoints() )
            M_pointPointers.reset();
    }

    //! return the geometric mapping if a mesh was set
    gm_ptrtype gm() const
    {
        if ( M_commonData != nullptr )
            return M_commonData->gm();
        return gm_ptrtype();
    }

    //! return the geometric mapping if a mesh was set
    gm1_ptrtype gm1() const
    {
        if ( M_commonData != nullptr )
            return M_commonData->gm1();
        return gm1_ptrtype();
    }

    /**
     * \return the mesh to which this geometric entity belongs to
     */
    MeshBase<> const* mesh() const
    {
        if ( M_commonData != nullptr )
            return M_commonData->mesh();
        return nullptr;
    }

    /**
     * \return true if points have been inserted in elements, false
     * otherwise
     */
    bool hasPoints() const
    {
        if ( this->canResolveCompactPoints() )
            return true;
        auto const& pointers = this->pointPointers();
        return std::ranges::find( pointers, nullptr ) ==
               pointers.end();
    }
    //! @return true if the element has at least a point with marker1 active
    bool hasPointWithMarker() const
    {
        for ( uint16_type point = 0;
              point < this->nPointsRuntime(); ++point )
            if ( this->point( point ).hasMarker() )
                return true;
        return false;
    }

    /**
     * \return the number of points in convex
     */
    uint16_type nPoints() const
    {
        return numPoints;
    }

    /**
     * the number of neighbors is equal to the number of
     * faces. Neighbor identifiers are stored inline and invalid identifiers
     * denote boundary faces.
     *
     *\return the number of neighbors
     */
    uint16_type nNeighbors() const
    {
        return M_neighborsInitialized ? numNeighbors : 0;
    }

    /**
     * @brief Return the neighbor identifier across one local face.
     *
     * @param n local face number
     * @return neighbor element identifier, or the invalid identifier for a
     * boundary face
     */
    size_type neighbor( uint16_type n ) const
    {
        CHECK( M_neighborsInitialized )
            << "neighbor connectivity has not been initialized";
        CHECK_LT( n, numNeighbors )
            << "invalid local neighbor index";
        return M_neighbors[n];
    }

    /**
     * @brief Test whether a valid neighbor is registered across a local face.
     * @param n local face number
     * @return true when a valid neighboring element identifier is present
     */
    bool hasNeighbor( uint16_type n ) const
    {
        if ( !M_neighborsInitialized || n >= numNeighbors )
            return false;
        return M_neighbors[n] != invalid_v<size_type>;
    }

    /**
     * @brief Set the neighbor identifier across one local face.
     * @param n local face number
     * @param neigh_id neighboring element identifier
     */
    void setNeighbor( uint16_type n, size_type neigh_id )
    {
        CHECK_LT( n, numNeighbors )
            << "invalid local neighbor index";
        if ( !M_neighborsInitialized )
        {
            M_neighbors.fill( invalid_v<size_type> );
            M_neighborsInitialized = true;
        }
        M_neighbors[n] = neigh_id;
    }

    bool isNeighbor( self_type const& G ) const
    {
        for ( uint16_type i = 0; i < this->nNeighbors(); ++i )
            if ( this->neighbor( i ) == G.id() ) return true;

        return false;
    }

    /**
     * \return the barycenter of the element
     */
    node_type barycenter() const
    {
        auto M = glas::average( this->G() );
        return ublas::column( M, 0 );
    }

    /**
     * \return the barycenter at the faces of the element
     */
    node_type faceBarycenter( uint16_type f ) const
    {
        constexpr int nPtsInFace = GEOSHAPE::topological_face_type::numPoints;

        node_type n( nRealDim );
        em_node_type<value_type> en( n.data().begin(), n.size() );
        en.setZero();
        for ( uint16_type p = 0; p < nPtsInFace; ++p )
        {
            // get pt id in element  from local pt id  in face
            int ptid = this->fToP( f, p );
            auto const& node = this->point( ptid ).node();
            em_node_type<value_type> emnode( const_cast<typename point_type::value_type /*double*/*>( node.data().begin() ), node.size() );
            en += emnode;
        }
        en /= nPtsInFace;
        return n;
    }

    /**
     * \return the barycenters at the faces of the element
     */
    FEELPP_DEPRECATED matrix_node_type faceBarycenters() const
    {
        matrix_node_type n;
        return n;
    }

    /**
     * \return permutation
     */
    permutation_type permutation( uint16_type /*f*/ ) const
    {
        return permutation_type();
    }

    //! \return identity permutation of entity co-dim 1 (topological faces)
    virtual typename super::template PermutationSubEntity<1> permutation( uint16_type /*i*/, mpl::int_<1> ) const
    {
        return typename super::template PermutationSubEntity<1>();
    }

    //!
    //! @return points set with constness
    //!
    point_pointer_container_type const& points() const
    {
        return this->pointPointers();
    }

    /**
     * @brief Return mutable compatibility point pointers.
     *
     * Pointer replacement through this container must be followed by
     * `synchronizePointConnectivity()` before index-based consumers run.
     *
     * @return mutable point-pointer container
     */
    point_pointer_container_type& points()
    {
        return this->pointPointers();
    }

    /** @brief Return whether the legacy point-pointer cache is materialized. */
    [[nodiscard]] bool hasPointPointerCache() const noexcept
    {
        return static_cast<bool>( M_pointPointers );
    }

    /** @brief Return retained capacity bytes of the legacy pointer cache. */
    [[nodiscard]] size_type pointPointerCacheCapacityBytes() const noexcept
    {
        return M_pointPointers
                   ? static_cast<size_type>(
                         M_pointPointers.capacity() *
                         sizeof( point_type* ) )
                   : 0;
    }

    /**
     * @brief Return compact point identifiers in local mathematical order.
     *
     * The returned span is owned by the mesh entity storage after insertion.
     * Detached entities expose their lazy construction buffer.
     */
    [[nodiscard]] std::span<index_type const> pointIds() const
    {
        return this->connectivityIds();
    }

    /**
     * @brief Rebuild compact point connectivity from compatibility pointers.
     *
     * This is used at importer/container boundaries while pointer-based
     * consumers are migrated to the compact view.
     */
    void synchronizePointConnectivity()
    {
        if ( !M_pointPointers )
            return;
        std::vector<index_type> pointIds;
        pointIds.reserve( M_pointPointers.size() );
        for ( auto const* point : M_pointPointers )
            pointIds.push_back(
                point ? point->id() : invalid_v<index_type> );
        this->setConnectivity( pointIds );
    }

    /**
     * It returns the reference to an point object (possibly derived from
     * Geo0D)
     */
    point_type& point( uint16_type i )
    {
        if ( this->canResolveCompactPoints() )
            return M_commonData->point(
                this->pointIds()[i] );
        auto& pointers = this->pointPointers();
        CHECK_LT( i, pointers.size() );
        CHECK( pointers[i] )
            << "invalid detached geometric point";
        return *pointers[i];
    }

    /**
     * It returns the reference to an point object (possibly derived from
     * Geo0D)
     */
    point_type const& point( uint16_type i ) const
    {
        if ( this->canResolveCompactPoints() )
            return std::as_const( *M_commonData ).point(
                this->pointIds()[i] );
        auto const& pointers = this->pointPointers();
        CHECK_LT( i, pointers.size() );
        CHECK( pointers[i] )
            << "invalid detached geometric point";
        return *pointers[i];
    }
    /**
     * It returns the pointer to an point object (possibly derived from
     * Geo0D)
     */
    point_type* pointPtr( uint16_type i )
    {
        return std::addressof( this->point( i ) );
    }

    /**
     * It returns the pointer to an point object (possibly derived from
     * Geo0D)
     */
    point_type const* pointPtr( uint16_type i ) const
    {
        return std::addressof( this->point( i ) );
    }

    /**
     *  The ith point (starting from the end)
     *
     *  It returns the reference to an point object (possibly derived from
     *  Geo0D). It starts from the last point, yet it follows the rule: vertices
     *  first. It may be used to access the points of a Geometry Element in a
     *  reverse way (i.e. with the opposite GeoElement permutation)
     */
    point_type& reversepoint( uint16_type const i )
    {
        return this->point(
            Feel::detail::ReversePoint<GEOSHAPE>().operate( i ) );
    }

    /**
     *  The ith point (starting from the end)
     *
     *  It returns the reference to an point object (possibly derived from
     *  Geo0D). It starts from the last point, yet it follows the rule: vertices
     *  first. It may be used to access the points of a Geometry Element in a
     *  reverse way (i.e. with the opposite GeoElement permutation)
     */
    point_type const& reversepoint( uint16_type const i ) const
    {
        return this->point(
            Feel::detail::ReversePoint<GEOSHAPE>().operate( i ) );
    }

    /**
     * Inserts a point.  Uses point references
     * put point
     */
    void setPoint( uint16_type i, point_type& p );

    /**
     * show information about the geoND
     *
     * @param verbose true if verbose mode, false otherwise
     * @param c output stream
     *
     * @return the output stream
     */
    std::ostream& showMe( bool verbose = false, std::ostream& c = std::cout ) const;

    /**
     * Swap Points
     *
     *  This is a member function to be used ONLY by routines for checking or
     *  amending meshes. You must give the local id
     *
     * @param pt1 1st point to swap with 2nd point
     * @param pt2 2nd point to swap with 1st point
     */
    void swapPoints( const uint16_type& pt1, const uint16_type& pt2 );

    /**
     *  Exchange Points
     *
     *  Exchanges points according to a list of old2new local id numbering !
     *  old2new[i] is the new local id of a point whose old local id was ! i+1
     *  (remember the numbering from 1 of the id's!. This is a member function
     *  to be used ONLY by routines for checking or amending meshes. You must
     *  give uint16_type (which start
     */
    void exchangePoints( const uint16_type otn[numPoints] );

    /**
     * matrix of geometric nodes
     * retrieve the matrix of geometric nodes (Dim x NumPoints) the
     * matrix is column oriented, the column i contains the coordinate
     * of the i-th geometric node of the element
     *
     * \return the matrix of geometric nodes
     */
    matrix_node_type G() const
    {
        matrix_node_type G( nRealDim, numPoints );
        for ( uint16_type i = 0; i < numPoints; ++i )
            ublas::column( G, i ) = this->point( i ).node();

        return G;
    }

    /**
     * matrix of vertices geometric nodes
     * retrieve the matrix of geometric nodes (Dim x NumPoints) the
     * matrix is column oriented, the column i contains the coordinate
     * of the i-th geometric node of the element
     *
     * \return the matrix of vertices geometric nodes
     */
    //matrix_node_type  vertices() const { return ublas::subrange( M_G, 0, nRealDim, 0, numVertices ); }
    matrix_node_type vertices() const
    {
        matrix_node_type G( nRealDim, numVertices );
        for ( uint16_type i = 0; i < numVertices; ++i )
            ublas::column( G, i ) = this->point( i ).node();
        return G;
        // return ublas::subrange( M_G, 0, nRealDim, 0, numVertices );
    }

    //! update G with all points of the element
    template <typename TheNodesType>
    void updateG( TheNodesType & G,
                  std::enable_if_t<is_eigen_matrix_v<TheNodesType> >* = nullptr ) const
    {
        DCHECK( nRealDim == G.rows() && G.cols() >= numPoints ) << "G is not compatible";
        for ( uint16_type i = 0; i < numPoints; ++i )
            G.col( i ) = em_fixed_size_cmatrix_t<nRealDim,1,value_type>(
                this->point( i ).node().data().begin() );
    }

    //! update G with all vertices of the element
    template <typename TheNodesType>
    void updateVertices( TheNodesType & G,
                         std::enable_if_t<is_eigen_matrix_v<TheNodesType> >* = nullptr ) const
    {
        DCHECK( nRealDim == G.rows() && G.cols() >= numVertices ) << "G is not compatible";
        for ( uint16_type i = 0; i < numVertices; ++i )
            G.col( i ) = em_fixed_size_cmatrix_t<nRealDim,1,value_type>(
                this->point( i ).node().data().begin() );
    }

    point_iterator beginPoint()
    {
        return this->pointPointers().begin();
    }
    point_const_iterator beginPoint() const
    {
        return this->pointPointers().begin();
    }
    point_iterator endPoint()
    {
        return this->pointPointers().end();
    }
    point_const_iterator endPoint() const
    {
        return this->pointPointers().end();
    }

    /**
     * get the max length of the edges of the element
     *
     *
     * @return the max length of the edges of the element
     */
    value_type h() const
    {
        matrix_node_type nodesG = this->G();
        em_matrix_col_type<value_type> G( const_cast<value_type*>( nodesG.data().begin() ), nodesG.size1(), nodesG.size2() );
        value_type res = 0.;
        for ( uint16_type __e = 0; __e < numLocalEdges; ++__e )
        {
            int col1 = this->eToP( __e, 0 );
            int col2 = this->eToP( __e, 1 );
            value_type r = ( G.col( col1 ) - G.col( col2 ) ).norm();

            res = ( res > r ) ? res : r;
        }
        return res;
    }
    /**
     * @brief get the minimum edge length in the element
     * @return the minimum edge length in the element
     */
    value_type hMin() const
    {
        matrix_node_type nodesG = this->G();
        em_matrix_col_type<value_type> G( const_cast<value_type*>( nodesG.data().begin() ), nodesG.size1(), nodesG.size2() );

        value_type res = 1e10;
        for ( uint16_type __e = 0; __e < numLocalEdges; ++__e )
        {
            int col1 = this->eToP( __e, 0 );
            int col2 = this->eToP( __e, 1 );
            value_type r = ( G.col( col1 ) - G.col( col2 ) ).norm();
            res = ( res > r ) ? r : res;
        }
        return res;
    }
    /**
     * get the max length of the edge in the local face \c f
     *
     * @param f local id of the face
     *
     * @return the max length of the edges of the local face
     */
    value_type hFace( uint16_type f ) const
    {
        if ( nRealDim == 1 )
            return 1;

        constexpr int nEdges = GEOSHAPE::topological_face_type::numEdges;
        matrix_node_type nodesG = this->G();
        em_matrix_col_type<value_type> G( const_cast<value_type*>( nodesG.data().begin() ), nodesG.size1(), nodesG.size2() );

        value_type res = 0.;
        for ( uint16_type e = 0; e < nEdges; ++e )
        {
            // get edge id in face from local edge id
            int edg = this->f2e( f, ( nDim == 2 ) ? f : e );
            int col1 = this->eToP( edg, 0 );
            int col2 = this->eToP( edg, 1 );
            value_type r = ( G.col( col1 ) - G.col( col2 ) ).norm();
            res = ( res > r ) ? res : r;
        }
        return res;
    }

    value_type hEdge( uint16_type f ) const
    {
        int col1 = this->eToP( f, 0 );
        int col2 = this->eToP( f, 1 );
        auto const& node1 = this->point( col1 ).node();
        em_node_type<value_type> emnode1( const_cast<value_type*>( node1.data().begin() ), node1.size() );
        auto const& node2 = this->point( col2 ).node();
        em_node_type<value_type> emnode2( const_cast<value_type*>( node2.data().begin() ), node2.size() );
        return ( emnode1 - emnode2 ).norm();
    }

    struct tt
    {
        static uint16_type fToP( uint16_type const _localFace, uint16_type const _point )
        {
            return super::eToP( _localFace, _point );
        }
    };


    /**
     * \return the normals at the barycenter of the faces
     */
    matrix_node_type normals() const
    {
        matrix_node_type _normals( nRealDim, numTopologicalFaces );
        if ( nDim != nRealDim )
        {
            CHECK( false ) << "normal when nDim != nRealDim is not implemented";
            return _normals;
        }
        gm1_ptrtype gm1 = this->gm1();
        if ( !gm1.use_count() )
            gm1 = gm1_ptrtype( new gm1_type );

        auto const& baryOnRefFaces = gm1->referenceConvex().barycenterFaces();
        matrix_node_type baryOnRefFace( nRealDim, 1 );

        std::vector<std::map<uint16_type, matrix_node_type>> ctxPtsOnRefFaces( numTopologicalFaces );
        for ( uint16_type f = 0; f < numTopologicalFaces; ++f )
        {
            ublas::column( baryOnRefFace, 0 ) = ublas::column( baryOnRefFaces, f );
            for ( permutation_type __p( permutation_type::IDENTITY );
                  __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                ctxPtsOnRefFaces[f][__p.value()] = baryOnRefFace;
        }
        auto pcf = gm1->preComputeOnFaces( gm1, ctxPtsOnRefFaces );
        static const size_type gmc_context_v = /*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN;
        auto ctx = gm1->template context<gmc_context_v>( *this, pcf, 0 );
        em_matrix_col_type<value_type> ns( _normals.data().begin(), _normals.size1(), _normals.size2() );
        for ( uint16_type f = 0; f < numTopologicalFaces; ++f )
        {
            ctx->template update<gmc_context_v>( *this, f );
            ns.col( f ) = ctx->unitNormal( 0 );
            //ublas::column( _normals, f ) = ctx->unitNormal( 0 );
        }
        return _normals;
    }

    /**
     * \return the normal at the barycenter of the face \p f
     */
    node_type normal( uint16_type f ) const
    {
        if ( nDim != nRealDim )
        {
            CHECK( false ) << "normal when nDim != nRealDim is not implemented";
            return node_type( 0 );
        }

        gm1_ptrtype gm1 = this->gm1();
        if ( !gm1.use_count() )
            gm1 = gm1_ptrtype( new gm1_type );

        matrix_node_type baryOnFace( nRealDim, 1 );
        ublas::column( baryOnFace, 0 ) = gm1->referenceConvex().faceBarycenter( f );
        auto pcf = gm1->preComputeOnFaces( gm1, baryOnFace );
        static const size_type gmc_context_v =/*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN;
        auto ctx = gm1->template context<gmc_context_v>( *this, pcf, f );
        ctx->template update<gmc_context_v>( *this, f );

        node_type n( nRealDim );
        em_node_type<value_type> en( n.data().begin(), n.size() );
        en = ctx->unitNormal( 0 );
        return n;
    }

    /**
     * Get the local id of the point in the element
     *
     * @param _localFace local id of a face in the element
     * @param _point local id of a point in the face
     *
     * @return the local id of the point in the element
     */
    static uint16_type fToP( uint16_type const _localFace, uint16_type const _point )
    {
        if constexpr ( super::nDim != 2 )
            return super::fToP( _localFace, _point );
        else
            return tt::fToP( _localFace, _point );
    }

    /**
     * get the number of opposite points per face
     *
     * @return the number of opposite points per face
     */
    uint16_type nOppositePointsPerFace() const
    {
        return super::nbOppositePointsPerFace;
    }

    /**
     * faceToOppositePoint(i,j) = localId of jth opposite point to ith local face
     *
     * @return the localId of _point-th opposite point to _localFace-th local face
     */
    uint16_type faceToOppositePoint( uint16_type const _localFace, uint16_type const _point ) const
    {
        return super::faceToOppositePoint( _localFace, _point );
    }

    /**
     * Determines if the local numbering of a 2D element
     * is oriented anticlockwise
     */
    bool isAnticlockwiseOriented() const
    {

        // Calculate vectors originating from vertex zero
        ublas::matrix<T> orientation_matrix( nRealDim, nRealDim );

        for ( int i = 0; i < nRealDim; ++i )
        {
            matrix_node_type nodesG = this->G();
            ublas::row( orientation_matrix, i ) = ( ublas::column( nodesG, i + 1 ) -
                                                    ublas::column( nodesG, 0 ) );
        }

        LU<ublas::matrix<T>> lu( orientation_matrix );
        T sgn = lu.det();

        return ( sgn > 0 ) ? 1 : 0;
    }

    FEELPP_DEPRECATED
    void setPointCoordG( int i, ublas::vector<double> const& u )
    {
        // ublas::column( M_G, i ) = u;
    }
    void applyDisplacement( int i, ublas::vector<double> const& u )
    {
        // ublas::column( M_G, i ) += u;
        this->point( i ) += u;
    }
    FEELPP_DEPRECATED
    void applyDisplacementG( int i, ublas::vector<double> const& u )
    {
        // ublas::column( M_G, i ) += u;
    }

    //! \return the number of point element neighbors
    size_type numberOfPointElementNeighbors() const
    {
        return pointElementNeighborIds().size();
    }
    //! \return the set of ids of point element neighbors
    std::set<size_type> pointElementNeighborIds() const
    {
        std::set<size_type> n;
        for ( uint16_type __p = 0; __p < numPoints; ++__p )
        {
            auto const& point = this->point( __p );
            std::for_each( point.elements().begin(),
                           point.elements().end(),
                           [&n]( auto const& e ) { n.insert( e.first ); } );
        }
        return n;
    }

    value_type measure() const override
        {
            if constexpr ( UseMeasuresStorage )
                return this->measureMetadata().measure.value_or( 0. );
            else
            {
                gm_ptrtype gm = this->gm();
                if ( !gm.use_count() )
                    gm = gm_ptrtype( new gm_type );

                quad_meas_type thequad( meas_quad_order );
                auto pc = gm->preCompute( gm, thequad.points() );
                auto ctx = gm->template context<vm::JACOBIAN>( *this, pc );
                return this->computeMeasureImpl( thequad, ctx );
            }
        }
    value_type faceMeasure( uint16_type f ) const
        {
            if constexpr ( UseMeasuresStorage )
            {
                auto const& measures = this->measureMetadata().faceMeasures;
                return f < measures.size()
                           ? measures[f].value_or( 0. )
                           : value_type{ 0 };
            }
            else
            {
                gm_ptrtype gm = this->gm();
                if ( !gm.use_count() )
                    gm = gm_ptrtype( new gm_type );

                quad_meas_type thequad( meas_quad_order );
                auto pcf = gm->preComputeOnFaces( gm, thequad.allfpoints() );
                auto ctxf = gm->template context<vm::NORMAL | vm::KB | vm::JACOBIAN>( *this, pcf, f );
                return this->computeFaceMeasureImpl( thequad, ctxf, f );
            }
        }
    Eigen::Matrix<value_type,numTopologicalFaces,1> faceMeasures() const
        {
            Eigen::Matrix<value_type,numTopologicalFaces,1> m;
            for( int f = 0; f < numTopologicalFaces; ++f )
                m( f ) = faceMeasure( f );
            return m;
        }

    void setMeasurePointElementNeighbors( value_type meas )
    {
        if constexpr ( UseMeasuresStorage )
            this->mutableMeasureMetadata()
                .pointElementNeighborsMeasure = meas;
    }

    /**
     * @brief Return the cached point-element-neighbor patch measure.
     */
    [[nodiscard]] value_type measurePointElementNeighbors() const
    {
        if constexpr ( UseMeasuresStorage )
            return this->measureMetadata()
                .pointElementNeighborsMeasure.value_or( 0. );
        else
            return value_type{ 0 };
    }

    /**
     * @brief Cache the entity measure in mesh-owned sparse metadata.
     */
    void setMeasure( value_type meas ) const
        requires( UseMeasuresStorage )
    {
        this->mutableMeasureMetadata().measure = meas;
    }

    /**
     * @brief Cache a topological-face measure in sparse metadata.
     * @param f local topological-face index
     * @param meas face measure
     */
    void setFaceMeasure( uint16_type f, value_type meas ) const
        requires( UseMeasuresStorage )
    {
        CHECK_LT( f, numTopologicalFaces )
            << "invalid local face index " << f;
        auto& measures =
            this->mutableMeasureMetadata().faceMeasures;
        if ( measures.size() < numTopologicalFaces )
            measures.resize( numTopologicalFaces );
        measures[f] = meas;
    }

    /**
     * @return a quadrature to compute the measure of the element
     */
    static quad_meas_type imMeasure() { return quad_meas_type( meas_quad_order ); }

    /**
     * @return a quadrature to compute the measure of the order 1 approximation
     * of the element
     */
    static quad_meas1_type imMeasureOrder1() { return quad_meas1_type( 0 ); }

    //!
    //! check if an element is valid: check points are != nullptr
    //!
    bool isValid() const;

    void update();

    template <typename QuadType, typename CtxType, typename CtxFaceType>
    void updateWithCtx( QuadType const& thequad,
                        std::shared_ptr<CtxType>& ctx,
                        std::shared_ptr<CtxFaceType>& ctxf ) const;
    template <typename CtxType, typename CtxFaceType>
    void updateWithCtx1( quad_meas1_type const& thequad,
                         std::shared_ptr<CtxType>& ctx,
                         std::shared_ptr<CtxFaceType>& ctxf ) const;

  private:

    template <typename QuadType, typename CtxType>
        value_type computeMeasureImpl( QuadType const& thequad,
                                       std::shared_ptr<CtxType>& ctx ) const;
    template <typename QuadType, typename CtxType>
        value_type computeFaceMeasureImpl( QuadType const& thequad,
                                           std::shared_ptr<CtxType>& ctxf,
                                           uint16_type f ) const;

  private:
    /**
     * @brief Return whether compact connectivity can resolve every point.
     */
    [[nodiscard]] bool canResolveCompactPoints() const
    {
        if ( !M_commonData ||
             !M_commonData->hasPointResolver() )
            return false;
        auto const ids = this->pointIds();
        return ids.size() == this->nPointsRuntime() &&
               std::ranges::none_of(
                   ids,
                   []( index_type id )
                   {
                       return id == invalid_v<index_type>;
                   } );
    }

    /**
     * @brief Materialize the legacy point-pointer facade on demand.
     *
     * Compact connectivity remains canonical. This cache exists only for
     * source-compatible `points()` and point-iterator APIs.
     */
    void materializePointPointers() const
    {
        auto const pointCount =
            this->nPointsRuntime();
        M_pointPointers.resize( pointCount );
        M_pointPointers.materialize();

        if ( this->canResolveCompactPoints() )
        {
            auto const ids = this->pointIds();
            for ( size_type localPoint = 0;
                  localPoint < pointCount; ++localPoint )
                M_pointPointers[localPoint] =
                    std::addressof(
                        M_commonData->point(
                            ids[localPoint] ) );
        }
    }

    /** @brief Return the mutable lazy compatibility pointer cache. */
    [[nodiscard]] point_pointer_container_type&
    pointPointers()
    {
        this->materializePointPointers();
        return M_pointPointers;
    }

    /** @brief Return the immutable lazy compatibility pointer cache. */
    [[nodiscard]] point_pointer_container_type const&
    pointPointers() const
    {
        this->materializePointPointers();
        return M_pointPointers;
    }

    friend class boost::serialization::access;
    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
        DVLOG( 2 ) << "Serializing GeoND...\n";
        DVLOG( 2 ) << "  - base class...\n";
        ar& boost::serialization::base_object<super>( *this );
        DVLOG( 2 ) << "  - points...\n";
        auto& pointPointers = this->pointPointers();
        for ( uint16_type i = 0;
              i < this->nPointsRuntime(); ++i )
        {
            if ( Archive::is_loading::value )
            {
                // WARNING ! a memory leak is introduce here
                // we need to register in M_commonData for example : in the mesh if exist else in this struct
                if ( !pointPointers[i] )
                    pointPointers[i] = new point_type;
            }
            ar&*( pointPointers[i] );
        }
        if ( version >= 1 )
        {
            ar& M_neighborsInitialized;
            if ( M_neighborsInitialized )
            {
                for ( auto& neighborId : M_neighbors )
                    ar& neighborId;
            }
        }
        else if constexpr ( Archive::is_loading::value )
        {
            M_neighborsInitialized = false;
        }
        // DVLOG(2) << "  - G...\n";
        // ar & M_G;
    }

  private:
    /**
     * Lazy source-compatibility cache of geometric point facades.
     *
     * Stored mesh entities leave this empty and resolve through compact
     * connectivity. Detached importer entities materialize it temporarily.
     */
    mutable point_pointer_container_type M_pointPointers;

    //! Neighbor element identifiers in local-face order.
    std::array<size_type, numNeighbors> M_neighbors{};

    //! Whether the mesh topology builder initialized neighbor connectivity.
    bool M_neighborsInitialized = false;

    //! common data shared in a collection of multi-dimensional geometrical entity
    mutable GeoNDCommon<self_type> * M_commonData;

    //! @brief Runtime order storage (only occupies space when is_order_dynamic is true)
    [[no_unique_address]]
    std::conditional_t<is_order_dynamic, uint16_type, std::monostate> M_runtime_order{};
};

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
inline void
GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::setPoint( uint16_type i, point_type& p )
{
    if ( M_pointPointers )
    {
        CHECK_LT( i, M_pointPointers.size() );
        M_pointPointers[i] = std::addressof( p );
    }
    this->setConnectivityEntry(
        i, p.id(), this->nPointsRuntime() );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
std::ostream&
GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::showMe( bool verbose, std::ostream& out ) const
{
    out << "----- BEGIN OF GeoND data ---" << std::endl
        << std::endl;
    out << " GeoND object of shape " << Shape << std::endl;
    out << " Number of Vertices = " << numVertices << std::endl;
    out << "   Number of Points = " << numPoints << std::endl;
    out << "                 id = " << this->id() << std::endl;
    out << "                  G = " << this->G() /*M_G*/ << "\n";

    for ( int i = 0; i < numVertices; i++ )
    {
        out << "POINT id = " << i << std::endl;
        point( i ).showMe( verbose, out );
    }

    out << "----- END OF GeoND data ---" << std::endl
        << std::endl;
    return out;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::swapPoints( const uint16_type& pt1, const uint16_type& pt2 )
{
    auto connectivity = std::vector<index_type>(
        this->pointIds().begin(), this->pointIds().end() );
    CHECK_LT( pt1, connectivity.size() );
    CHECK_LT( pt2, connectivity.size() );
    std::swap( connectivity[pt1], connectivity[pt2] );
    this->setConnectivity( connectivity );
    M_pointPointers.reset();
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::exchangePoints( const uint16_type otn[numPoints] )
{
    auto const current = this->pointIds();
    std::vector<index_type> connectivity(
        current.begin(), current.end() );
    for ( uint16_type i = 0; i < numPoints; ++i )
        connectivity[i] = current[otn[i]];
    this->setConnectivity( connectivity );
    M_pointPointers.reset();
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
bool GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::isValid() const
{
    auto const valid = this->hasPoints();
    DCHECK( valid )
        << "element with invalid point : " << this->id();
    return valid;
}
template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::update()
{
    if constexpr ( !UseMeasuresStorage )
       return;

    gm_ptrtype gm = this->gm();
    if ( !gm.use_count() )
        gm = gm_ptrtype( new gm_type );

    quad_meas_type thequad( meas_quad_order );
    auto pc = gm->preCompute( gm, thequad.points() );
    auto pcf = gm->preComputeOnFaces( gm, thequad.allfpoints() );

    auto ctx = gm->template context<vm::JACOBIAN>( *this, pc );
    auto ctxf = gm->template context</*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN>( *this, pcf, 0 );
    this->updateWithCtx( thequad, ctx, ctxf );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
template <typename QuadType, typename CtxType>
typename GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::value_type
GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::computeMeasureImpl( QuadType const& thequad,
                                                                                    std::shared_ptr<CtxType>& ctx ) const
{
    value_type meas = 0.;
    if constexpr ( CtxType::is_linear )
    {
        meas = thequad.weightsSum() * ctx->J( 0 );
    }
    else
    {
        for ( int q = 0; q < thequad.nPoints(); ++q )
            meas += thequad.weight( q ) * ctx->J( q );
    }
    return meas;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
template <typename QuadType, typename CtxFaceType>
typename GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::value_type
GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::computeFaceMeasureImpl( QuadType const& thequad,
                                                                                        std::shared_ptr<CtxFaceType>& ctxf,
                                                                                        uint16_type f ) const
{
    value_type measFace = 0.;
    if constexpr ( nDim == nRealDim )
    {
        for ( int q = 0; q < thequad.nPointsOnFace( f ); ++q )
            measFace += thequad.weight( f, q ) * ctxf->J( q ) * ctxf->normalNorm( q );
    }
    return measFace;
}


template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
template <typename QuadType, typename CtxType, typename CtxFaceType>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::updateWithCtx( QuadType const& thequad,
                                                        std::shared_ptr<CtxType>& ctx,
                                                        std::shared_ptr<CtxFaceType>& ctxf ) const
{
    if constexpr ( UseMeasuresStorage )
    {
        this->setMeasure( this->computeMeasureImpl( thequad,ctx ) );

        if constexpr ( nDim == nRealDim )
        {
            for ( uint16_type f = 0; f < numTopologicalFaces; ++f )
            {
                ctxf->template update</*vm::POINT|*/ vm::NORMAL | vm::KB | vm::JACOBIAN>( dynamic_cast<typename CtxFaceType::element_type const&>( *this ), f );
                this->setFaceMeasure( f, computeFaceMeasureImpl( thequad,ctxf, f ) );
            }
        }
    }
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
template <typename CtxType, typename CtxFaceType>
void GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage>::updateWithCtx1( quad_meas1_type const& thequad,
                                                         std::shared_ptr<CtxType>& ctx,
                                                         std::shared_ptr<CtxFaceType>& ctxf ) const
{
    updateWithCtx( thequad, ctx, ctxf );
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
inline FEELPP_DEPRECATED
    DebugStream&
operator<<( DebugStream& __os, GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE,UseMeasuresStorage> const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n.showMe( true, __str );

        __os << __str.str() << "\n";
    }

    return __os;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
inline NdebugStream& FEELPP_DEPRECATED
operator<<( NdebugStream& __os, GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE, UseMeasuresStorage> const& __n )
{
    return __os;
}

template <uint16_type Dim, typename GEOSHAPE, typename T, typename IndexT, typename POINTTYPE, bool UseMeasuresStorage>
inline std::ostream&
operator<<( std::ostream& __os, GeoND<Dim, GEOSHAPE, T, IndexT, POINTTYPE,UseMeasuresStorage> const& __n )
{
    return __n.showMe( true, __os );
}

} // namespace Feel

namespace boost::serialization
{
/**
 * @brief Serialization version for geometric entities.
 *
 * Version 1 adds inline neighbor identifiers while retaining the ability to
 * read version-0 archives, whose adjacency is reconstructed by mesh update.
 */
template <Feel::uint16_type Dim, typename GeoShape, typename T,
          typename IndexT, typename PointType, bool UseMeasuresStorage>
struct version<
    Feel::GeoND<Dim, GeoShape, T, IndexT, PointType,
                UseMeasuresStorage>>
{
    using type = boost::mpl::int_<1>;
    using tag = boost::mpl::integral_c_tag;
    static constexpr int value = type::value;
};
} // namespace boost::serialization

#endif
