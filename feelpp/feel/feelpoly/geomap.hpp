/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2008-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file geomap.hpp
 * @brief Geometric mappings, pointwise mapping contexts, and inverse maps.
 */
#ifndef _GEOMAP_H
#define _GEOMAP_H

#include <boost/version.hpp>
#if ( BOOST_VERSION >= 103400 )
#include <boost/none.hpp>
#else
#include <boost/none_t.hpp>
#endif /* BOOST_VERSION >= 103400 */

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/StdVector>
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <Eigen/src/Core/util/ReenableStupidWarnings.h>

#include <boost/functional/hash.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>

// #include <boost/numeric/bindings/traits/traits.hpp>
// #include <boost/numeric/bindings/traits/ublas_vector.hpp>
// #include <boost/numeric/bindings/traits/ublas_matrix.hpp>
// #include <boost/numeric/bindings/atlas/cblas3.hpp>

#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/svd.hpp>

#include <feel/feelmesh/marker.hpp>
#include <feel/feelmesh/concepts.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/entitymarkers.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/context.hpp>
#include <feel/feelpoly/expansiontypes.hpp>
#include <feel/feelpoly/fekete.hpp>
#include <feel/feelpoly/geometricinverse.hpp>
#include <feel/feelpoly/geometryevaluation.hpp>
#include <feel/feelpoly/geometrynodallagrange.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/manifoldgeometry.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/piolatransform.hpp>

namespace Feel
{
/**
 * @enum GeomapStrategyType
 * @brief Select the geometric mapping used by an integration algorithm.
 */
enum class GeomapStrategyType
{
    /** Let the consumer select its optimized mapping strategy. */
    GEOMAP_OPT = 0,
    /** Force the affine order-one approximation of the geometry. */
    GEOMAP_O1 = 1,
    /** Use the full high-order geometry. */
    GEOMAP_HO = 2
};
namespace detail
{
/**
 * @brief Normalize a requested mapping strategy for an entity range.
 * @tparam RangeType mesh-entity range type
 * @param gs requested mapping strategy
 * @return the strategy applicable to the range geometry
 *
 * An affine mesh has no distinct high-order approximation, so its effective
 * strategy is normalized to `GEOMAP_HO`.
 */
template <typename RangeType>
GeomapStrategyType
geomapStrategy( RangeType const& /**/, GeomapStrategyType gs )
{
    constexpr uint16_type geoOrder = entity_range_t<RangeType>::nOrder;
    if ( geoOrder == 1 )
        return GeomapStrategyType::GEOMAP_HO;
    else
        return gs;
}
} // namespace detail

// namespace blas = boost::numeric::bindings::blas;
// namespace traits = boost::numeric::bindings::traits;

/** @brief Tag selecting inverse geometric-mapping operations. */
struct GeomapInverse
{
    // GeomapInverse
};

template <int O,
          template <uint16_type Dim> class PolySetType,
          typename ContinuityType,
          template <class, int, class> class Pts,
          uint16_type TheTag>
class Lagrange;

/**
 * @class GeoMap
 * @brief Polynomial mapping from a reference entity to a physical entity.
 * @author C. Prud'homme
 *
 * This class contains the geometrical transformation that maps the reference
 * element on the current element, and its values on integration points.
 * `Order` may be a positive compile-time order or `Dynamic`; the default is
 * the affine order-one mapping used by most meshes.
 *
 * @tparam Dim topological dimension of the reference entity
 * @tparam Order polynomial mapping order, defaulting to one
 * @tparam RealDim dimension of the ambient physical space
 * @tparam T scalar value type
 * @tparam Entity reference-entity family, such as `Simplex` or `Hypercube`
 * @tparam PP legacy scalar polynomial family customization point; the public
 *            compatibility façade currently requires nodal Lagrange semantics
 *
 */
template <int Dim,
          int Order = 1,
          int RealDim = Dim,
          typename T = double,
          template <int, int, int> class Entity = Simplex,
          template <int, template <uint16_type RDim> class PolySetType, typename ContinuityType,
                    template <class, int, class> class Pts, uint16_type>
          class PP = Lagrange>
class GeoMap
    : public PP<Order, Scalar, Continuous, PointSetEquiSpaced, 0>::template apply<Dim, RealDim /*Dim*/, T, Entity<Dim, Order, /*RealDim*/ Dim>>::result_type //,
                                                                                                                                                             // public boost::enable_shared_from_this<GeoMap<Dim, Order, RealDim, T, Entity, PP > >
// public PP<Order,Scalar, PointSetFekete>::template apply<Dim, T, Entity<Dim,Order,Dim> >::result_type
{
    static_assert( Order == Dynamic || Order >= 1,
                   "A geometric mapping polynomial order must be at least one" );

    /**
     * @brief Enforce the positive-order invariant of a geometric mapping.
     * @param order requested runtime polynomial order
     * @return the validated order tag
     * @throws std::invalid_argument if the requested order is zero
     */
    [[nodiscard]] static RuntimeOrder validatedGeometryOrder( RuntimeOrder order )
    {
        if ( order.value < 1 )
            throw std::invalid_argument( "GeoMap runtime polynomial order must be at least one" );
        return order;
    }

    // typedef typename PP<Order, Scalar, PointSetFekete>::template apply<Dim, T, Entity<Dim,Order,Dim> >::result_type super;
    typedef typename PP<Order, Scalar, Continuous, PointSetEquiSpaced, 0>::template apply<Dim, RealDim /*Dim*/, T, Entity<Dim, Order, /*RealDim*/ Dim>>::result_type super;

    // typedef boost::enable_shared_from_this<GeoMap<Dim, Order, RealDim, T, Entity, PP > > super_enable_this;

    static inline const uint16_type nRealDimCheck2d = mpl::if_<mpl::less_equal<mpl::int_<2>, mpl::int_<RealDim>>,
                                                               mpl::int_<RealDim>,
                                                               mpl::int_<Dim>>::type::value;

    typedef mpl::vector<boost::none_t, boost::none_t, boost::none_t,
                        GeoMap<1, Order, RealDim, T, Entity, PP>>
        geomap_edges_t;

    typedef mpl::vector<boost::none_t, boost::none_t,
                        GeoMap<1, Order, RealDim, T, Entity, PP>,
                        GeoMap<2, Order, nRealDimCheck2d /*RealDim*/, T, Entity, PP>>
        geomap_faces_t;

    typedef mpl::vector<GeoMap<1, Order, RealDim, T, Entity, PP>,
                        GeoMap<2, Order, nRealDimCheck2d /*RealDim*/, T, Entity, PP>,
                        GeoMap<3, Order, 3, T, Entity, PP>,
                        boost::none_t>
        geomap_elements_t;

  public:
    static constexpr uint16_type nDim = super::nDim;
    static constexpr uint16_type nRealDim = super::nRealDim;
    static constexpr uint16_type nDof = super::nDof;
    //! @brief True if Order is known at compile time (i.e., Order >= 0)
    static constexpr bool is_order_static = ( Order != Dynamic );
    //! @brief True if Order is determined at runtime (i.e., Order == Dynamic == -1)
    static constexpr bool is_order_dynamic = !is_order_static;
    //! @brief Template order parameter (use order() accessor for runtime value)
    static constexpr int nOrder_v = Order;
    //! @brief Static order (1 for dynamic as placeholder in template expressions)
    static constexpr uint16_type nOrder = ( Order >= 0 ) ? static_cast<uint16_type>( Order ) : 1;
    static constexpr uint16_type nNodes = super::nNodes;
    static constexpr fem::transformation_type trans = super::trans;
    // Dynamic-order geomaps are not linear at compile-time: linearity is decided at runtime.
    static constexpr bool is_linear = is_order_dynamic ? false : ( trans == fem::LINEAR );

    typedef typename super::value_type value_type;

    /**
     * @brief Explicit mathematical representation used by the compatibility façade.
     *
     * Existing `GeoMap` coefficient tables are physical interpolation nodes;
     * this alias makes that semantic contract visible to generic geometry code.
     */
    using representation_policy_type =
        geometry::NodalLagrange<Dim, Order, Entity, value_type>;
    using bezier_representation_policy_type =
        geometry::BernsteinBezier<Dim, Order, Entity, value_type>;
    /** Coefficient semantics of the current `GeoMap` compatibility façade. */
    /** @brief Return the coefficient semantics of the compatibility façade. */
    [[nodiscard]] static consteval geometry::CoefficientSemantics
    coefficientSemantics() noexcept
    {
        return representation_policy_type::coefficientSemantics();
    }

    typedef typename super::PreCompute precompute_type;
    typedef std::shared_ptr<precompute_type> precompute_ptrtype;

    typedef typename super::convex_type convex_type;
    static constexpr int reference_order =
        Order == Dynamic ? 1 : Order;
    typedef Reference<
        convex_type, nDim, reference_order,
        nDim /*nRealDim*/>
        reference_convex_type;

    typedef GeoMap<Dim, Order, RealDim, T, Entity, PP> self_type;
    typedef self_type geometric_mapping_type;
    typedef std::shared_ptr<geometric_mapping_type> geometric_mapping_ptrtype;

    typedef typename mpl::at<geomap_elements_t, mpl::int_<nDim>>::type element_gm_type;
    typedef std::shared_ptr<element_gm_type> element_gm_ptrtype;
    typedef element_gm_ptrtype gm_ptrtype;
    using gm_type = element_gm_type;

    typedef typename mpl::at<geomap_faces_t, mpl::int_<nDim>>::type face_gm_type;
    typedef std::shared_ptr<face_gm_type> face_gm_ptrtype;

    template <int N>
    struct face_gm
    {
        typedef typename mpl::at<geomap_faces_t, mpl::int_<N>>::type type;
    };

    typedef typename mpl::at<geomap_edges_t, mpl::int_<nDim>>::type edge_gm_type;
    typedef std::shared_ptr<edge_gm_type> edge_gm_ptrtype;

    template <int N>
    struct edge_gm
    {
        typedef typename mpl::at<geomap_edges_t, mpl::int_<N>>::type type;
    };

    typedef typename node<value_type>::type node_t_type;
    typedef typename matrix_node<value_type>::type matrix_node_t_type;

    typedef node_t_type normal_type;
    typedef ublas::vector<normal_type> normals_type;
    typedef typename normals_type::const_iterator normal_const_iterator;
    typedef node_t_type tangent_type;
    typedef ublas::vector<tangent_type> tangents_type;
    typedef typename tangents_type::const_iterator tangent_const_iterator;

    typedef typename ublas::vector<value_type> vector_type;
    typedef typename ublas::matrix<value_type, ublas::column_major> matrix_type;

    using eigen_vector_p_type = eigen_vector_type<Dim, value_type>;
    using vector_eigen_vector_p_type = vector_eigen_vector_type<Dim, value_type>;

    using eigen_vector_n_type = eigen_vector_type<RealDim, value_type>;
    using vector_eigen_vector_n_type = vector_eigen_vector_type<RealDim, value_type>;

    using hessian_basis_type = Eigen::Tensor<value_type, 3>;

    /** default constructor (static order only) */
    GeoMap()
        requires( is_order_static )
        : super(),
          M_is_cached( false ),
          _elementMap(),
          _boundaryMap(),
          M_g_linear( this->nbPoints(), nDim ),
          M_refconvex()
    {
        initGeometryRepresentation();
        initLinearGradientCache();
    }

    /** constructor with runtime order (dynamic order only) */
    explicit GeoMap( RuntimeOrder o )
        requires( is_order_dynamic )
        : super( validatedGeometryOrder( o ) ),
          M_is_cached( false ),
          _elementMap(),
          _boundaryMap(),
          M_g_linear( this->nbPoints(), nDim ),
          M_refconvex()
    {
        initGeometryRepresentation();
        initLinearGradientCache();
    }

    /** constructor with element and face maps (static order only) */
    GeoMap( element_gm_ptrtype const& e, face_gm_ptrtype const& f )
        requires( is_order_static )
        : super(),
          M_is_cached( false ),
          _elementMap( e ),
          _boundaryMap( f ),
          M_g_linear( this->nbPoints(), nDim ),
          M_refconvex()
    {
        initGeometryRepresentation();
        initLinearGradientCache();
    }
    /**
            destructor
            */
    ~GeoMap() override
    {
    }

    /**
            \return the dimension of the underlying element
            */
    constexpr uint16_type dim() const noexcept
    {
        return nDim;
    }

    constexpr uint16_type realDim() const noexcept
    {
        return nRealDim;
    }

    /**
     * @return the polynomial order of the geometric mapping
     * For static order, returns the template parameter.
     * For dynamic order, returns the runtime value.
     */
    [[nodiscard]] constexpr uint16_type order() const noexcept
    {
        if constexpr ( is_order_static )
            return nOrder;
        else if constexpr ( requires( super const& s ) { s.order(); } )
            return static_cast<uint16_type>( static_cast<super const&>( *this ).order() );
        else
            return nOrder;
    }

    /**
   \return true if the geometric mapping is linear, false otherwise
*/
    bool isLinear() const
    {
        if constexpr ( is_order_dynamic )
        {
            using static_order_one_type =
                GeoMap<Dim, 1, RealDim, T, Entity, PP>;
            return this->order() == 1 && static_order_one_type::is_linear;
        }
        else
            return is_linear;
    }

    /**
     *  \return the natural mapping for the element
     */
    const element_gm_ptrtype elementMap() const
    {
        return _elementMap;
    }

    /**
     *  \return the natural mapping for the boundary of the element
     */
    const face_gm_ptrtype& boundaryMap() const
    {
        return _boundaryMap;
    }

    /**
     * \return the i-th reference node associated with the geometric
     * mapping construction
     */
    ublas::vector<value_type> refNode( uint16_type i ) const
    {
        return ublas::column( this->points(), i );
    }

    /**
     * @brief Return the explicit nodal representation policy of this map.
     *
     * The policy is constructed from the same reference nodes as the legacy
     * finite-element base. It supplies representation-neutral evaluation and
     * exact nodal/Bezier conversion without duplicating `Context` kernels.
     *
     * @return nodal-Lagrange representation bound to this map's effective order
     */
    [[nodiscard]] representation_policy_type const& geometryRepresentation() const
    {
        return *M_geometryRepresentation;
    }

    /**
     * @brief Convert canonical stored coefficients to nodal interpolation data.
     * @tparam Matrix writable Eigen- or uBLAS-compatible matrix
     * @param coefficients component-by-coefficient table, modified in place
     * @param representation semantics of the stored table
     */
    template <typename Matrix>
    void convertStoredGeometryToNodal(
        Matrix& coefficients,
        geometry::GeometryRepresentationKind representation ) const
    {
        if ( representation ==
             geometry::GeometryRepresentationKind::NodalLagrange )
            return;
        this->validateRepresentation( representation );
        if ( this->matrixColumns( coefficients ) != this->nbPoints() )
            throw std::invalid_argument(
                "GeoMap geometry conversion has an invalid coefficient count" );

        auto canonical = this->copyToEigen( coefficients );
        decltype( canonical ) internal( canonical.rows(), canonical.cols() );
        for ( Eigen::Index coefficient = 0;
              coefficient < internal.cols(); ++coefficient )
            internal.col( coefficient ) = canonical.col(
                M_bernsteinToCanonical[coefficient] );
        this->copyFromEigen(
            M_geometryRepresentation->fromBezier( internal ),
            coefficients );
    }

    /**
     * @brief Convert nodal interpolation data to canonical stored coefficients.
     * @tparam Matrix writable Eigen- or uBLAS-compatible matrix
     * @param coefficients component-by-coefficient table, modified in place
     * @param representation requested storage semantics
     */
    template <typename Matrix>
    void convertNodalGeometryToStored(
        Matrix& coefficients,
        geometry::GeometryRepresentationKind representation ) const
    {
        if ( representation ==
             geometry::GeometryRepresentationKind::NodalLagrange )
            return;
        this->validateRepresentation( representation );
        if ( this->matrixColumns( coefficients ) != this->nbPoints() )
            throw std::invalid_argument(
                "GeoMap geometry conversion has an invalid coefficient count" );

        auto const internal = M_geometryRepresentation->toBezier(
            this->copyToEigen( coefficients ) );
        auto canonical = internal;
        for ( Eigen::Index coefficient = 0;
              coefficient < internal.cols(); ++coefficient )
            canonical.col( M_bernsteinToCanonical[coefficient] ) =
                internal.col( coefficient );
        this->copyFromEigen( canonical, coefficients );
    }

    /** @brief Resolve coefficient semantics from an entity or compact view. */
    template <typename GeometricEntity>
    [[nodiscard]] static geometry::GeometryRepresentationKind
    geometryRepresentationOf( GeometricEntity const& entity ) noexcept
    {
        if constexpr ( GeometricEntityViewType<GeometricEntity> )
            return entity.mesh().geometryRepresentation();
        else if constexpr ( requires { entity.mesh(); } )
        {
            auto const* mesh = entity.mesh();
            return mesh ? mesh->geometryRepresentation()
                        : geometry::GeometryRepresentationKind::NodalLagrange;
        }
        else
            return geometry::GeometryRepresentationKind::NodalLagrange;
    }

    /**
     * @return the reference convex
     */
    reference_convex_type const& referenceConvex() const
    {
        return M_refconvex;
    }

    /**
     * @return a positive number if the point is within the reference
     * convex
     */
    boost::tuple<bool, value_type> isIn( typename node<value_type>::type const& pt ) const
    {
        return M_refconvex.isIn( pt );
    }

    /**
     *  apply the geometric mapping to the point \c pt given the real
     *  geometric nodes stored in a NxNg matrix \c G
     */
    node_t_type transform( const node_t_type& ref_p,
                           const matrix_node_t_type& G ) const
    {
        matrix_node_t_type referencePoint( nDim, 1 );
        ublas::column( referencePoint, 0 ) = ref_p;
        auto const basisValues = super::evaluate( referencePoint );
        if ( basisValues.size1() != G.size2() )
            throw std::invalid_argument(
                "GeoMap::transform(): invalid geometric basis size" );

        typename node<value_type>::type real_p( G.size1() );
        real_p.clear();
        for ( size_type i = 0; i < G.size2(); ++i )
            real_p.plus_assign(
                basisValues( i, 0 ) * ublas::column( G, i ) );

        return real_p;
    }

    /**
     * apply the geometric mapping to the point index \c id_pt given the real
     * geometric nodes stored in a NxNg matrix \c G
     */
    node_t_type transform( uint16_type id_ref,
                           const matrix_node_t_type& G,
                           const precompute_type* pc ) const
    {
        // Initialize the real_p node with the appropriate size and clear it
        node_t_type real_p( G.size1() );
        real_p.clear();

        // Loop over all nodes
        const auto nPts = static_cast<size_type>( G.size2() );
        for ( size_type i = 0; i < nPts; ++i )
        {
            // Evaluate the transformation at the given point
            value_type phi_at_pt = pc->phi( id_ref, i );

            // Accumulate the result of phi_at_pt multiplied by the i-th column of G
            real_p.plus_assign( phi_at_pt * ublas::column( G, i ) );
        }

        return real_p;
    }

    /**
     * compute real coordinates from a matrix of ref coordinates
     */
    void transform( matrix_node_t_type const& G,
                    precompute_type const* pc,
                    matrix_type& x ) const
    {
        eigenMap( x ).noalias() = eigenMap( G ) * pc->phiEigen();
        // ublas::axpy_prod( G, pc->phi(), x, true );
    }

    /**
     *  compute the gradient of the transformation in the reference
     *  element
     *
     *  Compute the gradient at node \c x, pc is resized to
     *  [nbNodes() x dim()] if the transformation is linear, \c x is
     *  not used at all
     */
    void gradient( const node_t_type& __pt,
                   matrix_type& __g ) const
    {
        auto gradientView = eigenMap( __g );
        this->gradient( __pt, gradientView );
    }

    /**
     * @brief Evaluate the reference gradient into Eigen-compatible storage.
     *
     * This overload is the shared implementation used by the legacy uBLAS
     * facade and Eigen-owned inverse workspaces.
     *
     * @tparam Derived writable Eigen dense matrix expression
     * @param point reference point
     * @param gradientValues coefficient-by-reference-direction output
     */
    template <typename Derived>
    void gradient( node_t_type const& point,
                   Eigen::MatrixBase<Derived>& gradientValues ) const
    {
        if ( gradientValues.rows() != this->nbPoints() )
            throw std::invalid_argument(
                "GeoMap::gradient(): invalid coefficient count" );
        if ( gradientValues.cols() != dim() )
            throw std::invalid_argument(
                "GeoMap::gradient(): invalid reference dimension" );
        if ( this->isLinear() )
        {
            gradientValues = eigenMap( M_g_linear );
        }

        else
        {
            if ( point.size() != dim() )
                throw std::invalid_argument(
                    "GeoMap::gradient(): invalid point dimension" );

            matrix_node_t_type __pts( nDim, 1 );
            ublas::column( __pts, 0 ) = point;

            ublas::vector<ublas::matrix<value_type>> m = super::derivate( __pts );

            for ( uint16_type n = 0; n < nDim; ++n )
            {
                for ( size_type coefficient = 0;
                      coefficient < this->nbPoints(); ++coefficient )
                    gradientValues( coefficient, n ) =
                        m( n )( coefficient, 0 );
            }
        }
    }

    /**
   compute the gradient of the transformation in the reference
   element

   Compute the gradient at node \c x, pc is resized to
   [nbNodes() x dim()] if the transformation is linear, \c x is
   not used at all
*/
    void gradient( uint16_type __idref,
                   matrix_type& __g,
                   precompute_type const* __pc ) const
    {
        if ( !__pc )
            throw std::invalid_argument(
                "GeoMap::gradient(): missing reference precomputation" );

        if ( this->isLinear() )
        {
            __g = M_g_linear;
        }

        else
        {
            if ( __pc->dim() != dim() )
                throw std::invalid_argument(
                    "GeoMap::gradient(): incompatible precomputation dimension" );

            const auto nPts = static_cast<size_type>( __g.size1() );
            for ( size_type i = 0; i < nPts; ++i )
            {
                for ( uint16_type n = 0; n < nDim; ++n )
                {
                    __g( i, n ) = __pc->grad( i, 0, n, __idref );
                }
            }
        }
    }

    /**
     compute the gradient of the transformation in the reference
     element

     Compute the gradient at node \c x, pc is resized to
     [nbNodes() x dim()] if the transformation is linear, \c x is
     not used at all
     */
    void hessianBasisAtPoint( uint16_type __idref,
                              hessian_basis_type& _hessian,
                              precompute_type const* __pc ) const
    {
        DCHECK( __pc ) << "a PreCompute must be set first before using this function:" << __idref;

        const auto nPts = static_cast<size_type>( _hessian.dimension( 0 ) );
        for ( size_type i = 0; i < nPts; ++i )
        {
            _hessian.chip( i, 0 ) = __pc->hessian( i, __idref ).chip( 0, 2 );
        }
    }

    /**
     * get an estimate of the radius of the element defined by G
     *
     * @param G matrix of nodes defining the element
     * @return the estimate of the radius of the element defined by G
     */
    value_type radiusEstimate( matrix_node_t_type const& G ) const
    {
        size_type N = G.size1();
        size_type Nref = this->nbNodes();
        // size_type n = ( this->isLinear()) ? 1 : Nref;

        // compute K
        matrix_node_t_type __g( Nref, dim() );
        // xref is not defined at this point: it serves only as a dummy point
        // in the case of linear inversion (__g is a constant matrix in this case)
        this->gradient( this->refNode( 0 ), __g );

        if ( __g.size1() != G.size2() )
            throw std::invalid_argument(
                "GeoMap::radiusEstimate(): incompatible geometry coefficients" );

        matrix_node_t_type K( G.size1(), __g.size2() );
        ublas::axpy_prod( G, __g, K );

        SVD<matrix_node_t_type> __svd( K );
        value_type __max;
        value_type __min;
        __svd.conditionNumber( __max, __min );
        return __max * sqrt( value_type( N ) ) / value_type( N );
    }

    bool M_is_cached;
    void initLinearGradientCache()
    {
        if ( !this->isLinear() )
            return;

        matrix_node_t_type __dummy_pts( ublas::zero_matrix<value_type>( nDim, 1 ) );
        ublas::vector<ublas::matrix<value_type>> m = super::derivate( __dummy_pts );

        CHECK_EQ( M_g_linear.size2(), m.size() )
            << "invalid affine gradient dimension";
        CHECK_EQ( m( 0 ).size2(), 1 )
            << "invalid affine gradient point count";
        CHECK_EQ( M_g_linear.size1(), m( 0 ).size1() )
            << "invalid affine gradient coefficient count";

        const auto nPts = static_cast<size_type>( this->nbPoints() );
        for ( size_type i = 0; i < nPts; ++i )
            for ( uint16_type n = 0; n < nDim; ++n )
                M_g_linear( i, n ) = m( n )( i, 0 );
    }

    bool isCached() const
    {
        return M_is_cached;
    }

    std::unordered_map<uint32_type, value_type> M_J;
    using element_face_pair_t = std::pair<uint32_type, uint16_type>;
    std::unordered_map<element_face_pair_t, value_type, boost::hash<element_face_pair_t>> M_JFace, M_n_norm;
    using allocator_vector_n_t = Eigen::aligned_allocator<std::pair<const element_face_pair_t, eigen_vector_n_type>>;
    std::unordered_map<element_face_pair_t, eigen_vector_n_type, boost::hash<element_face_pair_t>, std::equal_to<element_face_pair_t>, allocator_vector_n_t> M_un_real;
    std::unordered_map<element_face_pair_t, int, boost::hash<element_face_pair_t>> M_permutation_element_face_neighbor;
    unordered_map_eigen_matrix_type<uint32_type, nRealDim, nDim, value_type> M_K;
    unordered_map_eigen_matrix_type<uint32_type, nRealDim, nDim, value_type> M_B;

    template <typename MeshType>
    void initCache( MeshType const* mesh )
    {
#if 0
        size_type nelts = mesh->numElements();
        LOG( INFO ) << "[Geomap] start caching J,K,B for " << nelts << " elements\n";

        M_J.reserve( nelts );
        M_K.clear();
        M_K.reserve( nelts );
        M_B.reserve( nelts );
        M_JFace.reserve( nelts*mesh->numLocalTopologicalFaces() );
        M_n_norm.reserve( nelts*mesh->numLocalTopologicalFaces() );
        M_un_real.reserve( nelts*mesh->numLocalTopologicalFaces() );
        M_permutation_element_face_neighbor.reserve( nelts*mesh->numLocalTopologicalFaces() );
        M_is_cached = true;
#endif
    }
    bool isCacheValid() const
    {
        if ( !( M_J.size() > 0 &&
                M_K.size() > 0 &&
                M_B.size() > 0 ) )
        {
            return false;
        }

        return true;
    }
    //! @return true if geomap data are cache for element @p e, false otherwise
    bool cachedK( int e ) const
    {
        return M_K.count( e ) > 0;
    }
    bool cachedB( int e ) const
    {
        return M_B.count( e ) > 0;
    }
    bool cachedJ( int e ) const
    {
        return M_J.count( e ) > 0;
    }
    //! return jacobian at element @p e
    value_type J( int e ) const
    {
        return M_J.at( e );
    }
    void setJacobian( int e, value_type v )
    {
        M_J[e] = v;
    }
    bool hasJacobian( int e ) const
    {
        return M_J.count( e );
    }
    value_type jacobianAtFace( int e, int f ) const
    {
        return M_JFace.at( { e, f } );
    }
    eigen_matrix_type<nDim, nRealDim, value_type> const& B( int e ) const
    {
        return M_B.at( e );
    }
    eigen_matrix_type<nDim, nRealDim, value_type> const& K( int e ) const
    {
        return M_K.at( e );
    }
    void setK( int e, eigen_matrix_type<nDim, nRealDim, value_type> const& K )
    {
        M_K[e] = K;
    }
    template <typename G = self_type>
    bool cacheK( int e,
                 eigen_matrix_type<nRealDim, nDim, value_type>& K )
    {
        if constexpr ( G::is_linear && ( nDim == nRealDim ) )
        {
            if ( !cachedK( e ) ) return false;
            K.noalias() = M_K.at( e );
            return true;
        }
        else
            return false;
    }
    template <typename G = self_type>
    bool cacheB( int e,
                 eigen_matrix_type<nRealDim, nDim, value_type>& B )
    {
        if constexpr ( G::is_linear && ( nDim == nRealDim ) )
        {
            if ( !cachedB( e ) ) return false;
            B.noalias() = M_B.at( e );
            return true;
        }
        else
            return false;
    }
    template <typename G = self_type>
    bool cacheJ( int e, value_type& J )
    {
        if constexpr ( G::is_linear && ( nDim == nRealDim ) )
        {
            if ( !cachedJ( e ) ) return false;
            J = M_J.at( e );
            return true;
        }
        else
            return false;
    }
    template <typename G = self_type>
    void updateCacheK( int e,
                       eigen_matrix_type<nRealDim, nDim, value_type> const& K )
    {
        if constexpr ( G::is_linear && ( nDim == nRealDim ) )
        {
            M_K[e] = K;
        }
    }
    template <typename G = self_type>
    void updateCacheB( int e,
                       eigen_matrix_type<nRealDim, nDim, value_type> const& B )
    {
        if constexpr ( G::is_linear && ( nDim == nRealDim ) )
        {
            M_B[e] = B;
        }
    }
    template <typename G = self_type>
    void updateCacheJ( int e,
                       value_type const& J )
    {
        if constexpr ( G::is_linear && ( nDim == nRealDim ) )
        {
            M_J[e] = J;
        }
    }
    void addJ( int e, value_type v )
    {
        DCHECK( e >= 0 && e < M_J.size() ) << "invalid element id " << e << "( " << M_J.size() << ") for geomap cache, are you using the proper mesh\n";
        M_J[e] = v;
    }
    bool hasJacobianAtFace( int e, int f ) const
    {
        return M_JFace.count( { e, f } );
    }
    void setJacobianAtFace( int e, int f, value_type v )
    {
        M_JFace[{ e, f }] = v;
    }
    value_type normalNormAtFace( int e, int f ) const
    {
        return M_n_norm.at( { e, f } );
    }
    bool hasNormalNormAtFace( int e, int f ) const
    {
        return M_n_norm.count( { e, f } );
    }
    void setNormalNormAtFace( int e, int f, value_type v )
    {
        M_n_norm[{ e, f }] = v;
    }
    eigen_vector_n_type const& unitNormalAtFace( int e, int f ) const
    {
        return M_un_real.at( { e, f } );
    }
    bool hasUnitNormalAtFace( int e, int f ) const
    {
        return M_un_real.count( { e, f } );
    }
    void setUnitNormalAtFace( int e, int f, eigen_vector_n_type const& v )
    {
        M_un_real[{ e, f }] = v;
    }
    int permutationWithNeighborFace( int e, int f ) const
    {
        return M_permutation_element_face_neighbor.at( { e, f } );
    }
    bool hasPermutationWithNeighborFace( int e, int f ) const
    {
        return M_permutation_element_face_neighbor.count( { e, f } );
    }
    void setPermutationWithNeighborFace( int e, int f, int v )
    {
        M_permutation_element_face_neighbor[{ e, f }] = v;
    }
    template <typename Derived>
    void setB( int e, const Eigen::MatrixBase<Derived>& B )
    {
        M_B[e].noalias() = B;
    }

    /**
     * @class Context
     * @brief Pointwise differential geometry of a mapped element or subentity.
     *
     * A context combines an element geometry with precomputed reference basis
     * data and evaluates the physical points and requested differential
     * quantities used by integration and finite-element transformations.
     *
     * @tparam ElementType geometric entity supplying physical coefficients
     * @tparam SubEntityCoDim codimension of the evaluated subentity
     */
    template <typename ElementType, int SubEntityCoDim = 0>
    class Context
    {
      public:
        static const int subEntityCoDim = SubEntityCoDim;
        static const int subEntityCoDimFix = SubEntityCoDim > 0 ? SubEntityCoDim : 1;
        // reference space dimension
        static inline const uint16_type PDim = ElementType::nDim;
        // real space dimension
        static constexpr uint16_type NDim = ElementType::nRealDim;
        static inline const uint16_type nDim = NDim;
        /** Codimension of the mapped entity in physical space. */
        static constexpr uint16_type manifold_codimension = NDim-PDim;
        // type of transformation (linear or not)
        static const fem::transformation_type trans = geometric_mapping_type::trans;
        static inline const bool is_linear = geometric_mapping_type::is_linear;

        static inline const bool condition = ( ( PDim == NDim ) || ( ( NDim >= 1 ) && ( PDim == NDim - 1 ) ) );
        // BOOST_MPL_ASSERT_MSG( condition, INVALID_DIM, (mpl::int_<NDim>, mpl::int_<PDim>, ElementType ) );
        typedef typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim>>,
                                  mpl::identity<GeoMap<Dim, Order, NDim, T, Entity, PP>>,
                                  typename mpl::if_<mpl::and_<mpl::greater_equal<mpl::int_<NDim>, mpl::int_<1>>,
                                                              mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim - 1>>>,
                                                    // typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-1> >,
                                                    //  mpl::identity<typename GeoMap<Dim, Order, T, Entity, PP >::template face_gm<NDim>::type>,
                                                    mpl::identity<typename GeoMap<NDim, Order, NDim, T, Entity, PP>::face_gm_type>,
                                                    typename mpl::if_<mpl::and_<mpl::equal_to<mpl::int_<NDim>, mpl::int_<3>>,
                                                                                mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim - 2>>>,
                                                                      // typename mpl::if_<mpl::equal_to<mpl::int_<PDim>, mpl::int_<NDim-1> >,
                                                                      //  mpl::identity<typename GeoMap<Dim, Order, T, Entity, PP >::template face_gm<NDim>::type>,
                                                                      mpl::identity<typename GeoMap<NDim, Order, NDim, T, Entity, PP>::edge_gm_type>,
                                                                      mpl::identity<boost::none_t>>::type>::type>::type::type gm_type;
        typedef std::shared_ptr<gm_type> gm_ptrtype;

        typedef typename gm_type::value_type value_type;

        typedef typename gm_type::precompute_ptrtype precompute_ptrtype;

        typedef Context<ElementType, SubEntityCoDim> gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;

        typedef node_t_type normal_type;
        typedef ublas::vector<normal_type> normals_type;
        typedef typename normals_type::const_iterator normal_const_iterator;
        typedef node_t_type tangent_type;
        typedef ublas::vector<tangent_type> tangents_type;
        typedef typename tangents_type::const_iterator tangent_const_iterator;

        typedef ElementType element_type;
        // typedef typename element_type::permutation_type permutation_type;
        typedef typename element_type::template PermutationSubEntity<subEntityCoDimFix> permutation_type;
        using mesh_marker_type = typename element_type::marker_type;

        using eigen_matrix_nx_type = eigen_matrix_type<NDim, Eigen::Dynamic, value_type>;
        using eigen_matrix_xn_type = eigen_matrix_type<Eigen::Dynamic, NDim, value_type>;
        using eigen_matrix_px_type = eigen_matrix_type<PDim, Eigen::Dynamic, value_type>;
        using eigen_matrix_xp_type = eigen_matrix_type<Eigen::Dynamic, PDim, value_type>;

        using eigen_matrix_nn_type = eigen_matrix_type<NDim, NDim, value_type>;
        static constexpr int facet_tangent_dimension = PDim > 0 ? PDim - 1 : 0;
        using eigen_matrix_nt_type =
            eigen_matrix_type<NDim, facet_tangent_dimension, value_type>;
        using eigen_matrix_pt_type =
            eigen_matrix_type<PDim, facet_tangent_dimension, value_type>;
        using vector_eigen_matrix_pp_type = vector_eigen_matrix_type<PDim, PDim, value_type>;
        using vector_eigen_matrix_nn_type = vector_eigen_matrix_type<NDim, NDim, value_type>;
        using eigen_matrix_np_type = eigen_matrix_type<NDim, PDim, value_type>;
        using eigen_matrix_pp_type = eigen_matrix_type<PDim, PDim, value_type>;
        using vector_eigen_matrix_np_type = vector_eigen_matrix_type<NDim, PDim, value_type>;
        using derivative_component_batch_type =
            eigen_matrix_type<NDim, Eigen::Dynamic, value_type>;
        using gradient_component_batches_type =
            std::array<derivative_component_batch_type, PDim>;
        using hessian_component_batches_type =
            std::array<derivative_component_batch_type, PDim * PDim>;

        using hessian_type = tensor3_fixed_size_t<NDim, PDim, PDim, value_type>;
        using vector_hessian_type = vector_tensor3_fixed_size_t<NDim, PDim, PDim, value_type>;

        /**
         * @brief Pointwise quality diagnostics for a geometric mapping.
         *
         * Diagnostics are computed only when explicitly requested and are not
         * stored in the hot geometric context. Singular values refer to the
         * physical Jacobian matrix, while `jacobian` is the positive
         * integration measure and `signed_jacobian` detects orientation for
         * full-dimensional maps.
         */
        struct mapping_quality_type
        {
            /** Positive integration measure, equal to `sqrt(det(K^T K))`. */
            value_type jacobian = std::numeric_limits<value_type>::quiet_NaN();
            /** Signed determinant for full-dimensional maps; measure otherwise. */
            value_type signed_jacobian = std::numeric_limits<value_type>::quiet_NaN();
            /** Smallest singular value of the physical Jacobian matrix. */
            value_type minimum_singular_value = std::numeric_limits<value_type>::quiet_NaN();
            /** Largest singular value of the physical Jacobian matrix. */
            value_type maximum_singular_value = std::numeric_limits<value_type>::quiet_NaN();
            /** Spectral condition number of the physical Jacobian matrix. */
            value_type condition_number = std::numeric_limits<value_type>::infinity();
            /** Whether all mapping inputs used by the diagnostic are finite. */
            bool finite = false;
            /** Whether the local map preserves orientation and positive measure. */
            bool orientation_preserving = false;
            /** Whether the local map is finite, full-rank, and orientation preserving. */
            bool valid = false;
        };

        static constexpr bool is_on_face = subEntityCoDim == 1;

        Context() = default;

        template <int CTX>
        Context( gm_ptrtype __gm,
                 element_type const& __e,
                 std::integral_constant<int, CTX> const& statctx,
                 precompute_ptrtype const& __pc = precompute_ptrtype(),
                 uint16_type __f = invalid_uint16_type_value,
                 size_type dynctx = 0 )
            : M_gm( __gm ),
              M_element( nullptr ),
              M_pc( __pc ),
              M_pc_faces(),
              M_npoints( ( M_pc ) ? M_pc->nPoints() : 0 ),
              M_computed_point_stride( __gm->isLinear() ? 0 : 1 ),
              M_evaluation_strategy(
                  geometry::geometry_evaluation_strategy::automatic ),
              M_G( NDim, __gm->nbPoints() ),
              M_g( M_G.cols(), PDim ),
              M_hessian_basis_at_pt( M_G.cols(), PDim, PDim ),
              M_CS(),
              M_CSi(),
              M_B3( boost::extents[NDim][NDim][PDim][PDim] ),
              M_id( __e.id() ),
              M_e_markers( __e.markers() ),
              M_elem_id_1( invalid_v<size_type> ),             // __e.ad_first() ),
              M_pos_in_elem_id_1( invalid_uint16_type_value ), //__e.pos_first() ),
              M_elem_id_2( invalid_v<size_type> ),             //__e.ad_second() ),
              M_pos_in_elem_id_2( invalid_uint16_type_value ), //__e.pos_second() ),
              M_face_id( __f ),
              M_h( 0 ),
              M_h_min( 0 ),
              M_h_face( 0 ),
              M_meas( 0 ),
              M_measface( 0 ),
              M_perm(),
              M_dynamic_context( dynctx ),
              M_context_updated( 0 )
        {
            this->setElement( __e );

            if constexpr ( subEntityCoDim == 0 )
            {
                if ( M_pc )
                {
                    resize<CTX, true>();
                    update<CTX>( __e );
                }
            }
        }

        template <int CTX>
        Context( gm_ptrtype __gm,
                 element_type const& __e,
                 std::integral_constant<int, CTX> const& statctx,
                 std::vector<std::map<permutation_type, precompute_ptrtype>>& __pc,
                 uint16_type __f,
                 size_type dynctx = 0 )
            : Context( __gm, __e, statctx, __pc.empty() ? precompute_ptrtype{} : __pc[__f][__e.permutation( __f, mpl::int_<subEntityCoDimFix>() )], __f, dynctx )
        {
            if ( !__pc.empty() )
                M_pc_faces = __pc;
            if ( M_pc )
            {
                resize<CTX, true>();
                update<CTX>( __e, M_face_id );
            }
        }

        Context( gmc_ptrtype& p )
            : Context( *p )
        {
        }

        Context( Context const& ) = default;
        Context( Context&& ) = default;

        /**
         * clone the context
         */
        gmc_ptrtype clone()
        {
            return std::make_shared<gmc_type>( *this );
        }

        /**
         * update information on this context
         *
         *  - update the coordinate of the real points
         *  - update the pseudo-inverse of the gradient of the transformation
         *
         *  compute \f$ K(x_{\mathrm{ref}}) = G \nabla_{\mathrm{ref}} \phi(x_{\mathrm{ref}}) \f$
         *  where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
         *
         *  compute \f$ B(x_{\mathrm{ref}}) = K ( K^T K )^{-1} \f$
         *  where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
         */

        template <size_type CTX>
        void update( element_type const& __e )
            requires( SubEntityCoDim == 0 )
        {
            this->setElement( __e );
            this->updateImpl<CTX>();
        }

        /**
         * @brief Update from a compact mathematical entity view.
         * @tparam CTX compile-time geometric context flags
         * @tparam View view satisfying the geometric entity contract
         * @param view entity view supplying compact connectivity and coordinates
         */
        template <size_type CTX, GeometricEntityViewType View>
        void update( View const& view )
            requires(
                SubEntityCoDim == 0 &&
                std::same_as<
                    typename std::remove_cvref_t<View>::entity_type,
                    element_type> )
        {
            this->setElement( view );
            this->updateImpl<CTX>();
        }

        template <size_type CTX>
        void update( element_type const& __e, uint16_type __f )
            requires( SubEntityCoDim != 0 )
        {
            this->setElement( __e );
            this->setSubEntity( __f );
            this->updateImpl<CTX>();
        }

        /**
         * @brief Update a mapped subentity from a compact cell view.
         * @tparam CTX compile-time geometric context flags
         * @tparam View view satisfying the geometric entity contract
         * @param view cell view supplying compact geometry
         * @param __f local subentity number
         */
        template <size_type CTX, GeometricEntityViewType View>
        void update( View const& view, uint16_type __f )
            requires(
                SubEntityCoDim != 0 &&
                std::same_as<
                    typename std::remove_cvref_t<View>::entity_type,
                    element_type> )
        {
            this->setElement( view );
            this->setSubEntity( __f );
            this->updateImpl<CTX>();
        }

        template <size_type CTX>
        void update( element_type const& __e, uint16_type __f, permutation_type __perm )
            requires( SubEntityCoDim != 0 )
        {
            this->setElement( __e );
            this->setSubEntity( __f, __perm );
            this->updateImpl<CTX>();
        }

        /**
         * @brief Update a permuted subentity from a compact cell view.
         * @tparam CTX compile-time geometric context flags
         * @tparam View view satisfying the geometric entity contract
         * @param view cell view supplying compact geometry
         * @param __f local subentity number
         * @param __perm subentity permutation
         */
        template <size_type CTX, GeometricEntityViewType View>
        void update( View const& view, uint16_type __f,
                     permutation_type __perm )
            requires(
                SubEntityCoDim != 0 &&
                std::same_as<
                    typename std::remove_cvref_t<View>::entity_type,
                    element_type> )
        {
            this->setElement( view );
            this->setSubEntity( __f, __perm );
            this->updateImpl<CTX>();
        }

        /**
         * update information on this context
         *
         *  - update the coordinate of the real points
         *  - update the pseudo-inverse of the gradient of the transformation
         *
         *  compute \f$ K(x_{\mathrm{ref}}) = G \nabla_{\mathrm{ref}} \phi(x_{\mathrm{ref}}) \f$
         *  where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
         *
         *  compute \f$ B(x_{\mathrm{ref}}) = K ( K^T K )^{-1} \f$
         *  where \f$G \f$ is the matrix representing the geometric nodes nDof x dim
         */

        template <size_type CTX>
        void update( element_type const& __e,
                     precompute_ptrtype const& __pc,
                     uint16_type __f = invalid_uint16_type_value,
                     bool updatePc = true )
        {
            M_pc = __pc;

            if ( M_npoints != M_pc->nPoints() )
            {
                M_npoints = M_pc.get()->nPoints();
                resize<CTX>();
            }

            if constexpr ( subEntityCoDim > 0 )
                update<CTX>( __e, __f /*, updatePc*/ );
            else
                update<CTX>( __e );
        }

        template <size_type CTX>
        void updateContext( size_type dynctx = 0 )
        {
            M_dynamic_context = M_dynamic_context | dynctx;

            if ( M_dynamic_context != 0 )
            {
                this->resize<CTX | vm::DYNAMIC, true>();
                this->updateImpl<CTX | vm::DYNAMIC>();
            }
            else
            {
                this->resize<CTX, true>();
                this->updateImpl<CTX>();
            }
        }

        FEELPP_STRONG_INLINE
        void setElement( element_type const& __e )
        {
            element_type const* newElement = std::addressof( __e );
            bool isSameElt = newElement == M_element; // ( dynamic_cast<void const*>( newElement ) == dynamic_cast<void const*>( M_element ) )
            if ( !isSameElt )
            {
                M_context_updated = 0;
                M_element = std::addressof( __e );
                const bool use_vertices = ( M_G.cols() == element_type::numVertices );
                if constexpr ( !is_reference_convex_v<element_type> )
                {
                    if ( use_vertices )
                        __e.updateVertices( M_G );
                    else
                        __e.updateG( M_G );
                }
                else
                {
                    auto const& elementNodes = use_vertices ? __e.vertices() : __e.G();
                    if constexpr ( std::is_same_v<value_type, typename element_type::value_type> )
                        M_G = em_cmatrix_col_type<value_type>( elementNodes.data().begin(), M_G.rows(), M_G.cols() );
                    else
                        M_G = em_cmatrix_col_type<typename element_type::value_type>( elementNodes.data().begin(), M_G.rows(), M_G.cols() ).template cast<value_type>();
                }

                if ( !use_vertices )
                    M_gm->convertStoredGeometryToNodal(
                        M_G, geometric_mapping_type::geometryRepresentationOf( __e ) );

                M_id = __e.id();
                M_e_markers = __e.markers();
            }
        }

        /**
         * @brief Bind a compact entity view and gather its geometric nodes.
         * @tparam View view satisfying the geometric entity contract
         * @param view cell view backed by compact mesh storage
         */
        template <GeometricEntityViewType View>
        FEELPP_STRONG_INLINE void setElement( View const& view )
            requires std::same_as<
                typename std::remove_cvref_t<View>::entity_type,
                element_type>
        {
            auto const& element = view.entity();
            M_context_updated = 0;
            M_element = std::addressof( element );
            bool const useVertices =
                M_G.cols() == element_type::numVertices;
            if ( useVertices )
                view.updateVertexCoordinates( M_G );
            else
                view.updateCoordinates( M_G );
            if ( !useVertices )
                M_gm->convertStoredGeometryToNodal(
                    M_G, geometric_mapping_type::geometryRepresentationOf( view ) );
            M_id = view.id();
            M_e_markers = view.markers();
        }

        void setSubEntity( uint16_type __f )
            requires( SubEntityCoDim != 0 )
        {
            M_face_id = __f;
            CHECK( M_face_id != invalid_v<uint16_type> ) << "sub entity id not defined";
            M_perm = this->element().permutation( M_face_id, mpl::int_<subEntityCoDim>() );
            M_pc = M_pc_faces[__f][M_perm];
            M_f_markers = entityMarkers<subEntityCoDim>( this->element(), __f );
        }
        void setSubEntity( uint16_type __f, permutation_type __perm )
            requires( SubEntityCoDim != 0 )
        {
            M_face_id = __f;
            CHECK( M_face_id != invalid_v<uint16_type> ) << "sub entity id not defined";
            M_perm = __perm;
            M_pc = M_pc_faces[__f][M_perm];
            M_f_markers = entityMarkers<subEntityCoDim>( this->element(), __f );
        }

        ~Context()
        {
        }

        /** @name Accessors
         */
        //@{

        /**
         \return the geometric mapping associated with the context
         */
        gm_ptrtype const& geometricMapping() const
        {
            return M_gm;
        }

        /**
         \return the dimension of the space of the real element
         */
        uint16_type N() const
        {
            return NDim;
        }

        /**
         \return the dimension of the space of the reference element
         */
        uint16_type P() const
        {
            return PDim;
        }

        /**
         * \return the real element to which this context is associated
         */
        element_type const& element() const
        {
            return *M_element;
        }

        element_type const& element_c() const
        {
            return *M_element;
        }

        //!
        //! @return true if geomap associated to a face, false otherwise
        //!
        bool isOnFace() const
        {
            return ( subEntityCoDim == 1 ) && ( M_face_id != invalid_uint16_type_value );
        }
        //!
        //! @return true if geomap associated to a face, edge or point, false otherwise
        //!
        bool isOnSubEntity() const
        {
            return ( subEntityCoDim > 0 ) && ( M_face_id != invalid_uint16_type_value );
            // return subEntityCoDim > 1 || ( (subEntityCoDim == 1) && (M_face_id != invalid_uint16_type_value) );
        }

        //!
        //! @return the number of points to transfer from reference to real element
        //!
        uint16_type nPoints() const
        {
            return M_npoints;
        }

        //!
        //! in the case of linear transformation, the jacobian is constant and
        //! we don't have to compute the information at all points, only at one
        //!
        //! @return the number of points at which the geometric mapping is really computed.
        //!
        [[nodiscard]] uint16_type nComputedPoints() const noexcept
        {
            if ( M_npoints == 0 )
                return 0;
            return static_cast<uint16_type>(
                1 + M_computed_point_stride * ( M_npoints - 1 ) );
        }

        /**
         * @brief Return whether the context stores one constant derivative.
         * @return true for compile-time or runtime affine geometry
         */
        [[nodiscard]] bool usesAffineStorage() const noexcept
        {
            return M_computed_point_stride == 0;
        }

        /**
         * @brief Select the geometry derivative evaluation strategy.
         *
         * The strategy affects only how mathematically identical local
         * derivatives are evaluated. Affine mappings always use their
         * constant one-point kernel.
         *
         * @param strategy automatic, pointwise, or batched evaluation
         */
        void setGeometryEvaluationStrategy(
            geometry::geometry_evaluation_strategy strategy ) noexcept
        {
            M_evaluation_strategy = strategy;
        }

        /** @brief Return the selected geometry derivative strategy. */
        [[nodiscard]] geometry::geometry_evaluation_strategy
        geometryEvaluationStrategy() const noexcept
        {
            return M_evaluation_strategy;
        }

        /**
         * @brief Return contiguous AoS Jacobian matrices owned by the context.
         * @return non-owning view over computed-point matrices
         */
        [[nodiscard]] std::span<eigen_matrix_np_type const>
        jacobianMatrices() const noexcept
        {
            return { M_K.data(), M_K.size() };
        }

        /**
         * @brief Return contiguous AoS inverse-transpose Jacobian matrices.
         * @return non-owning view over computed-point matrices
         */
        [[nodiscard]] std::span<eigen_matrix_np_type const>
        inverseJacobianTransposeMatrices() const noexcept
        {
            return { M_B.data(), M_B.size() };
        }

        /**
         * @brief Return contiguous integration Jacobians.
         * @return non-owning view over computed-point scalar values
         */
        [[nodiscard]] std::span<value_type const>
        integrationJacobians() const noexcept
        {
            return { M_J.data(), M_J.size() };
        }

        /**
         * \return the set of points in the reference convex
         */
        matrix_node_t_type const& xRefs() const
        {
            return M_pc->nodes();
        }

        /**
         * \return the q-th point in the reference convex
         */
        ublas::matrix_column<matrix_node_t_type const> xRef( int q ) const
        {
            return M_pc->node( q );
        }

        /**
         * \return the node in the real element
         */
        // matrix_node_t_type const& xReal() const
        matrix_type const& xReal() const
        {
            // BOOST_STATIC_ASSERT( vm::has_point<context>::value );
            return M_xrealq;
        }

        /**
         * \return the node in the real element
         */
        ublas::matrix_column<matrix_type const> xReal( int q ) const
        {
            // BOOST_STATIC_ASSERT( vm::has_point<context>::value );
            return ublas::column( M_xrealq, q );
        }

        /**
         * Get the jacobian of the transformation at the \p q -th point
         *
         * \param q index of the point where the jacobian is requested
         * \return the jacobian of the transformation
         */
        [[nodiscard]] value_type J( int q ) const
        {
            return M_J[this->computedPointIndex( q )];
        }

        /** @brief Return the positive integration Jacobian at point `q`. */
        [[nodiscard]] value_type integrationJacobian( int q ) const
        {
            return this->J( q );
        }
        /**
         * \return the matrix associated with the geometric nodes
         */
        // matrix_node_t_type const& G() const { return M_G; }
        // matrix_type const& G() const
        eigen_matrix_nx_type const& G() const
        {
            return M_G;
        }

        /**
         * Get the inverse of the transformation at the \p i -th point
         *
         * \param i the index of point where the pseudo-inverse is requested
         * \return the pseudo inverse of the transformation
         */
        [[nodiscard]] eigen_matrix_np_type const& B( int q ) const
        {
            return M_B[this->computedPointIndex( q )];
        }

        /** @brief Return the inverse Jacobian transpose (or pseudoinverse transpose). */
        [[nodiscard]] eigen_matrix_np_type const& inverseJacobianTranspose( int q ) const
        {
            return this->B( q );
        }
        [[nodiscard]] value_type const& B( int c1, int c2, int q ) const
        {
            return M_B[this->computedPointIndex( q )]( c1, c2 );
        }

        [[nodiscard]] eigen_matrix_np_type const& localBasis( int q ) const
        {
            return M_local_basis_real[this->computedPointIndex( q )];
        }
        /**
         * @brief Return the oriented orthonormal facet frame at point `q`.
         *
         * Column zero is the unit normal (or embedded conormal); remaining
         * columns form the oriented orthonormal tangent frame. This is a
         * convenience frame and does not carry Piola functional scaling.
         * Its orientation is cell-local: adjacent cells share the tangent
         * subspace but may use different tangent columns and opposite outward
         * normals according to their facet permutation.
         */
        [[nodiscard]] eigen_matrix_np_type const& orthonormalFrame( int q ) const
        {
            return this->localBasis( q );
        }

        /**
         * @brief Return the oriented orthonormal tangent frame at point `q`.
         *
         * Columns follow the cell-local reference-facet parametrization. Use
         * `tangentProjector(q)` when a permutation-invariant tangent subspace
         * is required across an interface.
         *
         * @return matrix with one tangent column for a line facet and two for
         *         a three-dimensional cell facet
         */
        [[nodiscard]] eigen_matrix_nt_type orientedTangentFrame( int q ) const
        {
            return this->orthonormalFrame( q ).rightCols(
                facet_tangent_dimension );
        }
        [[nodiscard]] typename eigen_matrix_np_type::ColXpr const& basisN( int q ) const
        {
            return M_local_basis_real[this->computedPointIndex( q )].col( 0 );
        }

        [[nodiscard]] value_type const& localBasis( int c1, int c2, int q ) const
        {
            return M_local_basis_real[this->computedPointIndex( q )]( c1, c2 );
        }
        [[nodiscard]] value_type const& basisN( int c1, int c2, int q ) const
        {
            return M_local_basis_real[this->computedPointIndex( q )]( c1, 0 );
        }
        /**
         * @brief Return the orthogonal projector onto the active tangent space.
         *
         * On a facet context this is the facet tangent projector. On an
         * element context it is `K (K^T K)^{-1} K^T`, hence the manifold
         * tangent projector for an embedded cell and the identity for a
         * full-dimensional cell.
         */
        [[nodiscard]] eigen_matrix_nn_type tangentProjector( int q ) const
        {
            if constexpr ( subEntityCoDim == 1 )
            {
                if ( !M_Ptangent.empty() )
                    return M_Ptangent[this->computedPointIndex( q )];
                eigen_matrix_nt_type const mappedTangents =
                    this->K( q )*
                    this->referenceFacetTangentFrame( M_face_id );
                auto const tangentMetric =
                    geometry::metricTensor( mappedTangents );
                return ( mappedTangents*tangentMetric.inverse()*
                         mappedTangents.transpose() ).eval();
            }
            else
                return geometry::tangentProjector( this->K( q ),
                                                   this->B( q ) );
        }
        /** @brief Return one component of the tangent projector. */
        [[nodiscard]] value_type tangentProjector( int c1, int c2,
                                                   int q ) const
        {
            if constexpr ( subEntityCoDim == 1 )
                return this->tangentProjector( q )( c1, c2 );
            else
                return this->tangentProjector( q )( c1, c2 );
        }

        /** @brief Return the first fundamental form `g = K^T K`. */
        [[nodiscard]] eigen_matrix_pp_type metricTensor( int q ) const
        {
            return geometry::metricTensor( this->K( q ) );
        }

        /** @brief Return the inverse first fundamental form `g^{-1}`. */
        [[nodiscard]] eigen_matrix_pp_type inverseMetricTensor( int q ) const
        {
            return geometry::inverseMetricTensor( this->B( q ) );
        }

        /** @brief Return the positive manifold Jacobian `sqrt(det(g))`. */
        [[nodiscard]] value_type manifoldJacobian( int q ) const
        {
            return geometry::manifoldJacobian( this->K( q ) );
        }

        /** @brief Return the orthogonal projector onto the normal space. */
        [[nodiscard]] eigen_matrix_nn_type normalProjector( int q ) const
        {
            return geometry::normalProjectorFromTangent(
                this->tangentProjector( q ) );
        }

        /** @deprecated Use tangentProjector(q). */
        [[nodiscard]] eigen_matrix_nn_type projectorTangent( int q ) const
        {
            return this->tangentProjector( q );
        }

        /** @deprecated Use tangentProjector(c1,c2,q). */
        [[nodiscard]] value_type projectorTangent( int c1, int c2,
                                                   int q ) const
        {
            return this->tangentProjector( c1, c2, q );
        }

        [[nodiscard]] eigen_matrix_np_type const& K( int q ) const
        {
            return M_K[this->computedPointIndex( q )];
        }

        /** @brief Return the physical Jacobian matrix at point `q`. */
        [[nodiscard]] eigen_matrix_np_type const& jacobianMatrix( int q ) const
        {
            return this->K( q );
        }
        [[nodiscard]] vector_eigen_matrix_np_type const& K() const
            requires( !is_linear )
        {
            return M_K;
        }
        [[nodiscard]] value_type const& K( int c1, int c2, int q ) const
        {
            return M_K[this->computedPointIndex( q )]( c1, c2 );
        }

        [[nodiscard]] hessian_type const& hessian( int q ) const
        {
            this->ensureMappingHessian( q );
            return M_hessian[this->computedPointIndex( q )];
        }

        /** Analytic d^3 F_component / (d xi_first d xi_second d xi_third). */
        [[nodiscard]] value_type thirdDerivative(
            int component, int first, int second, int third, int q ) const
        {
            if ( component < 0 || component >= NDim || first < 0 || first >= PDim ||
                 second < 0 || second >= PDim || third < 0 || third >= PDim )
                throw std::out_of_range( "Invalid geometry third derivative component" );
            this->ensureMappingThirdDerivatives();
            return M_third_derivative_component_batches[( first * PDim + second ) * PDim + third](
                component, this->computedPointIndex( q ) );
        }

        /**
         * the tensor of rank 4 for the transformation of 2nd
         * order derivatives from the reference element to the real
         * element. The tensor has the shape \f$[N,N,P,P]\f$.
         *
         * \return the tensor of rank 4
         */
        boost::multi_array<value_type, 4> const& B3() const
        {
            return M_B3;
        }

        /**
         * \return the barycenter of the reference nodes
         */
        node_t_type barycenterRef() const
        {
            node_t_type __barycenter( M_gm->dim() );
            __barycenter.clear();

            for ( uint16_type __c = 0; __c < M_gm->refNodes().size(); ++__c )
            {
                __barycenter += M_gm->refNode( __c );
            }

            __barycenter /= M_gm->refNodes().size();
            return __barycenter;
        }

        /**
         \return the barycenter of the geometric nodes
         */
        eigen_vector_type<NDim, value_type> barycenterReal() const
        {
            int nPts = M_G.cols();
            return ( 1. / nPts ) * M_G * eigen_vector_type<Eigen::Dynamic, value_type>::Ones( nPts );
        }
        /**
         * tell whether the point is on the surface of the convex
         * @return true if the point is on the surface, false otherwise
         */
        bool isOnConvexSurface() const
        {
            if constexpr ( is_linear )
            {
#if 0
                // x -x0 - K(0)\bar{x}
                return std::abs( ublas::norm_2( xReal() - ublas::column( M_G, 0 ) -
                                                ublas::prod( M_K, xRef() ) ) ) < 1e-10;
#else
                CHECK( false ) << "TODO";
                return false;
#endif
            }
            else
                return false;
        }

        node_t_type const& refNormal( int /*q*/ ) const
        {
            return M_gm->referenceConvex().normal( M_face_id );
        }

        /**
         * get the norm_2 of normal of the real element
         *
         * @return the norm_2 of the normal of the real element
         */
        [[nodiscard]] value_type normalNorm( int q ) const
        {
            return M_normal_norms[this->computedPointIndex( q )];
        }

        /**
         * get the normal of the real element
         *
         * @return the normal of the real element
         */
        [[nodiscard]] eigen_vector_n_type const& normal( int q ) const
        {
            return M_normals[this->computedPointIndex( q )];
        }

        //!
        //! @return the unit normal of the real element
        //!
        [[nodiscard]] vector_eigen_vector_n_type const& unitNormal() const
            requires( is_linear && subEntityCoDim == 1 )
        {
            return M_unit_normals;
        }

        //!
        //! @return the unit_normal at point @c q
        //!
        [[nodiscard]] decltype(auto) unitNormal( int q ) const
            requires( subEntityCoDim == 1 ||
                      ( subEntityCoDim == 0 && NDim == PDim+1 &&
                        NDim <= 3 ) )
        {
            if constexpr ( subEntityCoDim == 1 )
                return ( M_unit_normals[this->computedPointIndex( q )] );
            else
                return geometry::canonicalUnitNormal( this->K( q ) );
        }

        /**
         * @brief Reject a normal evaluation when a cell has no unique normal.
         *
         * Generic DSEL tensors compile the normal-symbol evaluation path even
         * when an expression contains no normal symbol. Keeping this overload
         * in the cell-context interface preserves that generic compilation
         * path while reporting an immediate mathematical error if a normal is
         * actually requested on a full-dimensional cell or an embedded cell
         * of codimension greater than one.
         *
         * @param q geometric evaluation-point index
         * @return unreachable zero vector used only to satisfy the return type
         */
        [[nodiscard]] eigen_vector_n_type const&
        unitNormal( int q ) const
            requires(
                subEntityCoDim == 0 &&
                !( NDim == PDim + 1 && NDim <= 3 ) )
        {
            CHECK( false )
                << "this geometric cell has no unique canonical unit normal; "
                   "evaluate on a codimension-one range or use the normal "
                   "frame/projector for higher-codimension geometry";
            static eigen_vector_n_type const unreachableNormal =
                eigen_vector_n_type::Zero();
            (void)q;
            return unreachableNormal;
        }

        [[nodiscard]] value_type unitNormal( int n, int q ) const
            requires( subEntityCoDim == 1 ||
                      ( subEntityCoDim == 0 && NDim == PDim+1 &&
                        NDim <= 3 ) )
        {
            return this->unitNormal( q )( n );
        }

        /**
         * @brief Reject a normal-component evaluation on a cell interior.
         * @param n ambient component
         * @param q geometric evaluation-point index
         * @return unreachable scalar used only to satisfy the return type
         */
        [[nodiscard]] value_type unitNormal( int n, int q ) const
            requires(
                subEntityCoDim == 0 &&
                !( NDim == PDim + 1 && NDim <= 3 ) )
        {
            return this->unitNormal( q )( n );
        }

        /**
         * @brief Return an orthonormal frame of the embedded normal space.
         *
         * The default frame is canonical only in codimension one. In higher
         * codimension it is a deterministic QR frame; algorithms should use
         * `normalProjector()` when their result must be frame independent.
         *
         * @tparam Policy normal-frame construction policy
         * @param q geometric evaluation-point index
         * @param policy explicit frame policy
         * @return normal vectors stored by columns
         */
        template <typename Policy = std::conditional_t<
                      manifold_codimension == 1,
                      geometry::canonical_normal_frame_policy,
                      geometry::qr_normal_frame_policy>>
        [[nodiscard]] auto normalFrame( int q, Policy policy = {} ) const
            requires( subEntityCoDim == 0 && PDim < NDim )
        {
            return geometry::normalFrame( this->K( q ), std::move( policy ) );
        }

        /**
         * @brief Return the codimension-one second fundamental form.
         *
         * The sign convention is `b_ij = n . d_ij F`.
         */
        [[nodiscard]] eigen_matrix_pp_type secondFundamentalForm( int q ) const
            requires( subEntityCoDim == 0 && NDim == PDim+1 &&
                      NDim <= 3 )
        {
            auto const normalValue = this->unitNormal( q );
            return geometry::secondFundamentalForm(
                this->K( q ), normalValue, this->hessian( q ) );
        }

        /** @brief Return the codimension-one Weingarten shape operator. */
        [[nodiscard]] eigen_matrix_pp_type shapeOperator( int q ) const
            requires( subEntityCoDim == 0 && NDim == PDim+1 &&
                      NDim <= 3 )
        {
            return geometry::shapeOperator( this->inverseMetricTensor( q ),
                                            this->secondFundamentalForm( q ) );
        }

        /** @brief Return averaged signed mean curvature `trace(S)/PDim`. */
        [[nodiscard]] value_type meanCurvature( int q ) const
            requires( subEntityCoDim == 0 && NDim == PDim+1 &&
                      NDim <= 3 )
        {
            return geometry::meanCurvature( this->shapeOperator( q ) );
        }

        /** @brief Return Gaussian curvature of an embedded surface in R3. */
        [[nodiscard]] value_type gaussianCurvature( int q ) const
            requires( subEntityCoDim == 0 && PDim == 2 && NDim == 3 )
        {
            return geometry::gaussianCurvature( this->shapeOperator( q ) );
        }

        /**
         * @brief Evaluate the scalar Laplace--Beltrami operator from reference derivatives.
         * @tparam GradientDerived Eigen reference-gradient expression
         * @tparam HessianDerived Eigen reference-Hessian expression
         * @param q geometric evaluation-point index
         * @param referenceGradient scalar gradient in reference coordinates
         * @param referenceHessian scalar Hessian in reference coordinates
         * @return intrinsic scalar Laplacian
         */
        template <typename GradientDerived, typename HessianDerived>
        [[nodiscard]] value_type laplaceBeltrami(
            int q,
            Eigen::MatrixBase<GradientDerived> const& referenceGradient,
            Eigen::MatrixBase<HessianDerived> const& referenceHessian ) const
            requires( subEntityCoDim == 0 && PDim < NDim )
        {
            return geometry::laplaceBeltrami(
                this->K( q ), this->inverseMetricTensor( q ),
                this->hessian( q ), referenceGradient, referenceHessian );
        }

        //!
        //! @return the scaled tangent to the current face
        //!
        [[nodiscard]] eigen_vector_n_type const& tangent() const
            requires( is_linear && subEntityCoDim == 1 && PDim == 2 )
        {
            return M_tangents[0];
        }

        /**
         * @brief scaled tangent at point

         * @param q index of the point at which the tangent is evaluated
         * @return the scaled tangent at the point
         */
        [[nodiscard]] eigen_vector_n_type const& tangent( int q ) const
            requires( subEntityCoDim == 1 && PDim == 2 )
        {
            return M_tangents[this->computedPointIndex( q )];
        }

        /**
         * @brief Return the mapped tangent of a one-dimensional cell.
         *
         * The reference line has one canonical coordinate direction, so the
         * first Jacobian column is its oriented physical tangent. This covers
         * both ordinary 1D meshes and curves embedded in a higher-dimensional
         * ambient space.
         *
         * @param q geometric evaluation-point index
         * @return oriented scaled tangent
         */
        [[nodiscard]] eigen_vector_n_type tangent( int q ) const
            requires( subEntityCoDim == 0 && PDim == 1 )
        {
            return this->K( q ).col( 0 );
        }

        [[nodiscard]] eigen_vector_n_type const& unitTangent( int q ) const
            requires( subEntityCoDim == 1 && PDim == 2 )
        {
            return M_unit_tangents[this->computedPointIndex( q )];
        }

        /**
         * @brief Return the oriented unit tangent of a one-dimensional cell.
         * @param q geometric evaluation-point index
         * @return normalized first Jacobian column
         */
        [[nodiscard]] eigen_vector_n_type unitTangent( int q ) const
            requires( subEntityCoDim == 0 && PDim == 1 )
        {
            auto mappedTangent = this->tangent( q );
            auto const tangentNorm = mappedTangent.norm();
            CHECK_GT( tangentNorm, value_type{} )
                << "degenerate one-dimensional geometric mapping";
            return mappedTangent / tangentNorm;
        }

        [[nodiscard]] value_type const& unitTangent( int c1, int c2, int q ) const
            requires( subEntityCoDim == 1 && PDim == 2 )
        {
            return M_unit_tangents[this->computedPointIndex( q )]( c1, c2 );
        }

        /**
         * get the id of the element
         *
         * @return the id of the element
         */
        size_type id() const
        {
            return M_id;
        }

        /*
         * \return the face id
         */
        uint16_type faceId() const
        {
            return M_face_id;
        }

        /**
         * \return true if the element is a face, false otherwise
         */
        bool elementIsAFace() const
        {
            return M_face_id != invalid_uint16_type_value;
        }

        /**
         * get the marker of the element
         *
         * @return the marker of the element
         */
        mesh_marker_type marker( uint16_type k ) const
        {
            auto itFindMarker = M_e_markers.find( k );
            if ( itFindMarker != M_e_markers.end() )
                return itFindMarker->second;
            else
                return mesh_marker_type{};
        }
        /**
         * get the marker of the element
         *
         * @return the marker of the element
         */
        mesh_marker_type marker() const
        {
            auto itFindMarker = M_e_markers.find( 1 );
            if ( itFindMarker != M_e_markers.end() )
                return itFindMarker->second;
            else
                return mesh_marker_type{};
        }

        /**
         * get the marker2 of the element
         *
         * @return the marker2 of the element
         */
        mesh_marker_type marker2() const
        {
            auto itFindMarker = M_e_markers.find( 2 );
            if ( itFindMarker != M_e_markers.end() )
                return itFindMarker->second;
            else
                return mesh_marker_type{};
        }

        /**
         * get the marker3 of the element
         *
         * @return the marker3 of the element
         */
        mesh_marker_type marker3() const
        {
            auto itFindMarker = M_e_markers.find( 3 );
            if ( itFindMarker != M_e_markers.end() )
                return itFindMarker->second;
            else
                return mesh_marker_type{};
        }

        /**
         * get the marker of the face of the element
         *
         * @return the marker of the face of the  element
         */
        mesh_marker_type entityMarker( uint16_type k = 1 ) const
        {
            if ( !isOnSubEntity() || !M_f_markers )
                return mesh_marker_type{};
            auto itFindMarker = M_f_markers->find( k );
            if ( itFindMarker != M_f_markers->end() )
                return itFindMarker->second;
            else
                return mesh_marker_type{};
        }

        /**
         * get the id of the first element containing the element
         *
         * @return the id of the first element containing the element
         */
        size_type id1() const
        {
            return M_elem_id_1;
        }

        /**
         * get the local id of the element in the first element containing the element
         *
         * @return the local id of element in the first element containing the element
         */
        uint16_type idIn1() const
        {
            return M_pos_in_elem_id_1;
        }

        /**
         * get the id of the second element containing the element
         *
         * @return the id of the second element containing the element
         */
        size_type id2() const
        {
            return M_elem_id_2;
        }

        /**
         * get the local id of the element in the second element containing the element
         *
         * @return the local id of element in the second element containing the element
         */
        uint16_type idIn2() const
        {
            return M_pos_in_elem_id_2;
        }

        /**
         * get an estimate of the radius of the current element
         *
         * @return the radius estimate of the current element
         */
        value_type radiusEstimate() const
        {
#if 0
            return M_gm->radiusEstimate( M_G );
#else
            CHECK( false ) << "TODO";
            return value_type( 0 );
#endif
        }

        /**
         * @return the max length of the edges of the element
         */
        value_type h() const
        {
            return M_h;
        }
        /**
         * @return the min length of the edges of the element
         */
        value_type hMin() const
        {
            return M_h_min;
        }
        /**
         * Get max length of the edge of the face of the element
         *
         * @return the max length of the edge of the face of the element
         */
        value_type hFace() const
        {
            return M_h_face;
        }

        /*
         * @return the measure of the element
         */
        value_type meas() const
        {
            return M_meas;
        }

        /*
         * @return the measure of the set of elements which share a vertex with \p M_elements including himself
         */
        value_type measurePointElementNeighbors() const
        {
            return M_element->measurePointElementNeighbors();
        }

        /*
         * @return the measure of the (current) face of the element
         */
        value_type measFace() const
        {
            return M_measface;
        }

        /**
         * \return the permutation associated with the face
         */
        permutation_type permutation() const
        {
            return permutation( mpl::bool_<( nDim >= 2 )>() );
        }

        /**
         * \return the precompute type
         */
        precompute_ptrtype const& pc() const
        {
            return M_pc;
        }

        /**
         * \return the precompute type for the faces
         */
        std::vector<std::map<permutation_type, precompute_ptrtype>> const& pcFaces() const
        {
            return M_pc_faces;
        }

        //! @return the dynamic context associated to the geomap
        size_type dynamicContext() const { return M_dynamic_context; }

        //@}

        /** @name  Mutators
         */
        //@{

        /**
     set some precomputed data on the reference element
     */
        void setPc( precompute_ptrtype const& __pc )
        {
            M_pc = __pc;
        }

        void setPcFaces( std::vector<std::map<permutation_type, precompute_ptrtype>> const& __pcfaces )
        {
            M_pc_faces = __pcfaces;
        }

        /**
         * @brief Compute an edge tangent at a geometric evaluation point.
         *
         * With `scaled == true`, this method returns the covariant-Piola test
         * direction `K(q) t_hat`, where `t_hat` is the oriented reference-edge
         * tangent.  Otherwise it returns the corresponding physical unit
         * tangent.  The pointwise form is required on non-affine cells.
         *
         * @param edgeId local edge index in the element
         * @param q geometric evaluation-point index
         * @param t output tangent
         * @param scaled preserve the reference-edge parametrization when true
         */
        void edgeTangent( int edgeId, int q, ublas::vector<value_type>& t,
                          bool scaled = false ) const
        {
            if ( scaled )
                this->piolaCovariantEdgeTangent( edgeId, q, t );
            else
                this->physicalUnitEdgeTangent( edgeId, q, t );
        }

        /**
         * @brief Return the physical unit tangent of an edge at a point.
         * @param edgeId local edge index in the element
         * @param q geometric evaluation-point index
         * @param tangent output unit tangent
         */
        void physicalUnitEdgeTangent( int edgeId, int q,
                                      ublas::vector<value_type>& tangent ) const
        {
            eigenMap( tangent ) = this->physicalUnitEdgeTangent( edgeId, q );
        }

        /**
         * @brief Return the physical unit tangent of an edge at a point.
         * @param edgeId local edge index
         * @param q geometric evaluation-point index
         * @return physical unit edge tangent
         */
        [[nodiscard]] eigen_vector_n_type
        physicalUnitEdgeTangent( int edgeId, int q ) const
        {
            eigen_vector_n_type tangent =
                this->piolaCovariantEdgeTangent( edgeId, q );
            tangent.normalize();
            return tangent;
        }

        /**
         * @brief Return the covariant-Piola edge test direction `K(q) t_hat`.
         * @param edgeId local edge index in the element
         * @param q geometric evaluation-point index
         * @param tangent output tangent carrying the reference parametrization
         */
        void covariantPiolaTangent( int edgeId, int q,
                                    ublas::vector<value_type>& tangent ) const
        {
            this->piolaCovariantEdgeTangent( edgeId, q, tangent );
        }

        /**
         * @brief Store the covariant Piola edge trace direction.
         * @param edgeId local edge index
         * @param q geometric evaluation-point index
         * @param tangent output direction `K(q) tHat`
         */
        void piolaCovariantEdgeTangent(
            int edgeId, int q, ublas::vector<value_type>& tangent ) const
        {
            eigenMap( tangent ) =
                this->piolaCovariantEdgeTangent( edgeId, q );
        }

        /**
         * @brief Return the covariant Piola edge trace direction.
         * @param edgeId local edge index
         * @param q geometric evaluation-point index
         * @return `K(q) tHat`, preserving reference-edge parametrization
         */
        [[nodiscard]] eigen_vector_n_type
        piolaCovariantEdgeTangent( int edgeId, int q ) const
        {
            return geometry::covariantPiolaTestDirection(
                this->K( q ), M_gm->referenceConvex().tangent( edgeId ) );
        }

        /**
         * @brief Compute an edge tangent using the first geometric point.
         *
         * This compatibility overload is exact for affine geometric maps.
         * Point-supported finite-element functionals on curved cells must use
         * the overload taking an explicit evaluation-point index.
         *
         * @param edgeId local edge index in the element
         * @param t output tangent
         * @param scaled preserve the reference-edge parametrization when true
         */
        void
        edgeTangent( int edgeId, ublas::vector<value_type>& t, bool scaled = false ) const
        {
            this->edgeTangent( edgeId, 0, t, scaled );
        }

        /**
         * @brief Compute a facet normal at a geometric evaluation point.
         *
         * With `scaled == true`, this method returns the contravariant-Piola
         * test direction
         *
         * `J(q) K(q)^{-T} n_hat |f_hat|`,
         *
         * where `n_hat` is the unit reference normal and `|f_hat|` is the
         * reference-facet measure.  Otherwise it returns the physical unit
         * normal.  The pointwise cofactor is essential on curved cells.
         *
         * @param faceId local facet index in the element
         * @param q geometric evaluation-point index
         * @param n output normal
         * @param scaled return the Piola-scaled normal when true
         */
        void faceNormal( int faceId, int q, ublas::vector<value_type>& n,
                         bool scaled = false ) const
        {
            if ( scaled )
                this->piolaContravariantFacetNormal( faceId, q, n );
            else
                this->physicalUnitFacetNormal( faceId, q, n );
        }

        /**
         * @brief Return the physical unit normal of a facet at a point.
         * @param faceId local facet index in the element
         * @param q geometric evaluation-point index
         * @param normal output unit normal
         */
        void physicalUnitFacetNormal( int faceId, int q,
                                      ublas::vector<value_type>& normal ) const
        {
            eigenMap( normal ) = this->physicalUnitFacetNormal( faceId, q );
        }

        /**
         * @brief Return the physical unit facet normal at a point.
         * @param faceId local facet index
         * @param q geometric evaluation-point index
         * @return unit normal or embedded conormal
         */
        [[nodiscard]] eigen_vector_n_type
        physicalUnitFacetNormal( int faceId, int q ) const
        {
            eigen_vector_n_type normal =
                this->B( q ) * M_gm->referenceConvex().normal( faceId );
            normal.normalize();
            return normal;
        }

        /**
         * @brief Return the contravariant-Piola facet test direction.
         *
         * The result is `J(q) K(q)^{-T} n_hat |f_hat|`, matching the reference
         * normal-flux functional without normalization or a cell-wide facet
         * measure approximation.
         *
         * @param faceId local facet index in the element
         * @param q geometric evaluation-point index
         * @param normal output scaled normal
         */
        void contravariantPiolaNormal( int faceId, int q,
                                       ublas::vector<value_type>& normal ) const
        {
            eigenMap( normal ) =
                this->piolaContravariantFacetNormal( faceId, q );
        }

        /**
         * @brief Store the contravariant Piola facet trace direction.
         * @param faceId local facet index
         * @param q geometric evaluation-point index
         * @param normal output direction `J K^{-T} nHat |fHat|`
         */
        void piolaContravariantFacetNormal(
            int faceId, int q, ublas::vector<value_type>& normal ) const
        {
            eigenMap( normal ) =
                this->piolaContravariantFacetNormal( faceId, q );
        }

        /**
         * @brief Return the contravariant Piola facet trace direction.
         * @param faceId local facet index
         * @param q geometric evaluation-point index
         * @return `J K^{-T} nHat |fHat|`
         */
        [[nodiscard]] eigen_vector_n_type
        piolaContravariantFacetNormal( int faceId, int q ) const
        {
            return geometry::contravariantPiolaTestDirection(
                this->B( q ), this->integrationJacobian( q ),
                M_gm->referenceConvex().normal( faceId ),
                this->referenceFacetMeasure( faceId ) );
        }

        /**
         * @brief Return the measure of a reference facet.
         * @param faceId local facet index
         * @return reference length in 2D, area in 3D, or one in 1D
         */
        [[nodiscard]] value_type referenceFacetMeasure( int faceId ) const
        {
            if constexpr ( PDim == 1 )
                return value_type( 1 );
            else if constexpr ( PDim == 2 )
                return this->referenceFacetTangentFrame( faceId ).col( 0 ).norm();
            else
            {
                auto const vertices =
                    M_gm->referenceConvex().faceVertices( faceId );
                auto const tangents =
                    this->referenceFacetTangentFrame( faceId );
                value_type const parallelogramMeasure =
                    tangents.col( 0 ).cross( tangents.col( 1 ) ).norm();

                // A simplex facet is a triangle, whereas the currently
                // supported tensor-product facet is a parallelogram.
                return vertices.size2() == 3
                           ? value_type( 0.5 ) * parallelogramMeasure
                           : parallelogramMeasure;
            }
        }

        /**
         * @brief Return a direct oriented reference-facet tangent frame.
         *
         * Columns preserve the canonical reference-vertex parametrization.
         * The last tangent is flipped when necessary so that
         * `[nHat,tHat_0,...]` is direct. This is a cell-local convention;
         * interelement conformity remains governed by the entity permutation
         * and finite-element orientation transform.
         *
         * @param faceId local facet index
         * @return reference tangent directions by columns
         */
        [[nodiscard]] eigen_matrix_pt_type
        referenceFacetTangentFrame( int faceId ) const
        {
            eigen_matrix_pt_type tangents;
            if constexpr ( facet_tangent_dimension > 0 )
            {
                auto const vertices =
                    M_gm->referenceConvex().faceVertices( faceId );
                for ( int component = 0; component < PDim; ++component )
                    tangents( component, 0 ) =
                        vertices( component, 1 ) - vertices( component, 0 );

                if constexpr ( facet_tangent_dimension == 2 )
                {
                    value_type maximumArea = value_type{};
                    tangents.col( 1 ).setZero();
                    for ( std::size_t vertex = 2; vertex < vertices.size2();
                          ++vertex )
                    {
                        eigen_vector_p_type candidate;
                        for ( int component = 0; component < PDim; ++component )
                            candidate( component ) =
                                vertices( component, vertex ) -
                                vertices( component, 0 );
                        value_type const area =
                            tangents.col( 0 ).cross( candidate ).squaredNorm();
                        if ( area > maximumArea )
                        {
                            maximumArea = area;
                            tangents.col( 1 ) = candidate;
                        }
                    }
                }

                eigen_matrix_pp_type orientation;
                orientation.col( 0 ) =
                    M_gm->referenceConvex().normal( faceId );
                orientation.rightCols( facet_tangent_dimension ) = tangents;
                if ( orientation.determinant() < value_type{} )
                    tangents.col( facet_tangent_dimension - 1 ) *=
                        value_type( -1 );
            }
            return tangents;
        }

        /**
         * @brief Return the pointwise mixed Piola trace frame of a facet.
         *
         * The first column is the H(div) normal-flux direction and remaining
         * columns are the H(curl) tangential directions. This frame preserves
         * reference functional scaling and is distinct from
         * `orthonormalFrame()`.
         *
         * @param faceId local facet index
         * @param q geometric evaluation-point index
         * @return mixed contravariant/covariant Piola frame
         */
        [[nodiscard]] eigen_matrix_np_type
        piolaTraceFrame( int faceId, int q ) const
        {
            return geometry::piolaTraceFrame(
                this->K( q ), this->B( q ), this->integrationJacobian( q ),
                M_gm->referenceConvex().normal( faceId ),
                this->referenceFacetTangentFrame( faceId ),
                this->referenceFacetMeasure( faceId ) );
        }

        /**
         * @brief Return the signed Jacobian determinant at a point.
         *
         * For full-dimensional maps this detects orientation reversal, unlike
         * `J(q)`, which is the positive integration measure.  Embedded maps
         * have no square determinant and return the positive metric Jacobian.
         *
         * @param q geometric evaluation-point index
         * @return signed determinant or positive metric Jacobian
         */
        [[nodiscard]] value_type signedJacobian( int q ) const
        {
            if constexpr ( NDim == PDim )
                return this->K( q ).determinant();
            else
                return this->J( q );
        }

        /**
         * @brief Compute pointwise rank, conditioning, orientation, and validity.
         * @param q geometric evaluation-point index
         * @param relativeTolerance relative singular-value and Jacobian tolerance
         * @return opt-in pointwise mapping diagnostics
         */
        [[nodiscard]] mapping_quality_type
        mappingQuality( int q,
                        value_type relativeTolerance =
                            value_type( 64 ) * std::numeric_limits<value_type>::epsilon() ) const
        {
            mapping_quality_type result;
            auto const& jacobianMatrix = this->K( q );
            result.jacobian = this->J( q );
            result.signed_jacobian = this->signedJacobian( q );
            result.finite = jacobianMatrix.allFinite() &&
                            std::isfinite( result.jacobian ) &&
                            std::isfinite( result.signed_jacobian );
            if ( !result.finite )
                return result;

            Eigen::JacobiSVD<eigen_matrix_np_type> svd( jacobianMatrix );
            auto const& singularValues = svd.singularValues();
            result.maximum_singular_value = singularValues.maxCoeff();
            result.minimum_singular_value = singularValues.minCoeff();
            if ( result.minimum_singular_value > value_type( 0 ) )
                result.condition_number = result.maximum_singular_value /
                                          result.minimum_singular_value;

            const value_type singularTolerance = relativeTolerance *
                                                 result.maximum_singular_value;
            const value_type jacobianScale = std::pow(
                result.maximum_singular_value, static_cast<int>( PDim ) );
            const value_type jacobianTolerance = relativeTolerance * jacobianScale;
            if constexpr ( NDim == PDim )
                result.orientation_preserving = result.signed_jacobian > jacobianTolerance;
            else
                result.orientation_preserving = result.jacobian > jacobianTolerance;
            result.valid = result.minimum_singular_value > singularTolerance &&
                           result.orientation_preserving;
            return result;
        }

        /** @brief Return the smallest singular value of the Jacobian at `q`. */
        [[nodiscard]] value_type minimumSingularValue( int q ) const
        {
            return this->mappingQuality( q ).minimum_singular_value;
        }

        /** @brief Return the Jacobian condition number at `q`. */
        [[nodiscard]] value_type conditionNumber( int q ) const
        {
            return this->mappingQuality( q ).condition_number;
        }

        /** @brief Return the minimum positive integration Jacobian over all points. */
        [[nodiscard]] value_type minimumJacobian() const
        {
            value_type result = std::numeric_limits<value_type>::infinity();
            for ( int q = 0; q < this->nPoints(); ++q )
                result = std::min( result, this->J( q ) );
            return result;
        }

        /** @brief Return the minimum signed Jacobian over all points. */
        [[nodiscard]] value_type minimumSignedJacobian() const
        {
            value_type result = std::numeric_limits<value_type>::infinity();
            for ( int q = 0; q < this->nPoints(); ++q )
                result = std::min( result, this->signedJacobian( q ) );
            return result;
        }

        /** @brief Return the maximum signed Jacobian over all points. */
        [[nodiscard]] value_type maximumSignedJacobian() const
        {
            value_type result = -std::numeric_limits<value_type>::infinity();
            for ( int q = 0; q < this->nPoints(); ++q )
                result = std::max( result, this->signedJacobian( q ) );
            return result;
        }

        /** @brief Return the maximum Jacobian condition number over all points. */
        [[nodiscard]] value_type maximumConditionNumber() const
        {
            value_type result = value_type( 0 );
            for ( int q = 0; q < this->nPoints(); ++q )
                result = std::max( result, this->conditionNumber( q ) );
            return result;
        }

        /**
         * @brief Check whether the geometric map is finite and non-degenerate.
         * @param q geometric evaluation-point index
         * @param relativeTolerance scale-aware singular-value and Jacobian tolerance
         * @return true if the Jacobian is finite, full rank, and orientation preserving
         */
        [[nodiscard]] bool isValid(
            int q,
            value_type relativeTolerance =
                value_type( 64 ) * std::numeric_limits<value_type>::epsilon() ) const
        {
            return this->mappingQuality( q, relativeTolerance ).valid;
        }

        /**
         * @brief Check every evaluation point with a scale-aware tolerance.
         * @param relativeTolerance relative singular-value and Jacobian tolerance
         * @return true if every point is finite, full rank, and orientation preserving
         */
        [[nodiscard]] bool isValid(
            value_type relativeTolerance =
                value_type( 64 ) * std::numeric_limits<value_type>::epsilon() ) const
        {
            for ( int q = 0; q < this->nPoints(); ++q )
                if ( !this->isValid( q, relativeTolerance ) )
                    return false;
            return this->nPoints() > 0;
        }

        /**
         * @brief Compute a facet normal using the first geometric point.
         *
         * This compatibility overload is exact for affine geometric maps.
         * Point-supported finite-element functionals on curved cells must use
         * the overload taking an explicit evaluation-point index.
         *
         * @param faceId local facet index in the element
         * @param n output normal
         * @param scaled return the Piola-scaled normal when true
         */
        void
        faceNormal( int faceId, ublas::vector<value_type>& n, bool scaled = false ) const
        {
            this->faceNormal( faceId, 0, n, scaled );
        }

        /**
         * update the geomap through a neighbor matching face.
         *
         * Matching means that the points coordinates from both geomap are the same
         * up to a permutation. We compute here the permutation that ensures this
         * and update the geomap accordingly
         *
         * @return true if the permutation has been found and the geomap updated,
         * false otherwise.
         */
        template <int CTX, typename EltType, typename NeighborGeoType>
        bool updateFromNeighborMatchingFace( EltType const& elt, uint16_type face_in_elt,
                                             std::shared_ptr<NeighborGeoType> const& gmc )
            requires( SubEntityCoDim == 1 )
        {
            if ( M_gm->hasPermutationWithNeighborFace( elt.id(), face_in_elt ) == false )
            {
                auto [found_permutation, perm] = this->updateFromMatchingNodes<CTX>( elt, face_in_elt, gmc );
#if 0 // NOT WORK IN ALL CASES
                    if ( found_permutation )
                        M_gm->setPermutationWithNeighborFace( elt.id(), face_in_elt, perm.value() );
#endif
                return found_permutation;
            }
            this->update<CTX>( elt, face_in_elt, permutation_type( M_gm->permutationWithNeighborFace( elt.id(), face_in_elt ) ) );
            return true;
        }

        template <int CTX, typename EltType, typename MatchingGmcType>
        std::pair<bool, permutation_type> updateFromMatchingNodes(
            EltType const& elt, uint16_type subentity_in_elt,
            std::shared_ptr<MatchingGmcType> const& gmc )
            requires( SubEntityCoDim != 0 )
        {
            this->setElement( elt );
            auto gmcptr = gmc.get();
            for ( permutation_type __p( permutation_type::IDENTITY );
                  __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                // update only xReal in current geomap
                this->setSubEntity( subentity_in_elt, __p );
                this->updatePoints();

                bool check = true;
                for ( uint16_type i = 0; i < gmc->nPoints() && check; ++i )
                {
                    for ( uint16_type d = 0; ( d < NDim ); ++d )
                    {
                        check = check && ( std::abs( gmcptr->xReal( i )[d] - this->xReal( i )[d] ) < 1e-8 );
                    }
                }
                // if check update full gmc context with the good permutation
                if ( check )
                {
                    this->update<CTX>( elt, subentity_in_elt, __p );
                    return std::make_pair( true, __p );
                }
            }
            return std::make_pair( false, permutation_type( permutation_type::IDENTITY ) );
        }

        //@}
      private:
        /**
         * @brief Map a quadrature-point index to the compact geometry storage.
         *
         * The zero-or-one stride is selected when the context is constructed,
         * so runtime affine access is branch-free in pointwise consumers.
         */
        [[nodiscard]] FEELPP_STRONG_INLINE uint16_type
        computedPointIndex( int q ) const noexcept
        {
            if constexpr ( is_linear )
                return 0;
            else
                return static_cast<uint16_type>( q ) *
                       M_computed_point_stride;
        }

        /** @brief Decide once per update whether to use component batches. */
        [[nodiscard]] bool useBatchedGeometryEvaluation(
            uint16_type pointCount ) const noexcept
        {
            using enum geometry::geometry_evaluation_strategy;
            switch ( M_evaluation_strategy )
            {
            case pointwise:
                return false;
            case batched:
                return pointCount > 1;
            case automatic:
                return pointCount >= 4;
            }
            return false;
        }

        /**
         * \return the permutation associated with the face
         */
        permutation_type permutation( mpl::bool_<false> ) const
        {
            return M_perm;
        }

        /**
         * \return the permutation associated with the face
         */
        permutation_type permutation( mpl::bool_<true> ) const
        {
            CHECK( M_face_id == invalid_uint16_type_value ||
                   ( M_face_id != invalid_uint16_type_value &&
                     M_perm != permutation_type(
                         permutation_type::NO_PERMUTATION ) ) )
                << "invalid subentity permutation for face " << M_face_id;
            return M_perm;
        }

        template <size_type CTX>
        void updateImpl()
        {
            DCHECK( M_G.cols() == M_gm->nbPoints() ) << "Invalid number of points got " << M_G.cols() << " expected: " << M_gm->nbPoints();
            DCHECK( M_pc ) << "invalid precompute data structure";
            M_mapping_hessian_current = false;
            M_mapping_third_derivatives_current = false;

            if constexpr ( vm::has_point_v<CTX> )
            {
                this->updatePoints();
            }
            else if constexpr ( vm::has_dynamic_v<CTX> )
            {
                if ( hasPOINT( M_dynamic_context ) )
                    this->updatePoints();
            }
            if constexpr ( vm::has_measure_v<CTX> )
            {
                this->updateMeasures();
            }
            else if constexpr ( vm::has_dynamic_v<CTX> )
            {
                if ( hasMEASURE( M_dynamic_context ) )
                    this->updateMeasures();
            }

            if constexpr ( subEntityCoDim == 1 ) // this->isOnFace()
            {
                if constexpr ( vm::has_measure_v<CTX> || vm::has_tangent_v<CTX> )
                    this->updateFaceMeasures();
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasMEASURE( M_dynamic_context ) || hasTANGENT( M_dynamic_context ) )
                        this->updateFaceMeasures();
                }
            }

            updateJKBN<CTX>();
            M_context_updated = CTX | M_dynamic_context;
        }

        FEELPP_STRONG_INLINE
        void updateGradient( int q ) noexcept
        {
            M_gm->gradient( q, M_g, M_pc.get() );
        }

        /** Compute and cache geometry Hessians independently of context bits. */
        void ensureMappingHessian( int requestedPoint ) const
        {
            (void)requestedPoint;
            if ( M_mapping_hessian_current )
                return;

            const uint16_type computedPoints = nComputedPoints();
            M_hessian.resize( computedPoints );
            if ( this->usesAffineStorage() )
            {
                M_hessian[0].setZero();
                M_mapping_hessian_current = true;
                return;
            }

            const bool batched =
                this->useBatchedGeometryEvaluation( computedPoints );
            if ( batched )
                this->evaluateMappingHessianBatches( computedPoints );
            for ( uint16_type q = 0; q < computedPoints; ++q )
                this->updateMappingHessian( q, batched );
            M_mapping_hessian_current = true;
        }

        /** Third geometry derivatives are needed only by full Piola Hessians. */
        void ensureMappingThirdDerivatives() const
        {
            if ( M_mapping_third_derivatives_current )
                return;
            auto const pointCount = this->nComputedPoints();
            for ( int first = 0; first < PDim; ++first )
                for ( int second = 0; second < PDim; ++second )
                    for ( int third = 0; third < PDim; ++third )
                    {
                        auto& result = M_third_derivative_component_batches[( first * PDim + second ) * PDim + third];
                        result.resize( NDim, pointCount );
                        if ( this->usesAffineStorage() )
                            result.setZero();
                        else
                            result.noalias() = M_G * M_pc->thirdDerivativeEigen( first, second, third );
                    }
            M_mapping_third_derivatives_current = true;
        }

        FEELPP_STRONG_INLINE
        void updateMeasures() noexcept
        {
            M_h = M_element->h();
            M_h_min = M_element->hMin();
            M_meas = M_element->measure();
        }

        FEELPP_STRONG_INLINE
        void updateFaceMeasures() noexcept
        {
            M_measface = M_element->faceMeasure( M_face_id );
            M_h_face = M_element->hFace( M_face_id );
        }

        FEELPP_STRONG_INLINE
        void updatePoints() noexcept
        {
            auto Xreal = eigenMap( M_xrealq );
            // em_matrix_col_type<value_type> Pts( M_G.data().begin(), M_G.size1(), M_G.size2() );
            Xreal.noalias() = M_G * M_pc->phiEigen();
        }

        void updateJacobian( eigen_matrix_np_type const& K, value_type& J )
        {
            if constexpr ( dimension_v<ElementType> == real_dimension_v<ElementType> )
                J = math::abs( K.determinant() );
            else
                J = math::sqrt( math::abs( M_CSi.determinant() ) );
        }
        void updateB( eigen_matrix_np_type const& K, eigen_matrix_np_type& B )
        {
            if constexpr ( dimension_v<ElementType> == real_dimension_v<ElementType> )
            {
                M_CS.noalias() = K.inverse();
                B.noalias() = M_CS.transpose();
            }
            else
            {
                M_CS = M_CSi.inverse();
                B.noalias() = K * M_CS;
            }
        }
        //!
        //! update hessian information
        //!
        template <int CTX>
        void updateHessian( eigen_matrix_np_type const& B )
        {
            if constexpr ( vm::has_hessian_v<CTX> || vm::has_laplacian_v<CTX> )
            {
#if 0
                for ( uint16_type k = 0; k < NDim; ++k )
                    for ( uint16_type l = 0; l < NDim; ++l )
                        for ( uint16_type i = 0; i < PDim; ++i )
                            for ( uint16_type j = 0; j < PDim; ++j )
                                M_B3[k][l][i][j] = B( k, i ) * B( l, j );
#endif
            }
        }
        //!
        //! update normal in the real element
        //!

        FEELPP_STRONG_INLINE
        void updateNormals( eigen_matrix_np_type const& B,
                            eigen_vector_n_type& N, eigen_vector_n_type& unitN, value_type& Nnorm ) noexcept
        {
            N.noalias() = B * M_gm->referenceConvex().normal( M_face_id );
            Nnorm = N.norm();
            unitN = N / Nnorm;

            // M_gm->setNormalNormAtFace( M_id, M_face_id, M_n_norm );
            // M_gm->setUnitNormalAtFace( M_id, M_face_id, M_un_real );
        }

        FEELPP_STRONG_INLINE
        void updateTangents( eigen_matrix_np_type const& K, eigen_vector_n_type const& unitN,
                             eigen_vector_n_type& Ta, eigen_vector_n_type& unitT, value_type& Tnorm,
                             eigen_matrix_nn_type& P ) noexcept
        {
            auto const referenceTangents =
                this->referenceFacetTangentFrame( M_face_id );
            eigen_matrix_nt_type const mappedTangents = K * referenceTangents;
            auto const frame = geometry::orientedOrthonormalFrame(
                unitN, mappedTangents );
            auto const tangentFrame =
                frame.rightCols( facet_tangent_dimension );
            P = geometry::subspaceProjector( tangentFrame );

            if constexpr ( PDim == 2 )
            {
                Ta = mappedTangents.col( 0 );
                Tnorm = Ta.norm();
                unitT = tangentFrame.col( 0 );
            }
            else
            {
                Ta.setZero();
                unitT.setZero();
                Tnorm = value_type{};
            }
        }
        //!
        //! update normal in the real element
        //!
        FEELPP_STRONG_INLINE
        void updateLocalBasis( eigen_matrix_np_type const& K,
                               eigen_matrix_np_type const& B,
                               eigen_matrix_pp_type& local_basis_ref, eigen_matrix_np_type& local_basis_real ) noexcept
        {
#if 0
            if  ( !this->isOnFace() )
                return;//throw std::logic_error("Local basis defined only on faces ");
#endif
            local_basis_ref.col( 0 ) =
                M_gm->referenceConvex().normal( M_face_id );
            local_basis_ref.rightCols( facet_tangent_dimension ) =
                this->referenceFacetTangentFrame( M_face_id );

            eigen_vector_n_type const mappedNormal =
                B * local_basis_ref.col( 0 );
            eigen_matrix_nt_type const mappedTangents =
                K * local_basis_ref.rightCols( facet_tangent_dimension );
            local_basis_real = geometry::orientedOrthonormalFrame(
                mappedNormal, mappedTangents );
        }

        /**
         * @brief Store one geometry Hessian from the selected evaluation
         * layout.
         * @param q compact computed-point index
         * @param batched use component-major batch results when true
         */
        void evaluateMappingHessianBatches( uint16_type pointCount ) const
        {
            for ( auto& component : M_hessian_component_batches )
                component.resize( NDim, pointCount );
            for ( uint16_type first = 0; first < PDim; ++first )
            {
                for ( uint16_type second = 0; second < PDim; ++second )
                {
                    auto const component = first * PDim + second;
                    M_hessian_component_batches[component].noalias() =
                        M_G * M_pc->hessianEigen( first, second );
                }
            }
        }

        void updateMappingHessian( uint16_type q, bool batched ) const
        {
            if ( batched )
            {
                for ( uint16_type component = 0; component < NDim;
                      ++component )
                    for ( uint16_type first = 0; first < PDim; ++first )
                        for ( uint16_type second = 0; second < PDim;
                              ++second )
                            M_hessian[q]( component, first, second ) =
                                M_hessian_component_batches
                                    [first * PDim + second]( component, q );
                return;
            }

            tensor_map_t<2, value_type> geometryPoints(
                M_G.data(), M_G.rows(), M_G.cols() );
            Eigen::array<dimpair_t, 1> contraction = {
                { dimpair_t( 1, 0 ) }
            };
            M_gm->hessianBasisAtPoint(
                q, M_hessian_basis_at_pt, M_pc.get() );
            M_hessian[q] = geometryPoints.contract(
                M_hessian_basis_at_pt, contraction );
        }

        template <size_type CTX, bool ForceResizeJKBN = false>
        void resize()
        {
            if constexpr ( vm::has_point_v<CTX> )
            {
                M_xrealq.resize( NDim, nPoints() );
            }
            else if constexpr ( vm::has_dynamic_v<CTX> )
            {
                if ( hasPOINT( M_dynamic_context ) )
                    M_xrealq.resize( NDim, nPoints() );
            }
            if constexpr ( !is_linear_polynomial_v<gm_type> || ForceResizeJKBN )
                resizeJKBN<CTX>();
        }
        template <size_type CTX>
        void resizeJKBN()
        {
            constexpr bool IsOnFace = subEntityCoDim == 1;
            if constexpr ( !vm::has_dynamic_v<CTX> )
            {
                if constexpr ( !IsOnFace )
                {
                    if constexpr ( !(vm::has_jacobian_v<CTX> || vm::has_kb_v<CTX> || vm::has_mapping_hessian_v<CTX>))
                        return;
                }
                else
                {
                    if constexpr ( !(vm::has_jacobian_v<CTX> || vm::has_kb_v<CTX> || vm::has_mapping_hessian_v<CTX> || vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> || vm::has_local_basis_v<CTX>))
                        return;
                }
            }

            constexpr bool hasStaticB = IsOnFace ? vm::has_kb_v<CTX> || vm::has_mapping_hessian_v<CTX> || vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> || vm::has_local_basis_v<CTX> : vm::has_kb_v<CTX> || vm::has_mapping_hessian_v<CTX>;
            bool hasJ = vm::has_jacobian_v<CTX>;
            bool hasLocalBasis = IsOnFace && vm::has_local_basis_v<CTX>;
            bool hasTangent = IsOnFace && vm::has_tangent_v<CTX>;
            bool hasNormal = IsOnFace && ( vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> );
            bool hasB = hasStaticB; // vm::has_normal_v<CTX> || vm::has_tangent_v<CTX>;
            bool hasSecondDerivative = vm::has_mapping_hessian_v<CTX>;
            if constexpr ( vm::has_dynamic_v<CTX> )
            {
                if constexpr ( !vm::has_jacobian_v<CTX> )
                    hasJ = hasJACOBIAN( M_dynamic_context );
                if constexpr ( !vm::has_mapping_hessian_v<CTX> )
                    hasSecondDerivative = hasMAPPING_HESSIAN( M_dynamic_context );
                if constexpr ( IsOnFace )
                {
                    hasLocalBasis = hasLOCAL_BASIS( M_dynamic_context );
                    hasTangent = hasTANGENT( M_dynamic_context );
                    hasNormal = hasTangent || hasNORMAL( M_dynamic_context );
                    if constexpr ( !hasStaticB )
                        hasB = hasNormal || hasKB( M_dynamic_context );
                }
                else
                {
                    if constexpr ( !hasStaticB )
                        hasB = hasKB( M_dynamic_context );
                }

                if constexpr ( !( vm::has_jacobian_v<CTX> || hasStaticB ) )
                    if ( !( hasJ || hasB ) )
                        return;
            }

            const uint16_type _nComputedPoints = nComputedPoints();

            // K
            M_K.resize( _nComputedPoints );
            if ( _nComputedPoints > 1 )
                for ( auto& component : M_gradient_component_batches )
                    component.resize( NDim, _nComputedPoints );

            // jacobian
            if constexpr ( vm::has_jacobian_v<CTX> )
            {
                M_J.resize( _nComputedPoints );
            }
            else if constexpr ( vm::has_dynamic_v<CTX> )
            {
                if ( hasJ )
                    M_J.resize( _nComputedPoints );
            }

            // B
            if constexpr ( hasStaticB )
            {
                M_B.resize( _nComputedPoints );
            }
            else if constexpr ( vm::has_dynamic_v<CTX> )
            {
                if ( hasB )
                    M_B.resize( _nComputedPoints );
            }

            // Second derivatives are represented explicitly for every order.
            // Affine mappings keep one zero tensor so static and runtime G1
            // expose the same mathematical contract and compact storage.
            if constexpr ( vm::has_mapping_hessian_v<CTX> )
            {
                M_hessian.resize( _nComputedPoints );
                if ( _nComputedPoints > 1 )
                {
                    for ( auto& component : M_hessian_component_batches )
                        component.resize( NDim, _nComputedPoints );
                }
            }
            else if constexpr ( vm::has_dynamic_v<CTX> )
            {
                if ( hasSecondDerivative )
                {
                    M_hessian.resize( _nComputedPoints );
                    if ( _nComputedPoints > 1 )
                        for ( auto& component : M_hessian_component_batches )
                            component.resize( NDim, _nComputedPoints );
                }
            }

            if constexpr ( IsOnFace )
            {
                // normal
                if constexpr ( vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> )
                {
                    M_normals.resize( _nComputedPoints );
                    M_unit_normals.resize( _nComputedPoints );
                    M_normal_norms.resize( _nComputedPoints );
                }
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasNormal )
                    {
                        M_normals.resize( _nComputedPoints );
                        M_unit_normals.resize( _nComputedPoints );
                        M_normal_norms.resize( _nComputedPoints );
                    }
                }
                // tangent
                if constexpr ( vm::has_tangent_v<CTX> )
                {
                    M_tangents.resize( _nComputedPoints );
                    M_unit_tangents.resize( _nComputedPoints );
                    M_tangent_norms.resize( _nComputedPoints );
                    M_Ptangent.resize( _nComputedPoints );
                }
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasTangent )
                    {
                        M_tangents.resize( _nComputedPoints );
                        M_unit_tangents.resize( _nComputedPoints );
                        M_tangent_norms.resize( _nComputedPoints );
                        M_Ptangent.resize( _nComputedPoints );
                    }
                }
                // local basis
                if constexpr ( vm::has_local_basis_v<CTX> )
                {
                    M_local_basis_ref.resize( _nComputedPoints );
                    M_local_basis_real.resize( _nComputedPoints );
                }
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasLocalBasis )
                    {
                        M_local_basis_ref.resize( _nComputedPoints );
                        M_local_basis_real.resize( _nComputedPoints );
                    }
                }
            }
        }

        template <int CTX>
        void updateJKBN() noexcept
        {
            if constexpr ( is_linear )
            {
                this->template updateJKBNKernel<
                    CTX, geometry::AffineGeometryOrder>();
            }
            else if ( M_gm->isLinear() )
            {
                this->template updateJKBNKernel<
                    CTX, geometry::AffineGeometryOrder>();
            }
            else
            {
                this->template updateJKBNKernel<
                    CTX, geometry::HighOrderGeometry>();
            }
        }

        /**
         * @brief Evaluate requested geometry data with an order Strategy.
         * @tparam CTX compile-time context mask
         * @tparam OrderStrategy affine or high-order evaluation Strategy
         */
        template <int CTX,
                  geometry::GeometricOrderStrategy OrderStrategy>
        void updateJKBNKernel() noexcept
        {
            constexpr bool IsOnFace = subEntityCoDim == 1;
            if constexpr ( !vm::has_dynamic_v<CTX> )
            {
                if constexpr ( !IsOnFace )
                {
                    if constexpr ( !(vm::has_jacobian_v<CTX> || vm::has_kb_v<CTX> || vm::has_mapping_hessian_v<CTX>))
                        return;
                }
                else
                {
                    if constexpr ( !(vm::has_jacobian_v<CTX> || vm::has_kb_v<CTX> || vm::has_mapping_hessian_v<CTX> || vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> || vm::has_local_basis_v<CTX>))
                        return;
                }
            }

            constexpr bool hasStaticB = IsOnFace ? vm::has_kb_v<CTX> || vm::has_mapping_hessian_v<CTX> || vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> || vm::has_local_basis_v<CTX> : vm::has_kb_v<CTX> || vm::has_mapping_hessian_v<CTX>;
            bool hasJ = vm::has_jacobian_v<CTX>;
            bool hasLocalBasis = IsOnFace && vm::has_local_basis_v<CTX>;
            bool hasTangent = IsOnFace && vm::has_tangent_v<CTX>;
            bool hasNormal = IsOnFace && ( vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> );
            bool hasB = hasStaticB; // vm::has_normal_v<CTX> || vm::has_tangent_v<CTX>;
            bool hasSecondDerivative = vm::has_mapping_hessian_v<CTX>;
            if constexpr ( vm::has_dynamic_v<CTX> )
            {
                if constexpr ( !vm::has_jacobian_v<CTX> )
                    hasJ = hasJACOBIAN( M_dynamic_context );
                if constexpr ( !vm::has_mapping_hessian_v<CTX> )
                    hasSecondDerivative = hasMAPPING_HESSIAN( M_dynamic_context );
                if constexpr ( IsOnFace )
                {
                    hasLocalBasis = hasLOCAL_BASIS( M_dynamic_context );
                    hasTangent = hasTANGENT( M_dynamic_context );
                    hasNormal = hasTangent || hasNORMAL( M_dynamic_context );
                    if constexpr ( !hasStaticB )
                        hasB = hasNormal || hasKB( M_dynamic_context );
                }
                else
                {
                    if constexpr ( !hasStaticB )
                        hasB = hasKB( M_dynamic_context );
                }

                if constexpr ( !( vm::has_jacobian_v<CTX> || hasStaticB ) )
                    if ( !( hasJ || hasB ) )
                        return;
            }

            auto GradPhi = eigenMap( M_g );

            constexpr bool affine = OrderStrategy::is_affine;
            const uint16_type nPts = affine ? 1 : nComputedPoints();
            const bool useBatchedEvaluation =
                !affine && this->useBatchedGeometryEvaluation( nPts );

            if ( useBatchedEvaluation )
            {
                for ( uint16_type direction = 0; direction < PDim;
                      ++direction )
                {
                    M_gradient_component_batches[direction].noalias() =
                        M_G * M_pc->gradientEigen( direction );
                }

                if ( hasSecondDerivative )
                    this->evaluateMappingHessianBatches( nPts );
            }

            for ( uint16_type q = 0; q < nPts; ++q )
            {
                auto& K = M_K[q];
                if ( useBatchedEvaluation )
                {
                    for ( uint16_type direction = 0; direction < PDim;
                          ++direction )
                        K.col( direction ) =
                            M_gradient_component_batches[direction].col( q );
                }
                else
                {
                    updateGradient( q );
                    K.noalias() = M_G * GradPhi;
                }

                if constexpr ( dimension_v<ElementType> != real_dimension_v<ElementType> )
                {
                    // CS = K^T K
                    M_CSi.noalias() = K.transpose() * K;
                }

                // jacobian
                if constexpr ( vm::has_jacobian_v<CTX> )
                {
                    updateJacobian( K, M_J[q] );
                }
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasJ )
                        updateJacobian( K, M_J[q] );
                }

                // B
                if constexpr ( hasStaticB )
                {
                    updateB( K, M_B[q] );
                }
                else if constexpr ( vm::has_dynamic_v<CTX> )
                {
                    if ( hasB )
                        updateB( K, M_B[q] );
                }

                // second derivative
                if constexpr ( !affine )
                {
                    if constexpr ( vm::has_mapping_hessian_v<CTX> )
                    {
                        this->updateMappingHessian(
                            q, useBatchedEvaluation );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasSecondDerivative )
                            this->updateMappingHessian(
                                q, useBatchedEvaluation );
                    }
                }
                else if ( hasSecondDerivative )
                {
                    M_hessian[0].setZero();
                }

                if constexpr ( IsOnFace )
                {
                    // normal
                    if constexpr ( vm::has_normal_v<CTX> || vm::has_tangent_v<CTX> )
                    {
                        updateNormals( M_B[q], M_normals[q], M_unit_normals[q], M_normal_norms[q] );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasNormal )
                            updateNormals( M_B[q], M_normals[q], M_unit_normals[q], M_normal_norms[q] );
                    }
                    // tangent
                    if constexpr ( vm::has_tangent_v<CTX> )
                    {
                        updateTangents( K, M_unit_normals[q], M_tangents[q], M_unit_tangents[q], M_tangent_norms[q], M_Ptangent[q] );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasTangent )
                            updateTangents( K, M_unit_normals[q], M_tangents[q], M_unit_tangents[q], M_tangent_norms[q], M_Ptangent[q] );
                    }
                    // local basis
                    if constexpr ( vm::has_local_basis_v<CTX> )
                    {
                        updateLocalBasis( K, M_B[q], M_local_basis_ref[q], M_local_basis_real[q] );
                    }
                    else if constexpr ( vm::has_dynamic_v<CTX> )
                    {
                        if ( hasLocalBasis )
                            updateLocalBasis( K, M_B[q], M_local_basis_ref[q], M_local_basis_real[q] );
                    }
                }
            }
            if ( hasSecondDerivative )
                M_mapping_hessian_current = true;
        }

        friend class boost::serialization::access;

        template <class Archive>
        void serialize( Archive& ar, const unsigned int version )
        {
            if constexpr ( subEntityCoDim == 0 )
            {
                ar& BOOST_SERIALIZATION_NVP( M_dynamic_context );
                ar& BOOST_SERIALIZATION_NVP( M_context_updated );
                if constexpr ( Archive::is_saving::value )
                {
                    ar& boost::serialization::make_nvp( "xRefs", this->xRefs() );
                }
                else if constexpr ( Archive::is_loading::value )
                {
                    size_type saveDynamicContext = M_dynamic_context;
                    matrix_node_t_type _xRefs;
                    ar& boost::serialization::make_nvp( "xRefs", _xRefs );
                    M_pc = M_gm->preCompute( _xRefs );
                    M_dynamic_context = M_context_updated;
                    M_context_updated = 0;
                    M_npoints = M_pc->nPoints();
                    resize<vm::DYNAMIC, true>();
                    update<vm::DYNAMIC>( this->element() /*, __f, updatePc*/ );
                    M_dynamic_context = saveDynamicContext;
                }
            }
            else // TODO apply same procedure as above
            {

                ar& BOOST_SERIALIZATION_NVP( M_npoints );
                ar& BOOST_SERIALIZATION_NVP( M_J );
                ar& BOOST_SERIALIZATION_NVP( M_G );
                // ar& BOOST_SERIALIZATION_NVP( M_ref_normals );
                ar& BOOST_SERIALIZATION_NVP( M_normals );
                ar& BOOST_SERIALIZATION_NVP( M_unit_normals );
                ar& BOOST_SERIALIZATION_NVP( M_normal_norms );
                ar& BOOST_SERIALIZATION_NVP( M_tangents );
                ar& BOOST_SERIALIZATION_NVP( M_unit_tangents );
                ar& BOOST_SERIALIZATION_NVP( M_tangent_norms );
                ar& BOOST_SERIALIZATION_NVP( M_xrealq );
                // ar& BOOST_SERIALIZATION_NVP( M_g_linear );
                ar& BOOST_SERIALIZATION_NVP( M_g );
                ar& BOOST_SERIALIZATION_NVP( M_K );
                ar& BOOST_SERIALIZATION_NVP( M_CS );
                ar& BOOST_SERIALIZATION_NVP( M_CSi );
                ar& BOOST_SERIALIZATION_NVP( M_B );
                ar& BOOST_SERIALIZATION_NVP( M_Ptangent );
                ar& BOOST_SERIALIZATION_NVP( M_B3 ); //
                ar& BOOST_SERIALIZATION_NVP( M_id );
                ar& BOOST_SERIALIZATION_NVP( M_e_markers );
                ar& BOOST_SERIALIZATION_NVP( M_elem_id_1 );
                ar& BOOST_SERIALIZATION_NVP( M_pos_in_elem_id_1 );
                ar& BOOST_SERIALIZATION_NVP( M_elem_id_2 );
                ar& BOOST_SERIALIZATION_NVP( M_pos_in_elem_id_2 );
                ar& BOOST_SERIALIZATION_NVP( M_face_id );
                ar& BOOST_SERIALIZATION_NVP( M_h );
                ar& BOOST_SERIALIZATION_NVP( M_h_min );
                ar& BOOST_SERIALIZATION_NVP( M_h_face );
                ar& BOOST_SERIALIZATION_NVP( M_meas );
                ar& BOOST_SERIALIZATION_NVP( M_measface );
                ar& BOOST_SERIALIZATION_NVP( M_perm );
            }
        }

      private:
        gm_ptrtype M_gm;

        element_type const* M_element;
        // std::shared_ptr<element_type const> M_element_c;
        // element_type M_element_c;

        precompute_ptrtype M_pc;
        std::vector<std::map<permutation_type, precompute_ptrtype>> M_pc_faces;
        uint16_type M_npoints;
        uint16_type M_computed_point_stride;
        geometry::geometry_evaluation_strategy M_evaluation_strategy;

        std::vector<value_type> M_J;

        // matrix_type M_G;
        eigen_matrix_nx_type M_G;
        // vector_eigen_vector_p_type M_ref_normals;
        vector_eigen_vector_n_type M_normals, M_unit_normals;
        std::vector<value_type> M_normal_norms;
        vector_eigen_vector_n_type M_tangents, M_unit_tangents;
        std::vector<value_type> M_tangent_norms;
        vector_eigen_matrix_np_type M_local_basis_real;
        vector_eigen_matrix_pp_type M_local_basis_ref;

        matrix_type M_xrealq;

        // matrix_type M_g_linear;
        matrix_type M_g;
        mutable hessian_basis_type M_hessian_basis_at_pt;

        vector_eigen_matrix_np_type M_K;
        gradient_component_batches_type M_gradient_component_batches;
        eigen_matrix_pp_type M_CS;
        eigen_matrix_pp_type M_CSi;
        vector_eigen_matrix_np_type M_B;
        mutable vector_hessian_type M_hessian;
        mutable hessian_component_batches_type M_hessian_component_batches;
        mutable bool M_mapping_hessian_current = false;
        mutable std::array<derivative_component_batch_type, PDim * PDim * PDim> M_third_derivative_component_batches;
        mutable bool M_mapping_third_derivatives_current = false;
        vector_eigen_matrix_nn_type M_Ptangent;

        boost::multi_array<value_type, 4> M_B3;

        size_type M_id;
        std::map<uint16_type, mesh_marker_type> M_e_markers;
        std::optional<std::map<uint16_type, mesh_marker_type>> M_f_markers;
        size_type M_elem_id_1;
        uint16_type M_pos_in_elem_id_1;
        size_type M_elem_id_2;
        uint16_type M_pos_in_elem_id_2;

        uint16_type M_face_id;

        value_type M_h;
        value_type M_h_min;
        value_type M_h_face;
        value_type M_meas;
        value_type M_measface;

        permutation_type M_perm;

        size_type M_dynamic_context;
        size_type M_context_updated;

    }; // Context

    template <size_type Context_v, int SubEntityCoDim_v = 0, typename ElementType_t>
    std::shared_ptr<Context<ElementType_t, SubEntityCoDim_v>>
    context( geometric_mapping_ptrtype gm, ElementType_t const& e, precompute_ptrtype const& pc, size_type dynctx = 0 )
    {
        return std::make_shared<Context<ElementType_t, SubEntityCoDim_v>>( gm,
                                                                           e, std::integral_constant<int, Context_v>{},
                                                                           pc, invalid_uint16_type_value, dynctx );
    }

    template <size_type Context_v, int SubEntityCoDim_v = 0, typename ElementType_t>
    std::shared_ptr<Context<ElementType_t, SubEntityCoDim_v>>
    context( ElementType_t const& e, precompute_ptrtype const& pc, size_type dynctx = 0 )
        requires( !GeometricEntityViewType<ElementType_t> )
    {
        return std::make_shared<Context<ElementType_t, SubEntityCoDim_v>>(
            // super_enable_this::shared_from_this(),
            std::dynamic_pointer_cast<GeoMap<Dim, Order, RealDim, T, Entity, PP>>( this->shared_from_this() ),
            e, std::integral_constant<int, Context_v>{},
            pc, invalid_uint16_type_value, dynctx );
    }

    /**
     * @brief Create a geometric context from a compact entity view.
     * @tparam Context_v compile-time geometric context flags
     * @tparam SubEntityCoDim_v mapped subentity codimension
     * @tparam View geometric entity view type
     * @param view compact mathematical entity view
     * @param pc precomputed reference-basis values
     * @param dynctx runtime geometric context flags
     * @return initialized geometric context
     */
    template <size_type Context_v, int SubEntityCoDim_v = 0,
              GeometricEntityViewType View>
    auto context( View const& view,
                  precompute_ptrtype const& pc,
                  size_type dynctx = 0 )
    {
        using entity_type =
            typename std::remove_cvref_t<View>::entity_type;
        auto result =
            std::make_shared<Context<
                entity_type, SubEntityCoDim_v>>(
                std::dynamic_pointer_cast<
                    GeoMap<Dim, Order, RealDim, T, Entity, PP>>(
                    this->shared_from_this() ),
                view.entity(),
                std::integral_constant<int, Context_v>{},
                pc, invalid_uint16_type_value, dynctx );
        result->template update<Context_v>( view );
        return result;
    }

    template <size_type Context_v, int SubEntityCoDim_v = 1, typename ElementType_t>
    std::shared_ptr<Context<ElementType_t, SubEntityCoDim_v>>
    context( geometric_mapping_ptrtype gm,
             ElementType_t const& e,
             std::vector<std::map<typename ElementType_t::permutation_type, precompute_ptrtype>>& pc,
             uint16_type f,
             size_type dynctx = 0 )
    {
        return std::make_shared<Context<ElementType_t, SubEntityCoDim_v>>( gm,
                                                                           e, std::integral_constant<int, Context_v>{},
                                                                           pc,
                                                                           f, dynctx );
    }
    template <size_type Context_v, int SubEntityCoDim_v = 1, typename ElementType_t>
    std::shared_ptr<Context<ElementType_t, SubEntityCoDim_v>>
    context( ElementType_t const& e,
             std::vector<std::map<typename ElementType_t::template PermutationSubEntity<SubEntityCoDim_v>, precompute_ptrtype>>& pc,
             uint16_type f, size_type dynctx = 0 )
        requires( !GeometricEntityViewType<ElementType_t> )
    {
        return std::make_shared<Context<ElementType_t, SubEntityCoDim_v>>(
            // super_enable_this::shared_from_this(),
            std::dynamic_pointer_cast<GeoMap<Dim, Order, RealDim, T, Entity, PP>>( this->shared_from_this() ),
            e, std::integral_constant<int, Context_v>{},
            pc,
            f, dynctx );
    }

    /**
     * @brief Create a mapped-subentity context from a compact cell view.
     * @tparam Context_v compile-time geometric context flags
     * @tparam SubEntityCoDim_v mapped subentity codimension
     * @tparam View geometric cell-view type
     * @param view compact mathematical cell view
     * @param pc permutation-indexed subentity precomputes
     * @param f local subentity number
     * @param dynctx runtime geometric context flags
     * @return initialized mapped-subentity context
     */
    template <size_type Context_v, int SubEntityCoDim_v = 1,
              GeometricEntityViewType View>
    auto context(
        View const& view,
        std::vector<std::map<
            typename std::remove_cvref_t<View>::entity_type::
                template PermutationSubEntity<SubEntityCoDim_v>,
            precompute_ptrtype>>& pc,
        uint16_type f, size_type dynctx = 0 )
    {
        using entity_type =
            typename std::remove_cvref_t<View>::entity_type;
        auto result =
            std::make_shared<Context<
                entity_type, SubEntityCoDim_v>>(
                std::dynamic_pointer_cast<
                    GeoMap<Dim, Order, RealDim, T, Entity, PP>>(
                    this->shared_from_this() ),
                view.entity(),
                std::integral_constant<int, Context_v>{},
                pc, f, dynctx );
        result->template update<Context_v>( view, f );
        return result;
    }

    /**
     * @class Inverse
     * @brief Inverse and closest-point facade for one mapped geometric entity.
     *
     * Affine maps use their cached pseudoinverse. Nonlinear full-dimensional
     * maps use Newton steps, while embedded maps use Gauss--Newton closest-point
     * steps. All paths expose the same structured result and retain the legacy
     * `setXReal()`, `xRef()`, and `isIn()` accessors.
     */
    class Inverse
    {
      public:
        using geometric_mapping_type = GeoMap;
        using geometric_mapping_ptrtype =
            std::shared_ptr<geometric_mapping_type>;
        using node_t_type = typename geometric_mapping_type::node_t_type;
        using node_type = node_t_type;
        using matrix_node_t_type =
            typename geometric_mapping_type::matrix_node_t_type;
        using matrix_node_type = matrix_node_t_type;
        using eigen_reference_type =
            Eigen::Matrix<value_type, Dim, 1>;
        using eigen_physical_type =
            Eigen::Matrix<value_type, RealDim, 1>;
        using eigen_jacobian_type =
            Eigen::Matrix<value_type, RealDim, Dim>;
        using eigen_basis_gradient_type =
            Eigen::Matrix<value_type, Eigen::Dynamic, Dim, Eigen::ColMajor>;
        using eigen_physical_points_type =
            Eigen::Matrix<value_type, RealDim, Eigen::Dynamic, Eigen::ColMajor>;
        using inverse_options_type = geometry::InverseOptions<value_type>;
        using inverse_result_type =
            geometry::InverseResult<value_type, Dim, RealDim>;
        using inverse_results_type = std::vector<inverse_result_type>;

        static constexpr fem::transformation_type trans =
            geometric_mapping_type::trans;

        /**
         * @brief Construct an inverse facade from all geometric coefficients.
         * @tparam GeoElem geometric entity type
         * @param gm geometric mapping
         * @param ge geometric entity supplying `G()`
         * @param worldComm retained compatibility communicator; no solver is built
         */
        template <typename GeoElem>
        Inverse( geometric_mapping_ptrtype gm, GeoElem const& ge,
                 [[maybe_unused]] worldcomm_ptr_t const& worldComm =
                     Environment::worldCommSeqPtr() )
            requires( !GeometricEntityViewType<GeoElem> )
            : M_gm( std::move( gm ) ),
              M_xref( Dim ),
              M_xreal( RealDim ),
              M_G( ge.G() ),
              M_K( RealDim, Dim ),
              M_B( RealDim, Dim ),
              M_g( M_gm->nbPoints(), Dim )
        {
            M_gm->convertStoredGeometryToNodal(
                M_G, geometric_mapping_type::geometryRepresentationOf( ge ) );
            this->init();
        }

        /**
         * @brief Construct an inverse facade from a compact entity view.
         * @tparam View geometric entity view type
         * @param gm geometric mapping
         * @param view entity view supplying compact coordinates
         * @param worldComm retained compatibility communicator; no solver is built
         */
        template <GeometricEntityViewType View>
        Inverse( geometric_mapping_ptrtype gm, View const& view,
                 [[maybe_unused]] worldcomm_ptr_t const& worldComm =
                     Environment::worldCommSeqPtr() )
            : M_gm( std::move( gm ) ),
              M_xref( Dim ),
              M_xreal( RealDim ),
              M_G( RealDim, M_gm->nbPoints() ),
              M_K( RealDim, Dim ),
              M_B( RealDim, Dim ),
              M_g( M_gm->nbPoints(), Dim )
        {
            view.updateCoordinates( M_G );
            M_gm->convertStoredGeometryToNodal(
                M_G, geometric_mapping_type::geometryRepresentationOf( view ) );
            this->init();
        }

        /**
         * @brief Construct the affine compatibility facade from entity vertices.
         * @tparam GeoElem geometric entity type
         * @param gm order-one geometric mapping
         * @param ge geometric entity supplying `vertices()`
         * @param worldComm retained compatibility communicator; no solver is built
         */
        template <typename GeoElem>
        Inverse( geometric_mapping_ptrtype gm, GeoElem const& ge,
                 geometry::AffineGeometryOrder,
                 [[maybe_unused]] worldcomm_ptr_t const& worldComm =
                     Environment::worldCommSeqPtr() )
            requires( !GeometricEntityViewType<GeoElem> )
            : M_gm( std::move( gm ) ),
              M_xref( Dim ),
              M_xreal( RealDim ),
              M_G( ge.vertices() ),
              M_K( RealDim, Dim ),
              M_B( RealDim, Dim ),
              M_g( M_gm->nbPoints(), Dim )
        {
            this->init();
        }

        /**
         * @brief Compatibility constructor for the legacy MPL affine tag.
         * @tparam GeoElem geometric entity type
         */
        template <typename GeoElem>
        Inverse( geometric_mapping_ptrtype gm, GeoElem const& ge,
                 mpl::int_<1>,
                 worldcomm_ptr_t const& worldComm =
                     Environment::worldCommSeqPtr() )
            requires( !GeometricEntityViewType<GeoElem> )
            : Inverse( std::move( gm ), ge,
                       geometry::AffineGeometryOrder{}, worldComm )
        {}

        /**
         * @brief Construct an affine inverse facade from compact vertices.
         * @tparam View geometric entity view type
         * @param gm order-one geometric mapping
         * @param view entity view supplying compact vertex coordinates
         * @param worldComm retained compatibility communicator; no solver is built
         */
        template <GeometricEntityViewType View>
        Inverse( geometric_mapping_ptrtype gm, View const& view,
                 geometry::AffineGeometryOrder,
                 [[maybe_unused]] worldcomm_ptr_t const& worldComm =
                     Environment::worldCommSeqPtr() )
            : M_gm( std::move( gm ) ),
              M_xref( Dim ),
              M_xreal( RealDim ),
              M_G( RealDim, M_gm->nbPoints() ),
              M_K( RealDim, Dim ),
              M_B( RealDim, Dim ),
              M_g( M_gm->nbPoints(), Dim )
        {
            view.updateVertexCoordinates( M_G );
            this->init();
        }

        /**
         * @brief Compatibility view constructor for the legacy MPL affine tag.
         * @tparam View geometric entity view type
         */
        template <GeometricEntityViewType View>
        Inverse( geometric_mapping_ptrtype gm, View const& view,
                 mpl::int_<1>,
                 worldcomm_ptr_t const& worldComm =
                     Environment::worldCommSeqPtr() )
            : Inverse( std::move( gm ), view,
                       geometry::AffineGeometryOrder{}, worldComm )
        {}

        /** @brief Replace the geometric coefficients from an entity. */
        template <typename GeoElem>
        void update( GeoElem const& ge )
            requires( !GeometricEntityViewType<GeoElem> )
        {
            M_G = ge.G();
            M_gm->convertStoredGeometryToNodal(
                M_G, geometric_mapping_type::geometryRepresentationOf( ge ) );
            this->init();
        }

        /**
         * @brief Replace geometric coefficients from a compact entity view.
         * @tparam View geometric entity view type
         * @param view entity view supplying compact coordinates
         */
        template <GeometricEntityViewType View>
        void update( View const& view )
        {
            view.updateCoordinates( M_G );
            M_gm->convertStoredGeometryToNodal(
                M_G, geometric_mapping_type::geometryRepresentationOf( view ) );
            this->init();
        }

        /** @brief Replace the order-one geometric coefficients by vertices. */
        template <typename GeoElem>
        void update( GeoElem const& ge, geometry::AffineGeometryOrder )
            requires( !GeometricEntityViewType<GeoElem> )
        {
            M_G = ge.vertices();
            this->init();
        }

        /** @brief Compatibility update using the legacy MPL affine tag. */
        template <typename GeoElem>
        void update( GeoElem const& ge, mpl::int_<1> )
            requires( !GeometricEntityViewType<GeoElem> )
        {
            this->update( ge, geometry::AffineGeometryOrder{} );
        }

        /**
         * @brief Replace affine coefficients from compact entity vertices.
         * @tparam View geometric entity view type
         * @param view entity view supplying compact vertex coordinates
         */
        template <GeometricEntityViewType View>
        void update( View const& view, geometry::AffineGeometryOrder )
        {
            view.updateVertexCoordinates( M_G );
            this->init();
        }

        /** @brief Compatibility view update using the legacy MPL affine tag. */
        template <GeometricEntityViewType View>
        void update( View const& view, mpl::int_<1> )
        {
            this->update( view, geometry::AffineGeometryOrder{} );
        }

        /** @brief Return the associated geometric mapping. */
        [[nodiscard]] geometric_mapping_ptrtype const& geometricMapping() const
        {
            return M_gm;
        }

        /** @brief Return the physical-space dimension. */
        [[nodiscard]] static constexpr uint16_type N() noexcept
        {
            return RealDim;
        }

        /** @brief Return the reference-space dimension. */
        [[nodiscard]] static constexpr uint16_type P() noexcept
        {
            return Dim;
        }

        /** @brief Return the last recovered reference coordinate. */
        [[nodiscard]] node_t_type const& xRef() const { return M_xref; }
        /** @brief Return the last requested physical coordinate. */
        [[nodiscard]] node_t_type const& xReal() const { return M_xreal; }
        /** @brief Return the current geometric coefficient matrix. */
        [[nodiscard]] matrix_type const& G() const { return M_G; }
        /** @brief Return the Jacobian at the last recovered coordinate. */
        [[nodiscard]] matrix_type const& K() const { return M_K; }
        /** @brief Return the inverse-Jacobian transpose or metric pseudoinverse. */
        [[nodiscard]] matrix_type const& B() const { return M_B; }

        /**
         * @brief Return a read-only zero-copy Eigen view of geometric coefficients.
         * @return component-by-coefficient view tied to this inverse facade
         */
        [[nodiscard]] auto geometryCoefficientsView() const noexcept
        {
            return eigenMap( M_G );
        }

        /**
         * @brief Return a read-only zero-copy Eigen view of the Jacobian.
         * @return `RealDim`-by-`Dim` Jacobian view tied to this inverse facade
         */
        [[nodiscard]] auto jacobianView() const noexcept
        {
            return eigenMap( M_K );
        }

        /**
         * @brief Return a read-only zero-copy Eigen view of `K^{-T}` or its metric analogue.
         * @return `RealDim`-by-`Dim` inverse-Jacobian-transpose view
         */
        [[nodiscard]] auto inverseJacobianTransposeView() const noexcept
        {
            return eigenMap( M_B );
        }

        /**
         * @brief Select pointwise or batched affine inverse evaluation.
         * @param strategy evaluation Strategy; automatic selects batching
         */
        void setGeometryEvaluationStrategy(
            geometry::geometry_evaluation_strategy strategy ) noexcept
        {
            M_evaluation_strategy = strategy;
        }

        /** @brief Return the affine inverse evaluation Strategy. */
        [[nodiscard]] geometry::geometry_evaluation_strategy
        geometryEvaluationStrategy() const noexcept
        {
            return M_evaluation_strategy;
        }

        /** @brief Return the positive pointwise integration Jacobian. */
        [[nodiscard]] value_type J() const
        {
            eigen_jacobian_type const jacobian = this->jacobianView();
            if constexpr ( Dim == RealDim )
                return std::abs( jacobian.determinant() );
            else
                return std::sqrt(
                    ( jacobian.transpose() * jacobian ).determinant() );
        }

        /** @brief Return the reference-cell vertex barycenter. */
        [[nodiscard]] node_type barycenterRef() const
        {
            return this->toNode( this->referenceBarycenter() );
        }

        /** @brief Return the arithmetic barycenter of geometric coefficients. */
        [[nodiscard]] node_type barycenterReal() const
        {
            node_type result( RealDim );
            eigenMap( result ) =
                this->geometryCoefficientsView().rowwise().mean();
            return result;
        }

        /** @brief Return whether the last target belongs to the mapped cell. */
        [[nodiscard]] bool isIn() const noexcept { return M_is_in; }

        /** @brief Return the complete result of the last inverse operation. */
        [[nodiscard]] inverse_result_type const& lastResult() const noexcept
        {
            return M_last_result;
        }

        /**
         * @brief Return whether the last target lies on an embedded mapped cell.
         * @return true for an exact admissible embedded inverse
         */
        [[nodiscard]] bool isOnConvexSurface() const noexcept
        {
            if constexpr ( Dim == RealDim )
                return false;
            else
                return M_last_result.exact() &&
                       M_last_result.admissible();
        }

        /**
         * @brief Invert or project one physical point.
         * @param physicalPoint physical target
         * @param options numerical and reference-domain policy
         * @param initialGuess optional explicit reference initial guess
         * @return structured inverse result
         */
        [[nodiscard]] inverse_result_type inverse(
            node_type const& physicalPoint,
            inverse_options_type const& options = {},
            std::optional<eigen_reference_type> initialGuess = std::nullopt )
        {
            return this->inverse( this->physicalEigen( physicalPoint ),
                                  options, std::move( initialGuess ) );
        }

        /**
         * @brief Invert or project one Eigen physical point.
         * @param target physical target
         * @param options numerical and reference-domain policy
         * @param initialGuess optional explicit reference initial guess
         * @return structured inverse result
         */
        [[nodiscard]] inverse_result_type inverse(
            eigen_physical_type const& target,
            inverse_options_type const& options = {},
            std::optional<eigen_reference_type> initialGuess = std::nullopt )
        {
            eigen_reference_type guess = this->referenceBarycenter();
            bool usedWarmStart = false;
            if ( initialGuess )
                guess = *initialGuess;
            else if ( options.warm_start && M_last_result.converged() )
            {
                guess = M_last_result.reference_point;
                usedWarmStart = true;
            }

            inverse_result_type result;
            if ( M_gm->isLinear() )
                result = this->affineInverse( target, options );
            else
            {
                result = geometry::solveInverseMap<Dim, RealDim>(
                    target, guess,
                    [this]( eigen_reference_type const& point )
                    { return this->evaluate( point ); },
                    options,
                    [this]( eigen_reference_type const& point )
                    { return this->projectReference( point ); } );
            }
            result.used_warm_start = usedWarmStart;
            result.location = this->classifyReference(
                result.reference_point, options.reference_tolerance );
            this->storeResult( target, result, options );
            return result;
        }

        /**
         * @brief Return the closest point on the mapped closed cell.
         *
         * This is a constrained inverse for square mappings and a constrained
         * Gauss--Newton projection for embedded mappings. It reuses the same
         * solver and structured result as `inverse()`.
         *
         * @param physicalPoint physical target
         * @param options numerical policy; the domain policy is made constrained
         * @param initialGuess optional reference initial guess
         * @return reference coordinate, closest physical point, distance, and status
         */
        [[nodiscard]] inverse_result_type closestPoint(
            node_type const& physicalPoint,
            inverse_options_type options = {},
            std::optional<eigen_reference_type> initialGuess = std::nullopt )
        {
            options.domain_policy =
                geometry::InverseDomainPolicy::constrained;
            return this->inverse( physicalPoint, options,
                                  std::move( initialGuess ) );
        }

        /** @brief Eigen overload of the constrained closest-point operation. */
        [[nodiscard]] inverse_result_type closestPoint(
            eigen_physical_type const& physicalPoint,
            inverse_options_type options = {},
            std::optional<eigen_reference_type> initialGuess = std::nullopt )
        {
            options.domain_policy =
                geometry::InverseDomainPolicy::constrained;
            return this->inverse( physicalPoint, options,
                                  std::move( initialGuess ) );
        }

        /**
         * @brief Invert or project a batch of physical points.
         *
         * The first point starts from the reference barycenter. Subsequent
         * points reuse the preceding converged coordinate when warm starts are
         * enabled. No general nonlinear solver or per-point solver object is
         * constructed.
         *
         * @param physicalPoints physical targets stored by columns
         * @param options numerical and reference-domain policy
         * @return one structured result per target
         */
        [[nodiscard]] inverse_results_type inverse(
            matrix_node_t_type const& physicalPoints,
            inverse_options_type const& options = {} )
        {
            return this->inverse( eigenMap( physicalPoints ), options );
        }

        /**
         * @brief Invert or project an Eigen-compatible batch of physical points.
         * @tparam Derived Eigen dense expression storing targets by column
         * @param physicalPoints physical targets stored by columns
         * @param options numerical and reference-domain policy
         * @return one structured result per target
         */
        template <typename Derived>
        [[nodiscard]] inverse_results_type inverse(
            Eigen::MatrixBase<Derived> const& physicalPoints,
            inverse_options_type const& options = {} )
        {
            if ( physicalPoints.rows() != RealDim )
                throw std::invalid_argument(
                    "GeoMap::Inverse::inverse(): invalid physical point dimension" );
            if ( M_gm->isLinear() &&
                 M_evaluation_strategy !=
                     geometry::geometry_evaluation_strategy::pointwise &&
                 options.domain_policy ==
                     geometry::InverseDomainPolicy::unrestricted &&
                 M_jacobian_rank == Dim && physicalPoints.allFinite() )
                return this->affineInverseBatch( physicalPoints, options );
            inverse_results_type results;
            results.reserve( physicalPoints.cols() );
            std::optional<eigen_reference_type> guess =
                this->referenceBarycenter();
            for ( Eigen::Index point = 0; point < physicalPoints.cols();
                  ++point )
            {
                eigen_physical_type const target = physicalPoints.col( point );
                auto result = this->inverse( target, options, guess );
                result.used_warm_start = options.warm_start && point > 0;
                if ( options.warm_start && result.converged() )
                    guess = result.reference_point;
                else
                    guess = this->referenceBarycenter();
                results.push_back( std::move( result ) );
            }
            return results;
        }

        /**
         * @brief Project a uBLAS batch onto the mapped closed cell.
         * @param physicalPoints physical targets stored by columns
         * @param options numerical policy; the domain policy is made constrained
         * @return one structured closest-point result per target
         */
        [[nodiscard]] inverse_results_type closestPoint(
            matrix_node_t_type const& physicalPoints,
            inverse_options_type options = {} )
        {
            options.domain_policy =
                geometry::InverseDomainPolicy::constrained;
            return this->inverse( physicalPoints, options );
        }

        /**
         * @brief Project an Eigen-compatible batch onto the mapped closed cell.
         * @tparam Derived Eigen dense expression storing targets by column
         * @param physicalPoints physical targets stored by columns
         * @param options numerical policy; the domain policy is made constrained
         * @return one structured closest-point result per target
         */
        template <typename Derived>
        [[nodiscard]] inverse_results_type closestPoint(
            Eigen::MatrixBase<Derived> const& physicalPoints,
            inverse_options_type options = {} )
        {
            options.domain_policy =
                geometry::InverseDomainPolicy::constrained;
            return this->inverse( physicalPoints, options );
        }

        /**
         * @brief Legacy mutator applying the default unrestricted policy.
         * @param physicalPoint physical target
         */
        void setXReal( node_type const& physicalPoint )
        {
            static_cast<void>( this->inverse( physicalPoint ) );
        }

        /**
         * @brief Compatibility batch returning reference coordinates only.
         * @param physicalPoints physical targets stored by columns
         * @param allowExtrapolation retained compatibility flag
         * @return recovered reference coordinates stored by columns
         */
        matrix_node_t_type operator()(
            matrix_node_t_type const& physicalPoints,
            [[maybe_unused]] bool allowExtrapolation = false )
        {
            auto const results = this->inverse( physicalPoints );
            matrix_node_t_type referencePoints( Dim, physicalPoints.size2() );
            for ( std::size_t point = 0; point < results.size(); ++point )
                for ( int component = 0; component < Dim; ++component )
                    referencePoints( component, point ) =
                        results[point].reference_point( component );
            return referencePoints;
        }

      private:
        /** @brief Validate storage and initialize the pointwise linear algebra. */
        void init()
        {
            CHECK_EQ( M_G.size2(), M_gm->nbPoints() )
                << "invalid inverse geometry coefficient count";
            CHECK_EQ( M_G.size1(), RealDim )
                << "invalid inverse physical dimension";
            M_last_result = inverse_result_type{};
            M_is_in = false;
            this->updateAt( this->referenceBarycenter() );
        }

        /** @brief Convert a uBLAS physical node to its fixed Eigen value. */
        [[nodiscard]] eigen_physical_type
        physicalEigen( node_type const& point ) const
        {
            if ( point.size() != RealDim )
                throw std::invalid_argument(
                    "GeoMap::Inverse::inverse(): invalid physical point dimension" );
            return eigenMap( point );
        }

        /** @brief Convert a fixed Eigen reference value to a uBLAS node. */
        [[nodiscard]] node_type toNode( eigen_reference_type const& point ) const
        {
            node_type result( Dim );
            eigenMap( result ) = point;
            return result;
        }

        /** @brief Return the reference-cell vertex barycenter in Eigen form. */
        [[nodiscard]] eigen_reference_type referenceBarycenter() const
        {
            if constexpr ( convex_type::is_hypercube )
                return eigen_reference_type::Zero();
            else
                return eigen_reference_type::Constant(
                    value_type( 1 - Dim ) / value_type( Dim + 1 ) );
        }

        /** @brief Evaluate the physical map and Jacobian at a reference point. */
        [[nodiscard]] geometry::InverseMapEvaluation<value_type, Dim, RealDim>
        evaluate( eigen_reference_type const& point )
        {
            auto const referencePoint = this->toNode( point );
            auto const physicalPoint = M_gm->transform( referencePoint, M_G );
            M_gm->gradient( referencePoint, M_g );
            eigenMap( M_K ).noalias() =
                this->geometryCoefficientsView() * M_g;

            geometry::InverseMapEvaluation<value_type, Dim, RealDim> result;
            result.physical_point = eigenMap( physicalPoint );
            result.jacobian = this->jacobianView();
            return result;
        }

        /** @brief Update cached Jacobian and pseudoinverse at a point. */
        void updateAt( eigen_reference_type const& point )
        {
            auto const evaluation = this->evaluate( point );
            Eigen::ColPivHouseholderQR<eigen_jacobian_type> factorization(
                evaluation.jacobian );
            factorization.setThreshold( value_type( 1e-12 ) );
            M_jacobian_rank = factorization.rank();
            if ( M_jacobian_rank == Dim )
            {
                Eigen::Matrix<value_type, Dim, RealDim> const pseudoinverse =
                    factorization.solve(
                        Eigen::Matrix<value_type, RealDim, RealDim>::Identity() );
                eigenMap( M_B ) = pseudoinverse.transpose();
            }
            else
                M_B.clear();
            M_xref = this->toNode( point );
        }

        /** @brief Project a point onto the supported closed reference cell. */
        [[nodiscard]] eigen_reference_type
        projectReference( eigen_reference_type const& point ) const
        {
            if constexpr ( Dim == 0 )
                return point;
            else if constexpr ( convex_type::is_hypercube )
            {
                return point.cwiseMax( value_type( -1 ) )
                    .cwiseMin( value_type( 1 ) );
            }
            else
            {
                eigen_reference_type simplexCoordinate =
                    ( point.array() + value_type( 1 ) ) / value_type( 2 );
                simplexCoordinate =
                    simplexCoordinate.cwiseMax( value_type{} );
                if ( simplexCoordinate.sum() > value_type( 1 ) )
                {
                    std::array<value_type, Dim> sorted{};
                    for ( int component = 0; component < Dim; ++component )
                        sorted[component] = simplexCoordinate( component );
                    std::sort( sorted.begin(), sorted.end(), std::greater<>() );
                    value_type cumulative = value_type{};
                    value_type threshold = value_type{};
                    for ( int component = 0; component < Dim; ++component )
                    {
                        cumulative += sorted[component];
                        value_type const candidate =
                            ( cumulative - value_type( 1 ) ) /
                            value_type( component + 1 );
                        if ( sorted[component] > candidate )
                            threshold = candidate;
                    }
                    simplexCoordinate =
                        ( simplexCoordinate.array() - threshold )
                            .max( value_type{} );
                }
                return ( value_type( 2 ) * simplexCoordinate.array() -
                         value_type( 1 ) )
                    .matrix();
            }
        }

        /** @brief Return the minimum signed reference-cell margin. */
        [[nodiscard]] value_type
        referenceMargin( eigen_reference_type const& point ) const
        {
            if constexpr ( Dim == 0 )
                return value_type{};
            else if constexpr ( convex_type::is_hypercube )
                return value_type( 1 ) - point.cwiseAbs().maxCoeff();
            else
            {
                eigen_reference_type const simplexCoordinate =
                    ( point.array() + value_type( 1 ) ) / value_type( 2 );
                return std::min(
                    simplexCoordinate.minCoeff(),
                    value_type( 1 ) - simplexCoordinate.sum() );
            }
        }

        /** @brief Classify a recovered coordinate against the reference cell. */
        [[nodiscard]] geometry::InverseLocation classifyReference(
            eigen_reference_type const& point, value_type tolerance ) const
        {
            if ( !point.allFinite() )
                return geometry::InverseLocation::indeterminate;
            value_type const margin = this->referenceMargin( point );
            if ( margin < -tolerance )
                return geometry::InverseLocation::outside;
            if ( margin <= tolerance )
                return geometry::InverseLocation::boundary;
            return geometry::InverseLocation::inside;
        }

        /** @brief Apply the cached affine pseudoinverse and classify its residual. */
        template <typename Derived>
        [[nodiscard]] inverse_results_type affineInverseBatch(
            Eigen::MatrixBase<Derived> const& targets,
            inverse_options_type const& options )
        {
            using reference_points_type =
                Eigen::Matrix<value_type, Dim, Eigen::Dynamic>;
            using physical_points_type =
                Eigen::Matrix<value_type, RealDim, Eigen::Dynamic>;

            eigen_physical_type origin;
            for ( int component = 0; component < RealDim; ++component )
                origin( component ) = M_G( component, 0 );
            eigen_reference_type referenceOrigin;
            for ( int component = 0; component < Dim; ++component )
                referenceOrigin( component ) =
                    M_gm->refNode( 0 )( component );

            reference_points_type referencePoints =
                referenceOrigin.replicate( 1, targets.cols() );
            physical_points_type const centeredTargets =
                targets.colwise() - origin;
            referencePoints.noalias() +=
                this->inverseJacobianTransposeView().transpose() *
                centeredTargets;
            physical_points_type physicalPoints =
                origin.replicate( 1, targets.cols() );
            physicalPoints.noalias() +=
                this->jacobianView() *
                ( referencePoints.colwise() - referenceOrigin );

            inverse_results_type results;
            results.reserve( targets.cols() );
            for ( Eigen::Index point = 0; point < targets.cols(); ++point )
            {
                inverse_result_type result;
                result.reference_point = referencePoints.col( point );
                result.physical_point = physicalPoints.col( point );
                result.residual =
                    result.physical_point - targets.col( point );
                result.residual_norm = result.residual.norm();
                result.gradient_norm =
                    ( this->jacobianView().transpose() * result.residual )
                        .norm();
                result.jacobian_rank = M_jacobian_rank;
                result.used_warm_start = options.warm_start && point > 0;
                result.location = this->classifyReference(
                    result.reference_point,
                    options.reference_tolerance );
                value_type const tolerance =
                    options.absolute_residual_tolerance +
                    options.relative_residual_tolerance *
                        std::max( value_type( 1 ),
                                  targets.col( point ).norm() );
                result.status =
                    result.residual_norm <= tolerance
                        ? geometry::InverseStatus::converged_exact
                        : geometry::InverseStatus::converged_projected;
                results.push_back( std::move( result ) );
            }
            if ( !results.empty() )
            {
                eigen_physical_type const lastTarget =
                    targets.col( targets.cols() - 1 );
                this->storeResult(
                    lastTarget, results.back(), options );
            }
            return results;
        }

        /** @brief Apply the cached affine pseudoinverse and classify its residual. */
        [[nodiscard]] inverse_result_type affineInverse(
            eigen_physical_type const& target,
            inverse_options_type const& options )
        {
            inverse_result_type result;
            result.constrained =
                options.domain_policy ==
                geometry::InverseDomainPolicy::constrained;
            result.jacobian_rank = M_jacobian_rank;
            if ( !target.allFinite() )
            {
                result.status = geometry::InverseStatus::non_finite;
                return result;
            }
            if ( M_jacobian_rank < Dim )
            {
                result.status = geometry::InverseStatus::rank_deficient;
                return result;
            }

            eigen_physical_type origin;
            for ( int component = 0; component < RealDim; ++component )
                origin( component ) = M_G( component, 0 );
            eigen_reference_type referenceOrigin;
            for ( int component = 0; component < Dim; ++component )
                referenceOrigin( component ) = M_gm->refNode( 0 )( component );
            result.reference_point = referenceOrigin +
                                     this->inverseJacobianTransposeView().transpose() * ( target - origin );
            if ( result.constrained )
                result.reference_point =
                    this->projectReference( result.reference_point );

            auto const evaluation = this->evaluate( result.reference_point );
            result.physical_point = evaluation.physical_point;
            result.residual = result.physical_point - target;
            result.residual_norm = result.residual.norm();
            result.gradient_norm =
                ( evaluation.jacobian.transpose() * result.residual ).norm();
            value_type const tolerance =
                options.absolute_residual_tolerance +
                options.relative_residual_tolerance *
                    std::max( value_type( 1 ), target.norm() );
            result.status = result.residual_norm <= tolerance
                                ? geometry::InverseStatus::converged_exact
                                : geometry::InverseStatus::converged_projected;
            return result;
        }

        /** @brief Synchronize compatibility state from a structured result. */
        void storeResult( eigen_physical_type const& target,
                          inverse_result_type const& result,
                          inverse_options_type const& /**/ )
        {
            eigenMap( M_xreal ) = target;
            if ( result.reference_point.allFinite() )
            {
                M_xref = this->toNode( result.reference_point );
                if ( result.converged() )
                    this->updateAt( result.reference_point );
            }
            M_last_result = result;
            M_is_in = result.converged() && result.admissible() &&
                      ( Dim == RealDim || result.exact() );
        }

      private:
        geometric_mapping_ptrtype M_gm;
        node_type M_xref;
        node_type M_xreal;
        bool M_is_in = false;
        matrix_type M_G;
        matrix_type M_K;
        matrix_type M_B;
        eigen_basis_gradient_type M_g;
        int M_jacobian_rank = 0;
        geometry::geometry_evaluation_strategy M_evaluation_strategy =
            geometry::geometry_evaluation_strategy::automatic;
        inverse_result_type M_last_result;
    }; // Inverse

  private:
    template <typename Matrix>
    [[nodiscard]] static size_type matrixRows( Matrix const& matrix )
    {
        if constexpr ( requires { matrix.rows(); } )
            return static_cast<size_type>( matrix.rows() );
        else
            return static_cast<size_type>( matrix.size1() );
    }

    template <typename Matrix>
    [[nodiscard]] static size_type matrixColumns( Matrix const& matrix )
    {
        if constexpr ( requires { matrix.cols(); } )
            return static_cast<size_type>( matrix.cols() );
        else
            return static_cast<size_type>( matrix.size2() );
    }

    template <typename Matrix>
    [[nodiscard]] static Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>
    copyToEigen( Matrix const& matrix )
    {
        Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> result(
            matrixRows( matrix ), matrixColumns( matrix ) );
        for ( Eigen::Index row = 0; row < result.rows(); ++row )
            for ( Eigen::Index column = 0;
                  column < result.cols(); ++column )
                result( row, column ) = matrix( row, column );
        return result;
    }

    template <typename Derived, typename Matrix>
    static void copyFromEigen( Eigen::MatrixBase<Derived> const& source,
                               Matrix& destination )
    {
        if ( matrixRows( destination ) !=
                 static_cast<size_type>( source.rows() ) ||
             matrixColumns( destination ) !=
                 static_cast<size_type>( source.cols() ) )
            throw std::invalid_argument(
                "GeoMap geometry conversion matrix dimensions disagree" );
        for ( Eigen::Index row = 0; row < source.rows(); ++row )
            for ( Eigen::Index column = 0;
                  column < source.cols(); ++column )
                destination( row, column ) = source( row, column );
    }

    static void validateRepresentation(
        geometry::GeometryRepresentationKind representation )
    {
        if ( !geometry::isValidGeometryRepresentation( representation ) )
            throw std::invalid_argument(
                "Unsupported geometry coefficient representation" );
    }

    /** Build immutable nodal/Bezier conversion metadata once per GeoMap. */
    void initGeometryRepresentation()
    {
        typename representation_policy_type::points_type referenceNodes(
            nDim, this->points().size2() );
        for ( Eigen::Index coefficient = 0;
              coefficient < referenceNodes.cols(); ++coefficient )
            for ( int direction = 0; direction < nDim; ++direction )
                referenceNodes( direction, coefficient ) =
                    this->points()( direction, coefficient );

        if constexpr ( is_order_dynamic )
            M_geometryRepresentation.emplace(
                RuntimeOrder{ this->order() }, referenceNodes );
        else
            M_geometryRepresentation.emplace( referenceNodes );

        bezier_representation_policy_type bezier = [&]()
        {
            if constexpr ( is_order_dynamic )
                return bezier_representation_policy_type{
                    RuntimeOrder{ this->order() } };
            else
                return bezier_representation_policy_type{};
        }();
        auto const controlLocations =
            bezier.coefficientReferencePoints();
        M_bernsteinToCanonical.resize( controlLocations.cols() );
        std::vector<bool> matched(
            static_cast<size_type>( referenceNodes.cols() ), false );
        value_type const tolerance =
            value_type( 64 ) * std::numeric_limits<value_type>::epsilon();
        for ( Eigen::Index control = 0;
              control < controlLocations.cols(); ++control )
        {
            bool found = false;
            for ( Eigen::Index canonical = 0;
                  canonical < referenceNodes.cols(); ++canonical )
            {
                if ( matched[canonical] )
                    continue;
                if ( ( controlLocations.col( control ) -
                       referenceNodes.col( canonical ) )
                         .norm() <= tolerance )
                {
                    M_bernsteinToCanonical[control] = canonical;
                    matched[canonical] = true;
                    found = true;
                    break;
                }
            }
            if ( !found )
                throw std::logic_error(
                    "GeoMap could not match a Bezier coefficient to its canonical topology slot" );
        }
    }

    element_gm_ptrtype _elementMap;
    face_gm_ptrtype _boundaryMap;
    matrix_type M_g_linear;

    friend class Inverse;

    reference_convex_type M_refconvex;
    std::optional<representation_policy_type> M_geometryRepresentation;
    std::vector<Eigen::Index> M_bernsteinToCanonical;
};

#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/facilities/identity.hpp>
#include <boost/preprocessor/list/at.hpp>
#include <boost/preprocessor/list/cat.hpp>
#include <boost/preprocessor/list/for_each_product.hpp>
#include <boost/preprocessor/logical/and.hpp>
#include <boost/preprocessor/logical/or.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/tuple/eat.hpp>
#include <boost/preprocessor/tuple/to_list.hpp>

template <int Dim, int Order, int RealDim, template <int, int, int> class Entity = Simplex, typename T = double>
struct GT_Lagrange
{
};

/**
 * @brief Genuine runtime-order Lagrange geometric transformation.
 *
 * Static orders retain the generated `GT_Lagrange` specializations below.
 * The dynamic specialization directly reuses the same `GeoMap` mathematical
 * kernel and therefore does not dispatch to, or instantiate, a finite list of
 * static polynomial orders.
 *
 * @tparam Dim topological dimension
 * @tparam RealDim ambient dimension
 * @tparam Entity reference-cell family
 * @tparam T coordinate scalar type
 */
template <int Dim, int RealDim,
          template <int, int, int> class Entity, typename T>
struct GT_Lagrange<Dim, Dynamic, RealDim, Entity, T>
    : public GeoMap<Dim, Dynamic, RealDim, T, Entity, Lagrange>
{
    using super =
        GeoMap<Dim, Dynamic, RealDim, T, Entity, Lagrange>;
    using value_type = typename super::value_type;
    using element_gm_type = typename super::element_gm_type;
    using face_gm_type = typename super::face_gm_type;

    static constexpr uint16_type nDim = Dim;
    static constexpr int nOrder_v = Dynamic;
    static constexpr uint16_type nRealDim = RealDim;

    /** @brief Construct a true runtime-order mapping. */
    explicit GT_Lagrange( RuntimeOrder order )
        : super( order )
    {}
};

template <int Dim, int Order, int RealDim, typename T = double>
struct GT_QK
{
};

#/* List of dims. */
#define FEELPP_GEOMAP       \
    BOOST_PP_TUPLE_TO_LIST( \
        1,                  \
        (                   \
            Lagrange ) )    \
    /**/
#/* List of dims. */
#define FEELPP_DIMS         \
    BOOST_PP_TUPLE_TO_LIST( \
        4,                  \
        (                   \
            0, 1, 2, 3 ) )  \
    /**/
#/* List of real dims. */
#define FEELPP_REALDIMS     \
    BOOST_PP_TUPLE_TO_LIST( \
        3,                  \
        (                   \
            1, 2, 3 ) )

#/* List of real dims. */

#define FEELPP_NEWDIMS                    \
    BOOST_PP_TUPLE_TO_LIST(               \
        9,                                \
        (                                 \
            ( 0, 1 ), ( 0, 2 ), ( 0, 3 ), \
            ( 1, 1 ), ( 1, 2 ), ( 1, 3 ), \
            ( 2, 2 ), ( 2, 3 ),           \
            ( 3, 3 ) ) )                  \
    /**/
#/* List of orders. */
#define FEELPP_ORDERS         \
    BOOST_PP_TUPLE_TO_LIST(   \
        5,                    \
        (                     \
            1, 2, 3, 4, 5 ) ) \
    /**/
#
#define FEELPP_ENTITY BOOST_PP_TUPLE_TO_LIST( 2, ( Simplex, Hypercube ) )
/**/
#
#define FEELPP_GEN_GT( GEOM, LDIM, LORDER ) \
    "GT_" #GEOM "(" #LDIM "," #LORDER ")" /**/
#
#/* Generates code for all dim and order. */
#define FEELPP_GT_FACTORY_OP( _, GDO ) \
    FEELPP_GT_FACTORY GDO /**/
#
#
#define FEELPP_GT_DIM( T )         \
    BOOST_PP_TUPLE_ELEM( 2, 0, T ) \
/**/
#define FEELPP_GT_REALDIM( T )     \
    BOOST_PP_TUPLE_ELEM( 2, 1, T ) \
/**/
#
#
#define FEELPP_GT_FACTORY( GEOM, LDIMS, LORDER, ENTITY )                                                               \
    template <typename T>                                                                                              \
    struct BOOST_PP_CAT( GT_, GEOM )<FEELPP_GT_DIM( LDIMS ), LORDER, FEELPP_GT_REALDIM( LDIMS ), ENTITY, T>            \
        : public GeoMap<FEELPP_GT_DIM( LDIMS ), LORDER, FEELPP_GT_REALDIM( LDIMS ), T, ENTITY, GEOM>                   \
    {                                                                                                                  \
        static inline const uint16_type nDim = FEELPP_GT_DIM( LDIMS );                                                 \
        static inline const uint16_type order = LORDER;                                                                \
        static inline const uint16_type nRealDim = FEELPP_GT_REALDIM( LDIMS );                                         \
        static inline const uint16_type nRealDimCheck2d = mpl::if_<mpl::less_equal<mpl::int_<2>, mpl::int_<nRealDim>>, \
                                                                   mpl::int_<nRealDim>,                                \
                                                                   mpl::int_<nDim>>::type::value;                      \
                                                                                                                       \
        typedef mpl::vector<boost::none_t,                                                                             \
                            boost::none_t,                                                                             \
                            GeoMap<1, LORDER, nRealDim, T, ENTITY, GEOM>,                                              \
                            GeoMap<2, LORDER, nRealDimCheck2d, T, ENTITY, GEOM>>                                       \
            geomap_faces_t;                                                                                            \
        typedef mpl::vector<GeoMap<1, LORDER, nRealDim, T, ENTITY, GEOM>,                                              \
                            GeoMap<2, LORDER, nRealDimCheck2d, T, ENTITY, GEOM>,                                       \
                            GeoMap<3, LORDER, 3, T, ENTITY, GEOM>,                                                     \
                            boost::none_t>                                                                             \
            geomap_elements_t;                                                                                         \
        typedef typename type_traits<T>::value_type value_type;                                                        \
        typedef GeoMap<nDim, LORDER, nRealDim, T, ENTITY, GEOM> super;                                                 \
        typedef BOOST_PP_CAT( GT_, GEOM )<nDim - 1, LORDER, nRealDim, ENTITY, T> face_geo_type;                        \
                                                                                                                       \
        static inline const uint16_type nDof = super::nDof;                                                            \
        static inline const uint16_type nNodes = super::nNodes;                                                        \
        typedef typename mpl::at<geomap_elements_t, mpl::int_<nDim>>::type element_gm_type;                            \
        typedef typename mpl::at<geomap_faces_t, mpl::int_<nDim>>::type face_gm_type;                                  \
        template <int N>                                                                                               \
        struct face_gm                                                                                                 \
        {                                                                                                              \
            typedef typename mpl::at<geomap_faces_t, mpl::int_<N>>::type type;                                         \
        };                                                                                                             \
                                                                                                                       \
        BOOST_PP_CAT( GT_, GEOM )                                                                                      \
        ()                                                                                                             \
            : super( std::make_shared<element_gm_type>(), std::make_shared<face_gm_type>() )                           \
        {                                                                                                              \
        }                                                                                                              \
    };                                                                                                                 \
    /**/

template <typename Convex, typename T>
typename mpl::if_<Feel::is_simplex<Convex>,
                  mpl::identity<std::shared_ptr<GT_Lagrange<Convex::nDim, 1, Convex::nDim, Simplex, T>>>,
                  mpl::identity<std::shared_ptr<GT_Lagrange<Convex::nDim, 1, Convex::nDim, Hypercube, T>>>>::type::type
makeGeometricTransformation()
{
    using gmt = typename mpl::if_<is_simplex<Convex>,
                                  mpl::identity<GT_Lagrange<Convex::nDim, 1, Convex::nDim, Simplex, T>>,
                                  mpl::identity<GT_Lagrange<Convex::nDim, 1, Convex::nDim, Hypercube, T>>>::type::type;
    return std::make_shared<gmt>();
}

#if 0
#define FEELPP_GT_FACTORY( GEOM, LDIM, LORDER, LREALDIM, ENTITY ) \
    FEELPP_GT_FACTORY_BIS(GEOM,LDIM,LORDER,LREALDIM,ENTITY)),             \
                 BOOST_PP_EMPTY )                                       \
    /**/
#endif

// BOOST_PP_LIST_FOR_EACH_PRODUCT(FEELPP_GT_FACTORY_OP, 5, (FEELPP_GEOMAP, FEELPP_DIMS, FEELPP_ORDERS, FEELPP_REALDIMS, FEELPP_ENTITY))
BOOST_PP_LIST_FOR_EACH_PRODUCT( FEELPP_GT_FACTORY_OP, 4, ( FEELPP_GEOMAP, FEELPP_NEWDIMS, FEELPP_ORDERS, FEELPP_ENTITY ) )

#undef FEELPP_DIMS
#undef FEELPP_ORDERS
#undef FEELPP_REALDIMS
#undef FEELPP_NEWDIMS

template <typename Elem, template <int, int, int> class Entity = Simplex, typename T = double>
class RealToReference
{
  public:
    static inline const uint16_type nDim = Elem::nDim;
    static inline const uint16_type nRealDim = Elem::nRealDim;

    typedef T value_type;
    typedef typename matrix_node<T>::type points_type;
    typedef GT_Lagrange<nDim, 1, nRealDim, Entity, value_type> gm_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::Inverse inverse_gm_type;

    RealToReference( Elem const& elem )
        : M_gm( new gm_type ),
          M_igm( M_gm, elem )
    {
    }

    points_type operator()( points_type const& pts ) const
    {
        return M_igm( pts );
    }

    value_type J() const
    {
        return M_igm.J();
    }

  private:
    gm_ptrtype M_gm;
    inverse_gm_type M_igm;
};

} // namespace Feel
#endif
