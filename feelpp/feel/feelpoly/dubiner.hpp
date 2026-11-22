/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2007 Joseph Fourier University
 * SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file dubiner.hpp
 * @brief Compatibility facade for regular-coordinate simplex modal tabulation.
 */
#ifndef FEELPP_FEELPOLY_DUBINER_HPP
#define FEELPP_FEELPOLY_DUBINER_HPP 1
#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/pointset.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/polynomialcoefficients.hpp>
#include <feel/feelpoly/referenceoperators.hpp>
#include <feel/feelpoly/simplextabulation.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/warpblend.hpp>
#include <feel/feelpoly/expansiontypes.hpp>
#include <array>
#include <utility>

namespace Feel
{
namespace detail
{
/** Mathematical operator identity independent of a static facade's seed degree. */
template <int Dimension, int RealDimension, typename T, bool Normalized, typename Matrix>
struct DubinerOperatorFamily;
} // namespace detail
/** @brief Established simplex-basis facade, parameterized by its legacy storage policy. */
template<uint16_type Dim, uint16_type RealDim, uint16_type Degree,
         typename NormalizationPolicy, typename T, template<class> class StoragePolicy>
class Dubiner;

/** @brief Geometry and compatibility-storage types for a Dubiner expansion. */
template<uint16_type Dim, uint16_type RealDim, uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>, typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
struct DubinerTraits
{
    /** Fixed topological dimension. */
    static inline constexpr uint16_type nDim = Dim;
    /** Fixed embedding dimension. */
    static inline constexpr uint16_type nRealDim = RealDim;
    /** Compile-time seed expansion degree. */
    static inline constexpr uint16_type nOrder = Degree;
    /** Interior interpolation lattice order used only by coefficient operators. */
    static inline constexpr uint16_type nConvexOrderDiff = nDim + nOrder + 1;
    /** Whether the modal convention is L2-orthonormal. */
    static inline constexpr bool is_normalized = NormalizationPolicy::is_normalized;
    /** Modal scalar type. */
    using value_type = T;
    /** @brief Simplex and reference geometry at a selected geometry order. */
    template<uint16_type Order, typename V = value_type>
    struct Convex
    {
        /** Topological reference simplex. */
        using type = Simplex<nDim, Order, nDim>;
        /** Reference geometry, independent of sampled modal storage. */
        using reference_type = Reference<type, nDim, Order, nDim, V>;
    };
    /** @brief Rebind the modal scalar while retaining geometry and storage policy. */
    template<typename NewT>
    struct ChangeValueType
    {
        /** Rebound basis facade. */
        using type = Dubiner<Dim, RealDim, Degree, NormalizationPolicy, NewT, StoragePolicy>;
        /** Rebound facade traits. */
        using traits_type = DubinerTraits<Dim, RealDim, Degree, NormalizationPolicy, NewT, StoragePolicy>;
    };
    /** @brief Rebind the compile-time expansion seed degree. */
    template<uint16_type NewOrder>
    struct ChangeOrder
    {
        /** Rebound basis facade. */
        using type = Dubiner<Dim, RealDim, NewOrder, NormalizationPolicy, T, StoragePolicy>;
        /** Rebound facade traits. */
        using traits_type = DubinerTraits<Dim, RealDim, NewOrder, NormalizationPolicy, T, StoragePolicy>;
    };
    /** Existing reference simplex type. */
    using convex_type = typename Convex<nOrder>::type;
    /** Existing reference geometry type. */
    using reference_convex_type = typename Convex<nOrder>::reference_type;
    /** Auxiliary geometry used only by lazy coefficient differentiation. */
    using diff_convex_type = typename Convex<nConvexOrderDiff>::type;
    /** Reference geometry associated with that auxiliary lattice. */
    using diff_reference_convex_type = typename Convex<nConvexOrderDiff>::reference_type;
    /** Existing interior-node policy, constructed only by an operator request. */
    using diff_pointset_type = if_t<nDim == 2,
        PointSetWarpBlend<diff_convex_type, nConvexOrderDiff, value_type>,
        PointSetEquiSpaced<diff_convex_type, nConvexOrderDiff, value_type>>;
    /** Compatibility representation policy, separate from native numerical storage. */
    using storage_policy = StoragePolicy<value_type>;
    /** Owning compatibility result matrix. */
    using matrix_type = typename storage_policy::matrix_type;
    /** Compatibility direction collection. */
    using vector_matrix_type = typename storage_policy::vector_matrix_type;
    /** Compatibility point-batch matrix. */
    using matrix_node_type = typename storage_policy::matrix_node_type;
    /** Established point-batch argument representation. */
    using points_type = typename storage_policy::points_type;
    /** Established scalar-point argument representation. */
    using node_type = typename storage_policy::node_type;
};

/** @brief Established dispatch tag for simplex dimension and expansion degree. */
template<int D, int O>
struct DubinerTag
{
    /** Reference dimension. */
    static constexpr int Dim = D;
    /** Expansion degree. */
    static constexpr int Order = O;
};

/** @brief Degree-ordered Dubiner expansion with native regular-coordinate tabulation.
 * @tparam Dim Reference dimension; native kernels cover one through three.
 * @tparam RealDim Embedding dimension retained for existing FE descriptors.
 * @tparam Degree Compile-time seed degree, independently of a requested runtime order.
 * @tparam NormalizationPolicy Existing normalized or non-normalized convention.
 * @tparam T Scalar type.
 * @tparam StoragePolicy Representation used only at compatibility result boundaries.
 *
 * Values, gradients and Hessians use the same division-free native recurrence.
 * No auxiliary point set or factorization is constructed until d() is requested.
 * Compile-time coefficient operators are immutable function-local data; their
 * references survive copies and destruction of individual facade objects.
 * @see https://doi.org/10.1145/1644001.1644006
 */
template<uint16_type Dim, uint16_type RealDim, uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>, typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
class Dubiner
{
public:
    /** Geometry and compatibility-storage traits. */
    using traits_type = DubinerTraits<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy>;
    /** Fixed reference dimension. */
    static inline constexpr uint16_type nDim = traits_type::nDim;
    /** Fixed embedding dimension. */
    static inline constexpr uint16_type nRealDim = traits_type::nRealDim;
    /** Compile-time expansion seed degree. */
    static inline constexpr uint16_type nOrder = traits_type::nOrder;
    /** Auxiliary interior lattice order for coefficient operators. */
    static inline constexpr uint16_type nConvexOrderDiff = traits_type::nConvexOrderDiff;
    /** Whether each mode is L2-normalized. */
    static inline constexpr bool is_normalized = traits_type::is_normalized;
    /** Established affine transformation contract. */
    static inline constexpr bool isTransformationEquivalent = true;
    /** Modal coefficients have no inter-element continuity. */
    static inline constexpr bool isContinuous = false;
    /** Established descriptor classification retained for compatibility. */
    static inline constexpr bool is_product = true;
    /** Modal continuity marker. */
    using continuity_type = Discontinuous;
    /** This facade type. */
    using self_type = Dubiner;
    /** The facade is its own expansion basis. */
    using basis_type = self_type;

    /** Native tabulation kernel with degree independent of this facade's seed. */
    template <int OrderSpec>
    using TabulationKernelFor = SimplexDubinerBasis<nDim, OrderSpec, T, is_normalized>;
    /** Scalar type. */
    using value_type = typename traits_type::value_type;
    /** Existing convex descriptor. */
    using convex_type = typename traits_type::convex_type;
    /** Existing reference geometry. */
    using reference_convex_type = typename traits_type::reference_convex_type;
    /** Existing interior point policy, used lazily by d(). */
    using diff_pointset_type = typename traits_type::diff_pointset_type;
    /** Compatibility storage policy. */
    using storage_policy = typename traits_type::storage_policy;
    /** Owning compatibility matrix. */
    using matrix_type = typename traits_type::matrix_type;
    /** Compatibility direction collection. */
    using vector_matrix_type = typename traits_type::vector_matrix_type;
    /** Compatibility point matrix. */
    using matrix_node_type = typename traits_type::matrix_node_type;
    /** Existing point-batch argument type. */
    using points_type = typename traits_type::points_type;
    /** Existing single-point argument type. */
    using node_type = typename traits_type::node_type;

    /** @brief Construct a stateless descriptor without geometry nodes or operator setup. */
    constexpr Dubiner() = default;
    /** @brief Copy the stateless descriptor; cached operators have independent ownership. */
    Dubiner( Dubiner const& ) = default;
    /** @brief Move the stateless descriptor. */
    Dubiner( Dubiner&& ) = default;
    /** @brief Assign an equivalent descriptor without changing operator references. */
    Dubiner& operator=( Dubiner const& ) = default;
    /** @brief Move an equivalent descriptor. */
    Dubiner& operator=( Dubiner&& ) = default;
    /** @brief Destroy the descriptor; cached operators retain independent lifetimes. */
    ~Dubiner() = default;

    /** @brief Evaluate the complete expansion at a single reference point. */
    matrix_type operator()( node_type const& point ) const
    {
        if ( point.size() != nDim ) throw std::invalid_argument( "Dubiner point dimension mismatch" );
        points_type points( nDim, 1 );
        if constexpr ( EigenMatrix<points_type> ) points.col( 0 ) = point;
        else eigenMap( points ).col( 0 ) = eigenMap( point );
        return evaluate( points );
    }
    /** @brief Evaluate the complete expansion at a point batch. */
    matrix_type operator()( points_type const& points ) const { return evaluate( points ); }
    /** @brief Number of scalar modes at the compile-time seed degree. */
    [[nodiscard]] size_type size() const { return ::Feel::detail::SimplexCardinality::total( nDim, nOrder ); }
    /** @brief Compile-time seed expansion degree. */
    [[nodiscard]] uint16_type degree() const noexcept { return nOrder; }
    /** @brief Return this facade as its own expansion basis. */
    [[nodiscard]] self_type const& basis() const noexcept { return *this; }
    /** @brief Whether the expansion uses L2-normalized modes. */
    [[nodiscard]] bool isNormalized() const noexcept { return is_normalized; }
    /** @brief Established mathematical family name. */
    [[nodiscard]] std::string familyName() const { return "dubiner"; }
    /** @brief Identity coefficients in scalar mode space, independent of geometry nodes. */
    [[nodiscard]] matrix_type coeff() const
    {
        auto const count = size();
        if constexpr ( EigenMatrix<matrix_type> ) return matrix_type::Identity( count, count );
        else return ublas::identity_matrix<value_type>( count, count );
    }

    /** @brief Evaluate compatibility values at the compile-time degree. */
    [[nodiscard]] static matrix_type evaluate( points_type const& points )
    { return evaluateWithOrder<nOrder>( points, RuntimeOrder{ nOrder } ); }
    /** @brief Evaluate compatibility values at a requested runtime degree. */
    [[nodiscard]] static matrix_type evaluate( points_type const& points, uint16_type order )
    { return evaluateWithOrder<Dynamic>( points, RuntimeOrder{ order } ); }

    /** @brief Evaluate native values at the compile-time degree. */
    template<typename Points>
    [[nodiscard]] static auto values( Points const& points ) requires( nDim > 0 )
    { return values<nOrder>( points, RuntimeOrder{ nOrder } ); }
    /** @brief Evaluate native values with a runtime expansion axis. */
    template<typename Points>
    [[nodiscard]] static auto values( Points const& points, RuntimeOrder order ) requires( nDim > 0 )
    { return values<Dynamic>( points, order ); }
    /** @brief Preserve an explicit static/runtime degree independently of the facade seed. */
    template<int OrderSpec, typename Points>
    [[nodiscard]] static auto values( Points const& points, RuntimeOrder order ) requires( nDim > 0 )
    {
        return ::Feel::detail::withEigenMatrix<value_type>( points, [&]( auto const& coordinates )
        { return SimplexDubinerBasis<nDim, OrderSpec, value_type, is_normalized>{ order }.values( coordinates ); } );
    }
    /** @brief Evaluate native reference gradients at the compile-time degree. */
    template<typename Points>
    [[nodiscard]] static auto gradients( Points const& points ) requires( nDim > 0 )
    { return gradients<nOrder>( points, RuntimeOrder{ nOrder } ); }
    /** @brief Evaluate native reference gradients with runtime modal extents. */
    template<typename Points>
    [[nodiscard]] static auto gradients( Points const& points, RuntimeOrder order ) requires( nDim > 0 )
    { return gradients<Dynamic>( points, order ); }
    /** @brief Preserve explicit static/runtime order in the direct gradient recurrence. */
    template<int OrderSpec, typename Points>
    [[nodiscard]] static auto gradients( Points const& points, RuntimeOrder order ) requires( nDim > 0 )
    {
        return ::Feel::detail::withEigenMatrix<value_type>( points, [&]( auto const& coordinates )
        { return SimplexDubinerBasis<nDim, OrderSpec, value_type, is_normalized>{ order }.gradients( coordinates ); } );
    }
    /** @brief Evaluate native reference Hessians at the compile-time degree. */
    template<typename Points>
    [[nodiscard]] static auto hessians( Points const& points ) requires( nDim > 0 )
    { return hessians<nOrder>( points, RuntimeOrder{ nOrder } ); }
    /** @brief Evaluate native reference Hessians at a runtime degree. */
    template<typename Points>
    [[nodiscard]] static auto hessians( Points const& points, RuntimeOrder order ) requires( nDim > 0 )
    { return hessians<Dynamic>( points, order ); }
    /** @brief Evaluate actual Hessians directly, without coefficient differentiation matrices. */
    template<int OrderSpec, typename Points>
    [[nodiscard]] static auto hessians( Points const& points, RuntimeOrder order ) requires( nDim > 0 )
    {
        return ::Feel::detail::withEigenMatrix<value_type>( points, [&]( auto const& coordinates )
        { return SimplexDubinerBasis<nDim, OrderSpec, value_type, is_normalized>{ order }.hessians( coordinates ); } );
    }

    /** @brief Compatibility gradient collection evaluated by the regular recurrence. */
    template<typename Derived>
    [[nodiscard]] static vector_matrix_type derivate( ublas::matrix_expression<Derived> const& points )
    { return derivativesWithOrder<nOrder>( points(), RuntimeOrder{ nOrder } ); }
    /** @brief Compatibility gradients at the requested effective degree. */
    template<typename Derived>
    [[nodiscard]] static vector_matrix_type derivate( ublas::matrix_expression<Derived> const& points, uint16_type order )
    { return derivativesWithOrder<Dynamic>( points(), RuntimeOrder{ order } ); }

    /** @brief Borrow a thread-safe immutable coefficient operator, constructing it on first use. */
    [[nodiscard]] matrix_type const& d( uint16_type direction ) const
    {
        if ( direction >= nDim ) throw std::out_of_range( "invalid Dubiner derivative direction" );
        return differentiationOperators()[direction];
    }
    /** @brief Established spelling for borrowing a coefficient differentiation operator. */
    [[nodiscard]] matrix_type const& derivate( uint16_type direction ) const { return d( direction ); }
    /** @brief Own a coefficient operator at a requested degree.
     * Interpolation/factorization is confined to this explicitly requested operator
     * path. The owning compatibility result shares its construction with every direction.
     */
    [[nodiscard]] matrix_type d( uint16_type direction, uint16_type order ) const
    {
        if constexpr ( nDim == 0 ) throw std::out_of_range( "a point has no derivative directions" );
        else
        {
            if ( direction >= nDim ) throw std::out_of_range( "invalid Dubiner derivative direction" );
            if ( order == nOrder ) return d( direction );
            return referenceOperators( RuntimeOrder{ order } )->matrix( direction );
        }
    }

    /** Fixed-dimensional owned coefficient operators. */
    using OperatorData = ::Feel::detail::ReferenceDerivativeOperators<matrix_type, nDim>;

    /** @brief Retain every requested-order direction after one interpolation factorization. */
    [[nodiscard]] static std::shared_ptr<OperatorData const> referenceOperators( RuntimeOrder degree )
    {
        if constexpr ( nDim == 0 )
            return std::make_shared<OperatorData const>( OperatorData{ degree.value, {} } );
        else
        {
            degree = ::Feel::detail::validatedNodalOrder<::Feel::detail::SimplexCardinality>( nDim, degree );

            // Keep the established static interpolation nodes and their numerical policy.
            if ( degree.value == nOrder )
            {
                static auto const operators = []
                {
                    diff_pointset_type points( 1 );
                    return std::make_shared<OperatorData const>( interpolateAtPoints<nOrder>(
                        points.points(), RuntimeOrder{ nOrder } ) );
                }();
                return operators;
            }

            // Non-seed runtime orders share their existing equispaced interpolation policy.
            using Family = ::Feel::detail::DubinerOperatorFamily<nDim, nRealDim, value_type, is_normalized, matrix_type>;
            auto build = [degree]
            {
                using Cell = Simplex<nDim, 1, nDim>;
                PointSetEquiSpaced<Cell, Dynamic, value_type> points{ degree };
                return interpolateAtPoints<Dynamic>( points.points(), degree );
            };
            return ::Feel::detail::ReferenceOperatorCache<Family, OperatorData>::get( degree, build );
        }
    }

  private:
    /** @brief Use one numerical implementation for static and runtime point-selection policies. */
    template <int OrderSpec, typename Points>
    [[nodiscard]] static OperatorData interpolateAtPoints( Points const& points, RuntimeOrder degree )
    {
        auto values = self_type::values<OrderSpec>( points, degree );
        auto gradients = self_type::gradients<OrderSpec>( points, degree );
        return ::Feel::detail::interpolateReferenceOperators<nDim, matrix_type>(
            degree, values, gradients, Feel::type_traits<value_type>::epsilon(), "Dubiner differentiation" );
    }
    /** @brief Materialize only the requested compatibility value representation. */
    template<int OrderSpec>
    [[nodiscard]] static matrix_type evaluateWithOrder( points_type const& points, RuntimeOrder order )
    {
        if constexpr ( nDim == 0 )
        {
            return ::Feel::detail::withEigenMatrix<value_type>( points, []( auto const& coordinates )
            {
                return ::Feel::detail::materializePolynomialMatrix<matrix_type>(
                    Eigen::Matrix<value_type, 1, Eigen::Dynamic>::Ones( 1, coordinates.cols() ) );
            } );
        }
        else
        {
            auto result = values<OrderSpec>( points, order );
            return ::Feel::detail::materializePolynomialMatrix<matrix_type>( result.matrix() );
        }
    }
    /** @brief Adapt the direct gradient collection once at the compatibility boundary. */
    template<int OrderSpec, typename Points>
    [[nodiscard]] static vector_matrix_type derivativesWithOrder( Points const& points, RuntimeOrder order )
    {
        vector_matrix_type result( nDim == 0 ? 1 : nDim );
        if constexpr ( nDim == 0 )
            result[0] = ::Feel::detail::materializePolynomialMatrix<matrix_type>(
                Eigen::Matrix<value_type, 1, Eigen::Dynamic>::Zero( 1, points.size2() ) );
        else
        {
            auto native = gradients<OrderSpec>( points, order );
            for ( int d = 0; d < nDim; ++d )
                result[d] = ::Feel::detail::materializePolynomialMatrix<matrix_type>( native[d].matrix() );
        }
        return result;
    }
    /** @brief Construct every compile-time coefficient direction with one lazy factorization. */
    [[nodiscard]] static std::array<matrix_type, nDim> const& differentiationOperators()
    {
        static auto const operators = referenceOperators( RuntimeOrder{ nOrder } );
        return operators->matrices;
    }
};
} // namespace Feel
#endif // FEELPP_FEELPOLY_DUBINER_HPP
