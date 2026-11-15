/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2006 EPFL
 * SPDX-FileCopyrightText: 2007 Université Joseph Fourier Grenoble 1
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file legendre.hpp
 * @brief Tensor-product Legendre bases and requested-order evaluation.
 */
#ifndef __Legendre_H
#define __Legendre_H 1

#include <vector>
#include <array>

#include <boost/lambda/if.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/referenceoperators.hpp>
#include <feel/feelpoly/expansions.hpp>
#include <feel/feelpoly/legendretabulation.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/policy.hpp>
#include <feel/feelpoly/gausslobatto.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/expansiontypes.hpp>

namespace Feel
{
namespace detail
{
/** Mathematical identity shared across static facade seeds and requested runtime degrees. */
template <int Dimension, int RealDimension, typename T, bool Normalized, typename Matrix>
struct LegendreOperatorFamily;
} // namespace detail
template< class Convex,
          int Order,
          typename T >
class PointSetGaussLobatto;

template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy,
         typename T,
         template<class> class StoragePolicy>
class Legendre;


template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
struct LegendreTraits
{
    static const uint16_type nDim = Dim;
    static const uint16_type nRealDim = RealDim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrderDiff = nOrder+2;
    static inline const bool is_normalized = NormalizationPolicy::is_normalized;

    /** @name Typedefs
     */
    //@{

    /*
     * numerical type
     */
    typedef T value_type;

    template<uint16_type order, typename V = value_type>
    struct Convex
    {
        typedef Hypercube<nDim, order, /*nRealDim*/nDim> type;
        typedef Reference<Hypercube<nDim, order, /*nRealDim*/nDim>, nDim, order, nDim/*nRealDim*/, V>  reference_type;
    };

    template<typename NewT>
    struct ChangeValueType
    {
        typedef Legendre<Dim, RealDim, Degree, NormalizationPolicy, NewT, StoragePolicy> type;
        typedef LegendreTraits<Dim, RealDim, Degree, NormalizationPolicy, NewT, StoragePolicy> traits_type;
    };

    template<uint16_type NewOrder>
    struct ChangeOrder
    {
        typedef Legendre<Dim, RealDim, NewOrder, NormalizationPolicy, T, StoragePolicy> type;
        typedef LegendreTraits<Dim, RealDim, NewOrder, NormalizationPolicy, T, StoragePolicy> traits_type;
    };

    /*
     * Geometry where the polynomials are defined and constructed
     */
    typedef typename Convex<nOrder>::type convex_type;
    typedef typename Convex<nOrder>::reference_type reference_convex_type;

    typedef typename Convex<nConvexOrderDiff>::type diff_convex_type;
    typedef typename Convex<nConvexOrderDiff>::reference_type diff_reference_convex_type;

    typedef PointSetGaussLobatto<diff_convex_type, nConvexOrderDiff, value_type> diff_pointset_type;

    /*
     * storage policy
     */
    typedef StoragePolicy<value_type> storage_policy;
    typedef typename storage_policy::matrix_type matrix_type;
    typedef typename storage_policy::vector_matrix_type vector_matrix_type;
    typedef typename storage_policy::matrix_node_type matrix_node_type;
    typedef typename storage_policy::points_type points_type;
    typedef typename storage_policy::node_type node_type;
};

template<int D, int O>
struct LegendreTag
{
    static const int Dim = D;
    static const int Order = O;
};
/**
 * \class Legendre
 * \brief Legendre polynomial orthonormal basis
 *
 * This class represents the Legendre polynomials up to degree \c
 * Degree on a hypercube in dimension \c Dim.
 *
 *
 * The legendre polynomials in 1D, the segment \f$[-1;1]\f$ are defined
 * using Jacobi polynomials as follows:
 * \f$ \phi_i(x) = P_i^{0,0}(x) \f$
 *where \f$P_i^{0,0}(x)\f$ is the i-th Jacobi polynomial evaluated at
 * \f$x \in [-1;1]\f$ with weights \f$(0,0)\f$.
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 *
 * @see G.E. Karniadakis and S.J. Sherwin, ''Spectral/hp Element
 * Methods for CFD,'' Oxford University Press, March 1999.
 *
 */
template<uint16_type Dim,
         uint16_type RealDim,
         uint16_type Degree,
         typename NormalizationPolicy = Normalized<true>,
         typename T = double,
         template<class> class StoragePolicy = StorageUBlas>
class Legendre
{
public:

    typedef LegendreTraits<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy> traits_type;
    static const uint16_type nDim = Dim;
    static const uint16_type nRealDim = RealDim;
    static const uint16_type nOrder = Degree;
    static const uint16_type nConvexOrder = nOrder+2;
    static inline const bool is_normalized = NormalizationPolicy::is_normalized;
    static inline const bool isTransformationEquivalent = true;
    static inline const bool isContinuous = false;
    static inline const bool is_product = true;
    typedef Discontinuous continuity_type;

    /** @name Typedefs
     */
    //@{

    typedef Legendre<Dim, RealDim, Degree, NormalizationPolicy, T, StoragePolicy> self_type;

    /*
     * for now can be used only as a basis but we might be interested
     * to have then expressed in other basis like in the moment basis
     */
    typedef self_type basis_type;

    typedef typename traits_type::value_type value_type;

    /*
     * Geometry where the polynomials are defined and constructed
     */
    typedef typename traits_type::convex_type convex_type;
    typedef typename traits_type::reference_convex_type reference_convex_type;

    typedef typename traits_type::diff_pointset_type diff_pointset_type;

    /*
     * storage policy
     */
    typedef typename traits_type::storage_policy storage_policy;
    typedef typename traits_type::matrix_type matrix_type;
    typedef typename traits_type::vector_matrix_type vector_matrix_type;
    typedef typename traits_type::matrix_node_type matrix_node_type;
    typedef typename traits_type::points_type points_type;
    typedef typename traits_type::node_type node_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    Legendre()
        :
        M_refconvex(),
        M_pts( M_refconvex.makePoints( nDim, 0 ) )
    {
        // Values and sampled derivatives do not require coefficient operators.
    }
    Legendre( Legendre const & d )
        :
        M_refconvex(),
        M_pts( d.M_pts )
    {

    }

    ~Legendre()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type const& operator=( self_type const& d )
    {
        if ( this != &d )
        {
            M_pts = d.M_pts;
        }

        return *this;
    }

    //ublas::matrix_column<matrix_type const> operator()( node_type const& pt ) const
    matrix_type operator()( node_type const& pt ) const
    {
        points_type pts( pt.size(), 1 );
        ublas::column( pts, 0 ) = pt;
        return evaluate( pts );
    }

    matrix_type operator()( points_type const& pts ) const
    {
        return evaluate( pts );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * Number of polynomials in set
     */
    size_type size() const
    {
        return convex_type::polyDims( nOrder );
    }

    /**
     * \return the maximum degree of the Legendre polynomial to be
     * constructed
     */
    uint16_type degree() const
    {
        return nOrder;
    }

    /**
     * \return self as a basis
     */
    self_type const& basis() const
    {
        return *this;
    }

    /**
     * \return true if the Legendre polynomials are normalized, false
     * otherwise
     */
    bool isNormalized() const
    {
        return is_normalized;
    }

    /**
     * \return the \c familyName()
     */
    std::string familyName() const
    {
        return "legendre";
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Legendre polynomials is an orthonormal basis, the coefficients
     * of the polynomials of the basis are the canonical vectors and
     * represented by the identity matrix (lines are polynomials and
     * columns are the polynomial basis )
     *
     * This function is correct only if we use the Legendre polynomials
     * as a basis
     */
    matrix_type coeff() const
    {
        // Modal coefficients are the identity in mode space, independently of geometry vertices.
        auto const count = ::Feel::detail::HypercubeCardinality::total( nDim, nOrder );
        if constexpr ( EigenMatrix<matrix_type> ) return matrix_type::Identity( count, count );
        else return ublas::identity_matrix<value_type>( count, count );
    }


    /**
     * evaluate the Legendre polynomials at a set of points \p __pts
     *
     * \arg __x is a set of points
     */
    static matrix_type evaluate( points_type const& __pts )
    {
        return evaluateWithOrder<nOrder>( __pts, RuntimeOrder{ nOrder } );
    }

    /**
     * @brief Evaluate Legendre polynomials up to given order
     *
     * A different requested order uses the shared runtime kernel. Lower-order
     * tensor-product modes are not a contiguous prefix of the higher basis.
     *
     * @param __pts Points to evaluate at
     * @param order Requested polynomial order, below/equal/above the template order
     * @return Evaluation matrix
     */
    static matrix_type evaluate( points_type const& __pts, uint16_type order )
    {
        return evaluateWithOrder<Dynamic>( __pts, RuntimeOrder{ order } );
    }

    /** @brief Evaluate native modal values with the compile-time expansion degree. */
    template<typename Points>
    [[nodiscard]] static auto values( Points const& points )
    { return valuesWithOrder<nOrder>( points, RuntimeOrder{ nOrder } ); }

    /** @brief Evaluate native modal values while preserving runtime order semantics. */
    template<typename Points>
    [[nodiscard]] static auto values( Points const& points, RuntimeOrder order )
    { return valuesWithOrder<Dynamic>( points, order ); }

    /** @brief Evaluate an explicit static/runtime extent independently of the seed degree. */
    template<int OrderSpec, typename Points>
    [[nodiscard]] static auto values( Points const& points, RuntimeOrder order )
    { return valuesWithOrder<OrderSpec>( points, order ); }

    /** @brief Evaluate native reference gradients at the compile-time degree. */
    template<typename Points>
    [[nodiscard]] static auto gradients( Points const& points )
    { return gradientsWithOrder<nOrder>( points, RuntimeOrder{ nOrder } ); }

    /** @brief Evaluate native reference gradients at a runtime degree. */
    template<typename Points>
    [[nodiscard]] static auto gradients( Points const& points, RuntimeOrder order )
    { return gradientsWithOrder<Dynamic>( points, order ); }

    /** @brief Retain an explicit order specification in native reference gradients. */
    template<int OrderSpec, typename Points>
    [[nodiscard]] static auto gradients( Points const& points, RuntimeOrder order )
    { return gradientsWithOrder<OrderSpec>( points, order ); }

    /** Native kernel preserving effective degree independently of the facade seed. */
    template <int OrderSpec>
    using TabulationKernelFor = TensorLegendreBasis<nDim, OrderSpec, value_type, is_normalized>;

    /** @brief Evaluate direct native reference Hessians at the compile-time degree. */
    template <typename Points>
    [[nodiscard]] static auto hessians( Points const& points )
    {
        return hessians<nOrder>( points, RuntimeOrder{ nOrder } );
    }

    /** @brief Evaluate direct reference Hessians with a runtime degree axis. */
    template <typename Points>
    [[nodiscard]] static auto hessians( Points const& points, RuntimeOrder order )
    {
        return hessians<Dynamic>( points, order );
    }

    /** @brief Retain the requested degree and sample extents in direct reference Hessians. */
    template <int OrderSpec, typename Points>
    [[nodiscard]] static auto hessians( Points const& points, RuntimeOrder order )
    {
        return ::Feel::detail::withEigenMatrix<value_type>( points, [&]( auto const& coordinates )
                                                            { return TabulationKernelFor<OrderSpec>{ order }.hessians( coordinates ); } );
    }

    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts )
    {
        return derivativesWithOrder<nOrder>( __pts(), RuntimeOrder{ nOrder } );
    }

    /**
     * @brief Derivate Legendre polynomials up to given order
     *
     * Use the requested tensor-product ordering for every derivative direction.
     *
     * @param __pts Points to evaluate at
     * @param order Requested polynomial order, below/equal/above the template order
     * @return Vector of derivation matrices (one per spatial dimension)
     */
    template<typename AE>
    static vector_matrix_type derivate( ublas::matrix_expression<AE>  const& __pts, uint16_type order )
    {
        return derivativesWithOrder<Dynamic>( __pts(), RuntimeOrder{ order } );
    }

    /**
     * \brief derivatives of Legendre polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    static matrix_type const& d( uint16_type i )
    {
        return differentiationOperators().at( i );
    }

    /**
     * @brief Compute derivation matrix at runtime for given order
     *
     * Used for dynamic polynomial order support. Computes the derivation
     * matrix D_i such that d/dx_i f = D_i * f for polynomial coefficients f.
     *
     * @param i Derivative direction (0: x, 1: y, 2: z)
     * @param order Polynomial order
     * @return Derivation matrix of size (nBasis x nBasis) for the given order
     */
    static matrix_type d( uint16_type i, uint16_type order )
    {
        if ( i >= nDim ) throw std::out_of_range( "invalid Legendre derivative direction" );
        if ( order == nOrder ) return differentiationOperators().at( i );
        return referenceOperators( RuntimeOrder{ order } )->matrix( i );
    }

    /** Fixed-dimensional coefficient operators, shared by both order policies. */
    using OperatorData = ::Feel::detail::ReferenceDerivativeOperators<matrix_type, nDim>;

    /** @brief Retain all directions at one checked effective degree.
     * @note P_n' = sum_(k<n, n-k odd) (2k+1) P_k. Normalized coefficients
     * are sqrt((2n+1)(2k+1)); tensor directions preserve every other mode index.
     * This construction needs neither points nor interpolation/factorization.
     */
    [[nodiscard]] static std::shared_ptr<OperatorData const> referenceOperators( RuntimeOrder degree )
    {
        degree = ::Feel::detail::validatedNodalOrder<::Feel::detail::HypercubeCardinality>( nDim, degree );
        using Family = ::Feel::detail::LegendreOperatorFamily<nDim, nRealDim, value_type, is_normalized, matrix_type>;
        auto build = [degree]
        {
            OperatorData result{ degree.value, {} };
            Eigen::Index const radix = degree.value + 1;
            Eigen::Index const count = ::Feel::detail::HypercubeCardinality::total( nDim, degree.value );
            (void)::Feel::detail::matrixElementCount( count, count, sizeof( value_type ) );

            // Decode the established last-axis-fastest tensor ordering once per direction/row.
            Eigen::Index stride = 1;
            for ( int direction = nDim - 1; direction >= 0; --direction, stride *= radix )
            {
                result.matrices[direction].resize( count, count );
                decltype( auto ) matrix = ::Feel::detail::referenceOperatorView( result.matrices[direction] );
                matrix.setZero();
                for ( Eigen::Index row = 0; row < count; ++row )
                {
                    int const high = static_cast<int>( row / stride % radix );
                    for ( int low = high - 1; low >= 0; low -= 2 )
                    {
                        value_type coefficient = value_type( 2 * low + 1 );
                        if constexpr ( is_normalized )
                        {
                            using std::sqrt;
                            coefficient = sqrt( value_type( 2 * high + 1 ) * coefficient );
                        }
                        matrix( row, row - ( high - low ) * stride ) = coefficient;
                    }
                }
            }
            return result;
        };

        // Share complete immutable data without repeating construction across facade seeds.
        return ::Feel::detail::ReferenceOperatorCache<Family, OperatorData>::get( degree, build );
    }

    /**
     * \brief derivatives of Legendre polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    static matrix_type const& derivate( uint16_type i )
    {
        return differentiationOperators().at( i );
    }


    //@}

private:
    /** @brief Borrow contiguous legacy/Eigen coordinates or materialize an expression once. */
    template<typename Points, typename Operation>
    static auto withEigenPoints( Points const& points, Operation&& operation )
    {
        return ::Feel::detail::withEigenMatrix<value_type>( points, std::forward<Operation>( operation ) );
    }

    /** @brief Materialize one result at the selected compatibility-storage boundary. */
    template<typename Derived>
    static matrix_type materialize( Eigen::MatrixBase<Derived> const& expression )
    {
        matrix_type result( expression.rows(), expression.cols() );
        if constexpr ( EigenMatrix<matrix_type> ) result = expression;
        else eigenMap( result ) = expression;
        return result;
    }

    /** @brief Preserve the requested order specification without a compatibility copy. */
    template<int OrderSpec, typename Points>
    static auto valuesWithOrder( Points const& points, RuntimeOrder order )
    {
        return withEigenPoints( points, [&]( auto const& coordinates )
        {
            return TensorLegendreBasis<nDim, OrderSpec, value_type, is_normalized>{ order }.values( coordinates );
        } );
    }

    /** @brief Materialize native values only at the legacy evaluation boundary. */
    template<int OrderSpec>
    static matrix_type evaluateWithOrder( points_type const& points, RuntimeOrder order )
    {
        auto values = valuesWithOrder<OrderSpec>( points, order );
        return materialize( values.matrix() );
    }

    /** @brief Forward all reference directions to the shared native gradient kernel. */
    template<int OrderSpec, typename Points>
    static auto gradientsWithOrder( Points const& points, RuntimeOrder order )
    {
        return withEigenPoints( points, [&]( auto const& coordinates )
        {
            return TensorLegendreBasis<nDim, OrderSpec, value_type, is_normalized>{ order }.gradients( coordinates );
        } );
    }

    /** @brief Forward all reference directions to the same native gradient kernel. */
    template<int OrderSpec, typename Points>
    static vector_matrix_type derivativesWithOrder( Points const& points, RuntimeOrder order )
    {
        auto gradients = gradientsWithOrder<OrderSpec>( points, order );
        vector_matrix_type result( nDim );
        for ( int d = 0; d < nDim; ++d ) result[d] = materialize( gradients[d].matrix() );
        return result;
    }

    /** @brief Initialize fixed-dimensional operators once, with thread-safe static initialization. */
    static std::array<matrix_type, nDim> const& differentiationOperators()
    {
        static auto const operators = referenceOperators( RuntimeOrder{ nOrder } );
        return operators->matrices;
    }

    /** Reference geometry retained for legacy node access. */
    reference_convex_type M_refconvex;
    /** Owned legacy reference nodes. */
    points_type M_pts;

}; // class Dubiner



} // namespace Feel
#endif /* __Legendre_H */
