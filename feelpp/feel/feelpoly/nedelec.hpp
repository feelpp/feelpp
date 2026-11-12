/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2006-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file nedelec.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-11-25
 */
#ifndef FEELPP_NEDELEC_HPP
#define FEELPP_NEDELEC_HPP 1

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>

#include <boost/type_traits.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/lu.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/pointset.hpp>


#include <feel/feelpoly/dualbasis.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/orthonormalpolynomialset.hpp>
#include <feel/feelpoly/functionalset.hpp>
#include <feel/feelpoly/operations.hpp>
#include <feel/feelpoly/pointfunctionals.hpp>
#include <feel/feelpoly/momentfunctionals.hpp>
#include <feel/feelpoly/hcurlfunctionals.hpp>
#include <feel/feelpoly/hcurlinterpolation.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>
#include <feel/feeldiscr/doflayout.hpp>
#include <feel/feelpoly/fe.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/fefamilyshell.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

#include <feel/feelpoly/hcurlpolynomialset.hpp>

namespace Feel
{
namespace detail
{
/**
 * @brief Build the change of Nedelec triangular-face moment coordinates.
 *
 * The face test space is `[P_degree]^2`, ordered by tangent, total degree,
 * then exponent.  The returned matrix maps canonical functional values to
 * values in the permuted local tangent/monomial frame.  Construction is a
 * small reference-cell operation cached in the core library.
 *
 * @param vertexPermutation Local-to-canonical face vertex indices.
 * @param degree Degree of each scalar tangential test polynomial.
 * @return Cached immutable matrix acting on canonical face functional values.
 * @pre `vertexPermutation` is a permutation of `{0,1,2}`.
 */
[[nodiscard]] std::shared_ptr<DofTransform::MatrixStorage const>
nedelecTriangularFaceMomentTransform(
    std::array<uint16_type,3> const& vertexPermutation,
    uint16_type degree );

/**
 * @brief Build the NED1 face-moment action from a triangular-face permutation.
 * @tparam FacePermutation Triangular-face permutation type.
 * @param permutation Local-to-canonical face permutation.
 * @param degree Degree of each scalar tangential test polynomial.
 * @return Cached immutable face functional transformation matrix.
 * @pre `permutation.value()` is in `[IDENTITY,N_PERMUTATIONS)`.
 */
template<typename FacePermutation>
[[nodiscard]] std::shared_ptr<DofTransform::MatrixStorage const>
nedelecTriangularFaceMomentTransform( FacePermutation permutation, uint16_type degree )
{
    std::array<std::array<uint16_type,3>,7> const vertexPermutations{{
        {{0,1,2}}, {{0,1,2}}, {{2,1,0}}, {{0,2,1}},
        {{1,0,2}}, {{1,2,0}}, {{2,0,1}}
    }};
    auto const value = static_cast<uint16_type>( permutation.value() );
    CHECK( value >= FacePermutation::IDENTITY &&
           value < FacePermutation::N_PERMUTATIONS )
        << "invalid triangular face permutation " << value;
    return nedelecTriangularFaceMomentTransform(
        vertexPermutations[value], degree );
}

/**
 * @brief Build the triangular-face RT moment action used by NED2.
 *
 * The test basis is RT_{k-1} for public NED2 order `k`. Its pullback is the
 * contravariant two-dimensional Piola map, unlike NED1's componentwise
 * polynomial tangent basis.
 *
 * @param vertexPermutation Local-to-canonical face vertex indices.
 * @param publicOrder Public NED2 family order `k`.
 * @return Cached immutable matrix acting on canonical RT face moments.
 * @pre `vertexPermutation` is a permutation of `{0,1,2}` and `publicOrder > 0`.
 * @throws std::length_error If the reference sampling indices exceed local index storage.
 */
[[nodiscard]] std::shared_ptr<DofTransform::MatrixStorage const>
nedelecTriangularFaceRTMomentTransform(
    std::array<uint16_type,3> const& vertexPermutation,
    uint16_type publicOrder );

/**
 * @brief Return the explicit local-to-face-entity vertex permutation.
 *
 * Deriving the action from vertex identifiers avoids depending on the
 * historical direction convention of `facePermutation()`.
 *
 * @tparam ElementType Mesh element exposing face and point topology.
 * @param element Element containing the triangular face.
 * @param localFace Element-local face identifier.
 * @return Canonical face vertex index for each of its three local vertices.
 * @pre `localFace` identifies a triangular face of `element`.
 */
template<typename ElementType>
[[nodiscard]] std::array<uint16_type,3>
nedelecLocalToCanonicalFaceVertices( ElementType const& element,
                                     uint16_type localFace )
{
    std::array<uint16_type,3> permutation{};
    auto const& face = element.face( localFace );
    for ( uint16_type localVertex = 0; localVertex < 3; ++localVertex )
    {
        auto const pointId = element.point(
            element.fToP( localFace, localVertex ) ).id();
        uint16_type canonicalVertex = 0;
        for ( ; canonicalVertex < 3; ++canonicalVertex )
            if ( face.point( canonicalVertex ).id() == pointId )
                break;
        CHECK_LT( canonicalVertex, 3 )
            << "local face vertex is absent from its canonical face entity";
        permutation[localVertex] = canonicalVertex;
    }
    return permutation;
}
} // namespace detail

/**
 * @enum NedelecKind
 * @brief Mathematical Nedelec family variant.
 */
enum class NedelecKind
{
    NED1, /**< First-kind/trimmed H(curl) polynomial family. */
    NED2  /**< Second-kind/full H(curl) polynomial family. */
};

/**
 * @brief Whether a simplex Nedelec combination is implemented and validated.
 * @tparam Dim reference-simplex dimension
 * @tparam Order public Feel++ family order
 * @tparam Kind first- or second-kind family
 */
template<uint16_type Dim, uint16_type Order, NedelecKind Kind>
inline constexpr bool nedelecSimplexIsSupported =
    ( Kind == NedelecKind::NED1 ) ? ( Dim == 2 || ( Dim == 3 && Order <= 3 ) )
                                 : ( Dim == 2 || ( Dim == 3 && Order <= 3 ) );

/**
 * @brief Mathematical order and dimension traits for simplex Nedelec spaces.
 *
 * Dimension formulae describe the intended complete Ciarlet family. The
 * separate `isSupported` flag records the smaller envelope implemented and
 * validated today.
 *
 * @tparam Kind first- or second-kind Nedelec family
 * @tparam Dim reference-simplex dimension
 * @tparam Order public Feel++ family order
 */
template<NedelecKind Kind, uint16_type Dim, uint16_type Order>
struct NedelecTraits
{
    static_assert( Dim == 2 || Dim == 3, "Nedelec simplex traits support dimensions 2 and 3." );

    static constexpr uint16_type publicOrder = Order; /**< Feel++ family order. */
    static constexpr uint16_type polynomialDegree = static_cast<uint16_type>( Order + 1 ); /**< Degree of the ambient polynomial space. */
    /** Number of ordered tangential DoFs attached to one edge. */
    static constexpr uint16_type edgeDof = Kind == NedelecKind::NED1
                                                ? static_cast<uint16_type>( Order + 1 )
                                                : static_cast<uint16_type>( Order + 2 );
    /** Number of face-interior tangential DoFs attached to one triangular face. */
    static constexpr uint16_type faceDof = Dim == 3
                                                ? ( Kind == NedelecKind::NED1
                                                        ? static_cast<uint16_type>( Order * ( Order + 1 ) )
                                                        : static_cast<uint16_type>( Order * ( Order + 2 ) ) )
                                                : 0;
    /** Number of cell-interior DoFs. */
    static constexpr uint16_type cellDof = Kind == NedelecKind::NED1
                                                ? ( Dim == 2
                                                        ? static_cast<uint16_type>( Order * ( Order + 1 ) )
                                                        : ( Order > 1
                                                                ? static_cast<uint16_type>( ( Order + 1 ) * Order * ( Order - 1 ) / 2 )
                                                                : 0 ) )
                                                : ( Dim == 2
                                                        ? static_cast<uint16_type>( Order * ( Order + 2 ) )
                                                        : ( Order > 1
                                                                ? static_cast<uint16_type>( Order * ( Order - 1 ) * ( Order + 2 ) / 2 )
                                                                : 0 ) );
    /** Mathematical dimension of the complete reference polynomial space. */
    static constexpr uint16_type totalDof = Kind == NedelecKind::NED1
                                                 ? ( Dim == 2
                                                         ? static_cast<uint16_type>( ( Order + 1 ) * ( Order + 3 ) )
                                                         : static_cast<uint16_type>( ( Order + 1 ) * ( Order + 3 ) * ( Order + 4 ) / 2 ) )
                                                 : ( Dim == 2
                                                         ? static_cast<uint16_type>( ( Order + 2 ) * ( Order + 3 ) )
                                                         : static_cast<uint16_type>( ( Order + 2 ) * ( Order + 3 ) * ( Order + 4 ) / 2 ) );
    /** Whether this combination is implemented and validated in Feel++. */
    static constexpr bool isSupported = nedelecSimplexIsSupported<Dim, Order, Kind>;
};

namespace detail
{
/**
 * @brief Convert a public NED1 order to its ambient polynomial degree.
 * @param publicOrder Public NED1 family order `k`.
 * @return The checked degree `k+1`.
 * @throws std::length_error If the degree exceeds local index storage.
 */
[[nodiscard]] constexpr uint16_type nedelec1InternalOrder( uint16_type publicOrder )
{
    return checkedLocalIndexValue(
        static_cast<size_type>( publicOrder ) + 1 );
}

/**
 * @brief Compute the NED1 local dimension before narrowing to local index storage.
 * @param dim Simplex dimension, either 2 or 3.
 * @param publicOrder Public NED1 family order.
 * @return Number of element-local degrees of freedom as a `size_type`.
 */
[[nodiscard]] constexpr size_type
nedelec1TotalDofSize( uint16_type dim,
                       uint16_type publicOrder ) noexcept
{
    const size_type order = publicOrder;
    return dim == 2
               ? ( order + 1 ) * ( order + 3 )
               : ( order + 1 ) * ( order + 3 ) *
                     ( order + 4 ) / 2;
}

/**
 * @brief Compute the NED2 local dimension before narrowing to local index storage.
 * @param dim Simplex dimension, either 2 or 3.
 * @param publicOrder Public NED2 family order.
 * @return Dimension of the full vector polynomial space as a `size_type`.
 */
[[nodiscard]] constexpr size_type
nedelec2TotalDofSize( uint16_type dim,
                       uint16_type publicOrder ) noexcept
{
    const size_type order = publicOrder;
    return dim == 2
               ? ( order + 2 ) * ( order + 3 )
               : ( order + 2 ) * ( order + 3 ) *
                     ( order + 4 ) / 2;
}

/**
 * @brief Check that a runtime Nedelec order and local dimension fit index storage.
 * @tparam Kind Nedelec family whose dimension formula is checked.
 * @param dim Simplex dimension, either 2 or 3.
 * @param publicOrder Requested public family order.
 * @return The unchanged, validated order.
 * @throws std::length_error If the ambient degree or local dimension exceeds local index storage.
 */
template<NedelecKind Kind>
[[nodiscard]] constexpr RuntimeOrder
validatedNedelecSimplexOrder( uint16_type dim,
                              RuntimeOrder publicOrder )
{
    (void)checkedLocalIndexValue(
        static_cast<size_type>( publicOrder.value ) + 1 );
    if constexpr ( Kind == NedelecKind::NED1 )
        (void)checkedLocalIndexValue(
            nedelec1TotalDofSize( dim, publicOrder.value ) );
    else
        (void)checkedLocalIndexValue(
            nedelec2TotalDofSize( dim, publicOrder.value ) );
    return publicOrder;
}

/**
 * @brief Count tangential point evaluations on one NED1 edge.
 * @param publicOrder Public NED1 family order `k`.
 * @return The checked count `k+1`.
 * @throws std::length_error If the count exceeds local index storage.
 */
[[nodiscard]] constexpr uint16_type nedelec1EdgeDof( uint16_type publicOrder )
{
    return nedelec1InternalOrder( publicOrder );
}

/**
 * @brief Count NED1 moments in the interior of a triangular entity.
 * @param dim Simplex dimension; the count is the same for a 2D cell or a 3D face.
 * @param publicOrder Public NED1 family order `k`.
 * @return The checked count `k*(k+1)`.
 * @throws std::length_error If the count exceeds local index storage.
 */
[[nodiscard]] constexpr uint16_type nedelec1FaceDof( uint16_type dim, uint16_type publicOrder )
{
    const size_type order = publicOrder;
    return checkedLocalIndexValue( order * ( order + 1 ) );
}

/**
 * @brief Count NED1 moments in the interior of a tetrahedron.
 * @param dim Simplex dimension, either 2 or 3.
 * @param publicOrder Public NED1 family order `k`.
 * @return The checked count `(k+1)*k*(k-1)/2`, or zero for 2D or `k < 2`.
 * @throws std::length_error If the count exceeds local index storage.
 */
[[nodiscard]] constexpr uint16_type nedelec1VolumeDof( uint16_type dim, uint16_type publicOrder )
{
    const size_type order = publicOrder;
    return dim == 3 && publicOrder > 1
               ? checkedLocalIndexValue(
                     ( order + 1 ) * order * ( order - 1 ) / 2 )
               : 0;
}

/**
 * @brief Return the checked local dimension of a NED1 simplex element.
 * @param dim Simplex dimension, either 2 or 3.
 * @param publicOrder Public NED1 family order.
 * @return Number of element-local degrees of freedom.
 * @throws std::length_error If the dimension exceeds local index storage.
 */
[[nodiscard]] constexpr uint16_type nedelec1TotalDof( uint16_type dim, uint16_type publicOrder )
{
    return checkedLocalIndexValue(
        nedelec1TotalDofSize( dim, publicOrder ) );
}

/** @return ambient polynomial degree of NED2 public order `publicOrder`. */
[[nodiscard]] constexpr uint16_type nedelec2InternalOrder( uint16_type publicOrder )
{
    return checkedLocalIndexValue(
        static_cast<size_type>( publicOrder ) + 1 );
}

/** @return number of pointwise tangential DoFs on one NED2 edge. */
[[nodiscard]] constexpr uint16_type nedelec2EdgeDof( uint16_type publicOrder )
{
    return checkedLocalIndexValue(
        static_cast<size_type>( publicOrder ) + 2 );
}

/** @return number of NED2 face-interior DoFs on one triangular face. */
[[nodiscard]] constexpr uint16_type nedelec2FaceDof( uint16_type dim, uint16_type publicOrder )
{
    const size_type order = publicOrder;
    return dim == 3
               ? checkedLocalIndexValue( order * ( order + 2 ) )
               : 0;
}

/** @return number of NED2 cell-interior DoFs. */
[[nodiscard]] constexpr uint16_type nedelec2CellDof( uint16_type dim, uint16_type publicOrder )
{
    const size_type order = publicOrder;
    if ( dim == 2 )
        return checkedLocalIndexValue( order * ( order + 2 ) );
    return publicOrder > 1
               ? checkedLocalIndexValue(
                     order * ( order - 1 ) * ( order + 2 ) / 2 )
               : 0;
}

/** @return total local dimension of the full-vector NED2 polynomial space. */
[[nodiscard]] constexpr uint16_type nedelec2TotalDof( uint16_type dim, uint16_type publicOrder )
{
    return checkedLocalIndexValue(
        nedelec2TotalDofSize( dim, publicOrder ) );
}

template<typename P>
struct times_rotx
{
    typedef typename P::value_type value_type;
    typedef typename P::points_type points_type;
    /**
     * @brief Select a component of the two-dimensional rotated-coordinate product.
     * @param p Scalar polynomial to multiply by `(y,-x)`.
     * @param c Component index, either 0 or 1.
     */
    times_rotx ( P const& p, int c  )
        :
        M_p ( p ),
        M_c( c )
    {
        //VLOG(4) << "component : " << c << std::endl;
    }
    /**
     * @brief Evaluate the selected component of `(y,-x)*p`.
     * @param __pts Reference points stored as columns.
     * @return One product value per point.
     */
    typename ublas::vector<value_type> operator() ( points_type const& __pts ) const
    {
#if 0
        VLOG(4) << "times_rotx(pts) : " << __pts << std::endl;
        VLOG(4) << "times_rotx(pts) : " << M_p.evaluate( __pts ) << std::endl;
        VLOG(4) << "times_rotx(coeff) : " << M_p.coefficients() << std::endl;
#endif

        // __pts[c] * p( __pts )
        if ( M_c == 0 )
            return ublas::element_prod( ublas::row( __pts, 1 ),
                                        ublas::row( M_p.evaluate( __pts ), 0 ) );

        // c == 1
        return ublas::element_prod( -ublas::row( __pts, 0 ),
                                    ublas::row( M_p.evaluate( __pts ), 0 ) );
    }
    //P const& M_p;
    P M_p;
    int M_c;
};

template<typename P>
struct times_crossx_component
{
    using value_type = typename P::value_type;
    using points_type = typename P::points_type;
    /**
     * @brief Select a component of the cross product `(e_axis cross x)*p`.
     * @param p Scalar polynomial multiplying the cross product.
     * @param axis Cartesian unit-vector index, in `[0,3)`.
     * @param component Cartesian output component, in `[0,3)`.
     */
    times_crossx_component( P const& p, int axis, int component )
        : M_p( p ), M_axis( axis ), M_component( component ) {}

    /**
     * @brief Evaluate the selected cross-product component at reference points.
     * @param points Three-dimensional reference points stored as columns.
     * @return One component value per point, or zeros when axis equals component.
     */
    ublas::vector<value_type> operator()( points_type const& points ) const
    {
        ublas::vector<value_type> result( points.size2(), value_type( 0 ) );
        const int coordinate = 3-M_axis-M_component;
        if ( M_axis == M_component )
            return result;
        const value_type sign = ( ( M_axis+1 )%3 == coordinate ) ? value_type( 1 ) : value_type( -1 );
        result = sign*ublas::element_prod( ublas::row( points, coordinate ),
                                           ublas::row( M_p.evaluate( points ), 0 ) );
        return result;
    }
    P M_p;
    int M_axis;
    int M_component;
};

/**
 * @brief Construct the NED1 primal coefficient matrix for static or runtime storage.
 *
 * `AmbientVector` and `AmbientScalar` select only the polynomial container
 * storage policy. The mathematical construction is shared verbatim by the
 * compile-time and runtime finite elements.
 *
 * @tparam Dim Simplex dimension, either 2 or 3.
 * @tparam AmbientVector Ambient vector polynomial storage type.
 * @tparam AmbientScalar Ambient scalar polynomial storage type.
 * @tparam VectorPolynomialSet Resulting vector polynomial container type.
 * @tparam ScalarPolynomialSet Scalar container used to construct the enrichment.
 * @param ambientVector Vector polynomial basis through degree `publicOrder+1`.
 * @param ambientScalar Scalar polynomial basis through degree `publicOrder+1`.
 * @param publicOrder Public NED1 family order.
 * @return Coefficient matrix of an independent basis in the ambient coordinates.
 */
template<uint16_type Dim, typename AmbientVector, typename AmbientScalar,
         typename VectorPolynomialSet, typename ScalarPolynomialSet>
[[nodiscard]] auto
buildNedelec1PrimalCoefficients( AmbientVector const& ambientVector,
                                AmbientScalar const& ambientScalar,
                                uint16_type publicOrder )
{
    using value_type = typename AmbientVector::value_type;
    using convex_type = typename AmbientVector::convex_type;
    using scalar_polynomial_type = typename ScalarPolynomialSet::polynomial_type;
    const uint16_type internalOrder = nedelec1InternalOrder( publicOrder );
    const uint16_type dimPk = convex_type::polyDims( publicOrder );
    const uint16_type dimPkm1 = publicOrder == 0 ? 0 : convex_type::polyDims( publicOrder-1 );

    auto const vectorSource = ambientVector.polynomialsUpToDimension( dimPk );
    auto const scalarSource = ambientScalar.polynomialsUpToDimension( dimPk );
    VectorPolynomialSet pkVector = polynomialSetFrom<VectorPolynomialSet>(
        vectorSource, internalOrder );
    ScalarPolynomialSet pkScalar = polynomialSetFrom<ScalarPolynomialSet>(
        scalarSource, internalOrder );
    IMGeneral<Dim, value_type, Simplex> quadrature( static_cast<uint16_type>( 2*internalOrder ) );

    if constexpr ( Dim == 2 )
    {
        ublas::matrix<value_type> enrichmentCoefficients(
            Dim*( dimPk-dimPkm1 ), pkScalar.coeff().size2() );
        for ( int scalar = dimPkm1, candidate = 0; scalar < dimPk; ++scalar, ++candidate )
            for ( int component = 0; component < Dim; ++component )
            {
                times_rotx<scalar_polynomial_type> rotated( pkScalar.polynomial( scalar ), component );
                ublas::row( enrichmentCoefficients, candidate*Dim+component ) =
                    ublas::row( Feel::project( ambientScalar, rotated, quadrature ).coeff(), 0 );
            }
        VectorPolynomialSet enrichment(
            typename AmbientVector::basis_type(), enrichmentCoefficients, true );
        if constexpr ( requires( VectorPolynomialSet& p ) { p.setOrder( internalOrder ); } )
            enrichment.setOrder( internalOrder );
        auto result = unite( pkVector, enrichment );
        CHECK_EQ( result.polynomialDimension(), nedelec1TotalDof( Dim, publicOrder ) );
        return result.coeff();
    }
    else
    {
        const uint16_type homogeneousDimension = dimPk-dimPkm1;
        ublas::matrix<value_type> enrichmentCoefficients(
            Dim*Dim*homogeneousDimension, pkScalar.coeff().size2() );
        enrichmentCoefficients.clear();
        for ( uint16_type scalar = dimPkm1, candidate = 0; scalar < dimPk; ++scalar )
            for ( uint16_type axis = 0; axis < Dim; ++axis, ++candidate )
                for ( uint16_type component = 0; component < Dim; ++component )
                {
                    times_crossx_component<scalar_polynomial_type> crossComponent(
                        pkScalar.polynomial( scalar ), axis, component );
                    ublas::row( enrichmentCoefficients, candidate*Dim+component ) =
                        ublas::row( Feel::project( ambientScalar, crossComponent, quadrature ).coeff(), 0 );
                }
        VectorPolynomialSet enrichment(
            typename AmbientVector::basis_type(), enrichmentCoefficients, true );
        ublas::matrix<value_type> stacked( pkVector.coeff().size1()+enrichment.coeff().size1(),
                                           pkVector.coeff().size2() );
        ublas::project( stacked, ublas::range( 0, pkVector.coeff().size1() ),
                        ublas::range( 0, stacked.size2() ) ) = pkVector.coeff();
        ublas::project( stacked, ublas::range( pkVector.coeff().size1(), stacked.size1() ),
                        ublas::range( 0, stacked.size2() ) ) = enrichment.coeff();
        auto matrixCoefficients = VectorPolynomialSet::polyset_type::toMatrix( stacked );
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> candidates(
            matrixCoefficients.data().begin(), matrixCoefficients.size1(), matrixCoefficients.size2() );
        Eigen::JacobiSVD<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic>> svd(
            candidates, Eigen::ComputeThinV );
        const value_type tolerance = std::numeric_limits<value_type>::epsilon()*
                                     std::max( candidates.rows(), candidates.cols() )*svd.singularValues()( 0 );
        const int rank = ( svd.singularValues().array() > tolerance ).count();
        CHECK_EQ( rank, nedelec1TotalDof( Dim, publicOrder ) );
        ublas::matrix<value_type> independent( rank, matrixCoefficients.size2() );
        for ( int i = 0; i < rank; ++i )
            for ( int j = 0; j < independent.size2(); ++j )
                independent( i, j ) = svd.matrixV()( j, i );
        return VectorPolynomialSet::polyset_type::toType( independent );
    }
}

/** @brief Sampling block for one triangular NED1 face-moment entity. */
struct Nedelec1FaceSamplingBlock
{
    uint16_type face = 0; /**< Reference tetrahedron face identifier. */
    uint16_type firstDof = 0; /**< First element-local face-moment DoF. */
    uint16_type momentCount = 0; /**< Number of tangential face moments. */
    uint32_type elementPointOffset = 0; /**< First sample in the element point table. */
    uint32_type facetPointOffset = 0; /**< First sample in the facet point table. */
    uint32_type pointCount = 0; /**< Number of quadrature samples. */
};

/**
 * @brief Shared expression-sampling layout for a static or runtime NED1 dual.
 *
 * DoF offsets and expression-point offsets are deliberately distinct: a
 * moment functional generally needs several quadrature samples.
 */
struct Nedelec1SamplingLayout
{
    uint32_type edgePointCount = 0; /**< Total edge point-evaluation samples. */
    std::vector<Nedelec1FaceSamplingBlock> faces; /**< Triangular face blocks. */
    uint16_type cellFirstDof = 0; /**< First cell-moment local DoF. */
    uint16_type cellMomentCount = 0; /**< Number of cell moments. */
    uint32_type cellPointOffset = 0; /**< First cell quadrature sample. */
    uint32_type cellPointCount = 0; /**< Number of cell quadrature samples. */
    uint32_type totalPointCount = 0; /**< Complete element sampling cardinality. */
};

/**
 * @brief Build the ordered NED1 dual functional set for static or runtime order.
 * @tparam Space NED1 primal polynomial space type.
 * @tparam ReferenceConvex Reference simplex type.
 * @tparam AmbientVector Ambient vector polynomial storage type.
 * @tparam VectorPolynomialSet Container for the cell moment polynomials.
 * @tparam EdgePointsProvider Callable returning reference sampling points for an edge.
 * @param primal Primal space on which the functionals act.
 * @param referenceConvex Reference simplex providing entity topology and geometry.
 * @param publicOrder Public NED1 family order.
 * @param ambientVector Ambient vector basis used for cell moments.
 * @param[out] elementPoints Complete element expression-sampling table.
 * @param[out] facetPoints Sampling tables for each facet closure, sized by the caller.
 * @param[out] samplingLayout Offsets and cardinalities of the sampling blocks.
 * @param edgePointsProvider Provider of the ordered tangential edge evaluation points.
 * @return Functionals ordered by edge, triangular face, then cell interior.
 */
template<typename Space, typename ReferenceConvex, typename AmbientVector,
         typename VectorPolynomialSet, typename EdgePointsProvider>
[[nodiscard]] std::vector<Functional<Space>>
buildNedelec1Functionals( Space const& primal,
                          ReferenceConvex const& referenceConvex,
                          uint16_type publicOrder,
                          AmbientVector const& ambientVector,
                          typename Space::points_type& elementPoints,
                          std::vector<typename Space::points_type>& facetPoints,
                          Nedelec1SamplingLayout& samplingLayout,
                          EdgePointsProvider&& edgePointsProvider )
{
    static constexpr uint16_type Dim = Space::nDim;
    using value_type = typename Space::value_type;
    using convex_type = typename ReferenceConvex::super;
    const uint16_type internalOrder = nedelec1InternalOrder( publicOrder );
    const uint16_type edgeDof = nedelec1EdgeDof( publicOrder );
    const uint16_type faceDof = Dim == 3 ? nedelec1FaceDof( Dim, publicOrder ) : 0;
    using face_type = typename convex_type::topological_face_type;
    std::vector<Functional<Space>> functionals;
    functionals.reserve( nedelec1TotalDof( Dim, publicOrder ) );

    samplingLayout = {};
    samplingLayout.edgePointCount =
        static_cast<uint32_type>( convex_type::numEdges*edgeDof );
    elementPoints.resize( Dim, samplingLayout.edgePointCount, false );
    const uint32_type edgePointsPerFacet = Dim == 2
        ? edgeDof
        : static_cast<uint32_type>( face_type::numEdges*edgeDof );
    for ( auto& points : facetPoints )
        points.resize( Dim, edgePointsPerFacet, false );

    appendHCurlEdgeTangentialPointFunctionals(
        primal, referenceConvex, edgeDof, elementPoints, facetPoints,
        functionals, std::forward<EdgePointsProvider>( edgePointsProvider ) );
    if constexpr ( Dim == 3 )
    {
        const int facePolynomialDegree = static_cast<int>( publicOrder )-1;
        if ( facePolynomialDegree >= 0 )
        {
            IMGeneral<3, value_type, Simplex> faceQuadrature(
                hcurlTriangularFaceMomentQuadratureOrder( facePolynomialDegree ) );
            uint32_type pointOffset = samplingLayout.edgePointCount;
            samplingLayout.faces.reserve( convex_type::numGeometricFaces );
            for ( uint16_type face = 0; face < convex_type::numGeometricFaces; ++face )
            {
                auto const& points = faceQuadrature.fpoints( face );
                const uint32_type pointCount = faceQuadrature.nPointsOnFace( face );
                const uint32_type nextPointOffset = pointOffset+pointCount;
                elementPoints.resize( Dim, nextPointOffset, true );
                ublas::subrange( elementPoints, 0, Dim, pointOffset, nextPointOffset ) = points;

                facetPoints[face].resize( Dim, edgePointsPerFacet+pointCount, true );
                ublas::subrange( facetPoints[face], 0, Dim,
                                 edgePointsPerFacet, edgePointsPerFacet+pointCount ) = points;
                samplingLayout.faces.push_back( Nedelec1FaceSamplingBlock{
                    .face = face,
                    .firstDof = static_cast<uint16_type>(
                        convex_type::numEdges*edgeDof+face*faceDof ),
                    .momentCount = faceDof,
                    .elementPointOffset = pointOffset,
                    .facetPointOffset = edgePointsPerFacet,
                    .pointCount = pointCount } );
                pointOffset = nextPointOffset;
            }
        }
        appendHCurlTriangularFaceMoments(
            primal, referenceConvex, facePolynomialDegree, functionals );
    }

    const int momentDegree = Dim == 2 ? static_cast<int>( publicOrder )-1
                                      : static_cast<int>( publicOrder )-2;
    samplingLayout.cellFirstDof = static_cast<uint16_type>(
        convex_type::numEdges*edgeDof+
        ( Dim == 3 ? convex_type::numGeometricFaces*faceDof : 0 ) );
    samplingLayout.cellMomentCount = Dim == 2
        ? nedelec1FaceDof( Dim, publicOrder )
        : nedelec1VolumeDof( Dim, publicOrder );
    samplingLayout.cellPointOffset = elementPoints.size2();
    if ( momentDegree >= 0 )
    {
        const uint16_type dimMoment = convex_type::polyDims(
            static_cast<uint16_type>( momentDegree ) );
        auto const source = ambientVector.polynomialsUpToDimension( dimMoment );
        VectorPolynomialSet moments = polynomialSetFrom<VectorPolynomialSet>(
            source, internalOrder );

        IMGeneral<Dim, value_type, Simplex> quadrature(
            static_cast<uint16_type>( 2*( internalOrder+1 ) ) );
        samplingLayout.cellPointCount = quadrature.nPoints();
        const uint32_type requiredPointCount =
            samplingLayout.cellPointOffset+samplingLayout.cellPointCount;
        elementPoints.resize( Dim, requiredPointCount, true );
        ublas::subrange( elementPoints, 0, Dim,
                         samplingLayout.cellPointOffset, requiredPointCount ) = quadrature.points();

        using moment_type = functional::IntegralMoment<Space, VectorPolynomialSet>;
        for ( int i = 0; i < moments.polynomialDimension(); ++i )
            functionals.emplace_back( moment_type( primal, moments.polynomial( i ) ) );
    }
    samplingLayout.totalPointCount = elementPoints.size2();
    CHECK_EQ( functionals.size(), nedelec1TotalDof( Dim, publicOrder ) );
    return functionals;
}

/**
 * @brief Build the canonical point-variant NED2 dual for static/runtime order.
 *
 * FIAT's convention is followed exactly: `k+2` edge tangential evaluations
 * for Feel++ public order `k`, RT_{k-1} moments on triangular faces and 2D
 * cells, and RT_{k-2} moments in tetrahedral cells.
 *
 * @tparam Space NED2 primal polynomial space type.
 * @tparam ReferenceConvex Reference simplex type.
 * @tparam EdgePointsProvider Callable returning reference sampling points for an edge.
 * @param primal Primal space on which the functionals act.
 * @param referenceConvex Reference simplex providing entity topology and geometry.
 * @param publicOrder Public NED2 family order `k`.
 * @param[out] elementPoints Complete element expression-sampling table.
 * @param[out] facetPoints Sampling tables for each facet closure, sized by the caller.
 * @param[out] samplingLayout Offsets and cardinalities of the sampling blocks.
 * @param edgePointsProvider Provider of the ordered tangential edge evaluation points.
 * @return Functionals ordered by edge, triangular face, then cell interior.
 */
template<typename Space, typename ReferenceConvex, typename EdgePointsProvider>
[[nodiscard]] std::vector<Functional<Space>>
buildNedelec2Functionals( Space const& primal,
                          ReferenceConvex const& referenceConvex,
                          uint16_type publicOrder,
                          typename Space::points_type& elementPoints,
                          std::vector<typename Space::points_type>& facetPoints,
                          Nedelec1SamplingLayout& samplingLayout,
                          EdgePointsProvider&& edgePointsProvider )
{
    static constexpr uint16_type Dim = Space::nDim;
    using value_type = typename Space::value_type;
    using convex_type = typename ReferenceConvex::super;
    using face_type = typename convex_type::topological_face_type;
    using cell_moment_space_type = RaviartThomasDynamicPolynomialSet<Dim, value_type>;
    using face_moment_space_type = RaviartThomasDynamicPolynomialSet<2, value_type>;

    const uint16_type internalOrder = nedelec2InternalOrder( publicOrder );
    const uint16_type edgeDof = nedelec2EdgeDof( publicOrder );
    const uint16_type faceDof = nedelec2FaceDof( Dim, publicOrder );
    const uint16_type cellDof = nedelec2CellDof( Dim, publicOrder );
    std::vector<Functional<Space>> functionals;
    functionals.reserve( nedelec2TotalDof( Dim, publicOrder ) );

    samplingLayout = {};
    samplingLayout.edgePointCount =
        static_cast<uint32_type>( convex_type::numEdges*edgeDof );
    elementPoints.resize( Dim, samplingLayout.edgePointCount, false );
    const uint32_type edgePointsPerFacet = Dim == 2
        ? edgeDof
        : static_cast<uint32_type>( face_type::numEdges*edgeDof );
    for ( auto& points : facetPoints )
        points.resize( Dim, edgePointsPerFacet, false );

    appendHCurlEdgeTangentialPointFunctionals(
        primal, referenceConvex, edgeDof, elementPoints, facetPoints,
        functionals, std::forward<EdgePointsProvider>( edgePointsProvider ) );

    if constexpr ( Dim == 3 )
    {
        if ( publicOrder > 0 )
        {
            face_moment_space_type moments(
                RuntimeOrder{ static_cast<uint16_type>( publicOrder-1 ) } );
            const uint16_type quadratureOrder = static_cast<uint16_type>(
                2*( moments.publicOrder()+4 ) );
            IMGeneral<3, value_type, Simplex> quadrature( quadratureOrder );
            uint32_type pointOffset = samplingLayout.edgePointCount;
            samplingLayout.faces.reserve( convex_type::numGeometricFaces );
            for ( uint16_type face = 0; face < convex_type::numGeometricFaces; ++face )
            {
                auto const& points = quadrature.fpoints( face );
                const uint32_type pointCount = quadrature.nPointsOnFace( face );
                const uint32_type nextPointOffset = pointOffset+pointCount;
                elementPoints.resize( Dim, nextPointOffset, true );
                ublas::subrange( elementPoints, 0, Dim,
                                 pointOffset, nextPointOffset ) = points;
                facetPoints[face].resize( Dim, edgePointsPerFacet+pointCount, true );
                ublas::subrange( facetPoints[face], 0, Dim,
                                 edgePointsPerFacet, edgePointsPerFacet+pointCount ) = points;
                samplingLayout.faces.push_back( Nedelec1FaceSamplingBlock{
                    .face = face,
                    .firstDof = static_cast<uint16_type>(
                        convex_type::numEdges*edgeDof+face*faceDof ),
                    .momentCount = faceDof,
                    .elementPointOffset = pointOffset,
                    .facetPointOffset = edgePointsPerFacet,
                    .pointCount = pointCount } );
                pointOffset = nextPointOffset;
            }
            appendHCurlTriangularFaceRTMoments(
                primal, referenceConvex, moments, functionals );
        }
    }

    samplingLayout.cellFirstDof = static_cast<uint16_type>(
        convex_type::numEdges*edgeDof+
        ( Dim == 3 ? convex_type::numGeometricFaces*faceDof : 0 ) );
    samplingLayout.cellMomentCount = cellDof;
    samplingLayout.cellPointOffset = elementPoints.size2();
    const int cellMomentOrder = Dim == 2 ? static_cast<int>( publicOrder )-1
                                         : static_cast<int>( publicOrder )-2;
    if ( cellMomentOrder >= 0 )
    {
        cell_moment_space_type moments(
            RuntimeOrder{ static_cast<uint16_type>( cellMomentOrder ) } );
        IMGeneral<Dim, value_type, Simplex> quadrature(
            static_cast<uint16_type>( 2*( internalOrder+2 ) ) );
        samplingLayout.cellPointCount = quadrature.nPoints();
        const uint32_type requiredPointCount =
            samplingLayout.cellPointOffset+samplingLayout.cellPointCount;
        elementPoints.resize( Dim, requiredPointCount, true );
        ublas::subrange( elementPoints, 0, Dim,
                         samplingLayout.cellPointOffset, requiredPointCount ) = quadrature.points();
        appendHCurlCellMoments<Dim>(
            primal, moments, static_cast<uint16_type>( 2*( internalOrder+2 ) ),
            functionals );
    }
    samplingLayout.totalPointCount = elementPoints.size2();
    CHECK_EQ( functionals.size(), nedelec2TotalDof( Dim, publicOrder ) );
    return functionals;
}

#if 0
template< class T >
struct extract_all_poly_indices
{
    T start;
    extract_all_poly_indices( T comp, T N )
        :
        start( comp*N )
    { }

    T operator()()
    {
        return start++;
    }
};
#endif

}// detail

template<uint16_type N,
         uint16_type O,
         NedelecKind Kind = NedelecKind::NED1,
         typename T = double,
         uint16_type TheTAG = 0 >
class NedelecPolynomialSet{};

template<uint16_type O,
         typename T,
         uint16_type TheTAG>
class NedelecPolynomialSet<2,O,NedelecKind::NED1,T,TheTAG>
    :
    public Feel::detail::OrthonormalPolynomialSet<2, O+1, 2, Vectorial, T, TheTAG, Simplex>
{
    static const int N = 2;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> super;

public:
    static const uint16_type Om1 = (O==0)?0:O-1;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Vectorial, T, TheTAG, Simplex> Pk_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> Pkp1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, Om1, N, Vectorial, T, TheTAG, Simplex> Pkm1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Scalar, T, TheTAG, Simplex> Pk_s_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Scalar, T, TheTAG, Simplex> Pkp1_s_type;

    typedef PolynomialSet<typename super::basis_type,Vectorial> vectorial_polynomialset_type;
    typedef typename vectorial_polynomialset_type::polynomial_type vectorial_polynomial_type;
    typedef PolynomialSet<typename super::basis_type,Scalar> scalar_polynomialset_type;
    typedef typename scalar_polynomialset_type::polynomial_type scalar_polynomial_type;

    typedef NedelecPolynomialSet<N, O, NedelecKind::NED1,T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::convex_type convex_type;
    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const uint16_type nComponents = super::nComponents;
    static inline const bool is_product = false;
    /**
     * @brief Construct the primal basis for the compile-time Nedelec family and order.
     */
    NedelecPolynomialSet()
        :
        super()
    {
        Pkp1_v_type ambientVector;
        Pkp1_s_type ambientScalar;
        this->setCoefficient(
            Feel::detail::buildNedelec1PrimalCoefficients<N, Pkp1_v_type, Pkp1_s_type,
                                                          vectorial_polynomialset_type,
                                                          scalar_polynomialset_type>(
                ambientVector, ambientScalar, O ), true );
    }


};

template<uint16_type O,
         typename T,
         uint16_type TheTAG >
class NedelecPolynomialSet<3,O,NedelecKind::NED1,T,TheTAG>
    :
        public Feel::detail::OrthonormalPolynomialSet<3, O+1, 3, Vectorial, T, TheTAG, Simplex>
{
    static const int N = 3;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> super;

public:
    static const uint16_type Om1 = (O==0)?0:O-1;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Vectorial, T, TheTAG, Simplex> Pk_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> Pkp1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, Om1, N, Vectorial, T, TheTAG, Simplex> Pkm1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Scalar, T, TheTAG, Simplex> Pk_s_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Scalar, T, TheTAG, Simplex> Pkp1_s_type;

    typedef PolynomialSet<typename super::basis_type,Vectorial> vectorial_polynomialset_type;
    typedef typename vectorial_polynomialset_type::polynomial_type vectorial_polynomial_type;
    typedef PolynomialSet<typename super::basis_type,Scalar> scalar_polynomialset_type;
    typedef typename scalar_polynomialset_type::polynomial_type scalar_polynomial_type;

    typedef NedelecPolynomialSet<N, O, NedelecKind::NED1, T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::convex_type convex_type;
    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const uint16_type nComponents = super::nComponents;
    static inline const bool is_product = false;
    /**
     * @brief Construct the primal basis for the compile-time Nedelec family and order.
     */
    NedelecPolynomialSet()
        :
        super()
    {
        Pkp1_v_type ambientVector;
        Pkp1_s_type ambientScalar;
        this->setCoefficient(
            Feel::detail::buildNedelec1PrimalCoefficients<N, Pkp1_v_type, Pkp1_s_type,
                                                          vectorial_polynomialset_type,
                                                          scalar_polynomialset_type>(
                ambientVector, ambientScalar, O ), true );
    }


};

/** @brief Runtime-order first-kind Nedelec polynomial space on a simplex. */
template<uint16_type N, typename T = double, uint16_type TheTAG = 0>
class Nedelec1DynamicPolynomialSet
    : public Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Simplex>
{
    using super = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Simplex>;

public:
    using Pkp1_v_type = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Simplex>;
    using Pkp1_s_type = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Scalar, T, TheTAG, Simplex>;
    using vectorial_polynomialset_type = PolynomialSet<typename super::basis_type, Vectorial, Dynamic>;
    using scalar_polynomialset_type = PolynomialSet<typename super::basis_type, Scalar, Dynamic>;
    using scalar_polynomial_type = typename scalar_polynomialset_type::polynomial_type;
    using value_type = typename super::value_type;
    using convex_type = typename super::convex_type;

    static_assert( N == 2 || N == 3 );
    static inline const uint16_type nDim = N;
    static inline const int nOrder = Dynamic;
    static inline const uint16_type nComponents = N;
    static inline const bool is_product = false;

    /**
     * @brief Construct the lowest-order runtime Nedelec primal space.
     */
    Nedelec1DynamicPolynomialSet() : Nedelec1DynamicPolynomialSet( RuntimeOrder{ 0 } ) {}

    /**
     * @brief Construct a runtime Nedelec primal space of the requested order.
     * @param publicOrder Public Nedelec family order.
     * @throws std::length_error If the ambient degree or local dimension exceeds local index storage.
     */
    explicit Nedelec1DynamicPolynomialSet( RuntimeOrder publicOrder )
        : super( RuntimeOrder{
              Feel::detail::nedelec1InternalOrder(
                  Feel::detail::validatedNedelecSimplexOrder<
                      NedelecKind::NED1>( N, publicOrder ).value ) } )
    {
        build( publicOrder.value );
    }

    /**
     * @return Public NED1 family order.
     */
    [[nodiscard]] uint16_type publicOrder() const noexcept { return M_publicOrder; }
    /**
     * @return Ambient polynomial degree, equal to the public order plus one.
     */
    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        return Feel::detail::nedelec1InternalOrder( M_publicOrder );
    }

private:
    /**
     * @brief Build the independent NED1 polynomial coefficients in the ambient basis.
     * @param publicOrder Public NED1 family order.
     */
    void build( uint16_type publicOrder )
    {
        M_publicOrder = publicOrder;
        const uint16_type internal = internalOrder();
        Pkp1_v_type ambientVector( RuntimeOrder{ internal } );
        Pkp1_s_type ambientScalar( RuntimeOrder{ internal } );
        this->setCoefficient(
            Feel::detail::buildNedelec1PrimalCoefficients<N, Pkp1_v_type, Pkp1_s_type,
                                                          vectorial_polynomialset_type,
                                                          scalar_polynomialset_type>(
                ambientVector, ambientScalar, publicOrder ), true );
    }

    uint16_type M_publicOrder = 0;
};

/** @brief Runtime-order second-kind full vector polynomial space. */
template<uint16_type N, typename T = double, uint16_type TheTAG = 0>
class Nedelec2DynamicPolynomialSet
    : public Feel::detail::OrthonormalPolynomialSet<
          N, Dynamic, N, Vectorial, T, TheTAG, Simplex>
{
    using super = Feel::detail::OrthonormalPolynomialSet<
        N, Dynamic, N, Vectorial, T, TheTAG, Simplex>;

public:
    using value_type = typename super::value_type;
    using convex_type = typename super::convex_type;
    using points_type = typename super::points_type;
    using Pkp1_v_type = super;
    using vectorial_polynomialset_type = PolynomialSet<typename super::basis_type,
                                                        Vectorial, Dynamic>;
    static inline const uint16_type nDim = N;
    static inline const int nOrder = Dynamic;
    static inline const uint16_type nComponents = N;
    static inline const bool is_product = false;

    /**
     * @brief Construct the lowest-order runtime Nedelec primal space.
     */
    Nedelec2DynamicPolynomialSet()
        : Nedelec2DynamicPolynomialSet( RuntimeOrder{ 0 } )
    {}

    /**
     * @brief Construct a runtime Nedelec primal space of the requested order.
     * @param publicOrder Public Nedelec family order.
     * @throws std::length_error If the ambient degree or local dimension exceeds local index storage.
     */
    explicit Nedelec2DynamicPolynomialSet( RuntimeOrder publicOrder )
        : super( RuntimeOrder{ Feel::detail::nedelec2InternalOrder(
              Feel::detail::validatedNedelecSimplexOrder<
                  NedelecKind::NED2>( N, publicOrder ).value ) } ),
          M_publicOrder( publicOrder.value )
    {}

    /** @return public NED2 order. */
    [[nodiscard]] uint16_type publicOrder() const noexcept { return M_publicOrder; }
    /** @return ambient full-vector polynomial degree. */
    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        return Feel::detail::nedelec2InternalOrder( M_publicOrder );
    }

private:
    uint16_type M_publicOrder = 0;
};

template<uint16_type N,
         uint16_type O,
         typename T,
         uint16_type TheTAG >
class NedelecPolynomialSet<N,O,NedelecKind::NED2,T,TheTAG>
    :
        public Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex>
{
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> super;

public:
    static const uint16_type Om1 = (O==0)?0:O-1;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Vectorial, T, TheTAG, Simplex> Pk_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Simplex> Pkp1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, Om1, N, Vectorial, T, TheTAG, Simplex> Pkm1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Scalar, T, TheTAG, Simplex> Pk_s_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Scalar, T, TheTAG, Simplex> Pkp1_s_type;

    typedef PolynomialSet<typename super::basis_type,Vectorial> vectorial_polynomialset_type;
    typedef typename vectorial_polynomialset_type::polynomial_type vectorial_polynomial_type;
    typedef PolynomialSet<typename super::basis_type,Scalar> scalar_polynomialset_type;
    typedef typename scalar_polynomialset_type::polynomial_type scalar_polynomial_type;

    typedef NedelecPolynomialSet<N, O, NedelecKind::NED2,T,TheTAG> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::convex_type convex_type;
    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const uint16_type nComponents = super::nComponents;
    static inline const bool is_product = false;
    /**
     * @brief Construct the primal basis for the compile-time Nedelec family and order.
     */
    NedelecPolynomialSet()
        :
        super()
    {}
};

namespace fem
{

namespace detail
{

template<typename Basis,
         template<class, int, class> class PointSetType>
class NedelecDualFirstKind
    :
public DualBasis<Basis>
{
    typedef DualBasis<Basis> super;
public:

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder= super::nOrder;
    static inline const NedelecKind kind = NedelecKind::NED1;

    typedef typename super::primal_space_type primal_space_type;
    typedef typename primal_space_type::value_type value_type;
    typedef typename primal_space_type::points_type points_type;
    typedef typename primal_space_type::matrix_type matrix_type;
    typedef typename primal_space_type::template convex<nOrder+1>::type convex_type;
    typedef Reference<convex_type, nDim, nOrder+1, nDim, value_type> reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename convex_type::topological_face_type face_type;

    typedef typename primal_space_type::Pkp1_v_type Pkp1_v_type;
    typedef typename primal_space_type::vectorial_polynomialset_type vectorial_polynomialset_type;

    // point set type associated with the functionals
    typedef PointSetType<convex_type, nOrder, value_type> pointset_type;

    static inline const uint16_type nbPtsPerVertex = 0;
    static inline const uint16_type nbPtsPerEdge = reference_convex_type::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace2D = reference_convex_type::nbPtsPerFace;
    static const uint16_type nbPtsPerFace3D = nOrder*( nOrder-1 );
    static inline const uint16_type nbPtsPerFace = nDim == 2 ? nbPtsPerFace2D*nDim : nbPtsPerFace3D;
    static inline const uint16_type nbPtsPerVolume = nDim == 3 && nOrder > 2
        ? nOrder*( nOrder-1 )*( nOrder-2 )/2 : 0;
    static inline const uint16_type numPoints = ( reference_convex_type::numGeometricFaces*nbPtsPerFace+reference_convex_type::numEdges*nbPtsPerEdge+nbPtsPerVolume );

    /** Number of degrees of freedom per vertex */
    static inline const uint16_type nDofPerVertex = 0;

    /** Number of degrees of freedom per edge */
    static inline const uint16_type nDofPerEdge = nbPtsPerEdge;

    /** Number of degrees of freedom per face */
    static inline const uint16_type nDofPerFace = nbPtsPerFace;

    /** Number of degrees  of freedom per volume */
    static inline const uint16_type nDofPerVolume = nbPtsPerVolume;

    /** Total number of degrees of freedom (equal to refEle::nDof) */
    static inline const uint16_type nLocalDof = numPoints;

    static constexpr uint16_type nFacesInConvex =
        ( nDim == 1 ) ? reference_convex_type::numVertices :
        ( nDim == 2 ) ? reference_convex_type::numEdges :
                        reference_convex_type::numGeometricFaces;

    /**
     * @brief Construct the NED1 dual and its expression-sampling layout.
     * @param primal Primal space whose dimension must match the compile-time dual cardinality.
     */
    NedelecDualFirstKind( primal_space_type const& primal )
        :
        super( primal ),
        M_convex_ref(),
        M_pts( nDim, numPoints ),
        M_pts_per_face( nFacesInConvex ),
        M_samplingLayout(),
        M_fset( primal )
    {
        CHECK( primal.polynomialDimension() == nLocalDof )
            << "Invalid Nedelec first-kind primal/dual cardinality: primal dimension "
            << primal.polynomialDimension() << " differs from dual dimension " << nLocalDof
            << " in dimension " << nDim << " at internal order " << nOrder;
        VLOG(4) << "Nedelec finite element(dual): \n";
        VLOG(4) << " o- dim   = " << nDim << "\n";
        VLOG(4) << " o- order = " << nOrder << "\n";
        VLOG(4) << " o- kind = " << static_cast<int>(kind) << "\n";
        VLOG(4) << " o- numPoints      = " << numPoints << "\n";
        VLOG(4) << " o- nbPtsPerVertex = " << ( int )nbPtsPerVertex << "\n";
        VLOG(4) << " o- nbPtsPerEdge   = " << ( int )nbPtsPerEdge << "\n";
        VLOG(4) << " o- nbPtsPerFace   = " << ( int )nbPtsPerFace << "\n";
        VLOG(4) << " o- nbPtsPerVolume = " << ( int )nbPtsPerVolume << "\n";
        VLOG(4) << " o- nLocalDof      = " << nLocalDof << "\n";

        if ( nDim == 2 )
            CHECK( nLocalDof == nOrder*(nOrder+2) )
                << "Invalid finite element dimension. Expected : "
                << nOrder*(nOrder+2) << " Actual: " << nLocalDof;

        //LOG(INFO) << "Nedelec finite element(dual): \n";
        //LOG(INFO) << " o- dim   = " << nDim << "\n";
        //LOG(INFO) << " o- order = " << nOrder << "\n";
        //LOG(INFO) << " o- kind = " << static_cast<int>(kind) << "\n";
        //LOG(INFO) << " o- numPoints      = " << numPoints << "\n";
        //LOG(INFO) << " o- nbPtsPerVertex = " << ( int )nbPtsPerVertex << "\n";
        //LOG(INFO) << " o- nbPtsPerEdge   = " << ( int )nbPtsPerEdge << "\n";
        //LOG(INFO) << " o- nbPtsPerFace   = " << ( int )nbPtsPerFace << "\n";
        //LOG(INFO) << " o- nbPtsPerVolume = " << ( int )nbPtsPerVolume << "\n";
        //LOG(INFO) << " o- nLocalDof      = " << nLocalDof << "\n";

        size_type nbDofPerFace = ( nDim == 2 )
            ? face_type::numEdges*nbPtsPerEdge
            : face_type::numTopologicalFaces*nbPtsPerEdge + nbPtsPerFace;

        for(auto& m : M_pts_per_face)
            m.resize(nDim,nbDofPerFace);

        Pkp1_v_type ambientVector;
        auto fset = Feel::detail::buildNedelec1Functionals<
            primal_space_type, reference_convex_type, Pkp1_v_type,
            vectorial_polynomialset_type>(
                primal, M_convex_ref, static_cast<uint16_type>( nOrder-1 ),
                ambientVector, M_pts, M_pts_per_face, M_samplingLayout,
                [this]( uint16_type edge )
                { return M_convex_ref.makePoints( 1, edge ); } );

        VLOG(4) << "[Nedelec1stKind Dual] done 3, n fset = " << fset.size() << std::endl;
        M_fset.setFunctionalSet( fset );
        VLOG(4) << "[Nedelec1stKind DUAL matrix] mat = " << M_fset.rep() << "\n";
        VLOG(4) << "[Nedelec1stKind Dual] done 4\n";

    }

    /**
     * @return Complete element expression-sampling point table.
     */
    points_type const& points() const
    {
        return M_pts;
    }


    /**
     * @brief Evaluate the dual functionals on a primal polynomial basis.
     * @param pset Polynomial basis to evaluate.
     * @return Matrix of functional values on the basis.
     */
    matrix_type operator()( primal_space_type const& pset ) const
    {
        //VLOG(4) << "Nedelec1stKind matrix = " << M_fset( pset ) << std::endl;
        return M_fset( pset );
    }

    /**
     * @param f Element-local facet identifier.
     * @return Expression-sampling points in the facet closure.
     * @pre `f` indexes an existing facet.
     */
    points_type const& points( uint16_type f ) const
    {
        return M_pts_per_face[f];
    }
    /** @return expression-sampling blocks shared with the interpolant. */
    [[nodiscard]] Feel::detail::Nedelec1SamplingLayout const& samplingLayout() const noexcept
    {
        return M_samplingLayout;
    }
    /**
     * @param f Element-local facet identifier.
     * @param __i Column index in the facet sampling table.
     * @return Immutable view of the selected reference sample.
     * @pre Both indices are within the facet sampling table bounds.
     */
    ublas::matrix_column<points_type const> point( uint16_type f, uint32_type __i ) const
    {
        return ublas::column( M_pts_per_face[f], __i );
    }
    /**
     * @param f Element-local facet identifier.
     * @param __i Column index in the facet sampling table.
     * @return Mutable view of the selected reference sample.
     * @pre Both indices are within the facet sampling table bounds.
     */
    ublas::matrix_column<points_type> point( uint16_type f, uint32_type __i )
    {
        return ublas::column( M_pts_per_face[f], __i );
    }

private:
    /**
     * @brief Replace a facet expression-sampling table.
     * @param f Element-local facet identifier.
     * @param n Reference points stored as columns.
     * @pre `f` indexes an existing facet.
     */
    void setPoints( uint16_type f, points_type const& n )
    {
        M_pts_per_face[f].resize( n.size1(), n.size2(), false );
        M_pts_per_face[f] = n;
    }

private:
    reference_convex_type M_convex_ref;
    points_type M_pts;
    std::vector<points_type> M_pts_per_face;
    Feel::detail::Nedelec1SamplingLayout M_samplingLayout;
    FunctionalSet<primal_space_type> M_fset;


};

/**
 * @brief Runtime-order dual of the first-kind Nedelec simplex family.
 *
 * The mathematical functional construction is shared with
 * `NedelecDualFirstKind`; only storage extents and the equispaced edge-point
 * provider are selected at runtime.
 *
 * @tparam Basis runtime-order NED1 primal polynomial space
 * @tparam PointSetType point-set policy used for edge representative points
 */
template<typename Basis,
         template<class, int, class> class PointSetType>
class Nedelec1DynamicDual : public DualBasis<Basis>
{
    using super = DualBasis<Basis>;

public:
    static inline const uint16_type nDim = super::nDim;
    static inline const int nOrder = Dynamic;
    static inline const NedelecKind kind = NedelecKind::NED1;

    using primal_space_type = typename super::primal_space_type;
    using value_type = typename primal_space_type::value_type;
    using points_type = typename primal_space_type::points_type;
    using matrix_type = typename primal_space_type::matrix_type;
    using convex_type = Simplex<nDim, 1, nDim>;
    using reference_convex_type = Reference<convex_type, nDim, 1, nDim, value_type>;
    using face_type = typename convex_type::topological_face_type;
    using Pkp1_v_type = typename primal_space_type::Pkp1_v_type;
    using vectorial_polynomialset_type = typename primal_space_type::vectorial_polynomialset_type;
    using pointset_type = PointSetType<convex_type, Dynamic, value_type>;

    // Compile-time values describe the order-zero compatibility envelope.
    // Runtime consumers must use the accessors below for dynamic orders.
    static inline const uint16_type nbPtsPerVertex = 0;
    static inline const uint16_type nbPtsPerEdge = 1;
    static inline const uint16_type nbPtsPerFace = 0;
    static inline const uint16_type nbPtsPerVolume = 0;
    static inline const uint16_type numPoints = convex_type::numEdges;
    static inline const uint16_type nDofPerVertex = nbPtsPerVertex;
    static inline const uint16_type nDofPerEdge = nbPtsPerEdge;
    static inline const uint16_type nDofPerFace = nbPtsPerFace;
    static inline const uint16_type nDofPerVolume = nbPtsPerVolume;
    static inline const uint16_type nLocalDof = numPoints;

    /** @brief Construct the dual matching the public order of `primal`. */
    explicit Nedelec1DynamicDual( primal_space_type const& primal )
        : super( primal ),
          M_publicOrder( primal.publicOrder() ),
          M_internalOrder( Feel::detail::nedelec1InternalOrder( M_publicOrder ) ),
          M_convex_ref(),
          M_pts( nDim, localDof() ),
          M_pts_per_face( convex_type::numTopologicalFaces ),
          M_samplingLayout(),
          M_fset( primal )
    {
        const uint16_type facetClosureDof = nDim == 2
            ? dofPerEdge()
            : static_cast<uint16_type>( face_type::numTopologicalFaces*dofPerEdge() +
                                        dofPerFace() );
        for ( auto& points : M_pts_per_face )
            points.resize( nDim, facetClosureDof );

        // A simplex point set of degree internalOrder()+1 owns exactly
        // publicOrder()+1 interior points on each edge, matching the static
        // Reference<Simplex<...,nOrder+1>>::makePoints convention.
        pointset_type edgePointSet(
            RuntimeOrder{ static_cast<uint16_type>( M_internalOrder+1 ) } );
        Pkp1_v_type ambientVector( RuntimeOrder{ M_internalOrder } );
        auto functionals = Feel::detail::buildNedelec1Functionals<
            primal_space_type, reference_convex_type, Pkp1_v_type,
            vectorial_polynomialset_type>(
                primal, M_convex_ref, M_publicOrder, ambientVector,
                M_pts, M_pts_per_face, M_samplingLayout,
                [&edgePointSet, this]( uint16_type edge )
                {
                    auto points = edgePointSet.pointsBySubEntity( 1, edge );
                    CHECK_EQ( points.size2(), dofPerEdge() )
                        << "Invalid NED1 dynamic edge point count on edge " << edge;
                    return points;
                } );
        M_fset.setFunctionalSet( functionals );
    }

    /** @return public runtime polynomial order. */
    [[nodiscard]] uint16_type publicOrder() const noexcept { return M_publicOrder; }
    /** @return internal ambient polynomial order. */
    [[nodiscard]] uint16_type internalOrder() const noexcept { return M_internalOrder; }
    /** @return number of degrees of freedom on one vertex. */
    [[nodiscard]] constexpr uint16_type dofPerVertex() const noexcept { return 0; }
    /** @return number of degrees of freedom on one edge. */
    [[nodiscard]] uint16_type dofPerEdge() const noexcept
    {
        return Feel::detail::nedelec1EdgeDof( M_publicOrder );
    }
    /** @return number of face-interior degrees of freedom. */
    [[nodiscard]] uint16_type dofPerFace() const noexcept
    {
        return Feel::detail::nedelec1FaceDof( nDim, M_publicOrder );
    }
    /** @return number of cell-interior degrees of freedom. */
    [[nodiscard]] uint16_type dofPerVolume() const noexcept
    {
        return nDim == 3 ? Feel::detail::nedelec1VolumeDof( nDim, M_publicOrder ) : 0;
    }
    /** @return total number of local degrees of freedom. */
    [[nodiscard]] uint16_type localDof() const noexcept
    {
        return Feel::detail::nedelec1TotalDof( nDim, M_publicOrder );
    }

    /** @return representative element point table. */
    [[nodiscard]] points_type const& points() const { return M_pts; }
    /** @return the dual matrix evaluated on `pset`. */
    [[nodiscard]] matrix_type operator()( primal_space_type const& pset ) const
    {
        return M_fset( pset );
    }
    /** @return representative closure point table for facet `f`. */
    [[nodiscard]] points_type const& points( uint16_type f ) const
    {
        return M_pts_per_face[f];
    }
    /** @return expression-sampling blocks shared with the interpolant. */
    [[nodiscard]] Feel::detail::Nedelec1SamplingLayout const& samplingLayout() const noexcept
    {
        return M_samplingLayout;
    }
    /** @return immutable representative point `i` on facet `f`. */
    [[nodiscard]] ublas::matrix_column<points_type const>
    point( uint16_type f, uint32_type i ) const
    {
        return ublas::column( M_pts_per_face[f], i );
    }
    /** @return mutable representative point `i` on facet `f`. */
    [[nodiscard]] ublas::matrix_column<points_type>
    point( uint16_type f, uint32_type i )
    {
        return ublas::column( M_pts_per_face[f], i );
    }

private:
    uint16_type M_publicOrder = 0;
    uint16_type M_internalOrder = 1;
    reference_convex_type M_convex_ref;
    points_type M_pts;
    std::vector<points_type> M_pts_per_face;
    Feel::detail::Nedelec1SamplingLayout M_samplingLayout;
    FunctionalSet<primal_space_type> M_fset;
};

/** @brief Runtime-order canonical dual for second-kind Nedelec simplices. */
template<typename Basis, template<class, int, class> class PointSetType>
class Nedelec2DynamicDual : public DualBasis<Basis>
{
    using super = DualBasis<Basis>;

public:
    static inline const uint16_type nDim = super::nDim;
    static inline const int nOrder = Dynamic;
    static inline const NedelecKind kind = NedelecKind::NED2;
    using primal_space_type = typename super::primal_space_type;
    using value_type = typename primal_space_type::value_type;
    using points_type = typename primal_space_type::points_type;
    using matrix_type = typename primal_space_type::matrix_type;
    using convex_type = Simplex<nDim, 1, nDim>;
    using reference_convex_type = Reference<convex_type, nDim, 1, nDim, value_type>;
    using face_type = typename convex_type::topological_face_type;
    using pointset_type = PointSetType<convex_type, Dynamic, value_type>;

    // Order-zero compatibility constants; runtime consumers use accessors.
    static inline const uint16_type nbPtsPerVertex = 0;
    static inline const uint16_type nbPtsPerEdge = 2;
    static inline const uint16_type nbPtsPerFace = 0;
    static inline const uint16_type nbPtsPerVolume = 0;
    static inline const uint16_type numPoints = convex_type::numEdges*nbPtsPerEdge;
    static inline const uint16_type nDofPerVertex = 0;
    static inline const uint16_type nDofPerEdge = nbPtsPerEdge;
    static inline const uint16_type nDofPerFace = 0;
    static inline const uint16_type nDofPerVolume = 0;
    static inline const uint16_type nLocalDof = numPoints;

    /**
     * @brief Construct the NED2 dual matching the runtime primal order.
     * @param primal Runtime NED2 primal space supplying the public order.
     */
    explicit Nedelec2DynamicDual( primal_space_type const& primal )
        : super( primal ),
          M_publicOrder( primal.publicOrder() ),
          M_internalOrder( Feel::detail::nedelec2InternalOrder( M_publicOrder ) ),
          M_referenceConvex(),
          M_points( nDim, 0 ),
          M_facetPoints( convex_type::numTopologicalFaces ),
          M_functionals( primal )
    {
        pointset_type edgePointSet(
            RuntimeOrder{ static_cast<uint16_type>( M_internalOrder+2 ) } );
        auto functionals = Feel::detail::buildNedelec2Functionals(
            primal, M_referenceConvex, M_publicOrder,
            M_points, M_facetPoints, M_samplingLayout,
            [&edgePointSet, this]( uint16_type edge )
            {
                auto points = edgePointSet.pointsBySubEntity( 1, edge );
                CHECK_EQ( points.size2(), dofPerEdge() );
                return points;
            } );
        M_functionals.setFunctionalSet( functionals );
    }

    /**
     * @return Public NED2 family order.
     */
    [[nodiscard]] uint16_type publicOrder() const noexcept { return M_publicOrder; }
    /**
     * @return Ambient polynomial degree.
     */
    [[nodiscard]] uint16_type internalOrder() const noexcept { return M_internalOrder; }
    /**
     * @return Zero: Nedelec functionals are not attached to vertices.
     */
    [[nodiscard]] constexpr uint16_type dofPerVertex() const noexcept { return 0; }
    /**
     * @return Number of tangential evaluation functionals on one edge.
     */
    [[nodiscard]] uint16_type dofPerEdge() const noexcept
    {
        return Feel::detail::nedelec2EdgeDof( M_publicOrder );
    }
    /**
     * @return Number of triangular interior moments, on the cell in 2D or a face in 3D.
     */
    [[nodiscard]] uint16_type dofPerFace() const noexcept
    {
        return nDim == 2
            ? Feel::detail::nedelec2CellDof( nDim, M_publicOrder )
            : Feel::detail::nedelec2FaceDof( nDim, M_publicOrder );
    }
    /**
     * @return Number of tetrahedral cell-interior moments, or zero in 2D.
     */
    [[nodiscard]] uint16_type dofPerVolume() const noexcept
    {
        return nDim == 3
            ? Feel::detail::nedelec2CellDof( nDim, M_publicOrder ) : 0;
    }
    /**
     * @return Total number of element-local degrees of freedom.
     */
    [[nodiscard]] uint16_type localDof() const noexcept
    {
        return Feel::detail::nedelec2TotalDof( nDim, M_publicOrder );
    }
    /**
     * @return Complete element expression-sampling point table.
     */
    [[nodiscard]] points_type const& points() const { return M_points; }
    /**
     * @param space Primal polynomial basis to evaluate.
     * @return Matrix of dual functional values on that basis.
     */
    [[nodiscard]] matrix_type operator()( primal_space_type const& space ) const
    {
        return M_functionals( space );
    }
    /**
     * @param face Element-local facet identifier.
     * @return Expression-sampling points in the facet closure.
     * @pre `face` indexes an existing facet.
     */
    [[nodiscard]] points_type const& points( uint16_type face ) const
    {
        return M_facetPoints[face];
    }
    /**
     * @return Sampling offsets and cardinalities used by the interpolant.
     */
    [[nodiscard]] Feel::detail::Nedelec1SamplingLayout const& samplingLayout() const noexcept
    {
        return M_samplingLayout;
    }
    /**
     * @param face Element-local facet identifier.
     * @param pointId Column index in the facet sampling table.
     * @return Immutable view of the selected reference sample.
     * @pre Both indices are within the facet sampling table bounds.
     */
    [[nodiscard]] ublas::matrix_column<points_type const>
    point( uint16_type face, uint32_type pointId ) const
    {
        return ublas::column( M_facetPoints[face], pointId );
    }
    /**
     * @param face Element-local facet identifier.
     * @param pointId Column index in the facet sampling table.
     * @return Mutable view of the selected reference sample.
     * @pre Both indices are within the facet sampling table bounds.
     */
    [[nodiscard]] ublas::matrix_column<points_type>
    point( uint16_type face, uint32_type pointId )
    {
        return ublas::column( M_facetPoints[face], pointId );
    }

private:
    uint16_type M_publicOrder = 0;
    uint16_type M_internalOrder = 1;
    reference_convex_type M_referenceConvex;
    points_type M_points;
    std::vector<points_type> M_facetPoints;
    Feel::detail::Nedelec1SamplingLayout M_samplingLayout;
    FunctionalSet<primal_space_type> M_functionals;
};

/**
 * @brief Canonical arbitrary-order NED2 dual on triangles/tetrahedra.
 *
 * Static and runtime-order elements both delegate their edge/face/cell
 * moments to `buildNedelec2Functionals`.
 */
template<typename Basis, template<class, int, class> class PointSetType>
class Nedelec2Dual : public DualBasis<Basis>
{
    using super = DualBasis<Basis>;

public:
    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const NedelecKind kind = NedelecKind::NED2;
    static constexpr uint16_type publicOrder = static_cast<uint16_type>( nOrder-1 );

    using primal_space_type = typename super::primal_space_type;
    using value_type = typename primal_space_type::value_type;
    using points_type = typename primal_space_type::points_type;
    using matrix_type = typename primal_space_type::matrix_type;
    using convex_type = typename primal_space_type::template convex<nOrder+2>::type;
    using reference_convex_type = Reference<convex_type, nDim, nOrder+2, nDim, value_type>;
    using node_type = typename reference_convex_type::node_type;
    using face_type = typename convex_type::topological_face_type;
    using pointset_type = PointSetType<convex_type, nOrder, value_type>;

    static inline const uint16_type nbPtsPerVertex = 0;
    static inline const uint16_type nbPtsPerEdge = Feel::detail::nedelec2EdgeDof( publicOrder );
    static inline const uint16_type nbPtsPerFace = nDim == 2
        ? Feel::detail::nedelec2CellDof( nDim, publicOrder )
        : Feel::detail::nedelec2FaceDof( nDim, publicOrder );
    static inline const uint16_type nbPtsPerVolume = nDim == 3
        ? Feel::detail::nedelec2CellDof( nDim, publicOrder ) : 0;
    static inline const uint16_type numPoints = Feel::detail::nedelec2TotalDof( nDim, publicOrder );
    static inline const uint16_type nDofPerVertex = 0;
    static inline const uint16_type nDofPerEdge = nbPtsPerEdge;
    static inline const uint16_type nDofPerFace = nbPtsPerFace;
    static inline const uint16_type nDofPerVolume = nbPtsPerVolume;
    static inline const uint16_type nLocalDof = Feel::detail::nedelec2TotalDof( nDim, publicOrder );
    static constexpr uint16_type nFacesInConvex = nDim == 2
        ? reference_convex_type::numEdges
        : reference_convex_type::numGeometricFaces;

    /**
     * @brief Construct the NED2 dual for the compile-time primal order.
     * @param primal Primal space whose dimension must match the compile-time dual cardinality.
     */
    explicit Nedelec2Dual( primal_space_type const& primal )
        : super( primal ),
          M_referenceConvex(),
          M_points( nDim, 0 ),
          M_facetPoints( nFacesInConvex ),
          M_functionals( primal )
    {
        CHECK_EQ( primal.polynomialDimension(), nLocalDof );
        auto functionals = Feel::detail::buildNedelec2Functionals(
            primal, M_referenceConvex, publicOrder,
            M_points, M_facetPoints, M_samplingLayout,
            [this]( uint16_type edge )
            {
                auto points = M_referenceConvex.makePoints( 1, edge );
                CHECK_EQ( points.size2(), nDofPerEdge );
                return points;
            } );
        M_functionals.setFunctionalSet( functionals );
    }

    /** @return complete expression-sampling point table. */
    [[nodiscard]] points_type const& points() const { return M_points; }
    /** @return the dual matrix evaluated on `space`. */
    [[nodiscard]] matrix_type operator()( primal_space_type const& space ) const
    {
        return M_functionals( space );
    }
    /** @return facet-closure expression-sampling points. */
    [[nodiscard]] points_type const& points( uint16_type face ) const
    {
        return M_facetPoints[face];
    }
    /** @return expression-sampling layout shared with interpolation. */
    [[nodiscard]] Feel::detail::Nedelec1SamplingLayout const& samplingLayout() const noexcept
    {
        return M_samplingLayout;
    }
    /**
     * @param face Element-local facet identifier.
     * @param pointId Column index in the facet sampling table.
     * @return Immutable view of the selected reference sample.
     * @pre Both indices are within the facet sampling table bounds.
     */
    [[nodiscard]] ublas::matrix_column<points_type const>
    point( uint16_type face, uint32_type pointId ) const
    {
        return ublas::column( M_facetPoints[face], pointId );
    }
    /**
     * @param face Element-local facet identifier.
     * @param pointId Column index in the facet sampling table.
     * @return Mutable view of the selected reference sample.
     * @pre Both indices are within the facet sampling table bounds.
     */
    [[nodiscard]] ublas::matrix_column<points_type>
    point( uint16_type face, uint32_type pointId )
    {
        return ublas::column( M_facetPoints[face], pointId );
    }

private:
    reference_convex_type M_referenceConvex;
    points_type M_points;
    std::vector<points_type> M_facetPoints;
    Feel::detail::Nedelec1SamplingLayout M_samplingLayout;
    FunctionalSet<primal_space_type> M_functionals;
};

}// detail

template<uint16_type N,
         uint16_type O,
         NedelecKind Kind = NedelecKind::NED1,
         typename T = double,
         uint16_type TheTAG = 0 >
struct NedelecBase
{
    using type = if_t<( Kind == NedelecKind::NED2 ),
                      FiniteElement<NedelecPolynomialSet<N, O, NedelecKind::NED2, T, TheTAG>,
                                    fem::detail::Nedelec2Dual, PointSetEquiSpaced>,
                      FiniteElement<NedelecPolynomialSet<N, O, NedelecKind::NED1, T>,
                                    fem::detail::NedelecDualFirstKind, PointSetEquiSpaced>>;
};

/**
 * \class Nedelec
 * \brief Nedelec Finite Element
 *
 * \f$ H(curl)\f$  conforming element
 * -\tparam  N topological dimension
 * -\tparam  O polynomial order
 * -\tparam  Kind Nedelec finite element : 1st or 2nd kind
 * -\tparam  T numerical type
 * -\tparam  TheTAG identifier for the finite element
 *
 * @author Christophe Prud'homme
 */
template<uint16_type N,
         uint16_type O,
         NedelecKind Kind = NedelecKind::NED1,
         typename T = double,
         uint16_type TheTAG = 0 >
class Nedelec
    :
    public HCurlPolynomialSet,
    public NedelecBase<N,O,Kind,T,TheTAG>::type,
    public std::enable_shared_from_this<Nedelec<N,O,Kind,T,TheTAG> >
{

public:
    static_assert( Kind != NedelecKind::NED1 || N == 2 || O <= 3,
                   "3D Nedelec first kind is currently validated through public order three." );
    static_assert( Kind != NedelecKind::NED2 || N == 2 || O <= 3,
                   "3D Nedelec second kind is currently validated through public order three." );

    typedef typename NedelecBase<N,O,Kind,T,TheTAG>::type super;

    BOOST_STATIC_ASSERT( N > 1 );

    /** @name Typedefs
     */
    //@{

    static inline const uint16_type nDim = N;
    //static inline const bool isTransformationEquivalent = false;
    static inline const bool isTransformationEquivalent = true;
    static inline const bool isContinuous = true;
    static inline const NedelecKind kind = Kind;
    typedef typename super::value_type value_type;
    typedef typename super::primal_space_type primal_space_type;
    typedef typename super::dual_space_type dual_space_type;
    typedef Continuous continuity_type;
    static const uint16_type TAG = TheTAG;

    //static const polynomial_transformation_type transformation = POLYNOMIAL_CONTEXT_NEEDS_1ST_PIOLA_TRANSFORMATION;

    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename super::polyset_type polyset_type;
    static inline const bool is_vectorial = polyset_type::is_vectorial;
    static inline const bool is_scalar = polyset_type::is_scalar;
    static inline const uint16_type nComponents = polyset_type::nComponents;
    static inline const bool is_product = false;


    typedef typename dual_space_type::convex_type convex_type;
    typedef typename dual_space_type::pointset_type pointset_type;
    typedef typename dual_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename reference_convex_type::points_type points_type;
    typedef typename convex_type::topological_face_type face_type;

    static inline const uint16_type nOrder =  dual_space_type::nOrder;
    static inline const uint16_type nbPtsPerVertex = reference_convex_type::nbPtsPerVertex;
    static inline const uint16_type nbPtsPerEdge = dual_space_type::nbPtsPerEdge;
    static inline const uint16_type nbPtsPerFace = dual_space_type::nbPtsPerFace;
    static inline const uint16_type nbPtsPerVolume = dual_space_type::nbPtsPerVolume;
    static inline const uint16_type numPoints = dual_space_type::numPoints;

    static inline const uint16_type nLocalDof = dual_space_type::nLocalDof;
    static inline const uint16_type nDofPerVertex = dual_space_type::nDofPerVertex;
    static inline const uint16_type nDofPerEdge = dual_space_type::nDofPerEdge;
    static inline const uint16_type nDofPerFace = dual_space_type::nDofPerFace;
    static inline const uint16_type nDofPerVolume = dual_space_type::nDofPerVolume;
    static inline const uint16_type nLocalFaceDof = ( face_type::numVertices * nDofPerVertex +
                                               face_type::numEdges * nDofPerEdge +
                                               face_type::numFaces * nDofPerFace );

    template<int subN>
    struct SubSpace
    {
        typedef Nedelec<N-1, O, kind, T, TheTAG> type;
    };

    struct SSpace
    {
        typedef Nedelec<N, O, kind, T, TheTAG> type;

    };

    template<uint16_type NewDim>
    struct ChangeDim
    {
        typedef Nedelec<NewDim, O, kind, T,  TheTAG> type;
    };

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * @brief Construct the nodal Nedelec element from its compile-time primal and dual spaces.
     */
    Nedelec()
        :
        super( dual_space_type( primal_space_type() ) ),
        M_refconvex()
    {
#if 1
        VLOG(4) << "[N] nPtsPerEdge = " << nbPtsPerEdge << "\n";
        VLOG(4) << "[N] nPtsPerFace = " << nbPtsPerFace << "\n";
        VLOG(4) << "[N] numPoints = " << numPoints << "\n";

        VLOG(4) << "[N] nDof = " << super::nDof << "\n";

        VLOG(4) << "[N] coeff : " << this->coeff() << "\n";
        VLOG(4) << "[N] pts : " << this->points() << "\n";
        VLOG(4) << "[N] eval at pts : " << this->evaluate( this->points() ) << "\n";
        VLOG(4) << "[N] is_product : " << is_product << "\n";
#endif
    }

    /**
     * @brief Copy a Nedelec element and initialize its reference simplex.
     * @param cr Element to copy.
     */
    Nedelec( Nedelec const & cr )
        :
        super( cr ),
        M_refconvex()
    {}

    /**
     * @brief Destroy the Nedelec element.
     */
    ~Nedelec() override
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @return Reference simplex associated with the Nedelec basis.
     */
    reference_convex_type const& referenceConvex() const
    {
        return M_refconvex;
    }

    /**
     * @return Finite-element family identifier, `nedelec`.
     */
    std::string familyName() const override
    {
        return "nedelec";
    }

    //! \return the component of a local dof
    uint16_type component( uint16_type localDofId ) const override
        {
            return 0;
        }

    //! \return a parent local dof id for each component (for example, the first component)
    uint16_type dofParent( uint16_type localDofId ) const override
        {
            return localDofId;
        }

    /**
     * @param localDofId Local functional index; the answer is independent of the index.
     * @return False: the vector functionals expose no scalar representative point.
     */
    bool dofHasRepresentativePoint( uint16_type localDofId ) const override
        {
            (void)localDofId;
            return false;
        }

    /**
     * @brief Classify a local functional by its topological attachment.
     * @param localDofId Element-local degree-of-freedom index.
     * @return InteriorMoment for cell DoFs, TangentialMoment for boundary DoFs, or Other for an invalid index.
     */
    uint16_type dofFunctionalKind( uint16_type localDofId ) const override
        {
            auto const attachment = this->dofAttachment( localDofId );
            if ( !attachment.isValid() )
                return static_cast<uint16_type>( DofFunctionalKind::Other );
            if ( static_cast<uint16_type>( attachment.entityDim ) == nDim )
                return static_cast<uint16_type>( DofFunctionalKind::InteriorMoment );
            return static_cast<uint16_type>( DofFunctionalKind::TangentialMoment );
        }

    /**
     * @brief Locate the mesh entity carrying a local degree of freedom.
     * @param localDofId Element-local degree-of-freedom index.
     * @return Entity dimension, identifier, and ordinal; invalid for an out-of-range index.
     */
    typename super::DofAttachment dofAttachment( uint16_type localDofId ) const override
        {
            return simplexDofAttachment<typename super::DofAttachment>(
                this->dofParent( localDofId ), nDim,
                reference_convex_type::numVertices,
                reference_convex_type::numEdges,
                reference_convex_type::numFaces,
                nDofPerVertex, nDofPerEdge, nDofPerFace, nDofPerVolume,
                this->dofType( localDofId ) );
        }

    /**
     * @brief Return the canonical-to-local orientation action for a degree of freedom.
     * @tparam ElementType Mesh element exposing entity permutations and vertex identifiers.
     * @param element Mesh element carrying the local functional.
     * @param localDofId Element-local degree-of-freedom index.
     * @return A triangular-face moment matrix row, or the generic entity orientation action.
     */
    template<typename ElementType>
    [[nodiscard]] DofTransform dofTransform( ElementType const& element, uint16_type localDofId ) const
    {
        auto const attachment = this->dofAttachment( localDofId );
        if constexpr ( nDim == 3 && Kind == NedelecKind::NED1 && nDofPerFace > 0 )
        {
            if ( attachment.isValid() && attachment.entityDim == 2 )
            {
                auto matrix = Feel::detail::nedelecTriangularFaceMomentTransform(
                    Feel::detail::nedelecLocalToCanonicalFaceVertices(
                        element, attachment.entityId ),
                    static_cast<uint16_type>( nOrder - 2 ) );
                return DofTransform::matrixOperator( std::move( matrix ), attachment.ordinal );
            }
        }
        else if constexpr ( nDim == 3 && Kind == NedelecKind::NED2 && nDofPerFace > 0 )
        {
            if ( attachment.isValid() && attachment.entityDim == 2 )
            {
                auto matrix = Feel::detail::nedelecTriangularFaceRTMomentTransform(
                    Feel::detail::nedelecLocalToCanonicalFaceVertices(
                        element, attachment.entityId ), O );
                return DofTransform::matrixOperator( std::move( matrix ), attachment.ordinal );
            }
        }
        return finiteElementEntityOrientationTransform( *this, element, localDofId );
    }

    //! \return the type of a local dof
    uint16_type dofType( uint16_type localDofId ) const override
        {
            return 1;
        }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{
    typedef Eigen::VectorXd local_interpolant_type;
    /**
     * @brief Allocate a zeroed element interpolation vector.
     * @param n Number of local DoF blocks to allocate.
     * @return Vector with `n` times the local DoF count entries.
     */
    local_interpolant_type
    localInterpolant( int n = 1 ) const
        {
            return local_interpolant_type::Zero( n*nLocalDof );
        }

    typedef Eigen::MatrixXd local_interpolants_type;
    /**
     * @brief Allocate a zeroed collection of local interpolation vectors.
     * @param p Number of columns.
     * @param n Number of local DoF blocks per column.
     * @return Matrix with `n` times the local DoF count rows and `p` columns.
     */
    local_interpolants_type
    localInterpolants( int p, int n = 1 ) const
        {
            return local_interpolants_type::Zero( n*nLocalDof, p );
        }

    /**
     * @brief Apply the element edge, face, and cell functionals to an expression.
     * @tparam ExprType Vector expression with a compatible geometric context.
     * @param expr Expression sampled at the dual sampling points.
     * @param[out] Ihloc Preallocated local DoF vector, cleared before interpolation.
     * @pre `Ihloc` has at least `nLocalDof` entries.
     */
    template<typename ExprType>
    void
    interpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            Ihloc.setZero();
            auto g = expr.geom();

            Feel::detail::interpolateHCurlEdgePointDofs(
                expr, Ihloc, convex_type::numEdges, nDofPerEdge, t,
                [&expr]( uint16_type edge, uint16_type point, auto& tangent )
                {
                    expr.geom()->piolaCovariantEdgeTangent( edge, point, tangent );
                } );
            auto const& sampling = this->dual().samplingLayout();
            if constexpr ( nDim == 3 && nDofPerFace > 0 )
            {
                if constexpr ( Kind == NedelecKind::NED1 )
                {
                    for ( auto const& block : sampling.faces )
                        Feel::detail::interpolateHCurlTriangularFaceMoments(
                            expr, Ihloc, M_refconvex, block.face,
                            static_cast<int>( nOrder )-2, block.firstDof,
                            block.elementPointOffset, block.momentCount );
                }
                else
                {
                    RaviartThomasDynamicPolynomialSet<2, value_type> moments(
                        RuntimeOrder{ static_cast<uint16_type>( O-1 ) } );
                    for ( auto const& block : sampling.faces )
                        Feel::detail::interpolateHCurlTriangularFaceRTMoments(
                            expr, Ihloc, M_refconvex, moments, block.face,
                            block.firstDof, block.elementPointOffset );
                }
            }
            constexpr uint16_type nInternalDof = nDim == 2 ? nDofPerFace : nDofPerVolume;
            if constexpr ( nInternalDof > 0 )
            {
                if constexpr ( Kind == NedelecKind::NED1 )
                {
                    constexpr uint16_type momentOrder = nDim == 2 ? nOrder-2 : nOrder-3;
                    Feel::detail::OrthonormalPolynomialSet<
                        nDim,momentOrder,nDim,Vectorial,value_type,0,Simplex> moments;
                    Feel::detail::interpolateHCurlCellMoments<nDim>(
                        expr, Ihloc, moments, nOrder,
                        sampling.cellFirstDof, sampling.cellPointOffset,
                        nInternalDof, nComponents );
                }
                else
                {
                    constexpr uint16_type momentOrder = nDim == 2 ? O-1 : O-2;
                    RaviartThomasDynamicPolynomialSet<nDim, value_type> moments(
                        RuntimeOrder{ momentOrder } );
                    Feel::detail::interpolateHCurlCellMoments<nDim>(
                        expr, Ihloc, moments, static_cast<uint16_type>( nOrder+1 ),
                        sampling.cellFirstDof, sampling.cellPointOffset,
                        nInternalDof, nComponents );
                }
            }
        }
    /**
     * @brief Allocate a zeroed facet-closure interpolation vector.
     * @return Vector containing all edge and face-interior DoFs in one facet closure.
     */
    local_interpolant_type
    faceLocalInterpolant() const
        {
            return local_interpolant_type::Zero( nLocalFaceDof, 1 );
        }
    /**
     * @brief Apply edge and face-interior functionals on the current facet closure.
     * @tparam ExprType Vector expression with a facet geometric context.
     * @param expr Expression sampled at the current facet sampling points.
     * @param[out] Ihloc Preallocated facet DoF vector, cleared before interpolation.
     * @pre `Ihloc` has at least `nLocalFaceDof` entries.
     */
    template<typename ExprType>
    void
    faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            Ihloc.setZero();
            auto g = expr.geom();

            Feel::detail::interpolateHCurlEdgePointDofs(
                expr, Ihloc, face_type::numEdges, nDofPerEdge, t,
                [&expr, &g]( uint16_type edge, uint16_type point, auto& tangent )
                {
                    const int elementEdge = nDim <= 2
                                                ? g->element().fToE( g->faceId(), g->faceId() )
                                                : g->element().fToE( g->faceId(), edge );
                    expr.geom()->piolaCovariantEdgeTangent( elementEdge, point, tangent );
                } );

            if constexpr ( nDim == 3 && nDofPerFace > 0 )
            {
                const uint16_type face = static_cast<uint16_type>( g->faceId() );
                auto const& blocks = this->dual().samplingLayout().faces;
                CHECK_LT( face, blocks.size() );
                auto const& block = blocks[face];
                if constexpr ( Kind == NedelecKind::NED1 )
                    Feel::detail::interpolateHCurlTriangularFaceMoments(
                        expr, Ihloc, M_refconvex, face,
                        static_cast<int>( nOrder )-2,
                        static_cast<uint16_type>( face_type::numEdges*nDofPerEdge ),
                        block.facetPointOffset, block.momentCount );
                else
                {
                    RaviartThomasDynamicPolynomialSet<2, value_type> moments(
                        RuntimeOrder{ static_cast<uint16_type>( O-1 ) } );
                    Feel::detail::interpolateHCurlTriangularFaceRTMoments(
                        expr, Ihloc, M_refconvex, moments, face,
                        static_cast<uint16_type>( face_type::numEdges*nDofPerEdge ),
                        block.facetPointOffset );
                }
            }

        }

    /**
     * @brief Apply all element functionals to each domain basis row.
     * @tparam ExprType Basis expression with an element geometric context.
     * @param expr Basis expression sampled at the dual sampling points.
     * @param[out] Ihloc Preallocated matrix with domain basis rows and `nLocalDof` columns.
     */
    template<typename ExprType>
        requires ( std::remove_cvref_t<ExprType>::gmc_type::subEntityCoDim == 0 )
    void
    interpolateBasisFunction( ExprType&& expr, local_interpolants_type& Ihloc ) const
    {
        Feel::detail::interpolateBasisFunctionRows( *this, expr, Ihloc );
    }

    /**
     * @brief Forward basis interpolation through the compatibility mapping overload.
     * @tparam ExprType Basis expression with a supported geometric context.
     * @param expr Basis expression passed to the matching context overload.
     * @param[out] Ihloc Preallocated interpolation matrix.
     * @param mapExprPointToDofPoint Compatibility argument; currently ignored.
     */
    template<typename ExprType>
    void
    interpolateBasisFunction( ExprType&& expr, local_interpolants_type & Ihloc, std::vector<uint16_type> const& mapExprPointToDofPoint ) const
        {
            //CHECK( false ) << "TODO";
            return this->interpolateBasisFunction( std::forward<ExprType>( expr ), Ihloc );
        }

    /**
     * @brief Interpolate oriented edge traces of basis rows on the current facet.
     * @tparam ExprType Basis expression with a facet geometric context.
     * @param expr Basis expression sampled at facet edge points.
     * @param[out] Ihloc Preallocated matrix indexed by domain basis row and element DoF.
     * @note This overload fills edge functionals only.
     */
    template<typename ExprType>
        requires ( std::remove_cvref_t<ExprType>::gmc_type::subEntityCoDim == 1 )
    void
    interpolateBasisFunction( ExprType&& expr, local_interpolants_type& Ihloc ) const
    {
        using shape = typename std::decay_t<ExprType>::shape;
        using expr_basis_t = typename std::decay_t<ExprType>::expr_type::test_basis;
        Ihloc.setZero();
        int nPoints = expr.nPoints();

        for( int e = 0; e < face_type::numEdges; ++e )
        {
            int edgeid_in_element;
            if( nDim <= 2 )
                edgeid_in_element = expr.geom()->element().fToE( expr.geom()->faceId(), expr.geom()->faceId());
            else
                edgeid_in_element = expr.geom()->element().fToE( expr.geom()->faceId(), e);

            int sign = 1;
#if 1
            if ( expr.geom()->element().edgePermutation( expr.geom()->faceId() ).value() == 2 )///*face_type::permutation_type*/ /*edge_permutation_type*/::REVERSE_PERMUTATION )
                sign = -1;
#endif


            for ( int l = 0; l < nDofPerEdge; ++l )
            {
                int q = e*nDofPerEdge+l;
                int q2 = edgeid_in_element*nDofPerEdge+l;
                expr.geom()->piolaCovariantEdgeTangent( edgeid_in_element, q, t );
                CHECK( q < nPoints) << "wrong point eval" << q << " vs " << nPoints;

                for( int i = 0; i < expr_basis_t::nLocalDof; ++i )
                {
                    int ncomp= ( expr_basis_t::is_product?expr_basis_t::nComponents1:1 );

                    for ( uint16_type c = 0; c < ncomp; ++c )
                    {
                        uint16_type I = expr_basis_t::nLocalDof*c + i;
                        for( int c1 = 0; c1 < shape::M; ++c1 )
                        {
                            Ihloc(I,q2) += sign*expr.evaliq( I, c1, 0, q )*t(c1);
                        }
                    }
                }
            }
        }

    }


    //@}



protected:
    reference_convex_type M_refconvex;
    mutable ublas::vector<value_type> t { nDim };
private:

};

namespace detail
{
template<uint16_type N, NedelecKind Kind, typename T, uint16_type Tag>
struct NedelecDynamicData;

template<uint16_type N, typename T, uint16_type Tag>
struct NedelecDynamicData<N, NedelecKind::NED1, T, Tag>
{
    using primal_type = Nedelec1DynamicPolynomialSet<N, T, Tag>;
    template<typename Basis, template<class, int, class> class PointSetType>
    using dual_type = Nedelec1DynamicDual<Basis, PointSetType>;
};

template<uint16_type N, typename T, uint16_type Tag>
struct NedelecDynamicData<N, NedelecKind::NED2, T, Tag>
{
    using primal_type = Nedelec2DynamicPolynomialSet<N, T, Tag>;
    template<typename Basis, template<class, int, class> class PointSetType>
    using dual_type = Nedelec2DynamicDual<Basis, PointSetType>;
};
} // namespace detail

/**
 * @brief True runtime-order Nedelec finite element on a simplex.
 *
 * The primal and dual delegate to the same order-parameterized builders as
 * the static family. This class owns only runtime FE metadata, layout, and
 * orientation dispatch.
 *
 * @tparam N reference-simplex dimension
 * @tparam Kind first- or second-kind family
 * @tparam T scalar type
 * @tparam TheTAG finite-element family tag
 */
template<uint16_type N, NedelecKind Kind = NedelecKind::NED1,
         typename T = double, uint16_type TheTAG = 0>
class NedelecDynamicSimplex
    : public HCurlPolynomialSet,
      public FiniteElement<
          typename detail::NedelecDynamicData<N, Kind, T, TheTAG>::primal_type,
          detail::NedelecDynamicData<N, Kind, T, TheTAG>::template dual_type,
                           PointSetEquiSpaced>,
      public std::enable_shared_from_this<NedelecDynamicSimplex<N, Kind, T, TheTAG>>
{
    using data_type = detail::NedelecDynamicData<N, Kind, T, TheTAG>;
    using super = FiniteElement<typename data_type::primal_type,
                                data_type::template dual_type,
                                PointSetEquiSpaced>;

public:
    static_assert( N == 2 || N == 3,
                   "Runtime Nedelec elements are implemented on triangles and tetrahedra." );

    static inline const uint16_type nDim = N;
    static inline const bool isTransformationEquivalent = true;
    static inline const bool isContinuous = true;
    static inline const NedelecKind kind = Kind;
    static inline const bool is_order_static = false;
    static inline const bool is_order_dynamic = true;
    static constexpr int nOrder_v = Dynamic;
    static constexpr int nOrder = Dynamic;
    static constexpr uint16_type TAG = TheTAG;

    using continuity_type = Continuous;
    using value_type = typename super::value_type;
    using primal_space_type = typename super::primal_space_type;
    using dual_space_type = typename super::dual_space_type;
    using polyset_type = typename super::polyset_type;
    using convex_type = typename dual_space_type::convex_type;
    using pointset_type = typename dual_space_type::pointset_type;
    using reference_convex_type = typename dual_space_type::reference_convex_type;
    using node_type = typename reference_convex_type::node_type;
    using points_type = typename reference_convex_type::points_type;
    using face_type = typename convex_type::topological_face_type;
    using Pkp1_v_type = typename primal_space_type::Pkp1_v_type;
    using vectorial_polynomialset_type = typename primal_space_type::vectorial_polynomialset_type;

    static inline const bool is_vectorial = polyset_type::is_vectorial;
    static inline const bool is_scalar = polyset_type::is_scalar;
    static inline const uint16_type nComponents = polyset_type::nComponents;
    static inline const bool is_product = false;

    // Order-zero compatibility constants. Runtime consumers use the virtual
    // layout API and the runtime accessors below.
    static inline const uint16_type nbPtsPerVertex = dual_space_type::nbPtsPerVertex;
    static inline const uint16_type nbPtsPerEdge = dual_space_type::nbPtsPerEdge;
    static inline const uint16_type nbPtsPerFace = dual_space_type::nbPtsPerFace;
    static inline const uint16_type nbPtsPerVolume = dual_space_type::nbPtsPerVolume;
    static inline const uint16_type numPoints = dual_space_type::numPoints;
    static inline const uint16_type nLocalDof = dual_space_type::nLocalDof;
    static inline const uint16_type nDofPerVertex = dual_space_type::nDofPerVertex;
    static inline const uint16_type nDofPerEdge = dual_space_type::nDofPerEdge;
    static inline const uint16_type nDofPerFace = dual_space_type::nDofPerFace;
    static inline const uint16_type nDofPerVolume = dual_space_type::nDofPerVolume;
    static inline const uint16_type nLocalFaceDof = face_type::numEdges*nDofPerEdge;

    /**
     * @brief Construct the lowest-order runtime Nedelec element.
     */
    NedelecDynamicSimplex() : NedelecDynamicSimplex( RuntimeOrder{ 0 } ) {}

    /**
     * @brief Construct a runtime Nedelec finite element of the requested order.
     * @param order Public Nedelec family order.
     * @throws std::length_error If the ambient degree or local dimension exceeds local index storage.
     */
    explicit NedelecDynamicSimplex( RuntimeOrder order )
        : super( dual_space_type( primal_space_type( order ) ) ),
          M_publicOrder( order.value ),
          M_refconvex()
    {}

    /**
     * @brief Copy the runtime element and initialize its reference simplex.
     * @param other Element whose order and finite-element basis are copied.
     */
    NedelecDynamicSimplex( NedelecDynamicSimplex const& other )
        : super( other ), M_publicOrder( other.M_publicOrder ), M_refconvex()
    {}

    /**
     * @brief Destroy the runtime Nedelec element.
     */
    ~NedelecDynamicSimplex() override = default;

    template<int subN>
    struct SubSpace
    {
        using type = NedelecDynamicSimplex<N-1, Kind, T, TheTAG>;
    };

    struct SSpace
    {
        using type = NedelecDynamicSimplex<N, Kind, T, TheTAG>;
    };

    template<uint16_type NewDim>
    struct ChangeDim
    {
        using type = NedelecDynamicSimplex<NewDim, Kind, T, TheTAG>;
    };

    /** @return reference simplex used for topology and orientation. */
    [[nodiscard]] reference_convex_type const& referenceConvex() const { return M_refconvex; }
    /** @return finite-element family name. */
    [[nodiscard]] std::string familyName() const override { return "nedelec"; }
    /** @return public Nedelec family order. */
    [[nodiscard]] uint16_type order() const noexcept { return M_publicOrder; }
    /** @return ambient polynomial degree. */
    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        if constexpr ( Kind == NedelecKind::NED1 )
            return Feel::detail::nedelec1InternalOrder( M_publicOrder );
        else
            return Feel::detail::nedelec2InternalOrder( M_publicOrder );
    }
    /** @return local cardinality. */
    [[nodiscard]] uint16_type localDof() const noexcept
    {
        if constexpr ( Kind == NedelecKind::NED1 )
            return Feel::detail::nedelec1TotalDof( nDim, M_publicOrder );
        else
            return Feel::detail::nedelec2TotalDof( nDim, M_publicOrder );
    }
    /** @return number of vertex DoFs. */
    [[nodiscard]] constexpr uint16_type dofPerVertex() const noexcept { return 0; }
    /** @return number of edge DoFs. */
    [[nodiscard]] uint16_type dofPerEdge() const noexcept
    {
        if constexpr ( Kind == NedelecKind::NED1 )
            return Feel::detail::nedelec1EdgeDof( M_publicOrder );
        else
            return Feel::detail::nedelec2EdgeDof( M_publicOrder );
    }
    /** @return number of face-interior DoFs. */
    [[nodiscard]] uint16_type dofPerFace() const noexcept
    {
        if constexpr ( Kind == NedelecKind::NED1 )
            return Feel::detail::nedelec1FaceDof( nDim, M_publicOrder );
        else if constexpr ( nDim == 2 )
            return Feel::detail::nedelec2CellDof( nDim, M_publicOrder );
        else
            return Feel::detail::nedelec2FaceDof( nDim, M_publicOrder );
    }
    /** @return number of volume-interior DoFs. */
    [[nodiscard]] uint16_type dofPerVolume() const noexcept
    {
        if constexpr ( nDim != 3 )
            return 0;
        else if constexpr ( Kind == NedelecKind::NED1 )
            return Feel::detail::nedelec1VolumeDof( nDim, M_publicOrder );
        else
            return Feel::detail::nedelec2CellDof( nDim, M_publicOrder );
    }
    /** @return number of DoFs in one facet closure. */
    [[nodiscard]] uint16_type localFacetDof() const noexcept
    {
        if constexpr ( nDim == 2 )
            return dofPerEdge();
        else
            return static_cast<uint16_type>( face_type::numEdges*dofPerEdge()+dofPerFace() );
    }

    /**
     * @return Total local dimension: Nedelec uses one vector basis rather than a component product.
     */
    uint16_type localDofPerComponent() const override { return localDof(); }
    /**
     * @param perComponent Compatibility flag; the vector-basis count is independent of it.
     * @return Total number of local degrees of freedom.
     */
    uint16_type localDofCount( bool perComponent = false ) const override
    {
        (void)perComponent;
        return localDof();
    }
    /**
     * @return Component index zero for every vector functional.
     */
    uint16_type component( uint16_type /*localDofId*/ ) const override { return 0; }
    /**
     * @param localDofId Element-local degree-of-freedom index.
     * @return The unchanged index, since the vector basis has no component-product parent.
     */
    uint16_type dofParent( uint16_type localDofId ) const override { return localDofId; }
    /**
     * @return False: the vector functionals expose no scalar representative point.
     */
    bool dofHasRepresentativePoint( uint16_type /*localDofId*/ ) const override { return false; }
    /**
     * @return Orientation-sensitive DoF type code, always 1.
     */
    uint16_type dofType( uint16_type /*localDofId*/ ) const override { return 1; }

    /**
     * @brief Classify a local functional by its topological attachment.
     * @param localDofId Element-local degree-of-freedom index.
     * @return InteriorMoment for cell DoFs, TangentialMoment for boundary DoFs, or Other for an invalid index.
     */
    uint16_type dofFunctionalKind( uint16_type localDofId ) const override
    {
        auto const attachment = this->dofAttachment( localDofId );
        if ( !attachment.isValid() )
            return static_cast<uint16_type>( DofFunctionalKind::Other );
        if ( static_cast<uint16_type>( attachment.entityDim ) == nDim )
            return static_cast<uint16_type>( DofFunctionalKind::InteriorMoment );
        return static_cast<uint16_type>( DofFunctionalKind::TangentialMoment );
    }

    /**
     * @brief Locate the mesh entity carrying a local degree of freedom.
     * @param localDofId Element-local degree-of-freedom index.
     * @return Entity dimension, identifier, and ordinal; invalid for an out-of-range index.
     */
    typename super::DofAttachment dofAttachment( uint16_type localDofId ) const override
    {
        return simplexDofAttachment<typename super::DofAttachment>(
            this->dofParent( localDofId ), nDim,
            reference_convex_type::numVertices,
            reference_convex_type::numEdges,
            reference_convex_type::numFaces,
            dofPerVertex(), dofPerEdge(), dofPerFace(), dofPerVolume(),
            this->dofType( localDofId ) );
    }

    /**
     * @brief Return the canonical-to-local orientation action for a degree of freedom.
     * @tparam ElementType Mesh element exposing entity permutations and vertex identifiers.
     * @param element Mesh element carrying the local functional.
     * @param localDofId Element-local degree-of-freedom index.
     * @return A triangular-face moment matrix row, or the generic entity orientation action.
     */
    template<typename ElementType>
    [[nodiscard]] DofTransform dofTransform( ElementType const& element,
                                             uint16_type localDofId ) const
    {
        auto const attachment = this->dofAttachment( localDofId );
        if constexpr ( nDim == 3 )
        {
            if ( attachment.isValid() && attachment.entityDim == 2 && dofPerFace() > 0 )
            {
                auto matrix = [&]()
                {
                    if constexpr ( Kind == NedelecKind::NED1 )
                        return Feel::detail::nedelecTriangularFaceMomentTransform(
                            Feel::detail::nedelecLocalToCanonicalFaceVertices(
                                element, attachment.entityId ),
                            static_cast<uint16_type>( M_publicOrder-1 ) );
                    else
                        return Feel::detail::nedelecTriangularFaceRTMomentTransform(
                            Feel::detail::nedelecLocalToCanonicalFaceVertices(
                                element, attachment.entityId ), M_publicOrder );
                }();
                return DofTransform::matrixOperator( std::move( matrix ), attachment.ordinal );
            }
        }
        return finiteElementEntityOrientationTransform( *this, element, localDofId );
    }

    using local_interpolant_type = Eigen::VectorXd;
    using local_interpolants_type = Eigen::MatrixXd;
    /**
     * @brief Allocate a zeroed element interpolation vector.
     * @param n Number of local DoF blocks to allocate.
     * @return Vector with `n` times the local DoF count entries.
     */
    [[nodiscard]] local_interpolant_type localInterpolant( int n = 1 ) const
    {
        return local_interpolant_type::Zero( n*localDof() );
    }
    /**
     * @brief Allocate a zeroed collection of local interpolation vectors.
     * @param p Number of columns.
     * @param n Number of local DoF blocks per column.
     * @return Matrix with `n` times the local DoF count rows and `p` columns.
     */
    [[nodiscard]] local_interpolants_type localInterpolants( int p, int n = 1 ) const
    {
        return local_interpolants_type::Zero( n*localDof(), p );
    }
    /**
     * @brief Allocate a zeroed facet-closure interpolation vector.
     * @return Vector containing all edge and face-interior DoFs in one facet closure.
     */
    [[nodiscard]] local_interpolant_type faceLocalInterpolant() const
    {
        return local_interpolant_type::Zero( localFacetDof(), 1 );
    }

    /**
     * @brief Interpolate an expression with runtime Nedelec edge/face/cell DoFs.
     * @tparam ExprType Vector expression with a compatible geometric context.
     * @param expr Expression sampled at the dual sampling points.
     * @param[out] localValues Preallocated local DoF vector, cleared before interpolation.
     * @pre `localValues` has at least `localDof()` entries.
     */
    template<typename ExprType>
    void interpolate( ExprType& expr, local_interpolant_type& localValues ) const
    {
        localValues.setZero();
        Feel::detail::interpolateHCurlEdgePointDofs(
            expr, localValues, convex_type::numEdges, dofPerEdge(), M_tangent,
            [&expr]( uint16_type edge, uint16_type point, auto& tangent )
            { expr.geom()->piolaCovariantEdgeTangent( edge, point, tangent ); } );

        auto const& sampling = this->dual().samplingLayout();
        if constexpr ( nDim == 3 )
        {
            if constexpr ( Kind == NedelecKind::NED1 )
            {
                for ( auto const& block : sampling.faces )
                    Feel::detail::interpolateHCurlTriangularFaceMoments(
                        expr, localValues, M_refconvex, block.face,
                        static_cast<int>( M_publicOrder )-1, block.firstDof,
                        block.elementPointOffset, block.momentCount );
            }
            else if ( M_publicOrder > 0 )
            {
                RaviartThomasDynamicPolynomialSet<2, value_type> moments(
                    RuntimeOrder{ static_cast<uint16_type>( M_publicOrder-1 ) } );
                for ( auto const& block : sampling.faces )
                    Feel::detail::interpolateHCurlTriangularFaceRTMoments(
                        expr, localValues, M_refconvex, moments, block.face,
                        block.firstDof, block.elementPointOffset );
            }
        }
        const int momentOrder = nDim == 2 ? static_cast<int>( M_publicOrder )-1
                                          : static_cast<int>( M_publicOrder )-2;
        if ( momentOrder >= 0 )
        {
            if constexpr ( Kind == NedelecKind::NED1 )
            {
                const uint16_type dimMoment = convex_type::polyDims(
                    static_cast<uint16_type>( momentOrder ) );
                Pkp1_v_type ambientVector( RuntimeOrder{ internalOrder() } );
                auto const source = ambientVector.polynomialsUpToDimension( dimMoment );
                auto moments = Feel::detail::polynomialSetFrom<vectorial_polynomialset_type>(
                    source, internalOrder() );
                Feel::detail::interpolateHCurlCellMoments<nDim>(
                    expr, localValues, moments, internalOrder(),
                    sampling.cellFirstDof, sampling.cellPointOffset,
                    sampling.cellMomentCount, nComponents );
            }
            else
            {
                RaviartThomasDynamicPolynomialSet<nDim, value_type> moments(
                    RuntimeOrder{ static_cast<uint16_type>( momentOrder ) } );
                Feel::detail::interpolateHCurlCellMoments<nDim>(
                    expr, localValues, moments,
                    static_cast<uint16_type>( internalOrder()+1 ),
                    sampling.cellFirstDof, sampling.cellPointOffset,
                    sampling.cellMomentCount, nComponents );
            }
        }
    }

    /**
     * @brief Apply all runtime Nedelec dual functionals to domain basis rows.
     *
     * Matrix assembly and ordinary interpolation share `interpolate()`, so
     * edge, face, and cell moments have one static/runtime implementation.
     *
     * @tparam ExprType Basis expression with an element geometric context.
     * @param expression Basis expression sampled at the dual sampling points.
     * @param[out] localMatrix Preallocated matrix with domain basis rows and `localDof()` columns.
     */
    template<typename ExprType>
        requires ( std::remove_cvref_t<ExprType>::gmc_type::subEntityCoDim == 0 )
    void interpolateBasisFunction( ExprType&& expression,
                                   local_interpolants_type& localMatrix ) const
    {
        Feel::detail::interpolateBasisFunctionRows(
            *this, expression, localMatrix );
    }

    /**
     * @brief Apply edge and face-interior functionals on the current facet closure.
     * @tparam ExprType Vector expression with a facet geometric context.
     * @param expr Expression sampled at the current facet sampling points.
     * @param[out] localValues Preallocated facet DoF vector, cleared before interpolation.
     * @pre `localValues` has at least `localFacetDof()` entries.
     */
    template<typename ExprType>
    void faceInterpolate( ExprType& expr, local_interpolant_type& localValues ) const
    {
        localValues.setZero();
        auto geometry = expr.geom();
        Feel::detail::interpolateHCurlEdgePointDofs(
            expr, localValues, face_type::numEdges, dofPerEdge(), M_tangent,
            [&expr, &geometry]( uint16_type edge, uint16_type point, auto& tangent )
            {
                const int elementEdge = nDim <= 2
                    ? geometry->element().fToE( geometry->faceId(), geometry->faceId() )
                    : geometry->element().fToE( geometry->faceId(), edge );
                expr.geom()->piolaCovariantEdgeTangent( elementEdge, point, tangent );
            } );
        if constexpr ( nDim == 3 )
        {
            if ( dofPerFace() > 0 )
            {
                const uint16_type face = static_cast<uint16_type>( geometry->faceId() );
                auto const& blocks = this->dual().samplingLayout().faces;
                CHECK_LT( face, blocks.size() );
                auto const& block = blocks[face];
                if constexpr ( Kind == NedelecKind::NED1 )
                    Feel::detail::interpolateHCurlTriangularFaceMoments(
                        expr, localValues, M_refconvex, face,
                        static_cast<int>( M_publicOrder )-1,
                        static_cast<uint16_type>( face_type::numEdges*dofPerEdge() ),
                        block.facetPointOffset, block.momentCount );
                else
                {
                    RaviartThomasDynamicPolynomialSet<2, value_type> moments(
                        RuntimeOrder{ static_cast<uint16_type>( M_publicOrder-1 ) } );
                    Feel::detail::interpolateHCurlTriangularFaceRTMoments(
                        expr, localValues, M_refconvex, moments, face,
                        static_cast<uint16_type>( face_type::numEdges*dofPerEdge() ),
                        block.facetPointOffset );
                }
            }
        }
    }

private:
    uint16_type M_publicOrder = 0;
    reference_convex_type M_refconvex;
    mutable ublas::vector<value_type> M_tangent { nDim };
};

namespace detail
{
template<int Order, uint16_type Dim, NedelecKind Kind, typename T, uint16_type Tag>
struct NedelecBackend
{
    using type = Nedelec<Dim, static_cast<uint16_type>( Order ), Kind, T, Tag>;
};

template<uint16_type Dim, NedelecKind Kind, typename T, uint16_type Tag>
struct NedelecBackend<Dynamic, Dim, Kind, T, Tag>
{
    using type = NedelecDynamicSimplex<Dim, Kind, T, Tag>;
};
} // namespace detail

/**
 * @brief Unified public Nedelec shell for static and runtime order.
 * @tparam Dim reference-simplex dimension
 * @tparam Order public static order or `Dynamic`
 * @tparam Kind first- or second-kind family
 * @tparam T scalar type
 * @tparam Tag compatibility family tag
 */
template<uint16_type Dim, int Order, NedelecKind Kind = NedelecKind::NED1,
         typename T = double, uint16_type Tag = 0>
class NedelecElement : public detail::NedelecBackend<Order, Dim, Kind, T, Tag>::type
{
    using super = typename detail::NedelecBackend<Order, Dim, Kind, T, Tag>::type;

public:
    using super::super;
    /**
     * @brief Construct the selected backend at its default public order.
     */
    NedelecElement() = default;

    /** @return public Nedelec family order, distinct from the ambient polynomial degree. */
    [[nodiscard]] uint16_type order() const noexcept
    {
        if constexpr ( Order == Dynamic )
            return super::order();
        else
            return static_cast<uint16_type>( Order );
    }
};

} // fem
template<int Order,
         NedelecKind Kind=NedelecKind::NED1,
         uint16_type TheTAG=0>
class Nedelec : public FiniteElementFamilyDescriptor<Order, TheTAG>
{
public:
    using descriptor_base = FiniteElementFamilyDescriptor<Order, TheTAG>;
    static constexpr bool is_order_static = descriptor_base::is_order_static;
    static constexpr bool is_order_dynamic = descriptor_base::is_order_dynamic;
    static constexpr int nOrder_v = Order;

    template<uint16_type N,
             uint16_type R = N,
             typename T = double,
             typename Convex=Simplex<N>>
    struct apply
    {
        static_assert( Convex::is_simplex,
                       "Nedelec hypercube support is not implemented in feelpoly; use simplex Nedelec or add a tensor-product H(curl) implementation first." );
        static_assert( Order == Dynamic || Kind != NedelecKind::NED1 || N == 2 || Order <= 3,
                       "3D Nedelec first kind is currently validated through public order three." );
        static_assert( Order == Dynamic || Kind != NedelecKind::NED2 || N == 2 || Order <= 3,
                       "3D Nedelec second kind is currently validated through public order three." );
        typedef fem::NedelecElement<N,Order,Kind,T,TheTAG> result_type;
        typedef result_type type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef Nedelec<Order,Kind,TheNewTAG> type;
    };

    //typedef Lagrange<Order,Scalar> component_basis_type;
    typedef Lagrange<( Order == Dynamic ? Dynamic : Order + 1 ),Scalar> component_basis_type;

    static inline const NedelecKind kind =  Kind;
};

} // Feel
#endif /* FEELPP_NEDELEC_HPP */
