/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file brezzidouglasmarini.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-02-18
 */
#ifndef __BrezziDouglasMarini_H
#define __BrezziDouglasMarini_H 1

#include <array>
#include <vector>

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/pointset.hpp>

#include <feel/feelpoly/dualbasis.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/functionalset.hpp>
#include <feel/feelpoly/operations.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>
#include <feel/feeldiscr/doflayout.hpp>
#include <feel/feelpoly/fe.hpp>
#include <feel/feelpoly/hdivpolynomialset.hpp>
#include <feel/feelpoly/hdivinterpolation.hpp>
#include <feel/feelpoly/hdivfunctionals.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/fefamilyshell.hpp>

namespace Feel
{
namespace detail
{
[[nodiscard]] constexpr uint16_type
brezziDouglasMariniSimplexInternalOrder( uint16_type publicOrder )
{
    return checkedLocalIndexValue(
        static_cast<size_type>( publicOrder ) + 1 );
}

[[nodiscard]] constexpr size_type
brezziDouglasMariniSimplexTotalDofSize(
    uint16_type dim, uint16_type publicOrder ) noexcept
{
    const size_type internalOrder =
        static_cast<size_type>( publicOrder ) + 1;
    return ( dim == 2 )
               ? ( internalOrder + 1 ) * ( internalOrder + 2 )
               : ( internalOrder + 1 ) * ( internalOrder + 2 ) *
                     ( internalOrder + 3 ) / 2;
}

[[nodiscard]] constexpr RuntimeOrder
validatedBrezziDouglasMariniSimplexOrder(
    uint16_type dim, RuntimeOrder publicOrder )
{
    (void)checkedLocalIndexValue(
        static_cast<size_type>( publicOrder.value ) + 1 );
    (void)checkedLocalIndexValue(
        brezziDouglasMariniSimplexTotalDofSize(
            dim, publicOrder.value ) );
    return publicOrder;
}

[[nodiscard]] constexpr uint16_type
brezziDouglasMariniSimplexFacetDof( uint16_type dim, uint16_type publicOrder )
{
    const size_type internalOrder =
        static_cast<size_type>( publicOrder ) + 1;
    return checkedLocalIndexValue(
        ( dim == 2 )
            ? internalOrder + 1
            : ( internalOrder + 1 ) * ( internalOrder + 2 ) / 2 );
}

[[nodiscard]] constexpr uint16_type
brezziDouglasMariniSimplexInteriorDof( uint16_type dim, uint16_type publicOrder )
{
    const size_type internalOrder =
        static_cast<size_type>( publicOrder ) + 1;
    if ( internalOrder <= 1 )
        return 0;
    return checkedLocalIndexValue(
        ( dim == 2 )
            ? internalOrder * internalOrder - 1
            : ( internalOrder + 1 ) * ( internalOrder + 2 ) *
                  ( internalOrder - 1 ) / 2 );
}

[[nodiscard]] constexpr uint16_type
brezziDouglasMariniSimplexTotalDof( uint16_type dim, uint16_type publicOrder )
{
    return checkedLocalIndexValue(
        brezziDouglasMariniSimplexTotalDofSize( dim, publicOrder ) );
}

} // namespace detail

template<uint16_type N,
         uint16_type O,
         typename T = double,
         template<int, int, int> class Convex = Simplex,
         uint16_type TheTAG = 0>
class BrezziDouglasMariniPolynomialSet
    :
    public Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Convex>
{
    using super = Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Convex>;

public:
    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const uint16_type nComponents = super::nComponents;
    static inline const bool is_product = false;

    using value_type = typename super::value_type;
    using convex_type = typename super::convex_type;
    using matrix_type = typename super::matrix_type;
    using points_type = typename super::points_type;
    using Pkp1_v_type = Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Convex>;
    using vectorial_polynomialset_type = PolynomialSet<typename super::basis_type, Vectorial>;

    BrezziDouglasMariniPolynomialSet()
        :
        super()
    {}
};

template<uint16_type N,
         typename T = double,
         template<int, int, int> class Convex = Simplex,
         uint16_type TheTAG = 0>
class BrezziDouglasMariniDynamicPolynomialSet
    :
    public Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Convex>
{
    using super = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Convex>;

public:
    static inline const uint16_type nDim = super::nDim;
    static inline const int nOrder = Dynamic;
    static inline const uint16_type nComponents = super::nComponents;
    static inline const bool is_product = false;

    using value_type = typename super::value_type;
    using convex_type = typename super::convex_type;
    using matrix_type = typename super::matrix_type;
    using points_type = typename super::points_type;
    using Pkp1_v_type = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Convex>;
    using vectorial_polynomialset_type = PolynomialSet<typename super::basis_type, Vectorial, Dynamic>;

    static_assert( convex_type::is_simplex,
                   "BrezziDouglasMariniDynamicPolynomialSet currently implements the simplex BDM space only." );

    BrezziDouglasMariniDynamicPolynomialSet()
        :
        BrezziDouglasMariniDynamicPolynomialSet( RuntimeOrder{ 0 } )
    {}

    explicit BrezziDouglasMariniDynamicPolynomialSet( RuntimeOrder publicOrder )
        :
        super( RuntimeOrder{
            Feel::detail::brezziDouglasMariniSimplexInternalOrder(
                Feel::detail::validatedBrezziDouglasMariniSimplexOrder(
                    N, publicOrder ).value ) } ),
        M_publicOrder( publicOrder.value )
    {
        CHECK_EQ( this->polynomialDimension(),
                  Feel::detail::brezziDouglasMariniSimplexTotalDof(
                      N, publicOrder.value ) );
    }

    [[nodiscard]] uint16_type publicOrder() const noexcept
    {
        return M_publicOrder;
    }

    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexInternalOrder( M_publicOrder );
    }

private:
    uint16_type M_publicOrder = 0;
};

namespace fem
{

namespace detail
{
template<typename Basis,
         template<class, int, class> class PointSetType>
class BrezziDouglasMariniDual
    :
    public DualBasis<Basis>
{
    using super = DualBasis<Basis>;

public:
    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;

    using primal_space_type = typename super::primal_space_type;
    using value_type = typename primal_space_type::value_type;
    using points_type = typename primal_space_type::points_type;
    using matrix_type = typename primal_space_type::matrix_type;
    using convex_type = typename primal_space_type::template convex<nDim+nOrder>::type;
    using reference_convex_type = Reference<convex_type, nDim, nDim+nOrder, nDim, value_type>;
    using node_type = typename reference_convex_type::node_type;
    using pointset_type = PointSetType<convex_type, nOrder, value_type>;
    using interior_polyset_type = Nedelec1DynamicPolynomialSet<nDim, value_type>;

    static inline const uint16_type nbPtsPerVertex = 0;
    static inline const uint16_type nbPtsPerEdge = ( nDim == 2 ) ? reference_convex_type::nbPtsPerEdge : 0;
    static inline const uint16_type nbPtsPerFace2d = ( nOrder > 1 ) ? ( nOrder*nOrder - 1 ) : 0;
    static inline const uint16_type nbPtsPerFace3d = ( nDim == 3 ) ? reference_convex_type::nbPtsPerFace : 0;
    static inline const uint16_type nbPtsPerFace = ( nDim == 2 ) ? nbPtsPerFace2d : nbPtsPerFace3d;
    static inline const uint16_type nbPtsPerVolume = ( nDim == 3 && nOrder > 1 )
                                                      ? static_cast<uint16_type>( ( nOrder+1 ) * ( nOrder+2 ) * ( nOrder-1 ) / 2 )
                                                      : 0;
    static inline const uint16_type numPoints2d = static_cast<uint16_type>( ( nOrder+1 ) * ( nOrder+2 ) );
    static inline const uint16_type numPoints3d = static_cast<uint16_type>( ( nOrder+1 ) * ( nOrder+2 ) * ( nOrder+3 ) / 2 );
    static inline const uint16_type numPoints = ( nDim == 2 ) ? numPoints2d : numPoints3d;

    static inline const uint16_type nDofPerVertex = 0;
    static inline const uint16_type nDofPerEdge = nbPtsPerEdge;
    static inline const uint16_type nDofPerFace = nbPtsPerFace;
    static inline const uint16_type nDofPerVolume = nbPtsPerVolume;
    static inline const uint16_type nLocalDof = numPoints;

    BrezziDouglasMariniDual( primal_space_type const& primal )
        :
        super( primal ),
        M_convex_ref(),
        M_pts( nDim, numPoints ),
        M_pts_per_face( convex_type::numTopologicalFaces ),
        M_fset( primal )
    {
        const uint16_type interiorDof = nDim == 2
            ? nDofPerFace : nDofPerVolume;
        auto functionals = Feel::detail::buildHDivSimplexFunctionals(
            primal, M_convex_ref,
            nDim == 2 ? nDofPerEdge : nDofPerFace,
            interiorDof, nLocalDof,
            static_cast<uint16_type>( 2*( nOrder+1 ) ),
            M_pts, M_pts_per_face,
            [this]( uint16_type facet )
            { return M_convex_ref.makePoints( nDim-1, facet ); },
            []()
            {
                return interior_polyset_type(
                    RuntimeOrder{ static_cast<uint16_type>( nOrder-2 ) } );
            },
            [&primal]( auto const& polynomial )
            {
                return Feel::detail::makeHDivOrthonormalIntegralMomentFunctional(
                    primal, polynomial );
            } );
        M_fset.setFunctionalSet( functionals );
    }

    points_type const& points() const
    {
        return M_pts;
    }

    matrix_type operator()( primal_space_type const& pset ) const
    {
        return M_fset( pset );
    }

    points_type const& points( uint16_type f ) const
    {
        return M_pts_per_face[f];
    }
    ublas::matrix_column<points_type const> point( uint16_type f, uint32_type __i ) const
    {
        return ublas::column( M_pts_per_face[f], __i );
    }
    ublas::matrix_column<points_type> point( uint16_type f, uint32_type __i )
    {
        return ublas::column( M_pts_per_face[f], __i );
    }

private:
    reference_convex_type M_convex_ref;
    points_type M_pts;
    std::vector<points_type> M_pts_per_face;
    FunctionalSet<primal_space_type> M_fset;
};

template<typename Basis,
         template<class, int, class> class PointSetType>
class BrezziDouglasMariniDynamicDual
    :
    public DualBasis<Basis>
{
    using super = DualBasis<Basis>;

public:
    static inline const uint16_type nDim = super::nDim;
    static inline const int nOrder = Dynamic;

    using primal_space_type = typename super::primal_space_type;
    using value_type = typename primal_space_type::value_type;
    using points_type = typename primal_space_type::points_type;
    using matrix_type = typename primal_space_type::matrix_type;
    using convex_type = Simplex<nDim, 1, nDim>;
    using reference_convex_type = Reference<convex_type, nDim, 1, nDim, value_type>;
    using node_type = typename reference_convex_type::node_type;
    using pointset_type = PointSetType<convex_type, Dynamic, value_type>;
    using interior_polyset_type = Nedelec1DynamicPolynomialSet<nDim, value_type>;

    static inline const uint16_type nbPtsPerVertex = 0;
    static inline const uint16_type nbPtsPerEdge = ( nDim == 2 ) ? 2 : 0;
    static inline const uint16_type nbPtsPerFace = ( nDim == 3 ) ? 3 : 0;
    static inline const uint16_type nbPtsPerVolume = 0;
    static inline const uint16_type numPoints = ( nDim == 2 ) ? 6 : 10;

    static inline const uint16_type nDofPerVertex = 0;
    static inline const uint16_type nDofPerEdge = nbPtsPerEdge;
    static inline const uint16_type nDofPerFace = nbPtsPerFace;
    static inline const uint16_type nDofPerVolume = nbPtsPerVolume;
    static inline const uint16_type nLocalDof = numPoints;

    BrezziDouglasMariniDynamicDual( primal_space_type const& primal )
        :
        super( primal ),
        M_publicOrder( primal.publicOrder() ),
        M_internalOrder( Feel::detail::brezziDouglasMariniSimplexInternalOrder( M_publicOrder ) ),
        M_convex_ref(),
        M_pts( nDim, Feel::detail::brezziDouglasMariniSimplexTotalDof( nDim, M_publicOrder ) ),
        M_pts_per_face( convex_type::numTopologicalFaces ),
        M_fset( primal )
    {
        const uint16_type nFacetDof = facetDof();
        pointset_type facetPoints( RuntimeOrder{ static_cast<uint16_type>( nDim + M_internalOrder ) } );
        auto functionals = Feel::detail::buildHDivSimplexFunctionals(
            primal, M_convex_ref, nFacetDof, interiorDof(), localDof(),
            static_cast<uint16_type>( 2*( M_internalOrder+1 ) ),
            M_pts, M_pts_per_face,
            [&facetPoints]( uint16_type facet )
            { return facetPoints.pointsBySubEntity( nDim-1, facet ); },
            [this]()
            {
                return interior_polyset_type( RuntimeOrder{
                    static_cast<uint16_type>( M_internalOrder-2 ) } );
            },
            [&primal]( auto const& polynomial )
            {
                return Feel::detail::makeHDivOrthonormalIntegralMomentFunctional(
                    primal, polynomial );
            } );
        M_fset.setFunctionalSet( functionals );
    }

    [[nodiscard]] uint16_type publicOrder() const noexcept
    {
        return M_publicOrder;
    }

    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        return M_internalOrder;
    }

    [[nodiscard]] uint16_type facetDof() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexFacetDof( nDim, M_publicOrder );
    }

    [[nodiscard]] uint16_type interiorDof() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexInteriorDof( nDim, M_publicOrder );
    }

    [[nodiscard]] uint16_type localDof() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexTotalDof( nDim, M_publicOrder );
    }

    points_type const& points() const
    {
        return M_pts;
    }

    matrix_type operator()( primal_space_type const& pset ) const
    {
        return M_fset( pset );
    }

    points_type const& points( uint16_type f ) const
    {
        return M_pts_per_face[f];
    }
    ublas::matrix_column<points_type const> point( uint16_type f, uint32_type __i ) const
    {
        return ublas::column( M_pts_per_face[f], __i );
    }
    ublas::matrix_column<points_type> point( uint16_type f, uint32_type __i )
    {
        return ublas::column( M_pts_per_face[f], __i );
    }

private:
    uint16_type M_publicOrder = 0;
    uint16_type M_internalOrder = 1;
    reference_convex_type M_convex_ref;
    points_type M_pts;
    std::vector<points_type> M_pts_per_face;
    FunctionalSet<primal_space_type> M_fset;
};
} // namespace detail

template<uint16_type N,
         uint16_type O,
         typename T = double,
         template<int, int, int> class Convex = Simplex,
         uint16_type TheTAG = 0>
class BrezziDouglasMarini
    :
    public FiniteElement<BrezziDouglasMariniPolynomialSet<N, O, T, Convex>,
                         fem::detail::BrezziDouglasMariniDual,
                         PointSetEquiSpaced>,
    public HDivPolynomialSet,
    public std::enable_shared_from_this<BrezziDouglasMarini<N, O, T, Convex, TheTAG>>
{
    using super = FiniteElement<BrezziDouglasMariniPolynomialSet<N, O, T, Convex>,
                                fem::detail::BrezziDouglasMariniDual,
                                PointSetEquiSpaced>;

public:
    BOOST_STATIC_ASSERT( N > 1 );

    static inline const uint16_type nDim = N;
    static inline const bool isTransformationEquivalent = true;
    static inline const bool isContinuous = true;
    typedef Continuous continuity_type;
    static const uint16_type TAG = TheTAG;

    using value_type = typename super::value_type;
    using primal_space_type = typename super::primal_space_type;
    using dual_space_type = typename super::dual_space_type;
    using polyset_type = typename super::polyset_type;

    static inline const bool is_vectorial = polyset_type::is_vectorial;
    static inline const bool is_scalar = polyset_type::is_scalar;
    static inline const uint16_type nComponents = polyset_type::nComponents;
    static inline const bool is_product = false;

    using convex_type = typename dual_space_type::convex_type;
    using pointset_type = typename dual_space_type::pointset_type;
    using reference_convex_type = typename dual_space_type::reference_convex_type;
    using node_type = typename reference_convex_type::node_type;
    using points_type = typename reference_convex_type::points_type;
    using face_type = typename convex_type::topological_face_type;

    static inline const uint16_type nOrder = dual_space_type::nOrder;
    static inline const uint16_type nbPtsPerVertex = 0;
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

    BrezziDouglasMarini()
        :
        super( dual_space_type( primal_space_type() ) ),
        M_refconvex()
    {
    }

    template<int subN>
    struct SubSpace
    {
        typedef BrezziDouglasMarini<N-1, O, T, Convex, TheTAG> type;
    };

    struct SSpace
    {
        typedef BrezziDouglasMarini<N, O, T, Convex, TheTAG> type;
    };

    template<uint16_type NewDim>
    struct ChangeDim
    {
        typedef BrezziDouglasMarini<NewDim, O, T, Convex, TheTAG> type;
    };

    BrezziDouglasMarini( BrezziDouglasMarini const& cr )
        :
        super( cr ),
        M_refconvex()
    {
    }
    ~BrezziDouglasMarini() override {}

    reference_convex_type const& referenceConvex() const
    {
        return M_refconvex;
    }

    std::string familyName() const override
    {
        return "brezzidouglasmarini";
    }

    uint16_type component( uint16_type /*localDofId*/ ) const override
    {
        return 0;
    }

    uint16_type dofParent( uint16_type localDofId ) const override
    {
        return localDofId;
    }

    bool dofHasRepresentativePoint( uint16_type localDofId ) const override
    {
        (void)localDofId;
        return false;
    }

    uint16_type dofFunctionalKind( uint16_type localDofId ) const override
    {
        auto const attachment = this->dofAttachment( localDofId );
        if ( !attachment.isValid() )
            return static_cast<uint16_type>( DofFunctionalKind::Other );
        if ( static_cast<uint16_type>( attachment.entityDim ) == nDim )
            return static_cast<uint16_type>( DofFunctionalKind::InteriorMoment );
        return static_cast<uint16_type>( DofFunctionalKind::NormalMoment );
    }

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

    template<typename ElementType>
    [[nodiscard]] DofTransform dofTransform( ElementType const& element, uint16_type localDofId ) const
    {
        return finiteElementEntityOrientationTransform( *this, element, localDofId );
    }

    uint16_type dofType( uint16_type /*localDofId*/ ) const override
    {
        // Keep legacy kind value for compatibility with existing RT/Nedelec paths.
        return 1;
    }

    typedef Eigen::VectorXd local_interpolant_type;
    local_interpolant_type localInterpolant( int n = 1 ) const
    {
        return local_interpolant_type::Zero( n*nLocalDof );
    }

    typedef Eigen::MatrixXd local_interpolants_type;
    local_interpolants_type localInterpolants( int p, int n = 1 ) const
    {
        return local_interpolants_type::Zero( n*nLocalDof, p );
    }

    template<typename ExprType>
    void interpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
    {
        Ihloc.setZero();

        const uint16_type nFacetDof = ( nDim == 2 ) ? nDofPerEdge : nDofPerFace;
        Feel::detail::interpolateHDivFacetDofs(
            expr, Ihloc, convex_type::numTopologicalFaces, nFacetDof, n );

        if constexpr ( nOrder > 1 )
        {
            static const uint16_type interiorOrder = nOrder - 2;
            typedef Nedelec1DynamicPolynomialSet<nDim, value_type, TheTAG> interior_polyset_type;
            interior_polyset_type interiorSpace( RuntimeOrder{ interiorOrder } );

            IM<nDim, 2*( nOrder + 1 ), value_type, Simplex> im;
            const uint16_type nInternalDof = ( nDim == 2 ) ? nDofPerFace : nDofPerVolume;
            const int firstInternalDof = ( nDim == 2 )
                                         ? convex_type::numTopologicalFaces * nDofPerEdge
                                         : convex_type::numTopologicalFaces * nDofPerFace;
            Feel::detail::interpolateHDivInteriorMomentDofs(
                expr, Ihloc, interiorSpace, im, firstInternalDof,
                nInternalDof, nComponents );
        }
    }


    /** Apply every static BDM dual functional to domain basis rows. */
    template<typename ExprType>
    void interpolateBasisFunction( ExprType&& expression,
                                   local_interpolants_type& localMatrix ) const
    {
        Feel::detail::interpolateBasisFunctionRows(
            *this, expression, localMatrix );
    }

    local_interpolant_type faceLocalInterpolant() const
    {
        return local_interpolant_type::Zero( nLocalFaceDof, 1 );
    }

    template<typename ExprType>
    void faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
    {
        Ihloc.setZero();

        const uint16_type nLocalDof = ( nDim == 2 ) ? nDofPerEdge : nDofPerFace;
        Feel::detail::interpolateHDivFacetDofs(
            expr, Ihloc, 1, nLocalDof, n );
    }

protected:
    reference_convex_type M_refconvex;
    mutable ublas::vector<value_type> n{ nDim };
};

template<uint16_type N,
         typename T = double,
         template<int, int, int> class Convex = Simplex,
         uint16_type TheTAG = 0>
class BrezziDouglasMariniDynamicSimplex
    :
    public FiniteElement<BrezziDouglasMariniDynamicPolynomialSet<N, T, Convex, TheTAG>,
                         fem::detail::BrezziDouglasMariniDynamicDual,
                         PointSetEquiSpaced>,
    public HDivPolynomialSet,
    public std::enable_shared_from_this<BrezziDouglasMariniDynamicSimplex<N, T, Convex, TheTAG>>
{
    using super = FiniteElement<BrezziDouglasMariniDynamicPolynomialSet<N, T, Convex, TheTAG>,
                                fem::detail::BrezziDouglasMariniDynamicDual,
                                PointSetEquiSpaced>;

public:
    BOOST_STATIC_ASSERT( N > 1 );

    static inline const uint16_type nDim = N;
    static inline const bool isTransformationEquivalent = true;
    static inline const bool isContinuous = true;
    using continuity_type = Continuous;
    static const uint16_type TAG = TheTAG;

    using value_type = typename super::value_type;
    using primal_space_type = typename super::primal_space_type;
    using dual_space_type = typename super::dual_space_type;
    using polyset_type = typename super::polyset_type;
    using interior_polyset_type = Nedelec1DynamicPolynomialSet<nDim, value_type, TheTAG>;

    static inline const bool is_order_static = false;
    static inline const bool is_order_dynamic = true;
    static constexpr int nOrder_v = Dynamic;

    static inline const bool is_vectorial = polyset_type::is_vectorial;
    static inline const bool is_scalar = polyset_type::is_scalar;
    static inline const uint16_type nComponents = polyset_type::nComponents;
    static inline const bool is_product = false;

    using convex_type = typename dual_space_type::convex_type;
    using pointset_type = typename dual_space_type::pointset_type;
    using reference_convex_type = typename dual_space_type::reference_convex_type;
    using node_type = typename reference_convex_type::node_type;
    using points_type = typename reference_convex_type::points_type;
    using face_type = typename convex_type::topological_face_type;

    static inline const int nOrder = Dynamic;
    static inline const uint16_type nbPtsPerVertex = 0;
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

    BrezziDouglasMariniDynamicSimplex()
        :
        BrezziDouglasMariniDynamicSimplex( RuntimeOrder{ 0 } )
    {}

    explicit BrezziDouglasMariniDynamicSimplex( RuntimeOrder order )
        :
        super( dual_space_type( primal_space_type( order ) ) ),
        M_publicOrder( order.value ),
        M_refconvex()
    {}

    template<int subN>
    struct SubSpace
    {
        using type = BrezziDouglasMariniDynamicSimplex<N-1, T, Convex, TheTAG>;
    };

    struct SSpace
    {
        using type = BrezziDouglasMariniDynamicSimplex<N, T, Convex, TheTAG>;
    };

    template<uint16_type NewDim>
    struct ChangeDim
    {
        using type = BrezziDouglasMariniDynamicSimplex<NewDim, T, Convex, TheTAG>;
    };

    BrezziDouglasMariniDynamicSimplex( BrezziDouglasMariniDynamicSimplex const& cr )
        :
        super( cr ),
        M_publicOrder( cr.M_publicOrder ),
        M_refconvex()
    {}
    ~BrezziDouglasMariniDynamicSimplex() override {}

    reference_convex_type const& referenceConvex() const
    {
        return M_refconvex;
    }

    std::string familyName() const override
    {
        return "brezzidouglasmarini";
    }

    [[nodiscard]] uint16_type order() const noexcept
    {
        return M_publicOrder;
    }

    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexInternalOrder( M_publicOrder );
    }

    [[nodiscard]] uint16_type localDof() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexTotalDof( nDim, M_publicOrder );
    }

    uint16_type localDofPerComponent() const override
    {
        return localDof();
    }

    uint16_type localDofCount( bool perComponent = false ) const override
    {
        (void)perComponent;
        return localDof();
    }

    [[nodiscard]] uint16_type dofPerVertex() const noexcept
    {
        return 0;
    }

    [[nodiscard]] uint16_type dofPerEdge() const noexcept
    {
        return ( nDim == 2 ) ? Feel::detail::brezziDouglasMariniSimplexFacetDof( nDim, M_publicOrder ) : 0;
    }

    [[nodiscard]] uint16_type dofPerFace() const noexcept
    {
        if constexpr ( nDim == 2 )
            return Feel::detail::brezziDouglasMariniSimplexInteriorDof( nDim, M_publicOrder );
        else
            return Feel::detail::brezziDouglasMariniSimplexFacetDof( nDim, M_publicOrder );
    }

    [[nodiscard]] uint16_type dofPerVolume() const noexcept
    {
        return ( nDim == 3 ) ? Feel::detail::brezziDouglasMariniSimplexInteriorDof( nDim, M_publicOrder ) : 0;
    }

    [[nodiscard]] uint16_type localFacetDof() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexFacetDof( nDim, M_publicOrder );
    }

    uint16_type component( uint16_type /*localDofId*/ ) const override
    {
        return 0;
    }

    uint16_type dofParent( uint16_type localDofId ) const override
    {
        return localDofId;
    }

    bool dofHasRepresentativePoint( uint16_type localDofId ) const override
    {
        (void)localDofId;
        return false;
    }

    uint16_type dofFunctionalKind( uint16_type localDofId ) const override
    {
        auto const attachment = this->dofAttachment( localDofId );
        if ( !attachment.isValid() )
            return static_cast<uint16_type>( DofFunctionalKind::Other );
        if ( static_cast<uint16_type>( attachment.entityDim ) == nDim )
            return static_cast<uint16_type>( DofFunctionalKind::InteriorMoment );
        return static_cast<uint16_type>( DofFunctionalKind::NormalMoment );
    }

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

    template<typename ElementType>
    [[nodiscard]] DofTransform dofTransform( ElementType const& element, uint16_type localDofId ) const
    {
        return finiteElementEntityOrientationTransform( *this, element, localDofId );
    }

    uint16_type dofType( uint16_type /*localDofId*/ ) const override
    {
        return 1;
    }

    using local_interpolant_type = Eigen::VectorXd;
    local_interpolant_type localInterpolant( int n = 1 ) const
    {
        return local_interpolant_type::Zero( n*localDof() );
    }

    using local_interpolants_type = Eigen::MatrixXd;
    local_interpolants_type localInterpolants( int p, int n = 1 ) const
    {
        return local_interpolants_type::Zero( n*localDof(), p );
    }

    template<typename ExprType>
    void interpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
    {
        Ihloc.setZero();

        const uint16_type nFacetDof = localFacetDof();
        Feel::detail::interpolateHDivFacetDofs(
            expr, Ihloc, convex_type::numTopologicalFaces, nFacetDof, n );

        if ( internalOrder() > 1 )
        {
            interior_polyset_type interiorSpace( RuntimeOrder{ static_cast<uint16_type>( internalOrder() - 2 ) } );
            IMGeneral<nDim, value_type, Simplex> im( static_cast<uint16_type>( 2*( internalOrder() + 1 ) ) );
            const uint16_type nInternalDof = ( nDim == 2 ) ? dofPerFace() : dofPerVolume();
            const int firstInternalDof = convex_type::numTopologicalFaces * nFacetDof;
            Feel::detail::interpolateHDivInteriorMomentDofs(
                expr, Ihloc, interiorSpace, im, firstInternalDof,
                nInternalDof, nComponents );
        }
    }


    /** Apply every runtime BDM dual functional to domain basis rows. */
    template<typename ExprType>
    void interpolateBasisFunction( ExprType&& expression,
                                   local_interpolants_type& localMatrix ) const
    {
        Feel::detail::interpolateBasisFunctionRows(
            *this, expression, localMatrix );
    }

    local_interpolant_type faceLocalInterpolant() const
    {
        return local_interpolant_type::Zero( localFacetDof(), 1 );
    }

    template<typename ExprType>
    void faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
    {
        Ihloc.setZero();

        const uint16_type nLocalDof = localFacetDof();
        Feel::detail::interpolateHDivFacetDofs(
            expr, Ihloc, 1, nLocalDof, n );
    }

protected:
    uint16_type M_publicOrder = 0;
    reference_convex_type M_refconvex;
    mutable ublas::vector<value_type> n{ nDim };
};

namespace detail
{
template<int Order, uint16_type Dim, typename T,
         template<int, int, int> class Convex, uint16_type Tag>
struct BrezziDouglasMariniBackend
{
    using type = BrezziDouglasMarini<Dim, static_cast<uint16_type>( Order ), T, Convex, Tag>;
};

template<uint16_type Dim, typename T,
         template<int, int, int> class Convex, uint16_type Tag>
struct BrezziDouglasMariniBackend<Dynamic, Dim, T, Convex, Tag>
{
    using type = BrezziDouglasMariniDynamicSimplex<Dim, T, Convex, Tag>;
};
} // namespace detail

/**
 * @brief Unified public BDM finite-element shell for static and runtime order.
 * @tparam Dim reference-cell dimension
 * @tparam Order public static order or `Dynamic`
 * @tparam T scalar type
 * @tparam Convex reference-cell family
 * @tparam Tag compatibility family tag
 */
template<uint16_type Dim, int Order, typename T = double,
         template<int, int, int> class Convex = Simplex, uint16_type Tag = 0>
class BrezziDouglasMariniElement :
    public detail::BrezziDouglasMariniBackend<Order, Dim, T, Convex, Tag>::type
{
    using super = typename detail::BrezziDouglasMariniBackend<Order, Dim, T, Convex, Tag>::type;

public:
    using super::super;
    BrezziDouglasMariniElement() = default;

    /** @return public BDM family order, distinct from the internal polynomial degree. */
    [[nodiscard]] uint16_type order() const noexcept
    {
        if constexpr ( Order == Dynamic )
            return super::order();
        else
            return static_cast<uint16_type>( Order );
    }
};

} // namespace fem

template<int Order,
         uint16_type TheTAG = 0>
class BrezziDouglasMarini : public FiniteElementFamilyDescriptor<Order, TheTAG>
{
public:
    using descriptor_base = FiniteElementFamilyDescriptor<Order, TheTAG>;
    static constexpr bool is_order_static = descriptor_base::is_order_static;
    static constexpr bool is_order_dynamic = descriptor_base::is_order_dynamic;
    static constexpr int nOrder_v = Order;

    template<uint16_type N,
             uint16_type R = N,
             typename T = double,
             typename Convex = Simplex<N>>
    struct apply
    {
        static_assert( Convex::is_simplex,
                       "BrezziDouglasMarini hypercube support is not implemented in feelpoly; use simplex BDM or add a tensor-product H(div) implementation first." );
        using result_type = fem::BrezziDouglasMariniElement<N, Order, T, Simplex, TheTAG>;
        using type = result_type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef BrezziDouglasMarini<Order, TheNewTAG> type;
    };

    typedef Lagrange<Order, Scalar> component_basis_type;

};

} // namespace Feel
#endif /* __BrezziDouglasMarini_H */
