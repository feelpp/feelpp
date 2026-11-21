/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2006-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file raviartthomas.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-01-14
 */
#ifndef __RaviartThomas_H
#define __RaviartThomas_H 1

#include <array>

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

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/lu.hpp>

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
#include <feel/feelpoly/meta.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/fefamilyshell.hpp>

namespace Feel
{
namespace detail
{
template<typename P>
struct times_x
{
    typedef typename P::value_type value_type;
    typedef typename P::points_type points_type;
    times_x ( P const& p, int c  )
        :
        M_p ( p ),
        M_c( c )
    {
        //VLOG(1) << "component : " << c << std::endl;
    }
    typename ublas::vector<value_type> operator() ( points_type const& __pts ) const
    {
#if 0
        VLOG(1) << "times_x(pts) : " << __pts << std::endl;
        VLOG(1) << "times_x(pts) : " << M_p.evaluate( __pts ) << std::endl;
        VLOG(1) << "times_x(coeff) : " << M_p.coefficients() << std::endl;
#endif
        // __pts[c] * p( __pts )
        return ublas::element_prod( ublas::row( __pts, M_c ),
                                    ublas::row( M_p.evaluate( __pts ), 0 ) );
    }
    //P const& M_p;
    P M_p;
    int M_c;
};

template<uint16_type Dim, uint16_type InternalOrder>
struct RaviartThomasSimplexDofLayout
{
    static_assert( Dim == 2 || Dim == 3,
                   "Raviart-Thomas simplex layout is implemented for 2D triangles and 3D tetrahedra." );
    static_assert( InternalOrder > 0,
                   "Feel++ RaviartThomas<O> uses internal polynomial order O+1; the internal order must be positive." );

    static constexpr uint16_type publicOrder = static_cast<uint16_type>( InternalOrder - 1 );
    static constexpr uint16_type facetDof = ( Dim == 2 )
                                                ? InternalOrder
                                                : static_cast<uint16_type>( InternalOrder * ( InternalOrder + 1 ) / 2 );
    static constexpr uint16_type totalDof = ( Dim == 2 )
                                                ? static_cast<uint16_type>( InternalOrder * ( InternalOrder + 2 ) )
                                                : static_cast<uint16_type>( InternalOrder * ( InternalOrder + 1 ) * ( InternalOrder + 3 ) / 2 );
    static constexpr uint16_type interiorDof = ( Dim == 2 )
                                                   ? static_cast<uint16_type>( InternalOrder * ( InternalOrder - 1 ) )
                                                   : static_cast<uint16_type>( InternalOrder * ( InternalOrder + 1 ) * ( InternalOrder - 1 ) / 2 );
};

[[nodiscard]] constexpr uint16_type
raviartThomasSimplexInternalOrder( uint16_type publicOrder )
{
    return checkedLocalIndexValue(
        static_cast<size_type>( publicOrder ) + 1 );
}

[[nodiscard]] constexpr size_type
raviartThomasSimplexTotalDofSize( uint16_type dim,
                                  uint16_type publicOrder ) noexcept
{
    const size_type internalOrder =
        static_cast<size_type>( publicOrder ) + 1;
    return ( dim == 2 )
               ? internalOrder * ( internalOrder + 2 )
               : internalOrder * ( internalOrder + 1 ) *
                     ( internalOrder + 3 ) / 2;
}

[[nodiscard]] constexpr RuntimeOrder
validatedRaviartThomasSimplexOrder( uint16_type dim,
                                    RuntimeOrder publicOrder )
{
    (void)checkedLocalIndexValue(
        static_cast<size_type>( publicOrder.value ) + 1 );
    (void)checkedLocalIndexValue(
        raviartThomasSimplexTotalDofSize( dim, publicOrder.value ) );
    return publicOrder;
}

[[nodiscard]] constexpr uint16_type
raviartThomasSimplexFacetDof( uint16_type dim, uint16_type publicOrder )
{
    const size_type internalOrder =
        static_cast<size_type>( publicOrder ) + 1;
    return checkedLocalIndexValue(
        ( dim == 2 )
            ? internalOrder
            : internalOrder * ( internalOrder + 1 ) / 2 );
}

[[nodiscard]] constexpr uint16_type
raviartThomasSimplexTotalDof( uint16_type dim, uint16_type publicOrder )
{
    return checkedLocalIndexValue(
        raviartThomasSimplexTotalDofSize( dim, publicOrder ) );
}

[[nodiscard]] constexpr uint16_type
raviartThomasSimplexInteriorDof( uint16_type dim, uint16_type publicOrder )
{
    const size_type internalOrder =
        static_cast<size_type>( publicOrder ) + 1;
    return checkedLocalIndexValue(
        ( dim == 2 )
            ? internalOrder * ( internalOrder - 1 )
            : internalOrder * ( internalOrder + 1 ) *
                  ( internalOrder - 1 ) / 2 );
}

/**
 * @brief Construct RT primal coefficients for static or runtime-order storage.
 *
 * The storage types are policy parameters; `(P_k)^d + x P_k^hom` is built by
 * this single mathematical implementation for both public FE forms.
 */
template<uint16_type Dim, typename AmbientVector, typename AmbientScalar,
         typename VectorPolynomialSet, typename ScalarPolynomialSet>
[[nodiscard]] auto
buildRaviartThomasPrimalCoefficients( AmbientVector const& ambientVector,
                                      AmbientScalar const& ambientScalar,
                                      uint16_type publicOrder )
{
    using value_type = typename AmbientVector::value_type;
    using convex_type = typename AmbientVector::convex_type;
    using scalar_polynomial_type = typename ScalarPolynomialSet::polynomial_type;
    const uint16_type internalOrder = raviartThomasSimplexInternalOrder( publicOrder );
    const uint16_type dimPk = convex_type::polyDims( publicOrder );
    const uint16_type dimPkm1 = publicOrder == 0
        ? 0 : convex_type::polyDims( publicOrder-1 );

    VectorPolynomialSet pkVector = polynomialSetFrom<VectorPolynomialSet>(
        ambientVector.polynomialsUpToDimension( dimPk ), internalOrder );
    ScalarPolynomialSet pkScalar = polynomialSetFrom<ScalarPolynomialSet>(
        ambientScalar.polynomialsUpToDimension( dimPk ), internalOrder );

    IMGeneral<Dim, value_type, Simplex> quadrature(
        static_cast<uint16_type>( 2*internalOrder ) );
    ublas::matrix<value_type> enrichmentCoefficients(
        Dim*( dimPk-dimPkm1 ), pkScalar.coeff().size2() );
    for ( int scalar = dimPkm1, candidate = 0;
          scalar < dimPk; ++scalar, ++candidate )
        for ( int component = 0; component < Dim; ++component )
        {
            times_x<scalar_polynomial_type> product(
                pkScalar.polynomial( scalar ), component );
            ublas::row( enrichmentCoefficients, candidate*Dim+component ) =
                ublas::row( Feel::project(
                    ambientScalar, product, quadrature ).coeff(), 0 );
        }

    VectorPolynomialSet enrichment(
        typename AmbientVector::basis_type(), enrichmentCoefficients, true );
    if constexpr ( requires( VectorPolynomialSet& p ) {
                       p.setOrder( internalOrder );
                   } )
        enrichment.setOrder( internalOrder );
    auto result = unite( pkVector, enrichment );
    CHECK_EQ( result.polynomialDimension(),
              raviartThomasSimplexTotalDof( Dim, publicOrder ) );
    return result.coeff();
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
         typename T = double,
         template<int, int, int> class Convex = Simplex,
         uint16_type TheTAG = 0>
class RaviartThomasPolynomialSet
    :
    public Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Convex>
{
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Convex> super;

public:
    static const uint16_type Om1 = (O==0)?0:O-1;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Vectorial, T, TheTAG, Convex> Pk_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Vectorial, T, TheTAG, Convex> Pkp1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, Om1, N, Vectorial, T, TheTAG, Convex> Pkm1_v_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O, N, Scalar, T, TheTAG, Convex> Pk_s_type;
    typedef Feel::detail::OrthonormalPolynomialSet<N, O+1, N, Scalar, T, TheTAG, Convex> Pkp1_s_type;

    typedef PolynomialSet<typename super::basis_type,Vectorial> vectorial_polynomialset_type;
    typedef typename vectorial_polynomialset_type::polynomial_type vectorial_polynomial_type;
    typedef PolynomialSet<typename super::basis_type,Scalar> scalar_polynomialset_type;
    typedef typename scalar_polynomialset_type::polynomial_type scalar_polynomial_type;

    typedef RaviartThomasPolynomialSet<N, O, T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::convex_type convex_type;
    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;

    static_assert( convex_type::is_simplex,
                   "RaviartThomasPolynomialSet currently implements the simplex RT space only." );

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder = super::nOrder;
    static inline const uint16_type nComponents = super::nComponents;
    static inline const bool is_product = false;
    RaviartThomasPolynomialSet()
        :
        super()
    {
        Pkp1_v_type ambientVector;
        Pkp1_s_type ambientScalar;
        this->setCoefficient(
            Feel::detail::buildRaviartThomasPrimalCoefficients<
                N, Pkp1_v_type, Pkp1_s_type,
                vectorial_polynomialset_type, scalar_polynomialset_type>(
                    ambientVector, ambientScalar, O ), true );
    }


};

template<uint16_type N,
         typename T = double,
         template<int, int, int> class Convex = Simplex,
         uint16_type TheTAG = 0>
class RaviartThomasDynamicPolynomialSet
    :
    public Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Convex>
{
    using super = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Convex>;

public:
    using Pkp1_v_type = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Convex>;
    using Pkp1_s_type = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Scalar, T, TheTAG, Convex>;
    using vectorial_polynomialset_type = PolynomialSet<typename super::basis_type, Vectorial, Dynamic>;
    using vectorial_polynomial_type = typename vectorial_polynomialset_type::polynomial_type;
    using scalar_polynomialset_type = PolynomialSet<typename super::basis_type, Scalar, Dynamic>;
    using scalar_polynomial_type = typename scalar_polynomialset_type::polynomial_type;

    using value_type = typename super::value_type;
    using convex_type = typename super::convex_type;
    using matrix_type = typename super::matrix_type;
    using points_type = typename super::points_type;

    static_assert( convex_type::is_simplex,
                   "RaviartThomasDynamicPolynomialSet currently implements the simplex RT space only." );

    static inline const uint16_type nDim = super::nDim;
    static inline const int nOrder = Dynamic;
    static inline const uint16_type nComponents = super::nComponents;
    static inline const bool is_product = false;

    RaviartThomasDynamicPolynomialSet()
        :
        RaviartThomasDynamicPolynomialSet( RuntimeOrder{ 0 } )
    {}

    explicit RaviartThomasDynamicPolynomialSet( RuntimeOrder publicOrder )
        :
        super( RuntimeOrder{
            Feel::detail::raviartThomasSimplexInternalOrder(
                Feel::detail::validatedRaviartThomasSimplexOrder(
                    N, publicOrder ).value ) } )
    {
        build( publicOrder.value );
    }

    [[nodiscard]] uint16_type publicOrder() const noexcept
    {
        return M_publicOrder;
    }

    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        return Feel::detail::raviartThomasSimplexInternalOrder( M_publicOrder );
    }

private:
    void build( uint16_type publicOrder )
    {
        M_publicOrder = publicOrder;
        const uint16_type internalOrder = Feel::detail::raviartThomasSimplexInternalOrder( publicOrder );
        Pkp1_v_type ambientVector( RuntimeOrder{ internalOrder } );
        Pkp1_s_type ambientScalar( RuntimeOrder{ internalOrder } );
        this->setCoefficient(
            Feel::detail::buildRaviartThomasPrimalCoefficients<
                N, Pkp1_v_type, Pkp1_s_type,
                vectorial_polynomialset_type, scalar_polynomialset_type>(
                    ambientVector, ambientScalar, publicOrder ), true );
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
class RaviartThomasDual
    :
public DualBasis<Basis>
{
    typedef DualBasis<Basis> super;
public:

    static inline const uint16_type nDim = super::nDim;
    static inline const uint16_type nOrder= super::nOrder;
    using dof_layout = Feel::detail::RaviartThomasSimplexDofLayout<nDim, nOrder>;

    typedef typename super::primal_space_type primal_space_type;
    typedef typename primal_space_type::value_type value_type;
    typedef typename primal_space_type::points_type points_type;
    typedef typename primal_space_type::matrix_type matrix_type;
    typedef typename primal_space_type::template convex<nDim+nOrder-1>::type convex_type;
    typedef Reference<convex_type, nDim, nDim+nOrder-1, nDim, value_type> reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;

    typedef typename primal_space_type::Pkp1_v_type Pkp1_v_type;
    typedef typename primal_space_type::vectorial_polynomialset_type vectorial_polynomialset_type;

    // point set type associated with the functionals
    typedef PointSetType<convex_type, nOrder, value_type> pointset_type;

    static inline const uint16_type nbPtsPerVertex = 0;
    static constexpr uint16_type nbPtsPerEdge = (nDim==2)?dof_layout::facetDof:0;
    static constexpr uint16_type nbPtsPerFace2d = (nDim==2)?dof_layout::interiorDof:0;
    static constexpr uint16_type nbPtsPerFace3d = (nDim==3)?dof_layout::facetDof:0;
    static constexpr uint16_type nbPtsPerFace =(nDim==3)?nbPtsPerFace3d:nbPtsPerFace2d;
    static constexpr uint16_type nbInteriorMoments3d = (nDim==3)?dof_layout::interiorDof:0;
    static constexpr uint16_type nbPtsPerVolume = (nDim==3)?nbInteriorMoments3d:0;
    static constexpr uint16_type numPoints2d = (nDim==2)?dof_layout::totalDof:0;
    static constexpr uint16_type numPoints3d = (nDim==3)?dof_layout::totalDof:0;
    static constexpr uint16_type numPoints = (nDim==2)?numPoints2d:numPoints3d;

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

    RaviartThomasDual( primal_space_type const& primal )
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
                Pkp1_v_type ambient;
                const uint16_type dimension = convex_type::polyDims( nOrder-2 );
                return Feel::detail::polynomialSetFrom<vectorial_polynomialset_type>(
                    ambient.polynomialsUpToDimension( dimension ), nOrder );
            },
            [&primal]( auto const& polynomial )
            {
                return Feel::detail::makeHDivIntegralMomentFunctional(
                    primal, polynomial, static_cast<uint16_type>( nOrder+1 ) );
            } );
        M_fset.setFunctionalSet( functionals );
    }

    /**
     * \return the point set
     */
    //pointset_type const& pointSet() const { return M_pset; }

    points_type const& points() const
    {
        return M_pts;
    }


    matrix_type operator()( primal_space_type const& pset ) const
    {
        //VLOG(1) << "RT matrix = " << M_fset( pset ) << std::endl;
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
    /**
     * set the pointset at face \c f using points \c n
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
    FunctionalSet<primal_space_type> M_fset;


};

template<typename Basis,
         template<class, int, class> class PointSetType>
class RaviartThomasDynamicDual
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
    using Pkp1_v_type = typename primal_space_type::Pkp1_v_type;
    using vectorial_polynomialset_type = typename primal_space_type::vectorial_polynomialset_type;
    using pointset_type = PointSetType<convex_type, Dynamic, value_type>;

    static inline const uint16_type nbPtsPerVertex = 0;
    static constexpr uint16_type nbPtsPerEdge = ( nDim == 2 ) ? 1 : 0;
    static constexpr uint16_type nbPtsPerFace2d = 0;
    static constexpr uint16_type nbPtsPerFace3d = ( nDim == 3 ) ? 1 : 0;
    static constexpr uint16_type nbPtsPerFace = ( nDim == 3 ) ? nbPtsPerFace3d : nbPtsPerFace2d;
    static constexpr uint16_type nbPtsPerVolume = 0;
    static constexpr uint16_type numPoints = convex_type::numTopologicalFaces;

    static inline const uint16_type nDofPerVertex = 0;
    static inline const uint16_type nDofPerEdge = nbPtsPerEdge;
    static inline const uint16_type nDofPerFace = nbPtsPerFace;
    static inline const uint16_type nDofPerVolume = nbPtsPerVolume;
    static inline const uint16_type nLocalDof = numPoints;

    RaviartThomasDynamicDual( primal_space_type const& primal )
        :
        super( primal ),
        M_publicOrder( primal.publicOrder() ),
        M_internalOrder( Feel::detail::raviartThomasSimplexInternalOrder( M_publicOrder ) ),
        M_convex_ref(),
        M_pts( nDim, Feel::detail::raviartThomasSimplexTotalDof( nDim, M_publicOrder ) ),
        M_pts_per_face( convex_type::numTopologicalFaces ),
        M_fset( primal )
    {
        const uint16_type nFacetDof = facetDof();
        pointset_type facetPoints( RuntimeOrder{ static_cast<uint16_type>( nDim + M_internalOrder - 1 ) } );
        auto functionals = Feel::detail::buildHDivSimplexFunctionals(
            primal, M_convex_ref, nFacetDof, interiorDof(), localDof(),
            static_cast<uint16_type>( 2*( M_internalOrder+1 ) ),
            M_pts, M_pts_per_face,
            [&facetPoints]( uint16_type facet )
            { return facetPoints.pointsBySubEntity( nDim-1, facet ); },
            [this]()
            {
                Pkp1_v_type ambient( RuntimeOrder{ M_internalOrder } );
                const uint16_type dimension = convex_type::polyDims(
                    static_cast<uint16_type>( M_internalOrder-2 ) );
                return Feel::detail::polynomialSetFrom<vectorial_polynomialset_type>(
                    ambient.polynomialsUpToDimension( dimension ), M_internalOrder );
            },
            [&primal, this]( auto const& polynomial )
            {
                return Feel::detail::makeHDivIntegralMomentFunctional(
                    primal, polynomial,
                    static_cast<uint16_type>( M_internalOrder+1 ) );
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
        return Feel::detail::raviartThomasSimplexFacetDof( nDim, M_publicOrder );
    }

    [[nodiscard]] uint16_type interiorDof() const noexcept
    {
        return Feel::detail::raviartThomasSimplexInteriorDof( nDim, M_publicOrder );
    }

    [[nodiscard]] uint16_type localDof() const noexcept
    {
        return Feel::detail::raviartThomasSimplexTotalDof( nDim, M_publicOrder );
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
}// detail


/**
 * \class RaviartThomas
 * \brief RaviartThomas Finite Element
 *
 * \f$ H(div)\f$  conforming element
 *
 * @author Christophe Prud'homme
 */
template<uint16_type N,
         uint16_type O,
         typename T = double,
         template<int, int, int> class Convex = Simplex,
         uint16_type TheTAG=0 >
class RaviartThomas
    :
public FiniteElement<RaviartThomasPolynomialSet<N, O, T, Convex>,
    fem::detail::RaviartThomasDual,
    PointSetEquiSpaced >,
public HDivPolynomialSet,
public std::enable_shared_from_this<RaviartThomas<N,O,T,Convex> >
{
    typedef FiniteElement<RaviartThomasPolynomialSet<N, O, T, Convex>,
            fem::detail::RaviartThomasDual,
            PointSetEquiSpaced > super;
public:

    BOOST_STATIC_ASSERT( N > 1 );

    /** @name Typedefs
     */
    //@{

    static inline const uint16_type nDim = N;
    //static inline const bool isTransformationEquivalent = false;
    static inline const bool isTransformationEquivalent = true;
    static inline const bool isContinuous = true;
    typedef Continuous continuity_type;
    static const uint16_type TAG = TheTAG;

    //static const polynomial_transformation_type transformation = POLYNOMIAL_CONTEXT_NEEDS_1ST_PIOLA_TRANSFORMATION;
    typedef typename super::value_type value_type;
    typedef typename super::primal_space_type primal_space_type;
    typedef typename super::dual_space_type dual_space_type;
    typedef typename primal_space_type::Pkp1_v_type Pkp1_v_type;
    typedef typename primal_space_type::vectorial_polynomialset_type vectorial_polynomialset_type;

    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename super::polyset_type polyset_type;
    static inline const bool is_vectorial = polyset_type::is_vectorial;
    static inline const bool is_scalar = polyset_type::is_scalar;
    static inline const uint16_type nComponents = polyset_type::nComponents;
    static inline const bool is_product = false;


    typedef typename dual_space_type::convex_type convex_type;
    static_assert( convex_type::is_simplex,
                   "RaviartThomas finite elements are implemented on simplices only; hypercube H(div) requires a tensor-product implementation." );
    typedef typename dual_space_type::pointset_type pointset_type;
    typedef typename dual_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename reference_convex_type::points_type points_type;
    typedef typename convex_type::topological_face_type face_type;

    static inline const uint16_type nOrder =  dual_space_type::nOrder;
    static inline const uint16_type nPublicOrder = O;
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
    //@}

    //@}

    /** @name Constructors, destructor
     */
    //@{

    RaviartThomas()
        :
        super( dual_space_type( primal_space_type() ) ),
        M_refconvex()
    {
#if 0
        VLOG(1) << "[RT] nPtsPerEdge = " << nbPtsPerEdge << "\n";
        VLOG(1) << "[RT] nPtsPerFace = " << nbPtsPerFace << "\n";
        VLOG(1) << "[RT] numPoints = " << numPoints << "\n";

        VLOG(1) << "[RT] nDof = " << super::nDof << "\n";

        VLOG(1) << "[RT] coeff : " << this->coeff() << "\n";
        VLOG(1) << "[RT] pts : " << this->points() << "\n";
        VLOG(1) << "[RT] eval at pts : " << this->evaluate( this->points() ) << "\n";
        VLOG(1) << "[RT] is_product : " << is_product << "\n";
#endif
    }

    template<int subN>
    struct SubSpace
    {
        typedef RaviartThomas<N-1, O, T, Convex, TheTAG> type;
    };

    struct SSpace
    {
        typedef  RaviartThomas<N, O, T, Convex, TheTAG> type;

    };

    template<uint16_type NewDim>
    struct ChangeDim
    {
        typedef  RaviartThomas<NewDim, O, T, Convex,  TheTAG> type;
    };

    RaviartThomas( RaviartThomas const & cr )
        :
        super( cr ),
        M_refconvex()
    {}
    ~RaviartThomas() override
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
     * \return the reference convex associated with the lagrange polynomials
     */
    reference_convex_type const& referenceConvex() const
    {
        return M_refconvex;
    }

    std::string familyName() const override
    {
        return "raviartthomas";
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

    //! \return the type of a local dof
    uint16_type dofType( uint16_type localDofId ) const override
        {
            // Keep legacy functional kind value for compatibility.
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
    local_interpolant_type
    localInterpolant( int n = 1 ) const
        {
            return local_interpolant_type::Zero( n*nLocalDof );
        }

    typedef Eigen::MatrixXd local_interpolants_type;
    local_interpolants_type
    localInterpolants( int p, int n = 1 ) const
        {
            return local_interpolants_type::Zero( n*nLocalDof, p );
        }

    template<typename ExprType>
    void
    interpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            Ihloc.setZero();

            const uint16_type nFacetDof = ( nDim == 2 ) ? nDofPerEdge : nDofPerFace;
            Feel::detail::interpolateHDivFacetDofs(
                expr, Ihloc, convex_type::numTopologicalFaces, nFacetDof, n );

            if constexpr ( nOrder > 1 )
            {
                static const uint16_type momentOrder = nOrder - 2;
                const uint16_type dimPm1 = convex_type::polyDims( momentOrder );
                Pkp1_v_type Pkp1;
                vectorial_polynomialset_type Pkm1( Pkp1.polynomialsUpToDimension( dimPm1 ) );

                const uint16_type nInternalDof = ( nDim == 2 ) ? nDofPerFace : nDofPerVolume;
                IM<nDim, 2*( nOrder + 1 ), value_type, Simplex> im;
                const int firstInternalDof = convex_type::numTopologicalFaces * nFacetDof;
                Feel::detail::interpolateHDivInteriorMomentDofs(
                    expr, Ihloc, Pkm1, im, firstInternalDof,
                    nInternalDof, nComponents );
            }
        }
    local_interpolant_type
    faceLocalInterpolant() const
        {
            return local_interpolant_type::Zero( nLocalFaceDof, 1 );
        }
    template<typename ExprType>
    void
    faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            Ihloc.setZero();

            const uint16_type nLocalDof = (nDim==2) ? nDofPerEdge : nDofPerFace;
            Feel::detail::interpolateHDivFacetDofs(
                expr, Ihloc, 1, nLocalDof, n );
        }

    using apply_curl_t = bool_c<true>;
    using apply_id_t = bool_c<false>;

    template<typename ExprType>
    void
    interpolateBasisFunction( ExprType& expr, local_interpolants_type& Ihloc ) const
    {
        Feel::detail::interpolateBasisFunctionRows( *this, expr, Ihloc );
    }

    template<typename ExprType>
    static auto
    //isomorphism( ExprType expr ) -> decltype( Feel::vf::detJ()*( trans( Feel::vf::JinvT() )*expr )*Feel::vf::Nref() )
    isomorphism( ExprType& expr ) -> decltype( expr )
    {
        using namespace Feel::vf;
        //return detJ()*( trans( JinvT() )*expr )*Nref();
        return expr;
    }
#if 0
    /**
     *
     * \return the value of the expression at the dof
     */
    template<typename ExprType, typename ContextType>
    std::vector<value_type>
    interpolate( std::shared_ptr<ContextType>& ctx, ExprType & expr )
    {
        using namespace Feel::vf;
        typedef std::shared_ptr<ContextType> gmc_ptrtype;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;

        std::vector<value_type> v( nLocalDof );
        auto I = integrate( face, _expr=_e1*N() );
        // First deal with the face dof
        for ( int face = 0; face < numTopologicalFaces; ++face )
        {
            // update the geomap at dof on face
            ctx->update( _face=face, _element=ctx->id() );

            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( ctx ) );
            expr.update( mapgmc, face );

            for ( int q = 0; q < nDofPerFace; ++q )
            {
                int ldof = nDofPerFace*face+i;
                v[ldof] = expr.evalq( 0,0,i );
            }
        }

        // evaluate expr \cdot n  on each face

        // evaluate moments of the expression
    }
#endif

#if 0
    template<typename GMContext, typename PC, typename Phi, typename GPhi, typename HPhi >
    static void transform( std::shared_ptr<GMContext> gmc,  std::shared_ptr<PC> const& pc,
                           Phi& phi_t,
                           GPhi& g_phi_t, const bool do_gradient,
                           HPhi& h_phi_t, const bool do_hessian

                         )
    {
        transform ( *gmc, *pc, phi_t, g_phi_t, do_gradient, h_phi_t, do_hessian );
    }
    template<typename GMContext, typename PC, typename Phi, typename GPhi, typename HPhi >
    static void transform( GMContext const& gmc,
                           PC const& pc,
                           Phi& phi_t,
                           GPhi& g_phi_t, const bool do_gradient,
                           HPhi& h_phi_t, const bool do_hessian

                         )
    {
        //phi_t = phi; return ;
        typename GMContext::gm_type::matrix_type const B = gmc.B( 0 );
        typename GMContext::gm_type::matrix_type const K = gmc.K( 0 );
        typename GMContext::gm_type::matrix_type JB( K/gmc.J( 0 ) );
#if 0
        VLOG(1) << "K= " << gmc.K( 0 ) << "\n";
        VLOG(1) << "B= " << B << "\n";
        VLOG(1) << "J= " << gmc.J( 0 ) << "\n";
        VLOG(1) << "JB= " << JB << "\n";
#endif
        std::fill( phi_t.data(), phi_t.data()+phi_t.num_elements(), value_type( 0 ) );

        if ( do_gradient )
        {
            //VLOG(1) << "compute gradient\n";
            std::fill( g_phi_t.data(), g_phi_t.data()+g_phi_t.num_elements(), value_type( 0 ) );
        }

        if ( do_hessian )
            std::fill( h_phi_t.data(), h_phi_t.data()+h_phi_t.num_elements(), value_type( 0 ) );

        const uint16_type Q = gmc.nPoints();//M_grad.size2();

        // transform
        for ( uint16_type i = 0; i < nLocalDof; ++i )
        {
            for ( uint16_type l = 0; l < nDim; ++l )
            {

                for ( uint16_type p = 0; p < nDim; ++p )
                {
                    for ( uint16_type q = 0; q < Q; ++q )
                    {
                        // \warning : here the transformation depends on the
                        // numbering of the degrees of freedom of the finite
                        // element
                        //phi_t[i][l][0][q] =  pc.phi(i,l,0,q);
                        phi_t[i][l][0][q] += JB( l, p ) * pc.phi( i,p,0,q );
                        //phi_t[i][l][0][q] = gmc.J( 0 ) * B( p, l ) * pc.phi(i,p,0,q);
                        //VLOG(1) << "pc[" << i << "][" << l << "][" << q << "]=" << pc.phi(i,l,0,q) << "\n";
                        //VLOG(1) << "phi_t[" << i << "][" << l << "][" << q << "]=" << phi_t[i][l][0][q] << "\n";
                    }
                }
            }

            //if ( do_gradient )
            {

                for ( uint16_type p = 0; p < nDim; ++p )
                {
                    for ( uint16_type r = 0; r < nDim; ++r )
                    {
                        for ( uint16_type s = 0; s < Q; ++s )
                        {
                            g_phi_t[i][p][r][s] = 0;

                            for ( uint16_type q = 0; q < nDim; ++q )
                            {
                                g_phi_t[i][p][r][s] += JB( p, q ) * pc.grad( i,q,r,s );
                                //g_phi_t[i][p][r][s] = pc.grad(i,p,r,s);

                                //VLOG(1) << "J G[" << i << "][" << q << "][" << r << "][" << s << "=" << JB( p, q ) * pc.grad(i,q,r,s) << "\n";
                            }
                        }

                        //VLOG(1) << "g_phi_t[" << i << "][" << p << "][" << r << "][" << 0 << "=" << g_phi_t[i][p][r][0] << "\n";

                    }
                }
            }

            if ( do_hessian )
            {
            }
        }
    }
#endif

    //@}

protected:
    reference_convex_type M_refconvex;
    mutable ublas::vector<value_type> n{ nDim }; //normal

private:

};

template<uint16_type N,
         typename T = double,
         template<int, int, int> class Convex = Simplex,
         uint16_type TheTAG=0 >
class RaviartThomasDynamicSimplex
    :
public FiniteElement<RaviartThomasDynamicPolynomialSet<N, T, Convex, TheTAG>,
    fem::detail::RaviartThomasDynamicDual,
    PointSetEquiSpaced >,
public HDivPolynomialSet,
public std::enable_shared_from_this<RaviartThomasDynamicSimplex<N,T,Convex,TheTAG> >
{
    using super = FiniteElement<RaviartThomasDynamicPolynomialSet<N, T, Convex, TheTAG>,
            fem::detail::RaviartThomasDynamicDual,
            PointSetEquiSpaced >;

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
    using Pkp1_v_type = typename primal_space_type::Pkp1_v_type;
    using vectorial_polynomialset_type = typename primal_space_type::vectorial_polynomialset_type;
    using polyset_type = typename super::polyset_type;

    static inline const bool is_order_static = false;
    static inline const bool is_order_dynamic = true;
    static constexpr int nOrder_v = Dynamic;

    static inline const bool is_vectorial = polyset_type::is_vectorial;
    static inline const bool is_scalar = polyset_type::is_scalar;
    static inline const uint16_type nComponents = polyset_type::nComponents;
    static inline const bool is_product = false;

    using convex_type = typename dual_space_type::convex_type;
    static_assert( convex_type::is_simplex,
                   "RaviartThomas dynamic finite elements are implemented on simplices only." );
    using pointset_type = typename dual_space_type::pointset_type;
    using reference_convex_type = typename dual_space_type::reference_convex_type;
    using node_type = typename reference_convex_type::node_type;
    using points_type = typename reference_convex_type::points_type;
    using face_type = typename convex_type::topological_face_type;

    static inline const int nOrder = Dynamic;
    static inline const uint16_type nPublicOrder = 0;
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

    RaviartThomasDynamicSimplex()
        :
        RaviartThomasDynamicSimplex( RuntimeOrder{ 0 } )
    {}

    explicit RaviartThomasDynamicSimplex( RuntimeOrder order )
        :
        super( dual_space_type( primal_space_type( order ) ) ),
        M_publicOrder( order.value ),
        M_refconvex()
    {}

    template<int subN>
    struct SubSpace
    {
        using type = RaviartThomasDynamicSimplex<N-1, T, Convex, TheTAG>;
    };

    struct SSpace
    {
        using type = RaviartThomasDynamicSimplex<N, T, Convex, TheTAG>;
    };

    template<uint16_type NewDim>
    struct ChangeDim
    {
        using type = RaviartThomasDynamicSimplex<NewDim, T, Convex, TheTAG>;
    };

    RaviartThomasDynamicSimplex( RaviartThomasDynamicSimplex const& cr )
        :
        super( cr ),
        M_publicOrder( cr.M_publicOrder ),
        M_refconvex()
    {}
    ~RaviartThomasDynamicSimplex() override {}

    reference_convex_type const& referenceConvex() const
    {
        return M_refconvex;
    }

    std::string familyName() const override
    {
        return "raviartthomas";
    }

    [[nodiscard]] uint16_type order() const noexcept
    {
        return M_publicOrder;
    }

    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        return Feel::detail::raviartThomasSimplexInternalOrder( M_publicOrder );
    }

    [[nodiscard]] uint16_type localDof() const noexcept
    {
        return Feel::detail::raviartThomasSimplexTotalDof( nDim, M_publicOrder );
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
        return ( nDim == 2 ) ? Feel::detail::raviartThomasSimplexFacetDof( nDim, M_publicOrder ) : 0;
    }

    [[nodiscard]] uint16_type dofPerFace() const noexcept
    {
        if constexpr ( nDim == 2 )
            return Feel::detail::raviartThomasSimplexInteriorDof( nDim, M_publicOrder );
        else
            return Feel::detail::raviartThomasSimplexFacetDof( nDim, M_publicOrder );
    }

    [[nodiscard]] uint16_type dofPerVolume() const noexcept
    {
        return ( nDim == 3 ) ? Feel::detail::raviartThomasSimplexInteriorDof( nDim, M_publicOrder ) : 0;
    }

    [[nodiscard]] uint16_type localFacetDof() const noexcept
    {
        return Feel::detail::raviartThomasSimplexFacetDof( nDim, M_publicOrder );
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
            const uint16_type momentOrder = static_cast<uint16_type>( internalOrder() - 2 );
            const uint16_type dimPm1 = convex_type::polyDims( momentOrder );
            Pkp1_v_type Pkp1( RuntimeOrder{ internalOrder() } );
            vectorial_polynomialset_type Pkm1 =
                Feel::detail::polynomialSetFrom<vectorial_polynomialset_type>( Pkp1.polynomialsUpToDimension( dimPm1 ),
                                                                               internalOrder() );

            const uint16_type nInternalDof = ( nDim == 2 ) ? dofPerFace() : dofPerVolume();
            IMGeneral<nDim, value_type, Simplex> im( static_cast<uint16_type>( 2*( internalOrder() + 1 ) ) );
            const int firstInternalDof = convex_type::numTopologicalFaces * nFacetDof;
            Feel::detail::interpolateHDivInteriorMomentDofs(
                expr, Ihloc, Pkm1, im, firstInternalDof,
                nInternalDof, nComponents );
        }
    }

    /** Apply all runtime RT dual functionals to domain test-basis rows. */
    template<typename ExprType>
        requires ( std::remove_cvref_t<ExprType>::gmc_type::subEntityCoDim == 0 )
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

    using apply_curl_t = bool_c<true>;
    using apply_id_t = bool_c<false>;

    template<typename ExprType>
    static auto isomorphism( ExprType& expr ) -> decltype( expr )
    {
        return expr;
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
struct RaviartThomasBackend
{
    using type = RaviartThomas<Dim, static_cast<uint16_type>( Order ), T, Convex, Tag>;
};

template<uint16_type Dim, typename T,
         template<int, int, int> class Convex, uint16_type Tag>
struct RaviartThomasBackend<Dynamic, Dim, T, Convex, Tag>
{
    using type = RaviartThomasDynamicSimplex<Dim, T, Convex, Tag>;
};
} // namespace detail

/**
 * @brief Unified public RT finite-element shell for static and runtime order.
 * @tparam Dim reference-cell dimension
 * @tparam Order public static order or `Dynamic`
 * @tparam T scalar type
 * @tparam Convex reference-cell family
 * @tparam Tag compatibility family tag
 */
template<uint16_type Dim, int Order, typename T = double,
         template<int, int, int> class Convex = Simplex, uint16_type Tag = 0>
class RaviartThomasElement :
    public detail::RaviartThomasBackend<Order, Dim, T, Convex, Tag>::type
{
    using super = typename detail::RaviartThomasBackend<Order, Dim, T, Convex, Tag>::type;

public:
    using super::super;
    RaviartThomasElement() = default;

    /** @return public RT family order, distinct from the internal polynomial degree. */
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
         uint16_type TheTAG=0 >
class RaviartThomas : public FiniteElementFamilyDescriptor<Order, TheTAG>
{
public:
    using descriptor_base = FiniteElementFamilyDescriptor<Order, TheTAG>;
    static constexpr bool is_order_static = descriptor_base::is_order_static;
    static constexpr bool is_order_dynamic = descriptor_base::is_order_dynamic;
    static constexpr int nOrder_v = Order;

    template<uint16_type N,
             uint16_type R = N,
             typename T = double,
             typename Convex = Simplex<N> >
    struct apply
    {
        static_assert( Convex::is_simplex,
                       "RaviartThomas hypercube support is not implemented in feelpoly; use simplex RT or add a tensor-product H(div) implementation first." );
        using result_type = fem::RaviartThomasElement<N, Order, T, Simplex, TheTAG>;
        using type = result_type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef RaviartThomas<Order,TheNewTAG> type;
    };

    typedef Lagrange<Order,Scalar> component_basis_type;

};

} // Feel
#endif /* __RaviartThomas_H */
