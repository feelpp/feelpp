/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-02-18

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
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

namespace Feel
{
namespace detail
{
[[nodiscard]] constexpr uint16_type
brezziDouglasMariniSimplexInternalOrder( uint16_type publicOrder ) noexcept
{
    return static_cast<uint16_type>( publicOrder + 1 );
}

[[nodiscard]] constexpr uint16_type
brezziDouglasMariniSimplexFacetDof( uint16_type dim, uint16_type publicOrder ) noexcept
{
    const auto internalOrder = brezziDouglasMariniSimplexInternalOrder( publicOrder );
    return ( dim == 2 )
               ? static_cast<uint16_type>( internalOrder + 1 )
               : static_cast<uint16_type>( ( internalOrder + 1 ) * ( internalOrder + 2 ) / 2 );
}

[[nodiscard]] constexpr uint16_type
brezziDouglasMariniSimplexInteriorDof( uint16_type dim, uint16_type publicOrder ) noexcept
{
    const auto internalOrder = brezziDouglasMariniSimplexInternalOrder( publicOrder );
    if ( internalOrder <= 1 )
        return 0;
    return ( dim == 2 )
               ? static_cast<uint16_type>( internalOrder * internalOrder - 1 )
               : static_cast<uint16_type>( ( internalOrder + 1 ) * ( internalOrder + 2 ) * ( internalOrder - 1 ) / 2 );
}

[[nodiscard]] constexpr uint16_type
brezziDouglasMariniSimplexTotalDof( uint16_type dim, uint16_type publicOrder ) noexcept
{
    const auto internalOrder = brezziDouglasMariniSimplexInternalOrder( publicOrder );
    return ( dim == 2 )
               ? static_cast<uint16_type>( ( internalOrder + 1 ) * ( internalOrder + 2 ) )
               : static_cast<uint16_type>( ( internalOrder + 1 ) * ( internalOrder + 2 ) * ( internalOrder + 3 ) / 2 );
}

template<typename DynamicPolynomialSet, typename SourcePolynomialSet>
[[nodiscard]] DynamicPolynomialSet
brezziDouglasMariniDynamicPolynomialSetFrom( SourcePolynomialSet const& source, uint16_type order )
{
    DynamicPolynomialSet result( source.coeff(), true );
    result.setOrder( order );
    return result;
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
    {
        // BDM_k primal space is the full vector polynomial space (P_k)^d.
        const uint16_type dimPk = convex_type::polyDims( nOrder );
        Pkp1_v_type Pkp1_v;
        vectorial_polynomialset_type Pk_v( Pkp1_v.polynomialsUpToDimension( dimPk ) );
        this->setCoefficient( Pk_v.coeff(), true );
    }
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
        super( RuntimeOrder{ Feel::detail::brezziDouglasMariniSimplexInternalOrder( publicOrder.value ) } )
    {
        build( publicOrder.value );
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
    void build( uint16_type publicOrder )
    {
        M_publicOrder = publicOrder;
        const uint16_type internalOrder = Feel::detail::brezziDouglasMariniSimplexInternalOrder( publicOrder );
        const uint16_type dimPk = convex_type::polyDims( internalOrder );

        Pkp1_v_type Pkp1_v( RuntimeOrder{ internalOrder } );
        vectorial_polynomialset_type Pk_v =
            Feel::detail::brezziDouglasMariniDynamicPolynomialSetFrom<vectorial_polynomialset_type>( Pkp1_v.polynomialsUpToDimension( dimPk ),
                                                                                                     internalOrder );
        CHECK( Pk_v.polynomialDimension() == Feel::detail::brezziDouglasMariniSimplexTotalDof( N, publicOrder ) )
            << "Invalid dynamic BDM polynomial dimension: got " << Pk_v.polynomialDimension()
            << " expected " << Feel::detail::brezziDouglasMariniSimplexTotalDof( N, publicOrder );
        this->setCoefficient( Pk_v.coeff(), true );
    }

private:
    uint16_type M_publicOrder = 0;
};

template<uint16_type N,
         typename T = double,
         uint16_type TheTAG = 0>
class BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet;

template<typename T,
         uint16_type TheTAG>
class BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet<2, T, TheTAG>
    :
    public Feel::detail::OrthonormalPolynomialSet<2, Dynamic, 2, Vectorial, T, TheTAG, Simplex>
{
    static constexpr uint16_type N = 2;
    using super = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Simplex>;

public:
    using value_type = typename super::value_type;
    using convex_type = typename super::convex_type;
    using matrix_type = typename super::matrix_type;
    using points_type = typename super::points_type;
    using Pkp1_v_type = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Simplex>;
    using Pkp1_s_type = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Scalar, T, TheTAG, Simplex>;
    using vectorial_polynomialset_type = PolynomialSet<typename super::basis_type, Vectorial, Dynamic>;
    using scalar_polynomialset_type = PolynomialSet<typename super::basis_type, Scalar, Dynamic>;
    using scalar_polynomial_type = typename scalar_polynomialset_type::polynomial_type;

    static inline const uint16_type nDim = super::nDim;
    static inline const int nOrder = Dynamic;
    static inline const uint16_type nComponents = super::nComponents;
    static inline const bool is_product = false;

    BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet()
        :
        BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet( RuntimeOrder{ 0 } )
    {}

    explicit BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet( RuntimeOrder publicOrder )
        :
        super( RuntimeOrder{ static_cast<uint16_type>( publicOrder.value + 1 ) } )
    {
        build( publicOrder.value );
    }

private:
    void build( uint16_type publicOrder )
    {
        const uint16_type internalOrder = static_cast<uint16_type>( publicOrder + 1 );
        const uint16_type dimPk = convex_type::polyDims( publicOrder );
        const uint16_type dimPkm1 = ( publicOrder == 0 ) ? 0 : convex_type::polyDims( publicOrder - 1 );

        Pkp1_v_type Pkp1_v( RuntimeOrder{ internalOrder } );
        vectorial_polynomialset_type Pk_v =
            Feel::detail::brezziDouglasMariniDynamicPolynomialSetFrom<vectorial_polynomialset_type>( Pkp1_v.polynomialsUpToDimension( dimPk ),
                                                                                                     internalOrder );

        Pkp1_s_type Pkp1( RuntimeOrder{ internalOrder } );
        scalar_polynomialset_type Pk =
            Feel::detail::brezziDouglasMariniDynamicPolynomialSetFrom<scalar_polynomialset_type>( Pkp1.polynomialsUpToDimension( dimPk ),
                                                                                                  internalOrder );

        IMGeneral<convex_type::nDim, value_type> im( static_cast<uint16_type>( 2*internalOrder ) );
        ublas::matrix<value_type> xPkc( nComponents*( dimPk - dimPkm1 ), Pk.coeff().size2() );
        for ( int l = dimPkm1, i = 0; l < dimPk; ++l, ++i )
        {
            for ( int j = 0; j < convex_type::nDim; ++j )
            {
                Feel::detail::times_rotx<scalar_polynomial_type> xp( Pk.polynomial( l ), j );
                ublas::row( xPkc, i*nComponents + j ) =
                    ublas::row( Feel::project( Pkp1, xp, im ).coeff(), 0 );
            }
        }

        vectorial_polynomialset_type xPk( typename super::basis_type(), xPkc, true );
        xPk.setOrder( internalOrder );
        auto nedSpace = unite( Pk_v, xPk );
        const uint16_type expectedDim = static_cast<uint16_type>( ( publicOrder + 1 ) * ( publicOrder + 3 ) );
        CHECK_GE( nedSpace.polynomialDimension(), expectedDim )
            << "Invalid dynamic 2D Nedelec first-kind dimension: got " << nedSpace.polynomialDimension()
            << " expected at least " << expectedDim;
        if ( nedSpace.polynomialDimension() > expectedDim )
        {
            matrix_type coeff( nComponents*expectedDim, nedSpace.coeff().size2() );
            ublas::subrange( coeff, 0, coeff.size1(), 0, coeff.size2() ) =
                ublas::subrange( nedSpace.coeff(), 0, coeff.size1(), 0, coeff.size2() );
            nedSpace = vectorial_polynomialset_type( typename super::basis_type(), coeff, true );
            nedSpace.setOrder( internalOrder );
        }
        this->setCoefficient( nedSpace.coeff(), true );
    }
};

template<typename T,
         uint16_type TheTAG>
class BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet<3, T, TheTAG>
    :
    public Feel::detail::OrthonormalPolynomialSet<3, Dynamic, 3, Vectorial, T, TheTAG, Simplex>
{
    static constexpr uint16_type N = 3;
    using super = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Simplex>;

public:
    using value_type = typename super::value_type;
    using convex_type = typename super::convex_type;
    using matrix_type = typename super::matrix_type;
    using points_type = typename super::points_type;
    using Pkp1_v_type = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Vectorial, T, TheTAG, Simplex>;
    using Pkp1_s_type = Feel::detail::OrthonormalPolynomialSet<N, Dynamic, N, Scalar, T, TheTAG, Simplex>;
    using vectorial_polynomialset_type = PolynomialSet<typename super::basis_type, Vectorial, Dynamic>;
    using scalar_polynomialset_type = PolynomialSet<typename super::basis_type, Scalar, Dynamic>;

    static inline const uint16_type nDim = super::nDim;
    static inline const int nOrder = Dynamic;
    static inline const uint16_type nComponents = super::nComponents;
    static inline const bool is_product = false;

    BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet()
        :
        BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet( RuntimeOrder{ 0 } )
    {}

    explicit BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet( RuntimeOrder publicOrder )
        :
        super( RuntimeOrder{ static_cast<uint16_type>( publicOrder.value + 1 ) } )
    {
        build( publicOrder.value );
    }

private:
    void build( uint16_type publicOrder )
    {
        const uint16_type internalOrder = static_cast<uint16_type>( publicOrder + 1 );
        const uint16_type dimPk = convex_type::polyDims( publicOrder );
        const uint16_type dimPkm1 = ( publicOrder == 0 ) ? 0 : convex_type::polyDims( publicOrder - 1 );

        Pkp1_v_type Pkp1_v( RuntimeOrder{ internalOrder } );
        vectorial_polynomialset_type Pk_v =
            Feel::detail::brezziDouglasMariniDynamicPolynomialSetFrom<vectorial_polynomialset_type>( Pkp1_v.polynomialsUpToDimension( dimPk ),
                                                                                                     internalOrder );
        vectorial_polynomialset_type Pke_v =
            Feel::detail::brezziDouglasMariniDynamicPolynomialSetFrom<vectorial_polynomialset_type>( Pkp1_v.polynomialsUpToDimension( dimPk ),
                                                                                                     internalOrder );

        Pkp1_s_type Pkp1( RuntimeOrder{ internalOrder } );

        IMGeneral<convex_type::nDim, value_type> im( static_cast<uint16_type>( 2*internalOrder ) );
        ublas::matrix<value_type> xPkcV( nComponents*nComponents*dimPk, Pke_v.coeff().size2() );
        Eigen::Map<Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> xPkc( xPkcV.data().begin(),
                                                                                                  xPkcV.size1(), xPkcV.size2() );
        xPkc.setZero();

        auto ePkp1 = Pkp1.evaluate( im.points() );
        auto ePkv = Pk_v.evaluate( im.points() );
        auto const& X = im.points();
        auto const& W = im.weights();

        for ( int c = 0; c < convex_type::nDim; ++c )
        {
            for ( int l = c*convex_type::nDim + dimPkm1; l < c*convex_type::nDim + dimPk; ++l )
            {
                for ( int j = 0; j < convex_type::nDim; ++j )
                {
                    for ( int q = 0; q < im.nPoints(); ++q )
                    {
                        auto b = W( q )*( X( ( j + 2 )%3, q )*ePkv( l + ( j + 1 )%3, q ) -
                                           X( ( j + 1 )%3, q )*ePkv( l + ( j + 2 )%3, q ) );
                        for ( int p = 0; p < xPkc.cols(); ++p )
                            xPkc( l + j, p ) += b * ePkp1( p, q );
                    }
                }
            }
        }

        vectorial_polynomialset_type xPk( typename super::basis_type(), xPkcV, true );
        xPk.setOrder( internalOrder );
        auto nedSpace = unite( Pke_v, xPk );
        const uint16_type expectedDim = static_cast<uint16_type>( ( publicOrder + 1 ) * ( publicOrder + 3 ) * ( publicOrder + 4 ) / 2 );
        CHECK_GE( nedSpace.polynomialDimension(), expectedDim )
            << "Invalid dynamic 3D Nedelec first-kind dimension: got " << nedSpace.polynomialDimension()
            << " expected at least " << expectedDim;
        if ( nedSpace.polynomialDimension() > expectedDim )
        {
            matrix_type coeff( nComponents*expectedDim, nedSpace.coeff().size2() );
            ublas::subrange( coeff, 0, coeff.size1(), 0, coeff.size2() ) =
                ublas::subrange( nedSpace.coeff(), 0, coeff.size1(), 0, coeff.size2() );
            nedSpace = vectorial_polynomialset_type( typename super::basis_type(), coeff, true );
            nedSpace.setOrder( internalOrder );
        }
        this->setCoefficient( nedSpace.coeff(), true );
    }
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
    using interior_polyset_type = BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet<nDim, value_type>;

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
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( nDim, numPoints ),
        M_pts_per_face( convex_type::numTopologicalFaces ),
        M_fset( primal )
    {
        for ( int p = 0, e = M_convex_ref.entityRange( nDim-1 ).begin();
              e < M_convex_ref.entityRange( nDim-1 ).end();
              ++e )
        {
            points_type Gt( M_convex_ref.makePoints( nDim-1, e ) );
            M_pts_per_face[e] = Gt;

            if ( Gt.size2() )
            {
                ublas::subrange( M_pts, 0, nDim, p, p + Gt.size2() ) = Gt;
                p += Gt.size2();
            }
        }

        typedef Functional<primal_space_type> functional_type;
        std::vector<functional_type> fset;

        Feel::detail::appendHDivFacetNormalPointFunctionals( primal, M_convex_ref, M_pts_per_face, fset );

        if constexpr ( nOrder > 1 )
        {
            interior_polyset_type interiorSet( RuntimeOrder{ static_cast<uint16_type>( nOrder - 2 ) } );
            const uint16_type nInteriorDof = ( nDim == 2 ) ? nDofPerFace : nDofPerVolume;
            CHECK( interiorSet.polynomialDimension() == nInteriorDof )
                << "Invalid BDM interior moments: got " << interiorSet.polynomialDimension()
                << " expected " << nInteriorDof;

            IM<nDim, 2*( nOrder + 1 ), value_type, Simplex> im;
            const uint16_type firstInternalDof = ( nDim == 2 )
                                                 ? reference_convex_type::numEdges * nDofPerEdge
                                                 : reference_convex_type::numTopologicalFaces * nDofPerFace;
            const uint16_type requiredPointCount = static_cast<uint16_type>( firstInternalDof + im.nPoints() );
            if ( M_pts.size2() < requiredPointCount )
                M_pts.resize( nDim, requiredPointCount, true );
            ublas::subrange( M_pts, 0, nDim, firstInternalDof, requiredPointCount ) = im.points();

            for ( int i = 0; i < interiorSet.polynomialDimension(); ++i )
            {
                fset.push_back( Feel::detail::makeHDivOrthonormalIntegralMomentFunctional(
                    primal, interiorSet.polynomial( i ) ) );
            }
        }

        M_fset.setFunctionalSet( fset );
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
    std::vector<std::vector<uint16_type>> M_eid;
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
    using interior_polyset_type = BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet<nDim, value_type>;

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
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( nDim, Feel::detail::brezziDouglasMariniSimplexTotalDof( nDim, M_publicOrder ) ),
        M_pts_per_face( convex_type::numTopologicalFaces ),
        M_fset( primal )
    {
        const uint16_type nFacetDof = runtimeFacetDof();
        pointset_type facetPoints( RuntimeOrder{ static_cast<uint16_type>( nDim + M_internalOrder ) } );

        for ( int p = 0, e = M_convex_ref.entityRange( nDim-1 ).begin();
              e < M_convex_ref.entityRange( nDim-1 ).end();
              ++e )
        {
            points_type Gt( facetPoints.pointsBySubEntity( nDim-1, e ) );
            CHECK( Gt.size2() == nFacetDof )
                << "Invalid BDM dynamic facet point count on facet " << e
                << ": got " << Gt.size2() << " expected " << nFacetDof;
            M_pts_per_face[e] = Gt;

            if ( Gt.size2() )
            {
                ublas::subrange( M_pts, 0, nDim, p, p+Gt.size2() ) = Gt;
                p += Gt.size2();
            }
        }

        using functional_type = Functional<primal_space_type>;
        std::vector<functional_type> fset;

        Feel::detail::appendHDivFacetNormalPointFunctionals( primal, M_convex_ref, M_pts_per_face, fset );

        if ( M_internalOrder > 1 )
        {
            interior_polyset_type interiorSet( RuntimeOrder{ static_cast<uint16_type>( M_internalOrder - 2 ) } );
            const uint16_type nInteriorDof = runtimeInteriorDof();
            CHECK( interiorSet.polynomialDimension() == nInteriorDof )
                << "Invalid BDM dynamic interior moments: got " << interiorSet.polynomialDimension()
                << " expected " << nInteriorDof;

            IMGeneral<nDim, value_type, Simplex> im( static_cast<uint16_type>( 2*( M_internalOrder + 1 ) ) );
            const uint16_type firstInternalDof = static_cast<uint16_type>( convex_type::numTopologicalFaces * nFacetDof );
            const uint16_type requiredPointCount = static_cast<uint16_type>( firstInternalDof + im.nPoints() );
            if ( M_pts.size2() < requiredPointCount )
                M_pts.resize( nDim, requiredPointCount, true );
            ublas::subrange( M_pts, 0, nDim, firstInternalDof, requiredPointCount ) = im.points();

            for ( int i = 0; i < interiorSet.polynomialDimension(); ++i )
                fset.push_back( Feel::detail::makeHDivOrthonormalIntegralMomentFunctional(
                    primal, interiorSet.polynomial( i ) ) );
        }

        CHECK( fset.size() == Feel::detail::brezziDouglasMariniSimplexTotalDof( nDim, M_publicOrder ) )
            << "Invalid BDM dynamic functional count: got " << fset.size()
            << " expected " << Feel::detail::brezziDouglasMariniSimplexTotalDof( nDim, M_publicOrder );
        M_fset.setFunctionalSet( fset );
    }

    [[nodiscard]] uint16_type publicOrder() const noexcept
    {
        return M_publicOrder;
    }

    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        return M_internalOrder;
    }

    [[nodiscard]] uint16_type runtimeFacetDof() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexFacetDof( nDim, M_publicOrder );
    }

    [[nodiscard]] uint16_type runtimeInteriorDof() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexInteriorDof( nDim, M_publicOrder );
    }

    [[nodiscard]] uint16_type runtimeLocalDof() const noexcept
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
    std::vector<std::vector<uint16_type>> M_eid;
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
        const uint16_type parentLocalDofId = this->dofParent( localDofId );

        const uint16_type nV = static_cast<uint16_type>( reference_convex_type::numVertices * nDofPerVertex );
        const uint16_type nE = static_cast<uint16_type>( reference_convex_type::numEdges * nDofPerEdge );
        const uint16_type nF = static_cast<uint16_type>( reference_convex_type::numFaces * nDofPerFace );

        if constexpr ( nDofPerVertex > 0 )
        {
            if ( parentLocalDofId < nV )
            {
                return typename super::DofAttachment{
                    .entityDim = 0,
                    .entityId = static_cast<uint16_type>( parentLocalDofId / nDofPerVertex ),
                    .ordinal = static_cast<uint16_type>( parentLocalDofId % nDofPerVertex ),
                    .kind = this->dofType( localDofId ) };
            }
        }

        const uint16_type parentAfterVertex = static_cast<uint16_type>( parentLocalDofId - nV );
        if constexpr ( nDofPerEdge > 0 )
        {
            if ( parentAfterVertex < nE )
            {
                return typename super::DofAttachment{
                    .entityDim = 1,
                    .entityId = static_cast<uint16_type>( parentAfterVertex / nDofPerEdge ),
                    .ordinal = static_cast<uint16_type>( parentAfterVertex % nDofPerEdge ),
                    .kind = this->dofType( localDofId ) };
            }
        }

        const uint16_type parentAfterEdge = static_cast<uint16_type>( parentAfterVertex - nE );
        if constexpr ( nDofPerFace > 0 )
        {
            if ( parentAfterEdge < nF )
            {
                return typename super::DofAttachment{
                    .entityDim = 2,
                    .entityId = static_cast<uint16_type>( parentAfterEdge / nDofPerFace ),
                    .ordinal = static_cast<uint16_type>( parentAfterEdge % nDofPerFace ),
                    .kind = this->dofType( localDofId ) };
            }
        }

        if constexpr ( nDofPerVolume > 0 )
        {
            const uint16_type parentAfterFace = static_cast<uint16_type>( parentAfterEdge - nF );
            return typename super::DofAttachment{
                .entityDim = 3,
                .entityId = 0,
                .ordinal = static_cast<uint16_type>( parentAfterFace % nDofPerVolume ),
                .kind = this->dofType( localDofId ) };
        }

        return typename super::DofAttachment{
            .entityDim = -1,
            .entityId = super::DofAttachment::invalid_id,
            .ordinal = super::DofAttachment::invalid_id,
            .kind = this->dofType( localDofId ) };
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
        auto g = expr.geom();

        for ( int f = 0; f < convex_type::numTopologicalFaces; ++f )
        {
            if ( g->faceId() == invalid_uint16_type_value )
                expr.geom()->faceNormal( f, n, true );
            else
                expr.geom()->faceNormal( g->faceId(), n, true );

            const uint16_type nFacetDof = ( nDim == 2 ) ? nDofPerEdge : nDofPerFace;
            for ( int l = 0; l < nFacetDof; ++l )
            {
                const int q = f * nFacetDof + l;
                for ( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                    Ihloc( q ) += expr.evalq( c1, 0, q ) * n( c1 );
            }
        }

        if constexpr ( nOrder > 1 )
        {
            static const uint16_type interiorOrder = nOrder - 2;
            typedef BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet<nDim, value_type, TheTAG> interior_polyset_type;
            interior_polyset_type interiorSpace( RuntimeOrder{ interiorOrder } );

            IM<nDim, 2*( nOrder + 1 ), value_type, Simplex> im;
            auto interiorAtQuadPts = interiorSpace.evaluate( im.points() );

            const uint16_type nInternalDof = ( nDim == 2 ) ? nDofPerFace : nDofPerVolume;
            CHECK( interiorSpace.polynomialDimension() == nInternalDof )
                << "Invalid BDM interior interpolation dimension: got "
                << interiorSpace.polynomialDimension() << " expected " << nInternalDof;

            const int firstInternalDof = ( nDim == 2 )
                                         ? convex_type::numTopologicalFaces * nDofPerEdge
                                         : convex_type::numTopologicalFaces * nDofPerFace;
            if constexpr ( requires { expr.nPoints(); } )
            {
                CHECK( expr.nPoints() >= firstInternalDof + im.nPoints() )
                    << "BDM interior moment interpolation requires expression values at appended quadrature points: got "
                    << expr.nPoints() << " points, need " << firstInternalDof + im.nPoints();
            }

            for ( int l = 0; l < nInternalDof; ++l )
            {
                const int dof = firstInternalDof + l;
            for ( int q = 0; q < im.nPoints(); ++q )
            {
                const int exprPoint = firstInternalDof + q;
                const value_type scal = Feel::detail::hdivInteriorMomentIntegrand(
                    expr, interiorAtQuadPts, l, q, exprPoint, nComponents );
                Ihloc( dof ) += im.weight( q ) * scal;
            }
        }
        }
    }

    local_interpolant_type faceLocalInterpolant() const
    {
        return local_interpolant_type::Zero( nLocalFaceDof, 1 );
    }

    template<typename ExprType>
    void faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
    {
        auto g = expr.geom();
        Ihloc.setZero();

        int f = 0;
        {
            if ( g->faceId() == invalid_uint16_type_value )
                expr.geom()->faceNormal( f, n, true );
            else
                expr.geom()->faceNormal( g->faceId(), n, true );

            auto nLocalDof = ( nDim == 2 ) ? nDofPerEdge : nDofPerFace;
            for ( int l = 0; l < nLocalDof; ++l )
            {
                int q = ( nDim == 2 ) ? f*nDofPerEdge + l : f*nDofPerFace + l;
                for ( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                    Ihloc( q ) += expr.evalq( c1, 0, q ) * n( c1 );
            }
        }
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
    using interior_polyset_type = BrezziDouglasMariniNedelecFirstKindDynamicPolynomialSet<nDim, value_type, TheTAG>;

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

    [[nodiscard]] uint16_type runtimeOrder() const noexcept
    {
        return order();
    }

    [[nodiscard]] uint16_type internalOrder() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexInternalOrder( M_publicOrder );
    }

    [[nodiscard]] uint16_type localDof() const noexcept
    {
        return runtimeLocalDof();
    }

    [[nodiscard]] uint16_type runtimeLocalDof() const noexcept
    {
        return Feel::detail::brezziDouglasMariniSimplexTotalDof( nDim, M_publicOrder );
    }

    uint16_type localDofPerComponent() const override
    {
        return runtimeLocalDof();
    }

    uint16_type localDofCount( bool perComponent = false ) const override
    {
        (void)perComponent;
        return runtimeLocalDof();
    }

    [[nodiscard]] uint16_type dofPerVertex() const noexcept
    {
        return 0;
    }

    [[nodiscard]] uint16_type runtimeDofPerVertex() const noexcept
    {
        return dofPerVertex();
    }

    [[nodiscard]] uint16_type dofPerEdge() const noexcept
    {
        return ( nDim == 2 ) ? Feel::detail::brezziDouglasMariniSimplexFacetDof( nDim, M_publicOrder ) : 0;
    }

    [[nodiscard]] uint16_type runtimeDofPerEdge() const noexcept
    {
        return dofPerEdge();
    }

    [[nodiscard]] uint16_type dofPerFace() const noexcept
    {
        if constexpr ( nDim == 2 )
            return Feel::detail::brezziDouglasMariniSimplexInteriorDof( nDim, M_publicOrder );
        else
            return Feel::detail::brezziDouglasMariniSimplexFacetDof( nDim, M_publicOrder );
    }

    [[nodiscard]] uint16_type runtimeDofPerFace() const noexcept
    {
        return dofPerFace();
    }

    [[nodiscard]] uint16_type dofPerVolume() const noexcept
    {
        return ( nDim == 3 ) ? Feel::detail::brezziDouglasMariniSimplexInteriorDof( nDim, M_publicOrder ) : 0;
    }

    [[nodiscard]] uint16_type runtimeDofPerVolume() const noexcept
    {
        return dofPerVolume();
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
        const uint16_type parentLocalDofId = this->dofParent( localDofId );

        const uint16_type nV = static_cast<uint16_type>( reference_convex_type::numVertices * dofPerVertex() );
        const uint16_type nE = static_cast<uint16_type>( reference_convex_type::numEdges * dofPerEdge() );
        const uint16_type nF = static_cast<uint16_type>( reference_convex_type::numFaces * dofPerFace() );

        if ( parentLocalDofId < nV && dofPerVertex() > 0 )
        {
            return typename super::DofAttachment{
                .entityDim = 0,
                .entityId = static_cast<uint16_type>( parentLocalDofId / dofPerVertex() ),
                .ordinal = static_cast<uint16_type>( parentLocalDofId % dofPerVertex() ),
                .kind = this->dofType( localDofId ) };
        }

        const uint16_type parentAfterVertex = static_cast<uint16_type>( parentLocalDofId - nV );
        if ( parentAfterVertex < nE && dofPerEdge() > 0 )
        {
            return typename super::DofAttachment{
                .entityDim = 1,
                .entityId = static_cast<uint16_type>( parentAfterVertex / dofPerEdge() ),
                .ordinal = static_cast<uint16_type>( parentAfterVertex % dofPerEdge() ),
                .kind = this->dofType( localDofId ) };
        }

        const uint16_type parentAfterEdge = static_cast<uint16_type>( parentAfterVertex - nE );
        if ( parentAfterEdge < nF && dofPerFace() > 0 )
        {
            return typename super::DofAttachment{
                .entityDim = 2,
                .entityId = static_cast<uint16_type>( parentAfterEdge / dofPerFace() ),
                .ordinal = static_cast<uint16_type>( parentAfterEdge % dofPerFace() ),
                .kind = this->dofType( localDofId ) };
        }

        if ( dofPerVolume() > 0 )
        {
            const uint16_type parentAfterFace = static_cast<uint16_type>( parentAfterEdge - nF );
            return typename super::DofAttachment{
                .entityDim = 3,
                .entityId = 0,
                .ordinal = static_cast<uint16_type>( parentAfterFace % dofPerVolume() ),
                .kind = this->dofType( localDofId ) };
        }

        return typename super::DofAttachment{
            .entityDim = -1,
            .entityId = super::DofAttachment::invalid_id,
            .ordinal = super::DofAttachment::invalid_id,
            .kind = this->dofType( localDofId ) };
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
        return local_interpolant_type::Zero( n*runtimeLocalDof() );
    }

    using local_interpolants_type = Eigen::MatrixXd;
    local_interpolants_type localInterpolants( int p, int n = 1 ) const
    {
        return local_interpolants_type::Zero( n*runtimeLocalDof(), p );
    }

    template<typename ExprType>
    void interpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
    {
        Ihloc.setZero();
        auto g = expr.geom();

        const uint16_type nFacetDof = localFacetDof();
        for ( int f = 0; f < convex_type::numTopologicalFaces; ++f )
        {
            if ( g->faceId() == invalid_uint16_type_value )
                expr.geom()->faceNormal( f, n, true );
            else
                expr.geom()->faceNormal( g->faceId(), n, true );

            for ( int l = 0; l < nFacetDof; ++l )
            {
                const int q = f * nFacetDof + l;
                for ( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                    Ihloc( q ) += expr.evalq( c1, 0, q ) * n( c1 );
            }
        }

        if ( internalOrder() > 1 )
        {
            interior_polyset_type interiorSpace( RuntimeOrder{ static_cast<uint16_type>( internalOrder() - 2 ) } );
            IMGeneral<nDim, value_type, Simplex> im( static_cast<uint16_type>( 2*( internalOrder() + 1 ) ) );
            auto interiorAtQuadPts = interiorSpace.evaluate( im.points() );

            const uint16_type nInternalDof = ( nDim == 2 ) ? dofPerFace() : dofPerVolume();
            CHECK( interiorSpace.polynomialDimension() == nInternalDof )
                << "Invalid BDM dynamic interior interpolation dimension: got "
                << interiorSpace.polynomialDimension() << " expected " << nInternalDof;

            const int firstInternalDof = convex_type::numTopologicalFaces * nFacetDof;
            if constexpr ( requires { expr.nPoints(); } )
            {
                CHECK( expr.nPoints() >= firstInternalDof + im.nPoints() )
                    << "BDM dynamic interior moment interpolation requires expression values at appended quadrature points: got "
                    << expr.nPoints() << " points, need " << firstInternalDof + im.nPoints();
            }

            for ( int l = 0; l < nInternalDof; ++l )
            {
                const int dof = firstInternalDof + l;
                for ( int q = 0; q < im.nPoints(); ++q )
                {
                    const int exprPoint = firstInternalDof + q;
                    const value_type scal = Feel::detail::hdivInteriorMomentIntegrand(
                        expr, interiorAtQuadPts, l, q, exprPoint, nComponents );
                    Ihloc( dof ) += im.weight( q ) * scal;
                }
            }
        }
    }

    local_interpolant_type faceLocalInterpolant() const
    {
        return local_interpolant_type::Zero( localFacetDof(), 1 );
    }

    template<typename ExprType>
    void faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
    {
        auto g = expr.geom();
        Ihloc.setZero();

        int f = 0;
        if ( g->faceId() == invalid_uint16_type_value )
            expr.geom()->faceNormal( f, n, true );
        else
            expr.geom()->faceNormal( g->faceId(), n, true );

        auto nLocalDof = localFacetDof();
        for ( int l = 0; l < nLocalDof; ++l )
        {
            int q = f*nLocalDof + l;
            for ( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                Ihloc( q ) += expr.evalq( c1, 0, q ) * n( c1 );
        }
    }

protected:
    uint16_type M_publicOrder = 0;
    reference_convex_type M_refconvex;
    mutable ublas::vector<value_type> n{ nDim };
};

} // namespace fem

template<int Order,
         uint16_type TheTAG = 0>
class BrezziDouglasMarini
{
public:
    static constexpr bool is_order_static = ( Order != Dynamic );
    static constexpr bool is_order_dynamic = !is_order_static;
    static constexpr int nOrder_v = Order;

    template<uint16_type N,
             uint16_type R = N,
             typename T = double,
             typename Convex = Simplex<N>>
    struct apply
    {
        static_assert( is_order_static,
                       "BrezziDouglasMarini<Dynamic> needs a dedicated FE-level runtime dispatch path; the static BDM factory cannot accept Dynamic as a uint16_type order." );
        static_assert( Convex::is_simplex,
                       "BrezziDouglasMarini hypercube support is not implemented in feelpoly; use simplex BDM or add a tensor-product H(div) implementation first." );
        using result_type = fem::BrezziDouglasMarini<N, static_cast<uint16_type>( Order ), T, Simplex, TheTAG>;
        using type = result_type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef BrezziDouglasMarini<Order, TheNewTAG> type;
    };

    typedef Lagrange<Order, Scalar> component_basis_type;

    static inline const uint16_type nOrder = is_order_static ? static_cast<uint16_type>( Order ) : 0;
    static const uint16_type TAG = TheTAG;
};

template<uint16_type TheTAG>
class BrezziDouglasMarini<Dynamic, TheTAG>
{
public:
    static constexpr bool is_order_static = false;
    static constexpr bool is_order_dynamic = true;
    static constexpr int nOrder_v = Dynamic;

    template<uint16_type N,
             uint16_type R = N,
             typename T = double,
             typename Convex = Simplex<N>>
    struct apply
    {
        static_assert( Convex::is_simplex,
                       "BrezziDouglasMarini hypercube support is not implemented in feelpoly; use simplex BDM or add a tensor-product H(div) implementation first." );
        using result_type = fem::BrezziDouglasMariniDynamicSimplex<N, T, Simplex, TheTAG>;
        using type = result_type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef BrezziDouglasMarini<Dynamic, TheNewTAG> type;
    };

    typedef Lagrange<Dynamic, Scalar> component_basis_type;

    static inline const uint16_type nOrder = 0;
    static const uint16_type TAG = TheTAG;
};

} // namespace Feel
#endif /* __BrezziDouglasMarini_H */
