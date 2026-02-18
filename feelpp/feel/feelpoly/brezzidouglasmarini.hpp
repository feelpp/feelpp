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
#include <feel/feelpoly/functionals.hpp>
#include <feel/feelpoly/functionals2.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>
#include <feel/feelpoly/fe.hpp>
#include <feel/feelpoly/hdivpolynomialset.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelpoly/meta.hpp>

namespace Feel
{
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
    using interior_polyset_type = NedelecPolynomialSet<nDim, (nOrder > 1 ? nOrder-2 : 0), NedelecKind::NED1, value_type>;

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

        std::array<value_type, convex_type::numTopologicalFaces> jacobianScaling{};
        if constexpr ( nDim == 2 )
            jacobianScaling = { value_type(2.8284271247461903), value_type(2.0), value_type(2.0) };
        else if constexpr ( nDim == 3 )
            jacobianScaling = { value_type(3.464101615137754), value_type(2.0), value_type(2.0), value_type(2.0) };

        for ( int e = M_convex_ref.entityRange( nDim-1 ).begin();
              e < M_convex_ref.entityRange( nDim-1 ).end();
              ++e )
        {
            typedef Feel::functional::DirectionalComponentPointsEvaluation<primal_space_type> dcpe_type;
            node_type dir( nDim );
            em_node_type<value_type> edir( dir.data().begin(), dir.size() );
            if constexpr ( ( nDim == 2 || nDim == 3 ) && convex_type::is_simplex )
                edir = M_convex_ref.normal( e ) * jacobianScaling[e];
            else
                edir = M_convex_ref.normal( e );
            dcpe_type __dcpe( primal, dir, M_pts_per_face[e] );
            std::copy( __dcpe.begin(), __dcpe.end(), std::back_inserter( fset ) );
        }

        if constexpr ( nOrder > 1 )
        {
            interior_polyset_type interiorSet;
            const uint16_type nInteriorDof = ( nDim == 2 ) ? nDofPerFace : nDofPerVolume;
            CHECK( interiorSet.polynomialDimension() == nInteriorDof )
                << "Invalid BDM interior moments: got " << interiorSet.polynomialDimension()
                << " expected " << nInteriorDof;
            for ( int i = 0; i < interiorSet.polynomialDimension(); ++i )
            {
                typedef functional::IntegralMoment<primal_space_type, interior_polyset_type> fim_type;
                fset.push_back( fim_type( primal, interiorSet.polynomial( i ) ) );
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
            typedef NedelecPolynomialSet<nDim, interiorOrder, NedelecKind::NED1, value_type, TheTAG> interior_polyset_type;
            interior_polyset_type interiorSpace;

            IM<nDim, 2*( nOrder + 1 ), value_type, Simplex> im;
            auto interiorAtQuadPts = interiorSpace.evaluate( im.points() );

            const uint16_type nInternalDof = ( nDim == 2 ) ? nDofPerFace : nDofPerVolume;
            CHECK( interiorSpace.polynomialDimension() == nInternalDof )
                << "Invalid BDM interior interpolation dimension: got "
                << interiorSpace.polynomialDimension() << " expected " << nInternalDof;

            const int firstInternalDof = ( nDim == 2 )
                                         ? convex_type::numTopologicalFaces * nDofPerEdge
                                         : convex_type::numTopologicalFaces * nDofPerFace;

            for ( int l = 0; l < nInternalDof; ++l )
            {
                const int dof = firstInternalDof + l;
                for ( int q = 0; q < im.nPoints(); ++q )
                {
                    value_type scal = 0.;
                    for ( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                        scal += expr.evalq( c1, 0, q ) * interiorAtQuadPts( nComponents * l + c1, q );
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

} // namespace fem

template<int Order,
         uint16_type TheTAG = 0>
class BrezziDouglasMarini
{
public:
    template<uint16_type N,
             uint16_type R = N,
             typename T = double,
             typename Convex = Simplex<N>>
    struct apply
    {
        using result_type = if_t<Convex::is_simplex,
                                 fem::BrezziDouglasMarini<N, Order, T, Simplex, TheTAG>,
                                 fem::BrezziDouglasMarini<N, Order, T, Hypercube, TheTAG>>;
        using type = result_type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef BrezziDouglasMarini<Order, TheNewTAG> type;
    };

    typedef Lagrange<Order, Scalar> component_basis_type;

    static const uint16_type TAG = TheTAG;
};

} // namespace Feel
#endif /* __BrezziDouglasMarini_H */
