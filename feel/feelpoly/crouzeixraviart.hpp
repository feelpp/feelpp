/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-01-12

  Copyright (C) 2006 EPFL
  Copyright (C) 2011 Universite Joseph Fourier

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
   \file crouzeixraviart.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-01-12
 */
#ifndef __CrouzeixRaviart_H
#define __CrouzeixRaviart_H 1

#include <boost/ptr_container/ptr_vector.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/lu.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/pointset.hpp>

#include <feel/feelpoly/fe.hpp>
#include <feel/feelpoly/dualbasis.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/functionalset.hpp>
#include <feel/feelpoly/moment.hpp>

namespace Feel
{
namespace fem
{
namespace detail
{
template<uint16_type N,
         template<uint16_type Dim> class PolySetType = Scalar,
         typename T = double>
class RannacherTurekPolynomialSet
    :
public MomentPolynomialSet<N, 2, N, PolySetType, T, Hypercube>
{
    typedef MomentPolynomialSet<N, 2, N, PolySetType, T, Hypercube> super;

public:

    typedef RannacherTurekPolynomialSet<N, PolySetType, T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::convex_type convex_type;
    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nComponents = super::nComponents;
    static const bool is_product = true;

    RannacherTurekPolynomialSet()
        :
        super()
    {
        Moment<2,2,Hypercube<2> > m;

        for ( int c = 0; c < nComponents; ++c )
        {
            // 1
            this->insert( m.template pick<PolySetType>( 0, c ).toSet( true ), c==0 );
            // x
            this->insert( m.template pick<PolySetType>( 1, c ).toSet( true ) );
            // y
            this->insert( m.template pick<PolySetType>( 3, c ).toSet( true ) );
            // x^2 - y^2
            Polynomial<Moment<2,2,Hypercube<2> >, PolySetType> p( m );
            p = m.template pick<PolySetType>( 2, c );
            p -= m.template pick<PolySetType>( 6, c );
            //std::cout << "x^2-y^2:" << p.coeff() << "\n";
            this->insert( p.toSet( true ) );
        }

    }


};

template<typename Basis, template<class, uint16_type, class> class PointSetType>
class CrouzeixRaviartDual
    :
public DualBasis<Basis>
{
    typedef DualBasis<Basis> super;
public:

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder= super::nOrder;

    typedef typename super::primal_space_type primal_space_type;
    typedef typename primal_space_type::value_type value_type;
    typedef typename primal_space_type::points_type points_type;
    typedef typename primal_space_type::matrix_type matrix_type;
    typedef typename primal_space_type::template convex<2>::type convex_type;
    typedef Reference<convex_type, nDim, 2, nDim, value_type> reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;

    // point set type associated with the functionals
    typedef PointSet<convex_type, value_type> pointset_type;
    typedef PointSetType<convex_type, 2, value_type> equispaced_pointset_type;

    static const uint16_type nVertices = reference_convex_type::numVertices;
    static const uint16_type nFaces = reference_convex_type::numFaces;
    static const uint16_type nGeometricFaces = reference_convex_type::numFaces;
    static const uint16_type nEdges = reference_convex_type::numEdges;
    static const uint16_type nNormals = reference_convex_type::numNormals;

    static const uint16_type nbPtsPerVertex = 0;
    static const uint16_type nbPtsPerEdge = mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<2> >,
                             mpl::int_<reference_convex_type::nbPtsPerEdge>,
                             mpl::int_<0> >::type::value;
    static const uint16_type nbPtsPerFace = mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<3> >,
                             mpl::int_<reference_convex_type::nbPtsPerFace>,
                             mpl::int_<0> >::type::value;
    static const uint16_type nbPtsPerVolume = 0;
    static const uint16_type numPoints = ( reference_convex_type::numGeometricFaces*nbPtsPerFace+
                                           reference_convex_type::numEdges*nbPtsPerEdge );

    /** Number of degrees of freedom per vertex */
    static const uint16_type nDofPerVertex = nbPtsPerVertex;

    /** Number of degrees of freedom per edge */
    static const uint16_type nDofPerEdge = nbPtsPerEdge;

    /** Number of degrees of freedom per face */
    static const uint16_type nDofPerFace = nbPtsPerFace;

    /** Number of degrees  of freedom per volume */
    static const uint16_type nDofPerVolume = nbPtsPerVolume;

    /** Total number of degrees of freedom (equal to refEle::nDof) */
    static const uint16_type nLocalDof = numPoints;

    static const uint16_type nFacesInConvex = mpl::if_< mpl::equal_to<mpl::int_<nDim>, mpl::int_<1> >,
                             mpl::int_<nVertices>,
                             typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
                             mpl::int_<nEdges>,
                             mpl::int_<nFaces> >::type >::type::value;


    CrouzeixRaviartDual( primal_space_type const& primal )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( nDim, numPoints ),
        M_points_face( nFacesInConvex ),
        M_fset( primal )
    {
#if 1
        std::cout << "Lagrange finite element: \n";
        std::cout << " o- dim   = " << nDim << "\n";
        std::cout << " o- order = " << nOrder << "\n";
        std::cout << " o- numPoints      = " << numPoints << "\n";
        std::cout << " o- nbPtsPerVertex = " << ( int )nbPtsPerVertex << "\n";
        std::cout << " o- nbPtsPerEdge   = " << ( int )nbPtsPerEdge << "\n";
        std::cout << " o- nbPtsPerFace   = " << ( int )nbPtsPerFace << "\n";
        std::cout << " o- nbPtsPerVolume = " << ( int )nbPtsPerVolume << "\n";
#endif
        equispaced_pointset_type epts;
        // in d-dimension, consider only the d-1 entity mid-points
        int d = M_convex_ref.topologicalDimension()-1;
        int p = 0;

        // loop on each entity forming the convex of topological
        // dimension d
        for ( int e = M_convex_ref.entityRange( d ).begin();
                e < M_convex_ref.entityRange( d ).end();
                ++e )
        {
            M_points_face[e] = epts.pointsBySubEntity( nDim-1, e, 0 );
            ublas::subrange( M_pts, 0, nDim, p, p+M_points_face[e].size2() ) = M_points_face[e];
            p+=M_points_face[e].size2();
        }

        //std::cout << "[CrouzeixRaviartDual] points= " << M_pts << "\n";
        setFset( primal, M_pts, mpl::bool_<primal_space_type::is_scalar>() );


    }

    points_type const& points() const
    {
        return M_pts;
    }
    points_type const& points( uint16_type f ) const
    {
        return M_points_face[f];
    }


    matrix_type operator()( primal_space_type const& pset ) const
    {
        return M_fset( pset );
    }
private:

    void setFset( primal_space_type const& primal, points_type const& __pts, mpl::bool_<true> )
    {
        M_fset.setFunctionalSet( functional::PointsEvaluation<primal_space_type>( primal,
                                  __pts ) );
    }

    void setFset( primal_space_type const& primal, points_type const& __pts, mpl::bool_<false> )
    {
        M_fset.setFunctionalSet( functional::ComponentsPointsEvaluation<primal_space_type>( primal,
                                  __pts ) );
    }


private:
    reference_convex_type M_convex_ref;
    std::vector<std::vector<uint16_type> > M_eid;
    points_type M_pts;
    std::vector<points_type> M_points_face;
    FunctionalSet<primal_space_type> M_fset;


};
}// detail


/*!
  \class CrouzeixRaviart
  \brief CrouzeixRaviart Finite Element

  @author Christophe Prud'homme
*/
template<uint16_type N,
         uint16_type RealDim,
         template<uint16_type Dim> class PolySetType,
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Convex = Simplex,
         uint16_type TheTAG=0 >
class CrouzeixRaviart
    :
public FiniteElement<typename mpl::if_<mpl::bool_<Convex<N,1,N>::is_simplex>,
                                       mpl::identity<Feel::detail::OrthonormalPolynomialSet<N, 1, RealDim, PolySetType, T, TheTAG, Convex> >,
                                       mpl::identity<fem::detail::RannacherTurekPolynomialSet<N, PolySetType, T> > >::type::type,
    detail::CrouzeixRaviartDual,
    PointSetEquiSpaced >
{
    typedef FiniteElement<typename mpl::if_<mpl::bool_<Convex<N,1,N>::is_simplex>,
                                            mpl::identity<Feel::detail::OrthonormalPolynomialSet<N, 1, RealDim, PolySetType, T, TheTAG, Convex> >,
                                            mpl::identity<fem::detail::RannacherTurekPolynomialSet<N, PolySetType, T> > >::type::type,
                          detail::CrouzeixRaviartDual,
                          PointSetEquiSpaced > super;
public:

    BOOST_STATIC_ASSERT( N > 1 );

    /** @name Typedefs
     */
    //@{

    static const uint16_type nDim = N;
    static const uint16_type nOrder =  super::nOrder;
    static const bool isTransformationEquivalent = true;
    static const bool isContinuous = true;

    typedef typename super::value_type value_type;
    typedef typename super::primal_space_type primal_space_type;
    typedef typename super::dual_space_type dual_space_type;
    typedef Continuous continuity_type;
    static const uint16_type TAG = TheTAG;

    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename super::polyset_type polyset_type;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;
    static const bool is_product = true;

    typedef CrouzeixRaviart<N, RealDim, Scalar, T, Convex> component_basis_type;

    typedef typename dual_space_type::convex_type convex_type;
    typedef typename dual_space_type::pointset_type pointset_type;
    typedef typename dual_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename reference_convex_type::points_type points_type;
    typedef typename convex_type::topological_face_type face_type;

    static const uint16_type nbPtsPerVertex = 0;
    static const uint16_type nbPtsPerEdge = mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<2> >,
                             mpl::int_<reference_convex_type::nbPtsPerEdge>,
                             mpl::int_<0> >::type::value;
    static const uint16_type nbPtsPerFace = mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<3> >,
                             mpl::int_<reference_convex_type::nbPtsPerFace>,
                             mpl::int_<0> >::type::value;
    static const uint16_type nbPtsPerVolume = 0;
    static const uint16_type numPoints = ( reference_convex_type::numGeometricFaces*nbPtsPerFace+
                                           reference_convex_type::numEdges*nbPtsPerEdge );

    static const uint16_type nLocalDof = dual_space_type::nLocalDof;
    static const uint16_type nDofPerVertex = dual_space_type::nDofPerVertex;
    static const uint16_type nDofPerEdge = dual_space_type::nDofPerEdge;
    static const uint16_type nDofPerFace = dual_space_type::nDofPerFace;
    static const uint16_type nDofPerVolume = dual_space_type::nDofPerVolume;
    static const uint16_type nLocalFaceDof = ( face_type::numVertices * nDofPerVertex +
                                               face_type::numEdges * nDofPerEdge +
                                               face_type::numFaces * nDofPerFace );

    struct SSpace
    {
        static constexpr uint16_type TheOrder = (nOrder > 1)?nOrder-1:0;
        typedef typename mpl::if_<mpl::less_equal<mpl::int_<nOrder>, mpl::int_<1> >,
                                  mpl::identity<CrouzeixRaviart<nDim, RealDim, PolySetType, T, Convex, TheTAG> >,
                                  mpl::identity<CrouzeixRaviart<nDim, RealDim, PolySetType, T, Convex, TheTAG> > >::type::type type;

    };

    //@}

    /** @name Constructors, destructor
     */
    //@{

    CrouzeixRaviart()
        :
        super( dual_space_type( primal_space_type() ) ),
        M_refconvex()
    {}
    CrouzeixRaviart( CrouzeixRaviart const & cr )
        :
        super( cr ),
        M_refconvex()
    {}
    ~CrouzeixRaviart()
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

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    typedef Eigen::MatrixXd local_interpolant_type;
    local_interpolant_type
    localInterpolant() const
        {
            return local_interpolant_type::Zero( nComponents*nLocalDof, 1 );
        }

    template<typename ExprType>
    void
    interpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            BOOST_MPL_ASSERT_MSG( nComponents1==ExprType::shape::M,
                                  INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                                  (mpl::int_<nComponents1>,mpl::int_<ExprType::shape::M>));
            BOOST_MPL_ASSERT_MSG( nComponents2==ExprType::shape::N,
                                  INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                                  (mpl::int_<nComponents2>,mpl::int_<ExprType::shape::N>));
            for( int q = 0; q < nLocalDof; ++q )
                for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                    for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                        Ihloc( (c1+nComponents1*c2)*nLocalDof+q ) = expr.evalq( c1, c2, q );

        }
    local_interpolant_type
    faceLocalInterpolant() const
        {
            return local_interpolant_type::Zero( nComponents*nLocalFaceDof, 1 );
        }
    template<typename ExprType>
    void
    faceInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            BOOST_MPL_ASSERT_MSG( nComponents1==ExprType::shape::M,
                                  INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                                  (mpl::int_<nComponents1>,mpl::int_<ExprType::shape::M>));
            BOOST_MPL_ASSERT_MSG( nComponents2==ExprType::shape::N,
                                  INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                                  (mpl::int_<nComponents2>,mpl::int_<ExprType::shape::N>));
            for( int q = 0; q < nLocalFaceDof; ++q )
                for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                    for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                        Ihloc( (c1+nComponents1*c2)*nLocalFaceDof+q ) = expr.evalq( c1, c2, q );

        }

    template<typename ExprType>
    void
    interpolateBasisFunction( ExprType& expr, local_interpolant_type& Ihloc ) const
    {
        BOOST_MPL_ASSERT_MSG( nComponents1==ExprType::shape::M,
                              INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                              (mpl::int_<nComponents1>,mpl::int_<ExprType::shape::M>));
        BOOST_MPL_ASSERT_MSG( nComponents2==ExprType::shape::N,
                              INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                              (mpl::int_<nComponents2>,mpl::int_<ExprType::shape::N>));

        //for ( int cc1 = 0; cc1 < nComponents1; ++cc1 )
        typedef typename ExprType::tensor_expr_type::expression_type::fe_type fe_expr_type;

        for( int q = 0; q <expr.geom()->nPoints(); ++q )
            for( int i = 0; i < fe_expr_type::nLocalDof; ++i )
                for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                    for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                    {
                        int ldof = (c1+fe_expr_type::nComponents1*c2)*fe_expr_type::nLocalDof + i;
                        int ldof2 = (fe_expr_type::is_product)? ldof : i;
                        Ihloc( ldof, q ) = expr.evaliq( /*i*/ldof2, c1, c2, q );
                    }
    }

    /**
     * \return the family name of the finite element
     */
    std::string familyName() const
    {
        return "CrouzeixRaviart";
    }

    template<typename ExprType>
    static auto
    isomorphism( ExprType& expr ) -> decltype( expr )
    {
        return expr;
        //return expr;
    }
    //@}



protected:
    reference_convex_type M_refconvex;
private:

};
template<uint16_type N,
         uint16_type RealDim,
         template<uint16_type Dim> class PolySetType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex,
         uint16_type TheTAG >
const uint16_type CrouzeixRaviart<N,RealDim,PolySetType,T,Convex,TheTAG>::nDim;
template<uint16_type N,
         uint16_type RealDim,
         template<uint16_type Dim> class PolySetType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex,
         uint16_type TheTAG >
const uint16_type CrouzeixRaviart<N,RealDim,PolySetType,T,Convex,TheTAG>::nOrder;

} // fem

template<uint16_type Order,
         template<uint16_type Dim> class PolySetType = Scalar,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced,
         uint16_type TheTAG=0 >
class CrouzeixRaviart
{
public:
    template<uint16_type N,
             uint16_type RealDim,
             typename T = double,
             typename Convex = Simplex<N> >
    struct apply
    {
        typedef typename mpl::if_<mpl::bool_<Convex::is_simplex>,
                mpl::identity<fem::CrouzeixRaviart<N,RealDim,PolySetType,T,Simplex,TheTAG> >,
                mpl::identity<fem::CrouzeixRaviart<N,RealDim,PolySetType,T,Hypercube,TheTAG> > >::type::type result_type;
        typedef result_type type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef CrouzeixRaviart<Order,PolySetType,Pts,TheNewTAG> type;
    };

    typedef CrouzeixRaviart<Order,Scalar,Pts> component_basis_type;

    static const uint16_type nOrder =  Order;
    static const uint16_type TAG = TheTAG;

};

template<uint16_type Order,
         template<uint16_type Dim> class PolySetType,
         template<class, uint16_type, class> class Pts,
         uint16_type TheTAG >
const uint16_type CrouzeixRaviart<Order,PolySetType,Pts,TheTAG>::nOrder;


} // Feel
#endif /* __CrouzeixRaviart_H */
