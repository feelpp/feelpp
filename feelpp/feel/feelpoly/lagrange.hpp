/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-18

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2008-2012 Université de Grenoble 1 (Joseph Fourier)

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
   \file lagrange.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-18
 */
#ifndef __lagrange_H
#define __lagrange_H 1

#include <boost/ptr_container/ptr_vector.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/lu.hpp>

#include <feel/feelmesh/refentity.hpp>
#include <feel/feelmesh/pointset.hpp>
#include <feel/feelpoly/equispaced.hpp>


#include <feel/feelpoly/dualbasis.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/functionalset.hpp>
#include <feel/feelpoly/functionals.hpp>
#include <feel/feelpoly/fe.hpp>
#include <feel/feelpoly/isp0continuous.hpp>





namespace Feel
{

namespace fem
{

/// \cond detail
namespace details
{
template<typename Basis, template<class, uint16_type, class> class PointSetType>
class LagrangeDual
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
    typedef typename primal_space_type::convex_type convex_type;
    typedef typename primal_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;

    // point set type associated with the functionals
    typedef PointSetType<convex_type, nOrder, value_type> pointset_type;

    template< template<class, uint16_type, class> class TestPointSetType >
    static const bool is_pointset_v = std::is_base_of_v<TestPointSetType<convex_type, nOrder, value_type>,pointset_type >;

    static const uint16_type numPoints = reference_convex_type::numPoints;
    static const uint16_type nbPtsPerVertex = reference_convex_type::nbPtsPerVertex;
    static const uint16_type nbPtsPerEdge = reference_convex_type::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace = reference_convex_type::nbPtsPerFace;
    static const uint16_type nbPtsPerVolume = reference_convex_type::nbPtsPerVolume;

#if 0
    /**
     * for Dim >= 3 : n edges = n(vertices) + n(faces) - 2
     * thanks to Euler formula
     */
    typedef mpl::vector_c<uint16_type, 0, 1, 3, ( 4 ) + ( 4 ) - 2> edges_t;
    typedef mpl::vector_c<uint16_type, 0, 0, 1, 4> geo_faces_t;
    typedef mpl::vector_c<uint16_type, 0, 2, 3, 4> faces_t;
    typedef mpl::vector_c<uint16_type, 0, 2, 3, 4> normals_t;
#endif

    static const uint16_type nVertices = reference_convex_type::numVertices;
    static const uint16_type nFaces = reference_convex_type::numFaces;
    static const uint16_type nGeometricFaces = reference_convex_type::numFaces;
    static const uint16_type nEdges = reference_convex_type::numEdges;
    static const uint16_type nNormals = reference_convex_type::numNormals;


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

    LagrangeDual( LagrangeDual const& d )
        :
        super( d ),
        M_convex_ref(),
        M_eid( d.M_eid ),
        M_pts( d.M_pts ),
        M_points_face( d.M_points_face ),
        M_fset( d.M_fset )
        {}

    LagrangeDual( LagrangeDual && d ) = default;
    LagrangeDual( primal_space_type const& primal )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( nDim, numPoints ),
        M_points_face( nFacesInConvex ),
        M_fset( primal )
    {
        DVLOG(2) << "Lagrange finite element: \n";
        DVLOG(2) << " o- dim   = " << nDim << "\n";
        DVLOG(2) << " o- order = " << nOrder << "\n";
        DVLOG(2) << " o- numPoints      = " << numPoints << "\n";
        DVLOG(2) << " o- nbPtsPerVertex = " << nbPtsPerVertex << "\n";
        DVLOG(2) << " o- nbPtsPerEdge   = " << nbPtsPerEdge << "\n";
        DVLOG(2) << " o- nbPtsPerFace   = " << nbPtsPerFace << "\n";
        DVLOG(2) << " o- nbPtsPerVolume = " << nbPtsPerVolume << "\n";

        M_pts = M_pset.points();

        if ( nOrder > 0 )
        {
            for ( uint16_type e = M_convex_ref.entityRange( nDim-1 ).begin();
                    e < M_convex_ref.entityRange( nDim-1 ).end();
                    ++e )
            {
                M_points_face[e] = M_pset.pointsBySubEntity( nDim-1, e, 1 );
                DVLOG(2) << "face " << e << " pts " <<  M_points_face[e] << "\n";
            }
        }

        setFset( primal, M_pts, mpl::bool_<primal_space_type::is_scalar>() );
    }

    LagrangeDual( primal_space_type const& primal, pointset_type const& pts )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( pts.points() ),
        M_points_face( nFacesInConvex ),
        M_fset( primal ),
        M_pset( pts )
    {
        DVLOG(2) << "Lagrange finite element: \n";
        DVLOG(2) << " o- dim   = " << nDim << "\n";
        DVLOG(2) << " o- order = " << nOrder << "\n";
        DVLOG(2) << " o- numPoints      = " << numPoints << "\n";
        DVLOG(2) << " o- nbPtsPerVertex = " << nbPtsPerVertex << "\n";
        DVLOG(2) << " o- nbPtsPerEdge   = " << nbPtsPerEdge << "\n";
        DVLOG(2) << " o- nbPtsPerFace   = " << nbPtsPerFace << "\n";
        DVLOG(2) << " o- nbPtsPerVolume = " << nbPtsPerVolume << "\n";

        if ( nOrder > 0 )
        {
            for ( uint16_type e = M_convex_ref.entityRange( nDim-1 ).begin();
                    e < M_convex_ref.entityRange( nDim-1 ).end();
                    ++e )
            {
                M_points_face[e] = M_pset.pointsBySubEntity( nDim-1, e, 1 );
                DVLOG(2) << "face " << e << " pts " <<  M_points_face[e] << "\n";
            }
        }

        setFset( primal, M_pts, mpl::bool_<primal_space_type::is_scalar>() );
    }

    ~LagrangeDual() = default;
    LagrangeDual& operator=( LagrangeDual const& ) = default;

    points_type const& points() const
    {
        return M_pts;
    }

    points_type const& points( uint16_type f ) const
    {
        return M_points_face[f];
    }
    ublas::matrix_column<points_type const> point( uint16_type f, uint32_type __i ) const
    {
        return ublas::column( M_points_face[f], __i );
    }
    ublas::matrix_column<points_type> point( uint16_type f, uint32_type __i )
    {
        return ublas::column( M_points_face[f], __i );
    }

#if 0
    std::vector<point_type> points( int topodim ) const
        {
            std::vector<point_type> pts( ;
            for ( uint16_type e = M_convex_ref.entityRange( nDim-1 ).begin();
                  e < M_convex_ref.entityRange( nDim-1 ).end();
                  ++e )
            {
                M_points_face[e] = M_pset.pointsBySubEntity( nDim-1, e, 1 );
            }
            return M_pset.pointsBySubEntity( topodim, edge, 1 );
        }
#endif
     points_type points( int topodim, int entity ) const
            {
                return M_pset.pointsBySubEntity( topodim, entity, 1 );
            }


    points_type edgePoints(int edge) const
        {
            return M_pset.pointsBySubEntity( 1, edge, 1 );
        }


    points_type vertexPoints(int vertex) const
        {
            return M_pset.pointsBySubEntity( 0, vertex, 1 );
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

    /**
     * set the pointset at face \c f using points \c n
     */
    void setPoints( uint16_type f, points_type const& n )
    {
        M_points_face[f].resize( n.size1(), n.size2(), false );
        M_points_face[f] = n;
    }

private:
    reference_convex_type M_convex_ref;
    std::vector<std::vector<uint16_type> > M_eid;
    points_type M_pts;
    std::vector<points_type> M_points_face;
    FunctionalSet<primal_space_type> M_fset;
    pointset_type M_pset;

};
}// details
/// \endcond detail

    class LagrangePolynomialSet {};
/**
 * \class Lagrange
 * \brief Lagrange polynomial set
 *
 * The \p Lagrange polynomial set is parametrized by
 *
 * -# dimension of the geometrical space
 * -# order of the Lagrange polynomials
 * -# the numerical type
 * -# the geometry it applies to (convexes such as simplices or product of simplices)
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type N,
         uint16_type RealDim,
         uint16_type O,
         template<uint16_type Dim> class PolySetType,
         typename ContinuityType = Continuous,
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Convex = Simplex,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced,
         uint16_type TheTAG = 0 >
class Lagrange
    :
    public LagrangePolynomialSet,
    public FiniteElement<Feel::detail::OrthonormalPolynomialSet<N, O, RealDim, PolySetType, T, TheTAG, Convex>, details::LagrangeDual, Pts >
{
    typedef FiniteElement<Feel::detail::OrthonormalPolynomialSet<N, O, RealDim, PolySetType, T, TheTAG, Convex>, details::LagrangeDual, Pts > super;
public:

    BOOST_STATIC_ASSERT( ( boost::is_same<PolySetType<N>, Scalar<N> >::value ||
                           boost::is_same<PolySetType<N>, Vectorial<N> >::value ||
                           boost::is_same<PolySetType<N>, Tensor2<N> >::value ||
                           boost::is_same<PolySetType<N>, Tensor2Symm<N> >::value ) );

    /** @name Typedefs
     */
    //@{

    static const uint16_type nDim = N;
    static const uint16_type nRealDim = RealDim;
    static const uint16_type nOrder =  O;
    static const bool isTransformationEquivalent = true;
    static const bool isContinuous = ContinuityType::is_continuous;
    typedef typename super::value_type value_type;
    typedef typename super::primal_space_type primal_space_type;
    typedef typename super::dual_space_type dual_space_type;
    typedef ContinuityType continuity_type;
    static const uint16_type TAG = TheTAG;

    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename super::polyset_type polyset_type;
    static constexpr bool is_symm_v  = Feel::is_symm_v<polyset_type>;
    using is_symm  = Feel::is_symm<polyset_type>;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_tensor2symm = is_tensor2 && is_symm_v;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;

    static const bool is_product = true;
    static constexpr int Nm2 = (N>2)?N-2:0;
    static constexpr int Nm1 = (N>0)?N-1:0;

    typedef Lagrange<N, RealDim, O, PolySetType, ContinuityType, T, Convex,  Pts, TheTAG> this_type;
    typedef Lagrange<N, RealDim, O, Scalar, continuity_type, T, Convex,  Pts, TheTAG> component_basis_type;


    typedef typename mpl::if_<mpl::less_equal<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<boost::none_t>,
            mpl::identity< Lagrange<Nm1, RealDim, O, Scalar, continuity_type, T, Convex,  Pts, TheTAG> > >::type::type face_basis_type;

    typedef std::shared_ptr<face_basis_type> face_basis_ptrtype;
    typedef typename mpl::if_<mpl::less_equal<mpl::int_<nDim>, mpl::int_<2> >,
                              mpl::identity<boost::none_t>,
                              mpl::identity< Lagrange<Nm2, RealDim, O, Scalar, continuity_type, T, Convex,  Pts, TheTAG> > >::type::type edge_basis_type;

    typedef std::shared_ptr<edge_basis_type> edge_basis_ptrtype;

    typedef typename dual_space_type::convex_type convex_type;
    typedef typename dual_space_type::pointset_type pointset_type;
    typedef typename dual_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename reference_convex_type::points_type points_type;
    typedef typename convex_type::topological_face_type face_type;
    typedef typename convex_type::edge_type edge_type;

    static const uint16_type numPoints = reference_convex_type::numPoints;
    static const uint16_type nbPtsPerVertex = reference_convex_type::nbPtsPerVertex;
    static const uint16_type nbPtsPerEdge = reference_convex_type::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace = reference_convex_type::nbPtsPerFace;
    static const uint16_type nbPtsPerVolume = reference_convex_type::nbPtsPerVolume;
    static const uint16_type nLocalDof = dual_space_type::nLocalDof;
    static const uint16_type nDofPerVertex = dual_space_type::nDofPerVertex;
    static const uint16_type nDofPerEdge = dual_space_type::nDofPerEdge;
    static const uint16_type nDofPerFace = dual_space_type::nDofPerFace;
    static const uint16_type nDofPerVolume = dual_space_type::nDofPerVolume;
    static const uint16_type nLocalFaceDof = ( face_type::numVertices * nDofPerVertex +
                                               face_type::numEdges * nDofPerEdge +
                                               face_type::numFaces * nDofPerFace );
    static const uint16_type nLocalEdgeDof = ( edge_type::numVertices * nDofPerVertex +
                                               edge_type::numEdges * nDofPerEdge);
    static const uint16_type nLocalVertexDof = nDofPerVertex;
    template<int subN>
    struct SubSpace
    {
        typedef Lagrange<N-1, RealDim, O, PolySetType, continuity_type, T, Convex,  Pts, TheTAG> type;
    };

    struct SSpace
    {
        static constexpr uint16_type TheOrder = (O > 1)?O-1:0;
        typedef typename mpl::if_<mpl::less_equal<mpl::int_<O>, mpl::int_<1> >,
                                  mpl::identity<Lagrange<N, RealDim, 0, PolySetType, Discontinuous, T, Convex,  Pts, TheTAG> >,
                                  mpl::identity<Lagrange<N, RealDim, TheOrder, PolySetType, continuity_type, T, Convex,  Pts, TheTAG> > >::type::type type;

    };

    template<uint16_type NewDim>
    struct ChangeDim
    {
        typedef Lagrange<NewDim, RealDim, O, PolySetType, continuity_type, T, Convex,  Pts, TheTAG> type;
    };

    static const bool isLagrangeP0Continuous = isP0Continuous<this_type>::result;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    Lagrange()
        :
        super( dual_space_type( primal_space_type() ) ),
        M_refconvex()
        // M_bdylag( new face_basis_type )
    {

        if ( is_tensor2symm )
        {
            M_unsymm2symm.resize( nComponents*nLocalDof );
            for ( uint16_type l = 0; l < nLocalDof; ++l )
            {
                for ( int c1 = 0; c1 < nComponents1; ++c1)
                {
                    for ( int c2 = c1+1; c2 < nComponents2; ++c2 )
                    {
                        const int k = Feel::detail::symmetricIndex(c1,c2,nComponents1);
                        M_unsymm2symm[ nLocalDof*(nComponents1*c1+c2) + l ] = nLocalDof*k+l;
                        M_unsymm2symm[ nLocalDof*(nComponents1*c2+c1) + l ] = nLocalDof*k+l;
                    }
                    const int k = Feel::detail::symmetricIndex(c1,c1,nComponents1);
                    M_unsymm2symm[ nLocalDof*(nComponents1*c1+c1) + l ] = nLocalDof*k+l;
                }
            }
        }

        // std::cout << "[LagrangeDual] points= " << M_pts << "\n";
    }

    ~Lagrange() override {}

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

    /**
     * \return the family name of the finite element
     */
    std::string familyName() const override
    {
        return "lagrange";
    }


    //! \return the component of a local dof
    uint16_type component( uint16_type localDofId ) const override
        {
            uint16_type comp = localDofId/nLocalDof;
            DCHECK( comp < nComponents ) << "invalid localDofId " << localDofId;
            return comp;
        }

    //! \return a parent local dof id for each component (for example, the first component)
    uint16_type dofParent( uint16_type localDofId ) const override
        {
            uint16_type ldofParent = localDofId % nLocalDof;
            return ldofParent;
        }

    //! \return the type of a local dof
    uint16_type dofType( uint16_type localDofId ) const override
        {
            return 1;
        }

    //! give an unsymmetric dof index i, provide the symmetric one
    uint16_type unsymmToSymm( uint16_type i ) const override
        {
            if ( !is_tensor2symm )
                return i;
            DCHECK( M_unsymm2symm.size() > i ) << "invalid size of unsymm2symm container";
            return M_unsymm2symm[i];
        }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //!
    //! build interpolant object on \p n finite elements
    //!
    using  local_interpolant_type = Eigen::VectorXd;
    local_interpolant_type
    localInterpolant( int n = 1 ) const
        {
            return local_interpolant_type::Zero( n*nComponents*nLocalDof );
        }

    //!
    //! build \p p  interpolants object on \p n finite elements
    //!
    using  local_interpolants_type = Eigen::MatrixXd;
    local_interpolants_type
    localInterpolants( int p, int n = 1 ) const
        {
            return local_interpolants_type::Zero( n*nComponents*nLocalDof, p );
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
                {
                    if ( is_symm_v )
                    {
                        for( int c2 = 0; c2 < c1; ++c2 )
                        {
                            Ihloc( (c2+nComponents2*c1)*nLocalDof+q ) = expr.evalq( c1, c2, q );
                            Ihloc( (c1+nComponents2*c2)*nLocalDof+q ) = Ihloc( (c2+nComponents2*c1)*nLocalDof+q );
                        }
                        // diagonal
                        Ihloc( (c1+nComponents2*c1)*nLocalDof+q ) = expr.evalq( c1, c1, q );
                    }
                    else
                    {
                        for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                            Ihloc( (c2+nComponents2*c1)*nLocalDof+q ) = expr.evalq( c1, c2, q );
                    }
                }
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
                {
                    if ( is_symm_v )
                    {
                        for( int c2 = 0; c2 < c1; ++c2 )
                        {
                            Ihloc( (c2+nComponents2*c1)*nLocalFaceDof+q ) = expr.evalq( c1, c2, q );
                            Ihloc( (c1+nComponents2*c2)*nLocalFaceDof+q ) = Ihloc( (c2+nComponents2*c1)*nLocalFaceDof+q );
                        }
                        Ihloc( (c1+nComponents2*c1)*nLocalFaceDof+q ) = expr.evalq( c1, c1, q );
                    }
                    else
                    {
                        for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                            Ihloc( (c2+nComponents2*c1)*nLocalFaceDof+q ) = expr.evalq( c1, c2, q );
                    }
                }
        }

    local_interpolant_type
    edgeLocalInterpolant() const
        {
            return local_interpolant_type::Zero( nComponents*nLocalEdgeDof, 1 );
        }
    template<typename ExprType>
    void
    edgeInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            BOOST_MPL_ASSERT_MSG( nComponents1==ExprType::shape::M,
                                  INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                                  (mpl::int_<nComponents1>,mpl::int_<ExprType::shape::M>));
            BOOST_MPL_ASSERT_MSG( nComponents2==ExprType::shape::N,
                                  INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                                  (mpl::int_<nComponents2>,mpl::int_<ExprType::shape::N>));
            for( int q = 0; q < nLocalEdgeDof; ++q )
                for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                {
                    if ( is_symm_v )
                    {
                        for( int c2 = 0; c2 < c1; ++c2 )
                        {
                            Ihloc( (c2+nComponents2*c1)*nLocalEdgeDof+q ) = expr.evalq( c1, c2, q );
                            Ihloc( (c1+nComponents2*c2)*nLocalEdgeDof+q ) = Ihloc( (c2+nComponents2*c1)*nLocalEdgeDof+q );
                        }
                        Ihloc( (c1+nComponents2*c1)*nLocalEdgeDof+q ) = expr.evalq( c1, c1, q );

                    }
                    else
                    {
                        for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                            Ihloc( (c2+nComponents2*c1)*nLocalEdgeDof+q ) = expr.evalq( c1, c2, q );
                    }
                }
        }
    local_interpolant_type
    vertexLocalInterpolant() const
        {
            return local_interpolant_type::Zero( nComponents*nLocalVertexDof, 1 );
        }
    template<typename ExprType>
    void
    vertexInterpolate( ExprType& expr, local_interpolant_type& Ihloc ) const
        {
            BOOST_MPL_ASSERT_MSG( nComponents1==ExprType::shape::M,
                                  INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                                  (mpl::int_<nComponents1>,mpl::int_<ExprType::shape::M>));
            BOOST_MPL_ASSERT_MSG( nComponents2==ExprType::shape::N,
                                  INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                                  (mpl::int_<nComponents2>,mpl::int_<ExprType::shape::N>));
            for( int q = 0; q < nLocalVertexDof; ++q )
                for( int c1 = 0; c1 < ExprType::shape::M; ++c1 )
                {
                    if ( is_symm_v )
                    {
                        for( int c2 = 0; c2 < c1; ++c2 )
                        {
                            Ihloc( (c2+nComponents2*c1)*nLocalVertexDof+q ) = expr.evalq( c1, c2, q );
                            Ihloc( (c1+nComponents2*c2)*nLocalVertexDof+q ) = Ihloc( (c2+nComponents2*c1)*nLocalVertexDof+q );
                        }
                        Ihloc( (c1+nComponents2*c1)*nLocalVertexDof+q ) = expr.evalq( c1, c1, q );
                    }
                    else
                    {
                        for( int c2 = 0; c2 < ExprType::shape::N; ++c2 )
                            Ihloc( (c2+nComponents2*c1)*nLocalVertexDof+q ) = expr.evalq( c1, c2, q );
                    }

                }
        }
    template<typename ExprType>
    void
    interpolateBasisFunction( ExprType&& expr, local_interpolants_type & Ihloc ) const
    {
        using shape = typename std::decay_t<ExprType>::shape;
        /*
        BOOST_MPL_ASSERT_MSG( nComponents1==shape::M,
                              INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                              (mpl::int_<nComponents1>,mpl::int_<shape::M>));
        BOOST_MPL_ASSERT_MSG( nComponents2==shape::N,
                              INCOMPATIBLE_NUMBER_OF_COMPONENTS,
                              (mpl::int_<nComponents2>,mpl::int_<shape::N>));
         */
        //for ( int cc1 = 0; cc1 < nComponents1; ++cc1 )
        using expr_basis_t = typename std::decay_t<ExprType>::expr_type::test_basis;

        for( int q = 0; q < nLocalDof; ++q )
        {
            for( int i = 0; i < expr_basis_t::nLocalDof; ++i )
            {
                int ncomp1= ( expr_basis_t::is_product?expr_basis_t::nComponents1:1 );

                for ( uint16_type c = 0; c < ncomp1; ++c )
                {
                    uint16_type I = expr_basis_t::nLocalDof*c + i;
                    for( int c1 = 0; c1 < shape::M; ++c1 )
                    {
                        if ( is_symm_v )
                        {
                            for( int c2 = 0; c2 < c1; ++c2 )
                            {
                                int ldof = (c2+nComponents2*c1)*nLocalDof + q;
                                Ihloc( I, ldof) = expr.evaliq( I, c1, c2, q );
                                int ldof2 = (c1+nComponents2*c2)*nLocalDof + q;
                                Ihloc( I, ldof2) = Ihloc( I, ldof);
                            }
                            int ldof = (c1+nComponents2*c1)*nLocalDof + q;
                            Ihloc( I, ldof) = expr.evaliq( I, c1, c1, q );
                        }
                        else
                        {
                            for( int c2 = 0; c2 < shape::N; ++c2 )
                            {
                                int ldof = (c2+nComponents2*c1)*nLocalDof + q;
                                Ihloc( I, ldof) = expr.evaliq( I, c1, c2, q );
                            }
                        }
                    }
                }
            }
        }
    }
    //@}

private:

    reference_convex_type M_refconvex;
    face_basis_ptrtype M_bdylag;
    std::vector<uint16_type> M_unsymm2symm;
};
template<uint16_type N,
         uint16_type RealDim,
         uint16_type O,
         template<uint16_type Dim> class PolySetType,
         typename ContinuityType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex,
         template<class, uint16_type, class> class Pts,
         uint16_type TheTAG >
const uint16_type Lagrange<N,RealDim,O,PolySetType,ContinuityType,T,Convex,Pts,TheTAG>::nDim;

template<uint16_type N,
         uint16_type RealDim,
         uint16_type O,
         template<uint16_type Dim> class PolySetType,
         typename ContinuityType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex,
         template<class, uint16_type, class> class Pts,
         uint16_type TheTAG >
const uint16_type Lagrange<N,RealDim,O,PolySetType,ContinuityType,T,Convex,Pts,TheTAG>::nOrder;

template<uint16_type N,
         uint16_type RealDim,
         uint16_type O,
         template<uint16_type Dim> class PolySetType,
         typename ContinuityType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex,
         template<class, uint16_type, class> class Pts,
         uint16_type TheTAG >
const uint16_type Lagrange<N,RealDim,O,PolySetType,ContinuityType,T,Convex,Pts,TheTAG>::nLocalDof;

} // namespace fem
template<uint16_type Order,
         template<uint16_type Dim> class PolySetType = Scalar,
         typename ContinuityType = Continuous,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced,
         uint16_type TheTAG=0 >
class Lagrange
{
public:
    template<uint16_type N,
             uint16_type RealDim,
             typename T = double,
             typename Convex = Simplex<N> >
    struct apply
    {
        typedef typename mpl::if_<mpl::bool_<Convex::is_simplex>,
                mpl::identity<fem::Lagrange<N,RealDim,Order,PolySetType,ContinuityType,T,Simplex, Pts, TheTAG > >,
                mpl::identity<fem::Lagrange<N,RealDim,Order,PolySetType,ContinuityType,T,Hypercube, Pts, TheTAG> > >::type::type result_type;
        typedef result_type type;
    };

    template<uint16_type TheNewTAG>
    struct ChangeTag
    {
        typedef Lagrange<Order,PolySetType,ContinuityType,Pts,TheNewTAG> type;
    };

    typedef Lagrange<Order,Scalar,ContinuityType,Pts,TheTAG> component_basis_type;

    static const uint16_type nOrder =  Order;
    static const uint16_type TAG = TheTAG;

};
template<uint16_type Order,
         template<uint16_type Dim> class PolySetType,
         typename ContinuityType,
         template<class, uint16_type, class> class Pts,
         uint16_type TheTAG>
const uint16_type Lagrange<Order,PolySetType,ContinuityType,Pts,TheTAG>::nOrder;


template<typename P>
using is_lagrange_polynomialset = std::is_base_of<fem::LagrangePolynomialSet,P>;
template<typename P>
constexpr bool is_lagrange_polynomialset_v = boost::is_base_of<fem::LagrangePolynomialSet,P>::value;

template<typename P>
constexpr bool is_lagrange_polynomialset_P0d_v = is_lagrange_polynomialset_v<P> && (P::nOrder == 0) && !P::continuity_type::is_continuous;

} // namespace Feel
#endif /* __lagrange_H */
