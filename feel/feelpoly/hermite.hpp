/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-03-04

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file hermite.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-03-04
 */
#ifndef __hermite_H
#define __hermite_H 1

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
namespace Feel
{

namespace fem
{

/// \cond detail
namespace details
{
template<typename Basis, template<class, uint16_type, class> class PointSetType>
class HermiteDual
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
    typedef PointSetType<convex_type, 3, value_type> pointset_type;


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

    static const uint16_type nbPtsPerVertex = ( nDim+1 ); //reference_convex_type::nbPtsPerVertex;
    static const uint16_type nbPtsPerEdge = 0;//reference_convex_type::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace = ( nDim==2 )?1:0; //reference_convex_type::nbPtsPerFace;
    static const uint16_type nbPtsPerVolume = 0;//reference_convex_type::nbPtsPerVolume;
    static const uint16_type numPoints = nVertices*nbPtsPerVertex+nGeometricFaces*nbPtsPerFace;


    /** Number of degrees of freedom per vertex */
    static const uint16_type nDofPerVertex = nbPtsPerVertex;

    /** Number of degrees of freedom per edge */
    static const uint16_type nDofPerEdge =  0;//(nDim+1)*nbPtsPerEdge;

    /** Number of degrees of freedom per face */
    static const uint16_type nDofPerFace =  nbPtsPerFace;

    /** Number of degrees  of freedom per volume */
    static const uint16_type nDofPerVolume =  0;//(nDim+1)*nbPtsPerVolume;

    /** Total number of degrees of freedom (equal to refEle::nDof) */
    static const uint16_type nLocalDof =  nDofPerVertex*nVertices + nDofPerFace*nGeometricFaces;

    static const uint16_type nFacesInConvex = mpl::if_< mpl::equal_to<mpl::int_<nDim>, mpl::int_<1> >,
                             mpl::int_<nVertices>,
                             typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<2> >,
                             mpl::int_<nEdges>,
                             mpl::int_<nFaces> >::type >::type::value;

    HermiteDual( primal_space_type const& primal )
        :
        super( primal ),
        M_convex_ref(),
        M_eid( M_convex_ref.topologicalDimension()+1 ),
        M_pts( nDim, nLocalDof ),
        M_points_face( nFacesInConvex ),
        M_fset( primal )
    {
        DVLOG(2) << "Hermite finite element: \n";
        DVLOG(2) << " o- dim   = " << nDim << "\n";
        DVLOG(2) << " o- order = " << nOrder << "\n";
        DVLOG(2) << " o- nVertices      = " << nVertices << "\n";
        DVLOG(2) << " o- nGeometricFaces= " << nGeometricFaces << "\n";
        DVLOG(2) << " o- numPoints      = " << numPoints << "\n";
        DVLOG(2) << " o- nbPtsPerVertex = " << nbPtsPerVertex << "\n";
        DVLOG(2) << " o- nbPtsPerEdge   = " << nbPtsPerEdge << "\n";
        DVLOG(2) << " o- nbPtsPerFace   = " << nbPtsPerFace << "\n";
        DVLOG(2) << " o- nbPtsPerVolume = " << nbPtsPerVolume << "\n";

        DVLOG(2) << " o- nbDofPerVertex = " << nDofPerVertex << "\n";
        DVLOG(2) << " o- nbDofPerEdge   = " << nDofPerEdge << "\n";
        DVLOG(2) << " o- nbDofPerFace   = " << nDofPerFace << "\n";
        DVLOG(2) << " o- nbDofPerVolume = " << nDofPerVolume << "\n";
        DVLOG(2) << " o- nLocalDof = " << nLocalDof << "\n";

        pointset_type pts;

        //M_pts = pts.points();

        int i = 0;

        for ( uint16_type e = M_convex_ref.entityRange( 0 ).begin();
                e < M_convex_ref.entityRange( 0 ).end();
                ++e )
        {
            M_points_face[e] = pts.pointsBySubEntity( 0, e, 1 );
            points_type _pts = pts.pointsBySubEntity( 0, e, 1 );

            for ( int j = 0; j < _pts.size2(); ++j, ++i )
            {
                ublas::column( M_pts, i ) = ublas::column( _pts, j );
                DVLOG(2) << "pts " << i << " = " <<  ublas::column( M_pts, i ) << "\n";
            }
        }

        for ( int j = 0; j < nDim; ++j )
        {
            ublas::subrange( M_pts, 0, nDim, i, i+nDim+1 ) = ublas::subrange( M_pts, 0, nDim, 0, nDim+1 );
            i+=nDim+1;
        }

        if ( nDim == 2 )
        {
            points_type _pts = pts.pointsBySubEntity( 2, 0, 0 );
            ublas::column( M_pts, i ) = ublas::column( _pts, 0 );
            DVLOG(2) << "pts " << i << " = " <<  ublas::column( M_pts, i ) << "\n";
        }

        //std::cout << "pts = " << M_pts << "\n";
        setFset( primal, M_pts, mpl::bool_<primal_space_type::is_scalar>() );
    }
    ~HermiteDual()
    {

    }
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

    matrix_type operator()( primal_space_type const& pset ) const
    {
        matrix_type m = M_fset( pset );
        //std::cout << "m=" << m << "\n";
        return m;
    }
private:

    void setFset( primal_space_type const& primal, points_type const& __pts, mpl::bool_<true> )
    {
        int nfs = 1+nDim+( ( nDim==2 )?1:0 );
        std::vector<std::vector<Functional<primal_space_type> > > pd( nfs );
        points_type pts1 = ublas::project( __pts,
                                           ublas::range( 0,nDim ),
                                           ublas::range( 0, nDim+1 ) );
        pd[0] = functional::PointsEvaluation<primal_space_type>( primal, pts1 );

        //std::cout << "pd[" << 0 << "].size()=" << pd[0].size() << "\n";
        for ( int d = 0; d < nDim; ++d )
        {
            pd[d+1] = functional::PointsDerivative<primal_space_type>( primal, d, pts1 );
            //std::cout << "pd[" << d+1 << "].size()=" << pd[d+1].size() << "\n";
        }

        if ( nDim == 2 )
        {
            points_type pts3 = ublas::project( __pts,
                                               ublas::range( 0,nDim ),
                                               ublas::range( ( nDim+1 )*( nDim+1 ), numPoints ) );
            //std::cout << "pts3 = " << pts3 << "\n";
            pd[3] = functional::PointsEvaluation<primal_space_type>( primal, pts3 );
            //std::cout << "pd[2]=" << pd[3][0].coeff() << "\n";
        }

        std::vector<Functional<primal_space_type> > fs( nLocalDof );
        typename std::vector<Functional<primal_space_type> >::iterator it = fs.begin();

        for ( int i = 0; i < nfs; ++i )
        {
            //std::cout << "pd[" << i << "].size()=" << pd[i].size() << "\n";
            it = std::copy( pd[i].begin(), pd[i].end(), it );
        }

        M_fset.setFunctionalSet( fs );
    }

    void setFset( primal_space_type const& primal, points_type const& __pts, mpl::bool_<false> )
    {
        //M_fset.setFunctionalSet( functional::ComponentsPointsEvaluation<primal_space_type>( primal,
        //__pts ) );
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


};
}// details
/// \endcond detail

/**
 * \class Hermite
 * \brief Hermite polynomial set
 *
 * The \p Hermite polynomial set is parametrized by
 *
 * -# dimension of the geometrical space
 * -# order of the Hermite polynomials
 * -# the numerical type
 * -# the geometry it applies to (convexes such as simplices or product of simplices)
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type N,
         uint16_type O,
         template<uint16_type Dim> class PolySetType,
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Convex = Simplex,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced >
class Hermite
    :
public FiniteElement<detail::OrthonormalPolynomialSet<N, N, O, PolySetType, T, Convex>, details::HermiteDual, Pts >
{
    typedef FiniteElement<detail::OrthonormalPolynomialSet<N, N, O, PolySetType, T, Convex>, details::HermiteDual, Pts > super;
public:

    BOOST_STATIC_ASSERT( ( boost::is_same<PolySetType<N>, Scalar<N> >::value ||
                           boost::is_same<PolySetType<N>, Vectorial<N> >::value ||
                           boost::is_same<PolySetType<N>, Tensor2<N> >::value ) );

    /** @name Typedefs
     */
    //@{

    static const uint16_type nDim = N;
    static const uint16_type nOrder =  O;
    static const bool isTransformationEquivalent = false;
    static const bool isContinuous = true;

    typedef typename super::value_type value_type;
    typedef typename super::primal_space_type primal_space_type;
    typedef typename super::dual_space_type dual_space_type;
    typedef Continuous continuity_type;
    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename super::polyset_type polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;


    typedef Hermite<N, O, Scalar, T, Convex,  Pts> component_basis_type;

    typedef typename mpl::if_<mpl::equal_to<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<boost::none_t>,
            mpl::identity< Hermite<N-1, O, Scalar, T, Convex,  Pts> > >::type::type face_basis_type;

    typedef boost::shared_ptr<face_basis_type> face_basis_ptrtype;

    typedef typename dual_space_type::convex_type convex_type;
    typedef typename dual_space_type::pointset_type pointset_type;
    typedef typename dual_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename reference_convex_type::points_type points_type;

    static const uint16_type numPoints = reference_convex_type::numPoints;
    static const uint16_type nbPtsPerVertex = reference_convex_type::nbPtsPerVertex;
    static const uint16_type nbPtsPerEdge = reference_convex_type::nbPtsPerEdge;
    static const uint16_type nbPtsPerFace = reference_convex_type::nbPtsPerFace;
    static const uint16_type nbPtsPerVolume = reference_convex_type::nbPtsPerVolume;

    static const uint16_type nLocalDof = super::nLocalDof;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Hermite()
        :
        super( dual_space_type( primal_space_type() ) ),
        M_refconvex()
        // M_bdylag( new face_basis_type )
    {



    }

    virtual ~Hermite() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the reference convex associated with the hermite polynomials
     */
    reference_convex_type const& referenceConvex() const
    {
        return M_refconvex;
    }

    /**
     * \return the family name of the finite element
     */
    std::string familyName() const
    {
        return "hermite";
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{
    template<typename GMContext, typename PC, typename Phi, typename GPhi, typename HPhi >
    static void transform( boost::shared_ptr<GMContext> gmc,  boost::shared_ptr<PC> const& pc,
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
        typename GMContext::gm_type::matrix_type const& B = gmc.B( 0 );
        std::fill( phi_t.data(), phi_t.data()+phi_t.num_elements(), value_type( 0 ) );

        if ( do_gradient )
            std::fill( g_phi_t.data(), g_phi_t.data()+g_phi_t.num_elements(), value_type( 0 ) );

        if ( do_hessian )
            std::fill( h_phi_t.data(), h_phi_t.data()+g_phi_t.num_elements(), value_type( 0 ) );

        const uint16_type Q = gmc.nPoints();//M_grad.size2();

        // transform
        for ( uint16_type i = numPoints; i < nLocalDof; i+=nDim )
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
                        phi_t[i+l][0][0][q] += B( l, p ) * pc.phi( i+p,0,0,q );

                        if ( do_gradient )
                        {
                            for ( uint16_type j = 0; j < nDim; ++j )
                            {
                                g_phi_t[i+l][0][p][q] += B( l, j ) * pc.grad( i+j,0,p,q );
                            }
                        }

                        if ( do_hessian )
                        {
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
};
template<uint16_type N,
         uint16_type O,
         template<uint16_type Dim> class PolySetType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex,
         template<class, uint16_type, class> class Pts >
const uint16_type Hermite<N,O,PolySetType,T,Convex,Pts>::nDim;

template<uint16_type N,
         uint16_type O,
         template<uint16_type Dim> class PolySetType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex,
         template<class, uint16_type, class> class Pts >
const uint16_type Hermite<N,O,PolySetType,T,Convex,Pts>::nOrder;

template<uint16_type N,
         uint16_type O,
         template<uint16_type Dim> class PolySetType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex,
         template<class, uint16_type, class> class Pts >
const uint16_type Hermite<N,O,PolySetType,T,Convex,Pts>::numPoints;

} // namespace fem

template<uint16_type Order,
         template<uint16_type Dim> class PolySetType = Scalar,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced>
class Hermite
{
public:
    template<uint16_type N,
             typename T = double,
             typename Convex = Simplex<N> >
    struct apply
    {
        typedef typename mpl::if_<mpl::bool_<Convex::is_simplex>,
                mpl::identity<fem::Hermite<N,Order,PolySetType,T,Simplex, Pts> >,
                mpl::identity<fem::Hermite<N,Order,PolySetType,T,Hypercube, Pts> > >::type::type result_type;
        typedef result_type type;
    };

    typedef Hermite<Order,Scalar,Pts> component_basis_type;

    static const uint16_type nOrder =  Order;

};

} // namespace Feel
#endif /* __hermite_H */
