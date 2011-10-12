/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-01-12

  Copyright (C) 2006 EPFL

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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

namespace Feel
{
namespace fem
{
namespace detail
{
template<typename Basis, template<class, uint16_type, class> class PointSetType>
class CrouzeixRaviartDual
    :
        public DualBasis<Basis>
{
    typedef DualBasis<Basis> super;
public:

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder= 1;

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
    _M_convex_ref(),
    _M_eid(_M_convex_ref.topologicalDimension()+1),
    _M_pts( nDim, numPoints ),
    _M_points_face( nFacesInConvex ),
    _M_fset( primal )
{
#if 1
    std::cout << "Lagrange finite element: \n";
    std::cout << " o- dim   = " << nDim << "\n";
    std::cout << " o- order = " << nOrder << "\n";
    std::cout << " o- numPoints      = " << numPoints << "\n";
    std::cout << " o- nbPtsPerVertex = " << (int)nbPtsPerVertex << "\n";
    std::cout << " o- nbPtsPerEdge   = " << (int)nbPtsPerEdge << "\n";
    std::cout << " o- nbPtsPerFace   = " << (int)nbPtsPerFace << "\n";
    std::cout << " o- nbPtsPerVolume = " << (int)nbPtsPerVolume << "\n";
#endif
    equispaced_pointset_type epts;
    // in d-dimension, consider only the d-1 entity mid-points
    int d = _M_convex_ref.topologicalDimension()-1;
    int p = 0;
    // loop on each entity forming the convex of topological
    // dimension d
    for ( int e = _M_convex_ref.entityRange( d ).begin();
          e < _M_convex_ref.entityRange( d ).end();
          ++e )
        {
            _M_points_face[e] = epts.pointsBySubEntity(nDim-1, e, 0);
            ublas::subrange( _M_pts, 0, nDim, p, p+_M_points_face[e].size2() ) = _M_points_face[e];
            p+=_M_points_face[e].size2();
        }
     std::cout << "[CrouzeixRaviartDual] points= " << _M_pts << "\n";
    setFset( primal, _M_pts, mpl::bool_<primal_space_type::is_scalar>() );


}

    points_type const& points() const { return _M_pts; }
    points_type const& points(uint16_type f ) const { return _M_points_face[f]; }


matrix_type operator()( primal_space_type const& pset ) const
{
    return _M_fset( pset );
}
private:

void setFset( primal_space_type const& primal, points_type const& __pts, mpl::bool_<true> )
{
    _M_fset.setFunctionalSet( functional::PointsEvaluation<primal_space_type>( primal,
                                                                               __pts ) );
}

void setFset( primal_space_type const& primal, points_type const& __pts, mpl::bool_<false> )
{
    _M_fset.setFunctionalSet( functional::ComponentsPointsEvaluation<primal_space_type>( primal,
                                                                                         __pts ) );
}


private:
reference_convex_type _M_convex_ref;
std::vector<std::vector<uint16_type> > _M_eid;
points_type _M_pts;
std::vector<points_type> _M_points_face;
FunctionalSet<primal_space_type> _M_fset;


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
         template<uint16_type, uint16_type, uint16_type> class Convex = Simplex>
class CrouzeixRaviart
    :
    public FiniteElement<Feel::detail::OrthonormalPolynomialSet<N, 1, RealDim, PolySetType, T, Convex>,
                         detail::CrouzeixRaviartDual,
                         PointSetEquiSpaced >
{
    typedef FiniteElement<Feel::detail::OrthonormalPolynomialSet<N, 1, RealDim, PolySetType, T, Convex>,
                          detail::CrouzeixRaviartDual,
                          PointSetEquiSpaced > super;
public:

    BOOST_STATIC_ASSERT( N > 1 );

    /** @name Typedefs
     */
    //@{

    static const uint16_type nDim = N;
    static const uint16_type nOrder =  1;
    static const bool isTransformationEquivalent = true;
    static const bool isContinuous = true;

    typedef typename super::value_type value_type;
    typedef typename super::primal_space_type primal_space_type;
    typedef typename super::dual_space_type dual_space_type;
    typedef Continuous continuity_type;

    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename super::polyset_type polyset_type;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const bool is_product = true;

    typedef CrouzeixRaviart<N, RealDim, Scalar, T, Convex> component_basis_type;

    typedef typename dual_space_type::convex_type convex_type;
    typedef typename dual_space_type::pointset_type pointset_type;
    typedef typename dual_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename reference_convex_type::points_type points_type;


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
    //@}

    /** @name Constructors, destructor
     */
    //@{

    CrouzeixRaviart()
        :
        super( dual_space_type( primal_space_type() ) ),
        _M_refconvex()
    {}
    CrouzeixRaviart( CrouzeixRaviart const & cr )
        :
        super( cr ),
        _M_refconvex()
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
    reference_convex_type const& referenceConvex() const { return _M_refconvex; }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * \return the family name of the finite element
     */
    std::string familyName() const { return "CrouzeixRaviart"; }

    template<typename ExprType>
    static auto
    isomorphism( ExprType& expr ) -> decltype( expr )
        {
            return expr;
            //return expr;
        }
    //@}



protected:
reference_convex_type _M_refconvex;
private:

};
template<uint16_type N,
         uint16_type RealDim,
         template<uint16_type Dim> class PolySetType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex>
const uint16_type CrouzeixRaviart<N,RealDim,PolySetType,T,Convex>::nDim;
template<uint16_type N,
         uint16_type RealDim,
         template<uint16_type Dim> class PolySetType,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex>
const uint16_type CrouzeixRaviart<N,RealDim,PolySetType,T,Convex>::nOrder;

} // fem

template<uint16_type Order,
         template<uint16_type Dim> class PolySetType = Scalar,
         template<class, uint16_type, class> class Pts = PointSetEquiSpaced>
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
                                  mpl::identity<fem::CrouzeixRaviart<N,RealDim,PolySetType,T,Simplex> >,
                                  mpl::identity<fem::CrouzeixRaviart<N,RealDim,PolySetType,T,Hypercube> > >::type::type result_type;
        typedef result_type type;
    };

    typedef CrouzeixRaviart<Order,Scalar,Pts> component_basis_type;
    static const uint16_type nOrder =  Order;
};

template<uint16_type Order,
         template<uint16_type Dim> class PolySetType,
         template<class, uint16_type, class> class Pts>
const uint16_type CrouzeixRaviart<Order,PolySetType,Pts>::nOrder;


} // Feel
#endif /* __CrouzeixRaviart_H */
