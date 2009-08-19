/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-01-14

  Copyright (C) 2006 EPFL
  Copyright (C) 2008, 2009 Université de Grenoble 1 (Joseph Fourier)

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
   \file raviartthomas.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-01-14
 */
#ifndef __RaviartThomas_H
#define __RaviartThomas_H 1

#include <boost/ptr_container/ptr_vector.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/traits.hpp>
#include <life/lifealg/lu.hpp>

#include <life/lifemesh/refentity.hpp>
#include <life/lifemesh/pointset.hpp>


#include <life/lifepoly/dualbasis.hpp>
#include <life/lifepoly/polynomialset.hpp>
#include <life/lifepoly/functionalset.hpp>
#include <life/lifepoly/operations.hpp>
#include <life/lifepoly/functionals.hpp>
#include <life/lifepoly/functionals2.hpp>
#include <life/lifepoly/quadpoint.hpp>
#include <life/lifepoly/fe.hpp>
namespace Life
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
        _M_p ( p ),
        _M_c( c )
    {
        std::cout << "component : " << c << std::endl;
    }
    typename ublas::vector<value_type> operator() ( points_type const& __pts ) const
    {
        std::cout << "times_x(pts) : " << __pts << std::endl;
        std::cout << "times_x(pts) : " << _M_p.evaluate(__pts) << std::endl;
        std::cout << "times_x(coeff) : " << _M_p.coefficients() << std::endl;

        // __pts[c] * p( __pts )
        return ublas::element_prod( ublas::row( __pts, _M_c ),
                                    ublas::row( _M_p.evaluate( __pts ), 0 ) );
    }
    //P const& _M_p;
    P _M_p;
    int _M_c;
};

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
}// detail

template<uint16_type N,
         uint16_type O,
         typename T = double,
         template<uint16_type, uint16_type, uint16_type> class Convex = Simplex>
class RaviartThomasPolynomialSet
    :
        public detail::OrthonormalPolynomialSet<N, O+1, Vectorial, T, Convex>
{
    typedef detail::OrthonormalPolynomialSet<N, O+1, Vectorial, T, Convex> super;

public:
    typedef detail::OrthonormalPolynomialSet<N, O, Vectorial, T, Convex> Pk_v_type;
    typedef detail::OrthonormalPolynomialSet<N, O+1, Vectorial, T, Convex> Pkp1_v_type;
    typedef detail::OrthonormalPolynomialSet<N, O-1, Vectorial, T, Convex> Pkm1_v_type;
    typedef detail::OrthonormalPolynomialSet<N, O, Scalar, T, Convex> Pk_s_type;
    typedef detail::OrthonormalPolynomialSet<N, O+1, Scalar, T, Convex> Pkp1_s_type;

    typedef PolynomialSet<typename super::basis_type,Vectorial> vectorial_polynomialset_type;
    typedef typename vectorial_polynomialset_type::polynomial_type vectorial_polynomial_type;
    typedef PolynomialSet<typename super::basis_type,Scalar> scalar_polynomialset_type;
    typedef typename scalar_polynomialset_type::polynomial_type scalar_polynomial_type;

    typedef RaviartThomasPolynomialSet<N, O, T> self_type;

    typedef typename super::value_type value_type;
    typedef typename super::convex_type convex_type;
    typedef typename super::matrix_type matrix_type;
    typedef typename super::points_type points_type;

    static const uint16_type nDim = super::nDim;
    static const uint16_type nOrder = super::nOrder;
    static const uint16_type nComponents = super::nComponents;

    RaviartThomasPolynomialSet()
        :
        super()
    {
        uint16_type dim_Pkp1 = convex_type::polyDims( nOrder );
        uint16_type dim_Pk = convex_type::polyDims( nOrder-1 );
        uint16_type dim_Pkm1 = (nOrder==1)?0:convex_type::polyDims( nOrder-2 );
        std::cout << "[RTPset] dim_Pkp1 = " << dim_Pkp1 << "\n";
        std::cout << "[RTPset] dim_Pk   = " << dim_Pk << "\n";
        std::cout << "[RTPset] dim_Pkm1 = " << dim_Pkm1 << "\n";

        // (P_k)^d
        Pkp1_v_type Pkp1_v;
        vectorial_polynomialset_type Pk_v( Pkp1_v.polynomialsUpToDimension( dim_Pk ) );
        std::cout << "[RTPset] Pk_v =" << Pk_v.coeff() << "\n";

        // P_k
        Pkp1_s_type Pkp1;
        scalar_polynomialset_type Pk ( Pkp1.polynomialsUpToDimension( dim_Pk ) );
        std::cout << "[RTPset] Pk =" << Pk.coeff() << "\n";
        std::cout << "[RTPset] Pk(0) =" << Pk.polynomial( 0 ).coefficients() << "\n";

        // x P_k \ P_{k-1}
        IMGeneral<convex_type::nDim, 2*nOrder,value_type> im;
        std::cout << "[RTPset] im.points() = " << im.points() << std::endl;
        ublas::matrix<value_type> xPkc( nComponents*(dim_Pk-dim_Pkm1),Pk.coeff().size2() );

        std::cout << "[RTPset] before xPkc = " << xPkc << "\n";
        for( int l = dim_Pkm1, i = 0; l < dim_Pk; ++l, ++i )
            {
                for( int j = 0; j < convex_type::nDim; ++j )
                    {
                        detail::times_x<scalar_polynomial_type> xp( Pk.polynomial( l ), j );
                        ublas::row(xPkc,i*nComponents+j)=
                            ublas::row( Life::project( Pkp1,
                                                       xp,
                                                       im ).coeff(), 0);
                    }
            }


        std::cout << "[RTPset] after xPkc = " << xPkc << "\n";
        vectorial_polynomialset_type xPk( typename super::basis_type(), xPkc, true );
        std::cout << "[RTPset] here 1\n";
        // (P_k)^d + x P_k
        std::cout << "[RTPset] RT Poly coeff = " << unite( Pk_v, xPk ).coeff() << "\n";
        this->setCoefficient( unite( Pk_v, xPk ).coeff(), true );
        std::cout << "[RTPset] here 2\n";
    }


};

namespace fem
{

namespace detail
{



template<typename Basis,
         template<class, uint16_type, class> class PointSetType>
class RaviartThomasDual
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
    typedef typename primal_space_type::template convex<nDim+nOrder-1>::type convex_type;
    typedef Reference<convex_type, nDim, nDim+nOrder-1, nDim, value_type> reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;

    typedef typename primal_space_type::Pkp1_v_type Pkp1_v_type;
    typedef typename primal_space_type::vectorial_polynomialset_type vectorial_polynomialset_type;

    // point set type associated with the functionals
    typedef PointSetType<convex_type, nOrder, value_type> pointset_type;

    static const uint16_type nbPtsPerVertex = 0;
    static const uint16_type nbPtsPerEdge = mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<2> >,mpl::int_<reference_convex_type::nbPtsPerEdge>,mpl::int_<0> >::type::value;
    static const uint16_type nbPtsPerFace =mpl::if_<mpl::equal_to<mpl::int_<nDim>,mpl::int_<3> >,mpl::int_<reference_convex_type::nbPtsPerFace>,mpl::int_<0> >::type::value;
    static const uint16_type nbPtsPerVolume = 0;
    static const uint16_type numPoints = ( reference_convex_type::numGeometricFaces*nbPtsPerFace+reference_convex_type::numEdges*nbPtsPerEdge );

    /** Number of degrees of freedom per vertex */
    static const uint16_type nDofPerVertex = 0;

    /** Number of degrees of freedom per edge */
    static const uint16_type nDofPerEdge = nbPtsPerEdge;

    /** Number of degrees of freedom per face */
    static const uint16_type nDofPerFace = nbPtsPerFace;

    /** Number of degrees  of freedom per volume */
    static const uint16_type nDofPerVolume = 0;

    /** Total number of degrees of freedom (equal to refEle::nDof) */
    static const uint16_type nLocalDof = numPoints;

    RaviartThomasDual( primal_space_type const& primal )
        :
        super( primal ),
        _M_convex_ref(),
        _M_eid(_M_convex_ref.topologicalDimension()+1),
        _M_pts( nDim, numPoints ),
        _M_fset( primal )
    {
#if 1
        std::cout << "Raviart-Thomas finite element(dual): \n";
        std::cout << " o- dim   = " << nDim << "\n";
        std::cout << " o- order = " << nOrder << "\n";
        std::cout << " o- numPoints      = " << numPoints << "\n";
        std::cout << " o- nbPtsPerVertex = " << (int)nbPtsPerVertex << "\n";
        std::cout << " o- nbPtsPerEdge   = " << (int)nbPtsPerEdge << "\n";
        std::cout << " o- nbPtsPerFace   = " << (int)nbPtsPerFace << "\n";
        std::cout << " o- nbPtsPerVolume = " << (int)nbPtsPerVolume << "\n";
        std::cout << " o- nLocalDof      = " << nLocalDof << "\n";
#endif
        std::vector<points_type> pts_per_face( convex_type::numTopologicalFaces );
        // loop on each entity forming the convex of topological
        // dimension nDim-1 ( the faces)
        for ( int p = 0, e = _M_convex_ref.entityRange( nDim-1 ).begin();
              e < _M_convex_ref.entityRange( nDim-1 ).end();
              ++e )
            {
                points_type Gt ( _M_convex_ref.makePoints( nDim-1, e ) );
                pts_per_face[e] =  Gt ;
                if ( Gt.size2() )
                    {
                        //Debug() << "Gt = " << Gt << "\n";
                        //Debug() << "p = " << p << "\n";
                        ublas::subrange( _M_pts, 0, nDim, p, p+Gt.size2() ) = Gt;
                        //for ( size_type j = 0; j < Gt.size2(); ++j )
                        //_M_eid[d].push_back( p+j );
                        p+=Gt.size2();
                    }
            }
        std::cout << "[RT Dual] done 1\n";
        // compute  \f$ \ell_e( U ) = (U * n[e]) (edge_pts(e)) \f$
        typedef Functional<primal_space_type> functional_type;
        std::vector<functional_type> fset;
        double j[3] = {2.8284271247461903,2.0,2.0};

        //for( int k = 0; k < nDim; ++k )
            {
                // loopover the each edge entities and add the correponding functionals
                for ( int e = _M_convex_ref.entityRange( nDim-1 ).begin();
                      e < _M_convex_ref.entityRange( nDim-1 ).end();
                      ++e )
                    {
                        typedef Life::functional::DirectionalComponentPointsEvaluation<primal_space_type> dcpe_type;

                        dcpe_type __dcpe( primal, 1, _M_convex_ref.normal(e)*j[e], pts_per_face[e] );
                        std::copy( __dcpe.begin(), __dcpe.end(), std::back_inserter( fset ) );
                    }
            }
        std::cout << "[RT Dual] done 2" << std::endl;
        if ( nOrder-1 > 0 )
            {
                // we need more equations : add interior moment
                // indeed the space is orthogonal to Pk-1
                uint16_type dim_Pkp1 = convex_type::polyDims( nOrder );
                uint16_type dim_Pk = convex_type::polyDims( nOrder-1 );
                uint16_type dim_Pm1 = convex_type::polyDims( nOrder-2 );

                Pkp1_v_type Pkp1;

                vectorial_polynomialset_type Pkm1 ( Pkp1.polynomialsUpToDimension( dim_Pm1 ) );
                for( int i = 0; i < Pkm1.polynomialDimension(); ++i )
                    {
                        typedef functional::IntegralMoment<primal_space_type, vectorial_polynomialset_type> fim_type;
                        fset.push_back( fim_type( primal, Pkm1.polynomial( i ) ) );
                    }
            }
        std::cout << "[RT Dual] done 3, n fset = " << fset.size() << std::endl;
        _M_fset.setFunctionalSet( fset );
        std::cout << "[RT Dual] done 4\n";

    }

    points_type const& points() const { return _M_pts; }


    matrix_type operator()( primal_space_type const& pset ) const
    {
        std::cout << "RT matrix = " << _M_fset( pset ) << std::endl;
        return _M_fset( pset );
    }

private:
    reference_convex_type _M_convex_ref;
    std::vector<std::vector<uint16_type> > _M_eid;
    points_type _M_pts;
    FunctionalSet<primal_space_type> _M_fset;


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
         template<uint16_type, uint16_type, uint16_type> class Convex = Simplex>
class RaviartThomas
    :
        public FiniteElement<RaviartThomasPolynomialSet<N, O, T, Convex>,
                             detail::RaviartThomasDual,
                             PointSetEquiSpaced >
{
    typedef FiniteElement<RaviartThomasPolynomialSet<N, O, T, Convex>,
                          detail::RaviartThomasDual,
                          PointSetEquiSpaced > super;
public:

    BOOST_STATIC_ASSERT( N > 1 );

    /** @name Typedefs
     */
    //@{

    static const uint16_type nDim = N;
    static const polynomial_transformation_type transformation = POLYNOMIAL_CONTEXT_NEEDS_1ST_PIOLA_TRANSFORMATION;
    typedef typename super::value_type value_type;
    typedef typename super::primal_space_type primal_space_type;
    typedef typename super::dual_space_type dual_space_type;

    /**
     * Polynomial Set type: scalar or vectorial
     */
    typedef typename super::polyset_type polyset_type;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;


    typedef typename dual_space_type::convex_type convex_type;
    typedef typename dual_space_type::pointset_type pointset_type;
    typedef typename dual_space_type::reference_convex_type reference_convex_type;
    typedef typename reference_convex_type::node_type node_type;
    typedef typename reference_convex_type::points_type points_type;


    static const uint16_type nOrder =  dual_space_type::nOrder;
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

    RaviartThomas()
        :
        super( dual_space_type( primal_space_type() ) ),
        _M_refconvex()
    {
        std::cout << "[RT] nPtsPerEdge = " << nbPtsPerEdge << "\n";
        std::cout << "[RT] nPtsPerFace = " << nbPtsPerFace << "\n";
        std::cout << "[RT] numPoints = " << numPoints << "\n";

        std::cout << "[RT] nDof = " << super::nDof << "\n";

        std::cout << "[RT] coeff : " << this->coeff() << "\n";
        std::cout << "[RT] pts : " << this->points() << "\n";
        std::cout << "[RT] eval at pts : " << this->evaluate( this->points() ) << "\n";
    }
    RaviartThomas( RaviartThomas const & cr )
        :
        super( cr ),
        _M_refconvex()
    {}
    ~RaviartThomas()
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

    std::string familyName() const { return "raviartthomas"; }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:
    reference_convex_type _M_refconvex;
private:

};
template<uint16_type N,
         uint16_type O,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex>
const uint16_type RaviartThomas<N,O,T,Convex>::nDim;
template<uint16_type N,
         uint16_type O,
         typename T,
         template<uint16_type, uint16_type, uint16_type> class Convex>
const uint16_type RaviartThomas<N,O,T,Convex>::nOrder;

} // fem
template<uint16_type Order>
class RaviartThomas
{
public:
    template<uint16_type N,
             typename T = double,
             typename Convex = Simplex<N> >
    struct apply
    {
        typedef typename mpl::if_<mpl::bool_<Convex::is_simplex>,
                                  mpl::identity<fem::RaviartThomas<N,Order,T,Simplex> >,
                                  mpl::identity<fem::RaviartThomas<N,Order,T,SimplexProduct> > >::type::type result_type;
        typedef result_type type;
    };

    typedef Lagrange<Order,Scalar> component_basis_type;
};

} // Life
#endif /* __RaviartThomas_H */

