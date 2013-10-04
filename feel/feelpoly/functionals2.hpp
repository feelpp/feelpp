/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-01-22

  Copyright (C) 2006 EPFL
  Copyright (C) 2009 Université de Grenoble 1 (Joseph Fourier)

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
   \file functionals2.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-01-22
 */
#ifndef __FEELPP_FUNCTIONALS2_HPP
#define __FEELPP_FUNCTIONALS2_HPP 1

#include <feel/feelpoly/functional.hpp>

#include <feel/feelpoly/quadpoint.hpp>
#include <feel/feelpoly/im.hpp>

namespace Feel
{

namespace functional
{
/**
 * \class IntegralMoment
 * \brief functional that returns \f$\ell_u (v) = \int_\Omega( u\, v )\f$
 *
 * if the basis functions in which $u$ and $v$ are orthonormal then
 * the integral is simply the inner product of the coefficients of $u$
 * and $v$ thanks to Parseval identity.
 *
 * \author Christophe Prud'homme
 */
template<typename Space, typename Poly = Space>
class IntegralMoment
    :
public Functional<Space>
{
    typedef Functional<Space> super;
public:

    typedef IntegralMoment<Space,Poly> self_type;
    typedef typename super::space_type space_type;
    typedef typename Poly::polynomial_type polynomial_type;
    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;

    //BOOST_STATIC_ASSERT( ( boost::is_same<typename space_type::basis_type, typename Poly::basis_type>::value ) );

    IntegralMoment()
        :
        super()
    {}

    IntegralMoment( space_type const& p, polynomial_type const& q )
        :
        super( p )
    {
#if 0
        std::cout << "p.coeff = " << p.coeff() << "\n";
        std::cout << "p.basis.coeff = " << p.basis().coeff() << "\n";
        std::cout << "q.coeff = " << q.coeff() << "\n";
#endif
        //ublas::matrix<value_type> m ( ublas::prod( q.coeff(), ublas::trans( p.basis().coeff() ) ) );
        typename space_type::Pkp1_v_type l;
        ublas::matrix<value_type> m ( ublas::prod( q.coeff(), ublas::trans( l.coeff() ) ) );

        //ublas::matrix<value_type> m ( ublas::prod( q.coeff(), p.coeff() ) );
        //std::cout << "[IntegralMoment] m = " << m << "\n";
        this->setCoefficient( m );
    }

    IntegralMoment( IntegralMoment const& im )
        :
        super( im )
    {}

    IntegralMoment& operator=( IntegralMoment const& im )
    {
        if ( this != &im )
        {
            super::operator=( im );
        }

        return *this;
    }

};

namespace detail
{
template<typename P1, typename P2>
struct prod
{
    typedef typename P1::value_type value_type;
    prod ( P1 const& p1, P2 const& p2 )
        :
        M_p1 ( p1 ),
        M_p2 ( p2 )
    {}
#if 0
    value_type operator() ( typename node<value_type>::type const& n )
    {
        return M_p1.evaluate( n )( 0,0 ) * M_p2.evaluate( n )( 0,0 );
    }
#endif
    typename node<value_type>::type operator() ( typename node<value_type>::type const& n )
    {
        return ublas::column( ublas::element_prod( M_p1.evaluate( n ), M_p2.evaluate( n ) ), 0 );
    }
    P1 M_p1;
    P2 M_p2;
};
} // detail
/**
 * \class IntegralMomentOnFace
 * \brief functional that returns \f$\ell_u (v) = \int_{\Gamma} ( u\, v )\f$ where \f$\Gamma \subset \partial \Omega\f$
 *
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
class IntegralMomentOnFace
    :
public Functional<Space>
{
    typedef Functional<Space> super;
public:

    typedef IntegralMomentOnFace<Space> self_type;
    typedef typename super::space_type space_type;
    typedef typename super::polynomial_type polynomial_type;
    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;

    const static uint16_type nComponents = space_type::nComponents;

    /**
     * Construct the functional
     *
     * \param p polynomial space on which we apply the functional
     * \param k index of the polynomial to integrate against
     * \param face face of the convex over which to integrate
     */
    IntegralMomentOnFace( space_type const& p, uint16_type k, IntegrationFaceEnum face = ALL_FACES )
        :
        super( p ),
        M_k ( k ),
        M_q()
    {
        ublas::matrix<value_type> __rep( nComponents, p.polynomialDimensionPerComponent() );
        typedef detail::prod<typename space_type::polynomial_type, polynomial_type> prod_fun;

        for ( uint16_type i = 0; i < p.polynomialDimensionPerComponent(); ++i )
        {
            //if ( p.is_scalar )
            //__rep( 0, i ) = M_q.integrate( face, prod_fun( p.polynomial(i), p.polynomial( k ) ) );
            //else
            typedef typename node<value_type>::type node_type;
            ublas::column( __rep, i ) = M_q.integrate( face, prod_fun( p.polynomial( i ), p.polynomial( k ) ) );
        }

        this->setCoefficient( __rep );
    }

    /**
     * Construct the functional
     *
     * \param p polynomial space on which we apply the functional
     * \param k index of the polynomial to integrate against
     * \param c component index
     * \param face face of the convex over which to integrate
     */
    IntegralMomentOnFace( space_type const& p, uint16_type k, uint16_type c, IntegrationFaceEnum face = ALL_FACES )
        :
        super( p ),
        M_k ( k ),
        M_q()
    {
        ublas::matrix<value_type> __rep( ublas::zero_matrix<value_type>( nComponents, p.polynomialDimensionPerComponent() ) );
        typedef detail::prod<typename space_type::polynomial_type::component_type,
                typename polynomial_type::component_type> prod_fun;

        int nc = p.polynomialDimensionPerComponent()*c;
        int ind_p2 = nc + k;

        for ( uint16_type i = 0; i < p.polynomialDimensionPerComponent(); ++i )
        {
            int ind_p1 = nc + i;
            typedef typename node<value_type>::type node_type;
            __rep( c, i ) = M_q.integrate( face, prod_fun( p.polynomial( ind_p1 )[c],
                                            p.polynomial( ind_p2 )[c] ) )( 0 );
        }

        this->setCoefficient( __rep );
    }

private:

    // disabled
    IntegralMomentOnFace();

private:

    // polynomial degree integrate against
    uint16_type M_k;

    // quadrature rule on the element and faces of the element
    IM<Space::nDim,2*Space::nOrder+1, value_type> M_q;
};


/**
 * \class IntegralMomentsOnFace
 * \brief functional that returns \f$\ell_u (v) = \int_{\Gamma} ( u\, v )\f$ where \f$\Gamma \subset \partial \Omega\f$
 *
 *
 * \author Christophe Prud'homme
 */
template<typename Space,typename BasisType>
class IntegralMomentsOnFace
    :
public std::vector<Functional<Space> >
{
    typedef std::vector<Functional<Space> > super;
public:

    typedef Functional<Space> functional_type;
    typedef IntegralMomentsOnFace<Space,BasisType> self_type;
    typedef BasisType basis_type;
    typedef Space space_type;
    typedef typename space_type::reference_convex_type reference_convex_type;
    typedef typename super::polynomial_type polynomial_type;
    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;

    const static uint16_type nComponents = space_type::nComponents;

    /**
     * Construct the functional
     *
     * \param p polynomial space on which we apply the functional
     * \param k index of the polynomial to integrate against
     * \param face face of the convex over which to integrate
     */
    IntegralMomentsOnFace( space_type const& p,
                           basis_type const& l,
                           IntegrationFaceEnum face = ALL_FACES )
        :
        super()
    {
        reference_convex_type ref_convex;
        typedef typename reference_convex_type::topological_face_type  element_type;
        element_type ref_convex_face = ref_convex.topologicalFace( face );

        typedef GeoMap<reference_convex_type::nDim,1,reference_convex_type::nDim> gm_type;
        typedef typename gm_type::face_gm_type::precompute_type face_pc_type;
        typedef typename gm_type::face_gm_type::precompute_ptrtype face_pc_ptrtype;
        gm_type __gm;
        IM<reference_convex_type::nDim-1,2*space_type::nOrder-1> __qr_face;
        face_pc_ptrtype __geopc( new face_pc_type( __gm->boundaryMap(),__qr_face.points() ) );

        DVLOG(2) << "[nc] ref_convex_face "  << face << "=" << ref_convex_face.points() << "\n";


        typename gm_type::template Context<vm::POINT,element_type> __c( __gm->boundaryMap(),
                ref_convex_face,
                __geopc );

        __c.update( ref_convex_face, __geopc );
        DVLOG(2) << "[nc] ref_convex_face "  << face << " xref" << __c.xRefs() << "\n";
        DVLOG(2) << "[nc] ref_convex_face "  << face << " xreal" << __c.xReal() << "\n";

        for ( uint16_type k = 0; k < l.polynomialDimensionPerComponent(); ++k )
        {
            ublas::matrix<value_type> __rep( nComponents, p.polynomialDimensionPerComponent() );

            for ( uint16_type i = 0; i < p.polynomialDimensionPerComponent(); ++i )
            {
                /*
                typedef typename node<value_type>::type node_type;
                double __res = 0;
                double __len = 0;
                for ( uint16_type __ip = 0; __ip < __qr_face.nPoints();++__ip )
                    {
                        __res += ( __qr_face.weight( __ip )*
                                   p.polynomial(i).evaluate( __c.xReal(__ip) )*
                                   l.polynomial(k).evaluate( __c.xRef( __ip ) ));

                        __len += __qr_face.weight( __ip );
                    }
                          DVLOG(2) << "[nc] length = " << __len << "\n";
                      DVLOG(2) << "[nc] res = " << __res << "\n";
                      ublas::column( __rep, i ) = __res/__len;
                      }
                this->push_back( functional_type( p,  __rep ) );
                */
            }

        }
    }

private:

    // disabled
    IntegralMomentsOnFace();

private:


};

/**
 * \class IntegralMomentOfDerivative
 * \brief functional that returns \f$\ell_u^i (v) = \int_\Omega( \frac{d u}{d x_i} \, v )\f$
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
class IntegralMomentOfDerivative
    :
public Functional<Space>
{
    typedef Functional<Space> super;
public:

    typedef IntegralMoment<Space> self_type;
    typedef typename super::space_type space_type;
    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;

    IntegralMomentOfDerivative()
        :
        super()
    {}
    template<typename P>
    IntegralMomentOfDerivative( space_type const& b, uint16_type i )
        :
        super( b, b.d( i ) )
    {}
};

/**
 * \class IntegralMomentOfDivergence
 * \brief functional that returns \f$\ell_u^i (v) = \int_\Omega( \frac{d u}{d x_i} \, v )\f$
 *
 * \author Christophe Prud'homme
 */
template<typename Space, typename Poly = Space>
class IntegralMomentOfDivergence
    :
public Functional<Space>
{
    typedef Functional<Space> super;
public:

    BOOST_STATIC_ASSERT( Space::is_vectorial );

    typedef IntegralMoment<Space> self_type;
    typedef typename Poly::polynomial_type polynomial_type;
    typedef typename super::space_type space_type;
    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;

    IntegralMomentOfDivergence()
        :
        super()
    {}

    IntegralMomentOfDivergence( space_type const& p, polynomial_type const& q )
        :
        super( p )
    {
        ublas::matrix<value_type> __rep( space_type::nComponents, p.polynomialDimensionPerComponent() );

        //std::cout << "q = " << q.coeff() << "\n";
        for ( int i = 0; i < space_type::nComponents; ++i )
        {

            std::cout << "p.d("<< i << ")  = " << p.basis().d( i ) << "\n"
                      << "prod = " << ublas::prod( q.coeff(), p.basis().d( i ) ) << "\n";

            ublas::row( __rep,i ) = ublas::row( ublas::prod( q.coeff(), p.basis().d( i ) ), 0 ) ;
        }

        //std::cout << "[IntegralMomentOfDivergence] rep = " << __rep << "\n";
        this->setCoefficient( __rep );
    }
};

} // functional

} // Feel

#endif // __FEELPP_FUNCTIONALS2_HPP
