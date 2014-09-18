/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-01-22

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
   \file functionals.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-01-22
 */
#ifndef __FEELPP_FUNCTIONALS_HPP
#define __FEELPP_FUNCTIONALS_HPP 1

#include <feel/feelpoly/functional.hpp>

namespace Feel
{

namespace functional
{
/**
 * \class PointEvaluation
 * \brief generate the functional associated with a point evaluation
 *
 * Given a function space \f$ E \f$, generate the functional \f$ \ell
 * : E \longrightarrow R \f$ associated with the evaluation
 * of the basis functions of the function space at a point in the
 * geometric space.
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
class PointEvaluation
    :
public Functional<Space>
{
    typedef Functional<Space> super;
public:

    typedef PointEvaluation<Space> self_type;
    typedef typename super::space_type space_type;
    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;

    PointEvaluation()
        :
        super()
    {}
    PointEvaluation( space_type const& b, node_type const& __pt )
        : super( b, ublas::trans( b.basis()( __pt ) ) )
    {
        //std::cout << "[PointEvaluation] eval = " << b.evaluate( __pt ) << "\n";
    }
};

/**
 * \class ComponentPointEvaluation
 * \brief generate the functional associated with a point evaluation associated with a component
 *
 * Given a function space \f$ E \f$, generate the functional \f$ \ell
 * : E \longrightarrow R \f$ associated with the evaluation
 * of the basis functions of the function space at a point in the
 * geometric space.
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
class ComponentPointEvaluation
    :
public Functional<Space>
{
    typedef Functional<Space> super;
public:

    typedef ComponentPointEvaluation<Space> self_type;
    typedef typename super::space_type space_type;
    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;

    ComponentPointEvaluation()
        :
        super()
    {}
    ComponentPointEvaluation( space_type const& b, uint16_type __c, node_type const& __pt )
        :
        super( b )
    {
        ublas::matrix<value_type> m( ublas::zero_matrix<value_type>( space_type::nComponents, b.polynomialDimensionPerComponent() ) );

#if 0
        std::cout << "[ComponentPointEvaluation] c = " << __c << "\n"
                  << "[ComponentPointEvaluation] d = " << b.polynomialDimensionPerComponent() << "\n"
                  << "[ComponentPointEvaluation] eval = " << b.evaluate( __pt ) << "\n"
                  << "[ComponentPointEvaluation] eval = " << b.basis()( __pt ) << "\n"
                  << "[ComponentPointEvaluation] bvals=" << ublas::project( b.evaluate( __pt ),
                          ublas::slice( __c, space_type::nComponents, b.polynomialDimensionPerComponent() ),
                          ublas::slice( 0, 0, 1 ) ) << "\n"
                  << "[ComponentPointEvaluation] m = " << ublas::project( m,
                          ublas::slice( __c, 0, 1 ),
                          ublas::slice( 0, 1, b.polynomialDimensionPerComponent() ) )
                  << "\n";
#endif
        ublas::row( m, __c ) = ublas::column( b.basis()( __pt ), 0 );

        this->setCoefficient( m );
        //std::cout << "[ComponentPointEvaluation] m = " << m << "\n";

    }
};

/**
 * \class PointDerivative
 * \brief generate the functional associated with a point derivative
 *
 * Given a function space \f$ E \f$, generate the functional \f$ \ell
 * : E \longrightarrow R \f$ associated with the derivation
 * of the basis functions of the function space at a point in the
 * geometric space.
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
class PointDerivative
    :
public Functional<Space>
{
    typedef Functional<Space> super;
public:

    typedef PointDerivative<Space> self_type;
    typedef typename super::space_type space_type;
    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;

    PointDerivative()
        :
        super()
    {}
    PointDerivative( space_type const& b, int i, node_type const& __pt )
        : super( b, ublas::trans( b.derivate( i ).evaluate(  __pt ) ) )
    {

    }
};

/**
 * \class PointsEvaluation
 * \brief generate the functionals associated with  point set
 *
 * Given a function space \f$ E \f$ and a set of points in the
 * geometric space \f$ \{x_i\}_{i=1...N} \f$, generate the set of
 * functionals
 * \f{eqnarray*}
 * \ell_i :& E \longrightarrow R, i = 1...N\\
 *         & f \longrightarrow f(x_i)
 * \f}
 * associated with the evaluation of the basis functions of the
 * function space at the set of points in the geometric space.
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
struct PointsEvaluation
        :
    public std::vector<Functional<Space> >
{
    typedef std::vector<Functional<Space> > super;
public:

    typedef PointsEvaluation<Space> self_type;
    typedef Space space_type;
    typedef typename space_type::points_type points_type;


    PointsEvaluation()
        :
        super()
    {}
    PointsEvaluation( space_type const& b, points_type const& __pts )
        : super()
    {
        for ( uint16_type c = 0; c < __pts.size2(); ++c )
        {
            //std::cout << "[PointsEvaluation] eval at point " << ublas::column( __pts, c)  << "\n";
            this->push_back( PointEvaluation<Space>( b, ublas::column( __pts, c ) ) );
        }

    }
}; // PointsEvaluation


/**
 * \class PointsDerivative
 * \brief generate the functionals associated with  point set
 *
 * Given a function space \f$ E \f$ and a set of points in the
 * geometric space \f$ \{y_i\}_{i=1...N} \f$, generate the set of
 * functionals
 * \f{eqnarray*}
 * \ell_i,j :& E \longrightarrow R, i = 1...N\\
 *         & f \longrightarrow \partial f(y_i)/ \partial x_j
 * \f}
 * associated with the evaluation of the basis functions of the
 * function space at the set of points in the geometric space.
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
struct PointsDerivative
        :
    public std::vector<Functional<Space> >
{
    typedef std::vector<Functional<Space> > super;
public:

    typedef PointsDerivative<Space> self_type;
    typedef Space space_type;
    typedef typename space_type::points_type points_type;


    PointsDerivative()
        :
        super()
    {}
    PointsDerivative( space_type const& b, int i, points_type const& __pts )
        : super()
    {
        for ( uint16_type c = 0; c < __pts.size2(); ++c )
        {
            //std::cout << "[PointsDerivative] eval at point " << ublas::column( __pts, c)  << "\n";
            this->push_back( PointDerivative<Space>( b, i, ublas::column( __pts, c ) ) );
        }

    }
}; // PointsDerivative

/**
 * \class PointsGradient
 * \brief generate the functionals associated with  point set
 *
 * Given a function space \f$ E \f$ and a set of points in the
 * geometric space \f$ \{y_i\}_{i=1...N} \f$, generate the set of
 * functionals
 * \f{eqnarray*}
 * \ell_i,j :& E \longrightarrow R, i = 1...N\\
 *         & f \longrightarrow \partial f(y_i)/ \partial x_j
 * \f}
 * associated with the evaluation of the basis functions of the
 * function space at the set of points in the geometric space.
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
struct PointsGradient
        :
    public std::vector<Functional<Space> >
{
    typedef std::vector<Functional<Space> > super;
public:

    typedef PointsGradient<Space> self_type;
    typedef Space space_type;
    typedef typename space_type::points_type points_type;


    PointsGradient()
        :
        super()
    {}
    PointsGradient( space_type const& b, points_type const& __pts )
        : super()
    {
        for ( uint16_type c = 0; c < __pts.size2(); ++c )
        {
            for ( int j = 0; j < __pts.size1(); ++j )
                this->push_back( PointDerivative<Space>( b, j, ublas::column( __pts, c ) ) );
        }

    }
}; // PointsGradient

/**
 *
 * \brief generate the functionals associated with  point set
 *
 * Given a function space \f$ E \f$ and a set of points in the
 * geometric space \f$ \{x_i\}_{i=1...N} \f$, generate the set of
 * functionals
 * \f{eqnarray*}
 * \ell_i :& E \longrightarrow R, i = 1...N\\
 *         & f \longrightarrow f(x_i)
 * \f}
 * associated with the evaluation of the basis functions of the
 * function space at the set of points in the geometric space.
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
struct ComponentsPointsEvaluation
        :
    public std::vector<Functional<Space> >
{
    typedef std::vector<Functional<Space> > super;
public:

    typedef ComponentsPointsEvaluation<Space> self_type;
    typedef Space space_type;
    typedef typename space_type::points_type points_type;


    BOOST_STATIC_ASSERT( space_type::is_vectorial ||
                         space_type::is_tensor2 );

    ComponentsPointsEvaluation()
        :
        super()
    {}
    ComponentsPointsEvaluation( space_type const& b, points_type const& __pts )
        : super()
    {

        for ( uint16_type c = 0; c < space_type::nComponents; ++c )
        {
            for ( uint16_type pt = 0; pt < __pts.size2(); ++pt )
            {
                //std::cout << "[ComponentsPointsEvaluation] eval at point " << ublas::column( __pts, c)  << "\n";
                this->push_back( ComponentPointEvaluation<Space>( b, c, ublas::column( __pts, pt ) ) );
            }
        }

    }
}; // ComponentsPointsEvaluation


/**
 * \class DirectionalComponentPointEvaluation
 * \brief functional associate with directional component point evaluation
 *
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
struct DirectionalComponentPointEvaluation
        :
    public Functional<Space>
{
    typedef Functional<Space> super;
public:

    typedef DirectionalComponentPointEvaluation<Space> self_type;
    typedef Space space_type;
    typedef typename space_type::points_type points_type;

    typedef typename super::value_type value_type;
    typedef typename node<value_type>::type node_type;


    BOOST_STATIC_ASSERT( space_type::is_vectorial ||
                         space_type::is_tensor2 );

    DirectionalComponentPointEvaluation( space_type const& b,
                                         node_type const& d,
                                         node_type const& __pt )
        : super( b )
    {
        ublas::matrix<value_type> m( d.size(), b.basis().size() );
        ublas::matrix<value_type> pts( __pt.size(), 1 );
        ublas::column( pts, 0 ) = __pt;
        for ( int i = 0; i < d.size(); ++i )
        {
            ublas::row( m, i ) = d( i ) * ublas::column( b.basis().evaluate( pts ), 0 );
        }
        this->setCoefficient( m );
    }
}; // DirectionalComponentPointEvaluation

/**
 * \class DirectionalComponentPointsEvaluation
 * \brief functional associate with directional component point evaluation
 *
 *
 * \author Christophe Prud'homme
 */
template<typename Space>
struct DirectionalComponentPointsEvaluation
        :
    public std::vector<Functional<Space> >
{
    typedef std::vector<Functional<Space> > super;
public:

    typedef DirectionalComponentPointsEvaluation<Space> self_type;
    typedef Functional<Space> functional_type;
    typedef Space space_type;
    typedef typename space_type::points_type points_type;

    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;


    BOOST_STATIC_ASSERT( space_type::is_vectorial ||
                         space_type::is_tensor2 );

    DirectionalComponentPointsEvaluation()
        :
        super()
    {}
    DirectionalComponentPointsEvaluation( space_type const& b,
                                          node_type const& d,
                                          points_type const& __pts )
        : super()
    {
        for ( int j = 0; j < __pts.size2(); ++j )
        {
            ublas::matrix<value_type> m( ublas::zero_matrix<value_type>( d.size(), b.basis().size() ) );

            for ( int i = 0; i < d.size(); ++i )
            {
                ublas::row( m, i ) = d( i ) * ublas::column( b.basis()( ublas::column( __pts, j ) ), 0 );
            }

            this->push_back( functional_type( b, m ) );
        }
    }
}; // DirectionalComponentPointsEvaluation

/**
 * \class DirectionalComponentPointsEvaluation
 * \brief functional associate with directional component point evaluation
 *
 *
 * \author Christophe Prud'homme
 */
    template<typename Space>
struct SecondDirectionalComponentPointEvaluation
        :
    public std::vector<Functional<Space> >
{
    typedef std::vector<Functional<Space> > super;
public:

    typedef SecondDirectionalComponentPointEvaluation<Space> self_type;
    typedef Functional<Space> functional_type;
    typedef Space space_type;
    typedef typename space_type::points_type points_type;

    typedef typename space_type::value_type value_type;
    typedef typename node<value_type>::type node_type;


    BOOST_STATIC_ASSERT( space_type::is_vectorial ||
                         space_type::is_tensor2 );

    SecondDirectionalComponentPointEvaluation()
        :
        super()
    {}
    SecondDirectionalComponentPointEvaluation( space_type const& b,
                                               node_type const& d,
                                               points_type const& __pts,
                                               ublas::vector<value_type> poly)
        : super()
    {
        ublas::matrix<value_type> m( ublas::zero_matrix<value_type>( d.size(), b.basis().size() ) );

        // Evaluate poly on points of __pts (2 points)
        double q1_1 = 0;
        double q1_2 = 0;
        for(int k=0; k<__pts.size1(); ++k)
            {
                q1_1 += poly(k)*__pts(k,0);
                q1_2 += poly(k)*__pts(k,1);
            }
        q1_1 += poly(poly.size() - 1);
        q1_2 += poly(poly.size() - 1);

        // approximate integral
        for ( int i = 0; i < d.size(); ++i )
            {
                ublas::row( m, i ) = 0.5*d( i )* q1_1 *ublas::column( b.basis()( ublas::column( __pts, 0 ) ), 0 )
                    - 0.5*d( i )* q1_2 *ublas::column( b.basis()( ublas::column( __pts, 1 ) ), 0 );
            }

        std::cout << "[SecondDirectionalComponentPointEvaluation] m  = " << m << "\n";
        //this->push_back( functional_type( b, m ) );
        this->push_back( functional_type( b, m ) );
    }
}; // SecondDirectionalComponentPointEvaluation

} // functional

} // Feel

#endif // __FEELPP_FUNCTIONALS_HPP
