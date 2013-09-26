/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-06

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2011 Universite de Grenoble

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
   \file polynomial.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-06
 */
#ifndef __Polynomial_H
#define __Polynomial_H 1

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>


#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/policy.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

template<typename, template<uint16_type> class PolySetType > class PolynomialSet;

/**
 * \class Polynomial
 * \brief polynomial class
 *
 * The polynomial is expressed in the basis from \p Poly. The
 * coefficients of the polynomial in this basis are represented by a
 * matrix whose lines are the polymomial components coefficients (1 if
 * \code is_scalar == true \endcode, \p nDim if \code is_vectorial ==
 * true \endcode and columns are the basis
 *
 * Evaluating the polynomial at a set of points(or just one point) is
 * then simply a matrix-matrix product.
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see
 */
template<typename Poly,
         template<uint16_type> class PolySetType = Scalar,
         typename Container =  typename Poly::basis_type::matrix_type>
class Polynomial
{
public:

    /** @name Constants
     */
    //@{

    static const uint16_type nDim = Poly::nDim;
    static const uint16_type nOrder = Poly::nOrder;

    //@}


    /** @name Typedefs
     */
    //@{

    typedef Polynomial<Poly, PolySetType> self_type;
    typedef typename Poly::value_type value_type;
    typedef typename Poly::basis_type basis_type;



    typedef PolySetType<nDim> polyset_type;
    static const bool is_tensor2 = polyset_type::is_tensor2;
    static const bool is_vectorial = polyset_type::is_vectorial;
    static const bool is_scalar = polyset_type::is_scalar;
    static const uint16_type nComponents = polyset_type::nComponents;
    static const uint16_type nComponents1 = polyset_type::nComponents1;
    static const uint16_type nComponents2 = polyset_type::nComponents2;

    typedef typename Component<polyset_type>::type component_type;
    typedef Polynomial<Poly,Scalar> scalar_component_type;

    typedef typename basis_type::points_type points_type;
    typedef typename basis_type::matrix_type matrix_type;
    typedef Container container_type;

    typedef typename node<value_type>::type node_type;

    BOOST_STATIC_ASSERT( ( boost::is_same<typename matrix_type::value_type, value_type>::value ) );
    BOOST_STATIC_ASSERT( ( boost::is_same<typename matrix_type::value_type, typename points_type::value_type>::value ) );

    //@}

    /** @name Constructors, destructor
     */
    //@

    /**
     * default constructor
     */
    Polynomial()
        :
        M_basis(),
        M_coeff( M_basis.coeff() )
    {
    }


    /**
     * constructor giving only the underlying basis
     * \param __poly polynomial whose we take the basis
     */
    Polynomial( Poly const& __poly )
        :
        M_basis( __poly.basis() ),
        M_coeff( M_basis.coeff() )
    {
    }

    /**
     * constructor giving the underlying basis and the coefficient of
     * the polynomial in the basis
     *
     * \param __poly polynomial whose we take the basis
     * \param __coeff coefficients of the polynomial in the basis
     */
    Polynomial( Poly const& __poly, container_type const& __coeff, bool __as_is = false )
        :
        M_basis( __poly.basis() ),
        M_coeff( M_basis.coeff() )
    {
        setCoefficient( __coeff, __as_is );
    }

    /**
     * constructor giving the underlying basis and the coefficient of
     * the polynomial in the basis
     *
     * \param __coeff coefficients of the polynomial in the basis
     */
    Polynomial( container_type const& __coeff, bool __as_is = false )
        :
        M_basis(),
        M_coeff( M_basis.coeff() )
    {
        setCoefficient( __coeff, __as_is );
    }

    /**
     * constructor giving the underlying basis and the coefficient of
     * the polynomial in the basis
     *
     * \param __poly polynomial whose we take the basis
     * \param __coeff coefficients of the polynomial in the basis
     */
    template<class AE>
    Polynomial( Poly const& __poly, ublas::matrix_expression<AE> const& __coeff, bool __as_is = false )
        :
        M_basis( __poly.basis() ),
        M_coeff( M_basis.coeff() )
    {
        setCoefficient( __coeff, __as_is );
    }

    Polynomial( Polynomial const & p )
        :
        M_basis( p.M_basis ),
        M_coeff( p.M_coeff )
    {}

    ~Polynomial()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type const& operator=( self_type const& __p )
    {
        if ( this != &__p )
        {
            M_basis = __p.M_basis;
            M_coeff = __p.M_coeff;
        }

        return *this;
    }

    self_type const& operator-=( self_type const& __p )
    {
        M_coeff -= __p.M_coeff;
        return *this;
    }

    self_type const& operator()( self_type const& __p ) const
    {
        if ( this != &__p )
        {
            M_basis = __p.M_basis;
            M_coeff = __p.M_coeff;
        }

        return *this;
    }


    /**
     * \brief extract the i-th component of a vectorial polynomial
     *
     * \return the i-th component of the polynomial
     */
    component_type operator[]( int i ) const
    {
        const int ncols = M_coeff.size2();
        return component_type( Poly(), ublas::project( M_coeff,
                               ublas::slice( nComponents*i+i, 1, 1  ),
                               ublas::slice( 0, 1, ncols ) ) );
    }

    /**
     * Evaluate polynomial at point \p __x
     *
     * \param __x the coordinate of the point
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type
    operator()( node_type const& __x ) const
    {
        return ublas::prod( M_coeff, M_basis( __x ) );
    }

    /**
     * Evaluate polynomial at points \p __pts
     *
     * \param __pts the matrix of coordinates of the points
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type operator()( points_type const& __pts ) const
    {
        return ublas::prod( M_coeff, M_basis( __pts ) );
    }

    //@}

    /** @name Accessors
     */
    //@{


    /**
     * \return \c true if the polynomial set is zero, \c false otherwise
     */
    bool isZero() const
    {
        return ublas::norm_2( M_coeff ) < Feel::type_traits<value_type>::epsilon();

    }
    /**
     * \return the dof
     */
    matrix_type const& coeff() const
    {
        return M_coeff;
    }

    /**
     * \return the dof
     */
    matrix_type const& coefficients() const
    {
        return M_coeff;
    }

    /**
     * \return the basis in which the polynomial is expressed
     */
    basis_type const& basis() const
    {
        return M_basis;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the coefficient of the polynomial in the basis.
     */
    void setCoefficient( matrix_type const& __c, bool __as_is = false )
    {
        //FEELPP_ASSERT( __c.size1() == nComponents*nComponents && __c.size2() == M_coeff.size2() )
        //    ( is_scalar )( is_vectorial )( is_tensor2 )( __c )( M_coeff ).error( "invalid polynomial coefficients" );
        if ( !__as_is )
        {
            M_coeff = ublas::prod( polyset_type::toMatrix( __c ), polyset_type::toMatrix( M_coeff ) );
            M_coeff = polyset_type::toType( M_coeff );
        }

        else
            M_coeff = __c;
    }


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Evaluate polynomial at point \p __x
     *
     * \param __x the coordinate of the point
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type
    evaluate( node_type const& __x ) const
    {
        return ublas::prod( M_coeff, M_basis( __x ) );
    }

    /**
     * Evaluate polynomial at points \p __pts
     *
     * \param __pts the matrix of coordinates of the points
     * \return the evaluation of the polynomial at point \p __x
     */
    matrix_type evaluate( points_type const& __pts ) const
    {
        return ublas::prod( M_coeff, M_basis( __pts ) );
    }

    template<typename AE>
    matrix_type derivate( uint16_type i, ublas::matrix_expression<AE> const& pts ) const
    {
        ublas::vector<matrix_type> der( M_basis.derivate( pts ) );
        matrix_type res( M_coeff.size1(), pts().size2() );
        ublas::axpy_prod( M_coeff, der[i], res );
        return res;
    }

    template<typename AE>
    matrix_type derivate( uint16_type i, uint16_type j, ublas::matrix_expression<AE> const& pts ) const
    {
        //std::cout << "[derivate2] M_coeff = " << M_coeff << "\n";
        matrix_type eval( M_basis.evaluate( pts ) );
        //matrix_type res( M_coeff.size1(), pts().size2() );
        //ublas::axpy_prod( M_coeff, der[i], res );
        matrix_type p1 = ublas::prod( M_coeff, M_basis.d( i ) );
        matrix_type p2 = ublas::prod( p1, M_basis.d( j ) );
        return ublas::prod( p2, eval );
    }


    /**
     * \brief differentiation matrix of Dubiner polynomials
     * the derivatives are computed at the nodes of the lattice
     *
     * \arg i index of the derivative (0 : x, 1 : y, 2 : z )
     */
    matrix_type const& d( uint16_type i ) const
    {
        return M_basis.d( i );
    }

    /**
     * derivative of the polynomial in the l-th direction
     *
     * \param l direction of derivation
     * \return the derivative of the polynomial
     */
    self_type derivative( uint16_type l ) const
    {
        return self_type( Poly(), ublas::prod( M_coeff, M_basis.d( l ) ), true );
    }

    /**
     * transform a \p Polynomial into a PolynomialSet with one polynomial
     *
     * \return a \p PolynomialSet
     */
    PolynomialSet<Poly,PolySetType> toSet( bool asis = false ) const
    {
        return PolynomialSet<Poly,PolySetType>( Poly(), M_coeff, asis );
    }
#if 0
    Polynomial<Poly, PolySetType> operator-( Polynomial<Poly, PolySetType> const& p ) const
    {
        matrix_type c = M_coeff-p.M_coeff;
        std::cout << "c=" << c << "\n";
        return Polynomial<Poly, PolySetType>( Poly(), c );
    }
#endif
    //@}


protected:

private:

    basis_type M_basis;
    container_type M_coeff;
};
template<typename Poly, template<uint16_type> class PolySetType>
Polynomial<Poly, PolySetType> operator-( Polynomial<Poly, PolySetType> const& p1,Polynomial<Poly, PolySetType> const& p2 )
{
    auto c = p1.coeff()-p2.coeff();
    //std::cout << "operator- c=" << c << "\n";
    return Polynomial<Poly, PolySetType>( Poly(), c );
}

template<typename Poly,template<uint16_type> class PolySetType, typename Container> const bool Polynomial<Poly,PolySetType,Container>::is_scalar;
template<typename Poly,template<uint16_type> class PolySetType, typename Container> const bool Polynomial<Poly,PolySetType, Container>::is_vectorial;
template<typename Poly,template<uint16_type> class PolySetType, typename Container> const bool Polynomial<Poly,PolySetType, Container>::is_tensor2;

}
#endif /* __Polynomial_H */
