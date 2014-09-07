/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-24

  Copyright (C) 2005,2006 EPFL

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
   \file operations.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-24
 */
#ifndef __operations_H
#define __operations_H 1


#include <feel/feelalg/svd.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <Eigen/SVD>

namespace Feel
{

/**
 * \brief integrate a polynomial over a convex
 *
 * since we use a L2 orthonormal basis, the integral of the polynomial
 * is equal to the first first coefficient of the polynomial
 */
template<
typename Poly,
         template<uint16_type> class PsetType
         >
typename Poly::value_type
integrate( Polynomial<Poly, PsetType> const& p )
{
    return p.coeff()( 0, 0 );
}

/**
 * \brief compute \f$\frac{\partial \cdot}{\partial x}\f$  of a  polynomial \p p
 *
 * \return the gradient of the scalar polynomial
 */
template<typename Poly, template<uint16_type> class Type>
Polynomial<Poly, Type>
dx( Polynomial<Poly, Type> const& p )
{
    return Polynomial<Poly, Type>( ublas::prod( p.coeff(), p.basis().d( 0 ) ), true );
}
/**
 * \brief compute \f$\frac{\partial \cdot}{\partial x}\f$  of a  polynomial set \p p
 *
 * \return the gradient of the scalar polynomial
 */
template<typename Poly, template<uint16_type> class Type>
PolynomialSet<Poly, Type>
dx( PolynomialSet<Poly, Type> const& p )
{
    return PolynomialSet<Poly, Type>( ublas::prod( p.coeff(), p.basis().d( 0 ) ), true );
}
/**
 * \brief compute \f$\frac{\partial \cdot}{\partial y}\f$  of a polynomial \p p
 *
 * \return the gradient of the scalar polynomial
 */
template<typename Poly, template<uint16_type> class Type>
Polynomial<Poly, Type>
dy( Polynomial<Poly, Type> const& p )
{
    BOOST_STATIC_ASSERT( ( Polynomial<Poly,Type>::nDim >= 2 ) );
    return Polynomial<Poly, Type>( Poly(), ublas::prod( p.coeff(), p.basis().d( 1 ) ), true );
}
/**
 * \brief compute \f$\frac{\partial \cdot}{\partial y}\f$  of a polynomial set \p p
 *
 * \return the gradient of the scalar polynomial
 */
template<typename Poly, template<uint16_type> class Type>
PolynomialSet<Poly, Type>
dy( PolynomialSet<Poly, Type> const& p )
{
    BOOST_STATIC_ASSERT( ( Polynomial<Poly,Type>::nDim >= 2 ) );
    return PolynomialSet<Poly, Type>( Poly(), ublas::prod( p.coeff(), p.basis().d( 1 ) ), true );
}
/**
 * \brief compute \f$\frac{\partial \cdot}{\partial z}\f$  of a  polynomial \p p
 *
 * \return the gradient of the scalar polynomial
 */
template<typename Poly, template<uint16_type> class Type>
Polynomial<Poly, Type>
dz( Polynomial<Poly, Type> const& p )
{
    BOOST_STATIC_ASSERT( ( Polynomial<Poly,Type>::nDim >= 3 ) );
    return Polynomial<Poly, Type>( Poly(), ublas::prod( p.coeff(), p.basis().d( 2 ) ), true );
}
/**
 * \brief compute \f$\frac{\partial \cdot}{\partial z}\f$  of a  polynomial set \p p
 *
 * \return the gradient of the scalar polynomial
 */
template<typename Poly, template<uint16_type> class Type>
PolynomialSet<Poly, Type>
dz( PolynomialSet<Poly, Type> const& p )
{
    BOOST_STATIC_ASSERT( ( PolynomialSet<Poly,Type>::nDim >= 3 ) );
    return PolynomialSet<Poly, Type>( Poly(), ublas::prod( p.coeff(), p.basis().d( 2 ) ), true );
}
/**
 * \brief compute the gradient of a scalar polynomial \p p
 *
 * \return the gradient of the scalar polynomial
 */
template<typename Poly>
Polynomial<Poly, Vectorial>
gradient( Polynomial<Poly, Scalar> const& p )
{
    const int nDim = Polynomial<Poly, Scalar>::nDim;

    typedef typename Poly::value_type value_type;
    const int ndof = p.coeff().size2();
    ublas::matrix<value_type> c ( nDim*nDim, ndof );
    c.clear();

    for ( int i = 0; i <nDim; ++i )
        ublas::row( c, nDim*i+i ) = ublas::row( ublas::prod( p.coeff(), p.basis().d( i ) ), 0 );

    return Polynomial<Poly, Vectorial>( Poly(), c, true );
}
/**
 * \brief compute the gradient of a vectorial polynomial \p p
 *
 * \return the gradient of the vectorial polynomial
 */
template<typename Poly>
Polynomial<Poly, Tensor2>
gradient( Polynomial<Poly, Vectorial> const& p )
{
    const int nDim = Polynomial<Poly, Vectorial>::nDim;
    const int nComponents = Polynomial<Poly, Tensor2>::nComponents;

    typedef typename Poly::value_type value_type;
    const int ndof = p.coeff().size2();
    ublas::matrix<value_type> c ( nComponents, ndof );
    c.clear();

    for ( int i = 0; i < nDim; ++i )
        for ( int j = 0; j < nDim; ++j )
            ublas::row( c, nComponents*i+nDim*j ) = ublas::row( ublas::prod( p[j].coeff(), p.basis().d( i ) ), 0 );

    return Polynomial<Poly, Tensor2>( Poly(), c, true );
}
/**
 * \brief compute the gradient of a scalar polynomial set
 *
 * \return the gradient of the scalar polynomial set
 */
template<typename Poly>
PolynomialSet<Poly, Vectorial>
gradient( PolynomialSet<Poly, Scalar> const& p )
{
    const int nDim = PolynomialSet<Poly, Scalar>::nDim;

    typedef typename Poly::value_type value_type;
    const int n1 = p.coeff().size1();
    const int n2 = p.coeff().size2();
    ublas::matrix<value_type> c ( nDim*nDim*n1, n2 );
    c.clear();

    for ( int i = 0; i <nDim; ++i )
    {
        ublas::project( c,
                        ublas::slice( nDim*n1*i+i, nDim, n1 ),
                        ublas::slice( 0, 1, n2 ) )  = ublas::prod( p.coeff(), p.basis().d( i ) );
    }

    return PolynomialSet<Poly, Vectorial>( Poly(), c, true );
}
/**
 * \brief compute the gradient of a vectorial polynomial set
 *
 * \return the gradient of the vectorial polynomial set
 */
template<typename Poly>
PolynomialSet<Poly, Tensor2>
gradient( PolynomialSet<Poly, Vectorial> const& p )
{
    const int nDim = PolynomialSet<Poly, Scalar>::nDim;
    const int nComponents = PolynomialSet<Poly, Tensor2>::nComponents;

    typedef typename Poly::value_type value_type;
    const int n1 = p.coeff().size1();
    const int n2 = p.coeff().size2();
    ublas::matrix<value_type> c ( nComponents*n1, n2 );
    c.clear();

    for ( int i = 0; i <nDim; ++i )
        for ( int j = 0; j < nDim; ++j )
        {
            ublas::project( c,
                            ublas::slice( nDim*n1*i+i, nDim, n1 ),
                            ublas::slice( 0, 1, n2 ) )  = ublas::prod( p.coeff(), p.basis().d( i ) );
        }

    return PolynomialSet<Poly, Tensor2>( Poly(), c, true );
}
/**
 * \brief compute the gradient of a vectorial polynomial set
 *
 * \return the gradient of the vectorial polynomial set
 */
template<typename Poly>
PolynomialSet<Poly, Tensor2>
hessian( PolynomialSet<Poly, Scalar> const& p )
{
    const int nDim = PolynomialSet<Poly, Scalar>::nDim;
    const int nComponents = PolynomialSet<Poly, Tensor2>::nComponents;

    typedef typename Poly::value_type value_type;
    const int n1 = p.coeff().size1();
    const int n2 = p.coeff().size2();
    ublas::matrix<value_type> c ( nComponents*nComponents*n1, n2 );
    c.clear();

    for ( int i = 0; i <nDim; ++i )
        for ( int j = 0; j < nDim; ++j )
        {
            ublas::project( c,
                            ublas::slice( nComponents*n1*( nDim*j+i ), nComponents, n1 ),
                            ublas::slice( 0, 1, n2 ) )  = ublas::prod( p.coeff(), p.d( i, j ) );
        }

    return PolynomialSet<Poly, Tensor2>( Poly(), c, true );
}

/**
 * \brief compute the divergence of a vectorial polynomial \p p
 *
 * \return the divergence of the vectorial polynomial
 */
template<typename Poly>
Polynomial<Poly, Scalar>
divergence( Polynomial<Poly, Vectorial> const& p )
{
    const int nDim = Polynomial<Poly, Vectorial>::nDim;

    typedef typename Poly::value_type value_type;

    ublas::matrix<value_type> c ( ublas::zero_matrix<value_type>( 1, p.coeff().size2() ) );

    for ( int i = 0; i <nDim; ++i )
    {
        ublas::noalias( c ) += ublas::prod( p[i].coeff(), p.basis().d( i ) );
    }

    return Polynomial<Poly, Scalar>( Poly(), c, true );
}
/**
 * \brief compute the divergence of a vectorial polynomial set \p p
 *
 * \return the divergence of the vectorial polynomial set
 */
template<typename Poly>
PolynomialSet<Poly, Scalar>
divergence( PolynomialSet<Poly, Vectorial> const& p )
{
    const int nComponents = PolynomialSet<Poly, Vectorial>::nComponents;

    typedef typename Poly::value_type value_type;

    ublas::matrix<value_type> c ( ublas::zero_matrix<value_type>( p.coeff().size1()/nComponents,
                                  p.coeff().size2() ) );

    for ( int i = 0; i <nComponents; ++i )
    {
        ublas::noalias( c ) += ublas::prod( p[i].coeff(), p.basis().d( i ) );
    }

    return PolynomialSet<Poly, Scalar>( Poly(), c, true );
}
/**
 * \brief compute the curl of a vectorial polynomial \p p
 *
 * \return the curl of the vectorial polynomial
 */
template<typename Poly>
Polynomial<Poly, Vectorial>
curl( Polynomial<Poly, Vectorial> const& p )
{
    const int nDim = Polynomial<Poly, Vectorial>::nDim;
    const int nComponents = Polynomial<Poly, Vectorial>::nComponents;
    BOOST_STATIC_ASSERT( ( Polynomial<Poly, Vectorial>::nDim == 3 ) );

    typedef typename Poly::value_type value_type;

    ublas::matrix<value_type> c ( p.coeff().size1(), p.coeff().size2() );
    c.clear();
    const uint16_type dim_p = p.coeff().size1()/nComponents;
    const uint16_type ndof = p.coeff().size2();

    ublas::project( c,
                    ublas::slice( 0, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = p[2].derivate( 1 ).coeff() - p[1].derivate( 2 ).coeff();
    ublas::project( c,
                    ublas::slice( dim_p+1, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = p[0].derivate( 2 ).coeff() - p[2].derivate( 0 ).coeff();
    ublas::project( c,
                    ublas::slice( 2*dim_p+2, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = p[1].derivate( 0 ).coeff() - p[0].derivate( 1 ).coeff();

    return Polynomial<Poly, Vectorial>( Poly(), c, true );
}
/**
 * \brief compute the curl of a vectorial polynomial set \p p
 *
 * \return the curl of the vectorial polynomial
 */
template<typename Poly>
PolynomialSet<Poly, Vectorial>
curl( PolynomialSet<Poly, Vectorial> const& p )
{
    const int nDim = PolynomialSet<Poly, Vectorial>::nDim;
    const int nComponents = PolynomialSet<Poly, Vectorial>::nComponents;
    BOOST_STATIC_ASSERT( ( PolynomialSet<Poly, Vectorial>::nDim == 3 ) );

    typedef typename Poly::value_type value_type;

    ublas::matrix<value_type> c ( p.coeff().size1(), p.coeff().size2() );
    c.clear();
    const uint16_type dim_p = p.coeff().size1()/nComponents;
    const uint16_type ndof = p.coeff().size2();

    ublas::project( c,
                    ublas::slice( 0, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = p[2].derivate( 1 ).coeff() - p[1].derivate( 2 ).coeff();
    ublas::project( c,
                    ublas::slice( dim_p+1, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = p[0].derivate( 2 ).coeff() - p[2].derivate( 0 ).coeff();
    ublas::project( c,
                    ublas::slice( 2*dim_p+2, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = p[1].derivate( 0 ).coeff() - p[0].derivate( 1 ).coeff();

    return PolynomialSet<Poly, Vectorial>( Poly(), c, true );
}
/**
 * \brief compute the transpose of the curl of a vectorial polynomial \p p
 * \note swap the signs of the curl matrix operator
 * \return the transpose of the curl of the vectorial polynomial
 */
template<typename Poly>
Polynomial<Poly, Vectorial>
curlt( Polynomial<Poly, Vectorial> const& p )
{
    const int nDim = Polynomial<Poly, Vectorial>::nDim;
    const int nComponents = Polynomial<Poly, Vectorial>::nComponents;
    BOOST_STATIC_ASSERT( ( Polynomial<Poly, Vectorial>::nDim == 3 ) );

    typedef typename Poly::value_type value_type;

    ublas::matrix<value_type> c ( p.coeff().size1(), p.coeff().size2() );
    const uint16_type dim_p = p.coeff().size1()/nComponents;
    const uint16_type ndof = p.coeff().size2();

    ublas::project( c,
                    ublas::slice( 0, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = - p[2].derivate( 1 ).coeff() + p[1].derivate( 2 ).coeff();
    ublas::project( c,
                    ublas::slice( dim_p+1, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = - p[0].derivate( 2 ).coeff() + p[2].derivate( 0 ).coeff();
    ublas::project( c,
                    ublas::slice( 2*dim_p+2, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = - p[1].derivate( 0 ).coeff() + p[0].derivate( 1 ).coeff();

    return Polynomial<Poly, Vectorial>( Poly(), c, true );
}
/**
 * \brief compute the transpose of the curl of a vectorial polynomial \p p
 * \note swap the signs of the curl matrix operator
 * \return the transpose of the curl of the vectorial polynomial
 */
template<typename Poly>
PolynomialSet<Poly, Vectorial>
curlt( PolynomialSet<Poly, Vectorial> const& p )
{
    const int nDim = PolynomialSet<Poly, Vectorial>::nDim;
    const int nComponents = PolynomialSet<Poly, Vectorial>::nComponents;
    BOOST_STATIC_ASSERT( ( PolynomialSet<Poly, Vectorial>::nDim == 3 ) );

    typedef typename Poly::value_type value_type;

    ublas::matrix<value_type> c ( p.coeff().size1(), p.coeff().size2() );
    const uint16_type dim_p = p.coeff().size1()/nComponents;
    const uint16_type ndof = p.coeff().size2();

    ublas::project( c,
                    ublas::slice( 0, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = - p[2].derivate( 1 ).coeff() + p[1].derivate( 2 ).coeff();
    ublas::project( c,
                    ublas::slice( dim_p+1, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = - p[0].derivate( 2 ).coeff() + p[2].derivate( 0 ).coeff();
    ublas::project( c,
                    ublas::slice( 2*dim_p+2, nComponents, dim_p/nComponents ),
                    ublas::slice( 0, 1, ndof ) ) = - p[1].derivate( 0 ).coeff() + p[0].derivate( 1 ).coeff();

    return PolynomialSet<Poly, Vectorial>( Poly(), c, true );
}
/**
 * \f$ L_2\f$ Projection over the associated convex of a function \p f
 * onto a scalar polynomial set.
 *
 * \return the polynomial with the same
 * basis as the polynomial set and with coeff \f$ \sum_q w_q f_q
 * \phi_q \f$
 */
template<typename Pset,
         typename Func,
         typename IM>
Polynomial<Pset, Scalar>
project( Pset const& pset, Func const& f, IM const& im )
{
    typedef typename Func::value_type value_type;

    // evaluate f at quad nodes
    ublas::vector<value_type> fq( f( im.points() ) );

    // evaluate pset at quad nodes
    ublas::matrix<value_type> psetq( pset.evaluate( im.points() ) );

    ublas::matrix<value_type> c( 1, psetq.size1() );
    ublas::row( c, 0 ) = ublas::prod( psetq, ublas::element_prod( fq, im.weights() ) );
    return Polynomial<Pset, Scalar>( pset, c, true );
}

/**
 * \brief union of two polynomial sets \p P1 and \p P2 : \f$ P_1 \oplus P_2\f$
 *
 * \return the union of the two sets by appending the coefficient and
 * computing the range of the resulting set
 */
template<typename P, template<uint16_type> class Type,
         template<class, template<uint16_type> class> class Poly1,
         template<class, template<uint16_type> class> class Poly2 >
PolynomialSet<P, Type>
unite( Poly1<P, Type> const& pset1,
       Poly2<P, Type> const& pset2 )
{
    typedef PolynomialSet<P, Type> res_type;
    typedef typename res_type::value_type value_type;

    FEELPP_ASSERT( pset1.coeff().size2() == pset2.coeff().size2() )
    ( pset1.coeff().size2() )( pset2.coeff().size2() ).error( "incompatible size" );

    ublas::matrix<value_type> M( pset1.coeff().size1()+pset2.coeff().size1(), pset1.coeff().size2() );
    project( M,
             ublas::range( 0, pset1.coeff().size1() ),
             ublas::range( 0,pset1.coeff().size2() ) ) = pset1.coeff();
    project( M,
             ublas::range( pset1.coeff().size1(), pset1.coeff().size1()+pset2.coeff().size1() ),
             ublas::range( 0,pset1.coeff().size2() ) ) = pset2.coeff();

    M = res_type::polyset_type::toMatrix( M );
#if 1
    //std::cout << "before SVD M = " << M << "\n";
    Eigen::MatrixXd A ( Eigen::MatrixXd::Zero( M.size1(), M.size2() ) );

    for ( int i = 0; i < M.size1(); ++i )
        for ( int j = 0; j < M.size2(); ++j )
        {
            A( i,j ) = M( i,j );
        }

    Eigen::JacobiSVD<Eigen::MatrixXd> svdOfA( A, Eigen::ComputeThinU | Eigen::ComputeThinV );
    //std::cout << "SVDofA.S() = " << svdOfA.singularValues() << "\n";
    //std::cout << "SVDofA.U() = " << svdOfA.matrixU() << "\n";
    //std::cout << "SVDofA.V() = " << svdOfA.matrixV() << "\n";
#endif
    //SVD<ublas::matrix<value_type> > svd( M );
    //std::cout << "S() = " << svd.S() << "\n";
#if 1

    ublas::matrix<value_type> m ( svdOfA.singularValues().size(), M.size2() );

    for ( int i = 0; i < m.size1(); ++i )
        for ( int j = 0; j < m.size2(); ++j )
        {
            m( i,j ) = svdOfA.matrixV()( j, i );
        }
    //VLOG(1) << "V=" << m << "\n";
#else
    ublas::matrix<value_type> m ( ublas::subrange( svd.V(), 0, svd.S().size(), 0, M.size2() ) );

    //std::cout << "m=" << m << "\n";
#endif
    return res_type( P(), res_type::polyset_type::toType( m ), true  );
}

}

#endif /* __operations_H */
