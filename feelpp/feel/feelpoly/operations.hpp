// SPDX-FileCopyrightText: 2005-2006 EPFL
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

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

/** @brief A polynomial object with an explicit expansion-order policy. */
template<class P>
concept PolynomialObject = requires( P const& p )
{
    P::order_spec;
    p.coeff();
    p.basis();
    p.order();
    p.derivative( uint16_type{} );
};

/** @brief Integrate the first component using the normalized constant mode. */
template<class Poly, template<uint16_type> class Shape, class Container, int OrderSpec>
auto integrate( Polynomial<Poly, Shape, Container, OrderSpec> const& p )
{
    return p.toSet().firstMoments()(0,0);
}

/** @brief Differentiate each polynomial with respect to x. */
template<PolynomialObject P>
auto dx( P const& p ) { return p.derivative(0); }

/** @brief Differentiate each polynomial with respect to y. */
template<PolynomialObject P>
auto dy( P const& p )
{
    static_assert(P::nDim >= 2);
    return p.derivative(1);
}

/** @brief Differentiate each polynomial with respect to z. */
template<PolynomialObject P>
auto dz( P const& p )
{
    static_assert(P::nDim >= 3);
    return p.derivative(2);
}

/** @brief Gradient or Jacobian, preserving function count and expansion order. */
template<PolynomialObject P>
auto gradient( P const& p ) { return p.gradient(); }

/** @brief Hessian of each scalar polynomial, in component/direction order. */
template<PolynomialObject P>
    requires(P::is_scalar)
auto hessian( P const& p ) { return p.gradient().gradient(); }

/** @brief Divergence of every vector polynomial, retaining the ambient expansion. */
template<class Poly, int OrderSpec>
auto divergence( PolynomialSet<Poly, Vectorial, OrderSpec> const& p )
{
    using result_type=PolynomialSet<Poly,Scalar,OrderSpec>;
    typename result_type::matrix_type coefficients(
        ublas::zero_matrix<typename result_type::value_type>(p.polynomialDimension(),p.coeff().size2()));
    for(uint16_type direction=0;direction<p.nDim;++direction)
        coefficients+=p[direction].derivative(direction).coeff();
    return result_type(p,coefficients,true);
}

/** @brief Divergence of one vector polynomial using the shared set operation. */
template<class Poly, class Container, int OrderSpec>
auto divergence( Polynomial<Poly, Vectorial, Container, OrderSpec> const& p )
{
    return divergence(p.toSet()).polynomial(0);
}

/** @brief Curl of every vector polynomial on a three-dimensional reference cell. */
template<class Poly, int OrderSpec>
auto curl( PolynomialSet<Poly, Vectorial, OrderSpec> const& p )
{
    static_assert(Poly::nDim==3);
    using result_type=PolynomialSet<Poly,Vectorial,OrderSpec>;
    typename result_type::matrix_type coefficients(p.coeff().size1(),p.coeff().size2());
    for(uint16_type i=0;i<3;++i)
    {
        auto const j=static_cast<uint16_type>((i+1)%3);
        auto const k=static_cast<uint16_type>((i+2)%3);
        ublas::project(coefficients,ublas::slice(i,3,p.polynomialDimension()),
                       ublas::slice(0,1,p.coeff().size2())) =
            p[k].derivative(j).coeff()-p[j].derivative(k).coeff();
    }
    return result_type(p,coefficients,true);
}

/** @brief Curl of one vector polynomial using the shared set operation. */
template<class Poly, class Container, int OrderSpec>
auto curl( Polynomial<Poly, Vectorial, Container, OrderSpec> const& p )
{
    return curl(p.toSet()).polynomial(0);
}

/** @brief Apply the transposed curl convention, the negative of curl. */
template<PolynomialObject P>
auto curlt( P const& p )
{
    auto result=curl(p);
    typename decltype(result)::matrix_type coefficients=-result.coeff();
    result.setCoefficient(coefficients,true);
    return result;
}

/**
 * \f$ L_2\f$ Projection over the associated convex of a function \p f
 * onto a scalar polynomial set.
 * @pre The supplied set is orthonormal for the integration rule.
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
    return Polynomial<Pset, Scalar>( pset, c );
}

/**
 * \brief union of two polynomial sets \p P1 and \p P2 : \f$ P_1 \oplus P_2\f$
 *
 * \return the union of the two sets by appending the coefficient and
 * computing the range of the resulting set
 */
template<typename P, template<uint16_type> class Type,
         template<class, template<uint16_type> class, int> class Poly1,
         template<class, template<uint16_type> class, int> class Poly2,
         int OrderSpec1, int OrderSpec2>
PolynomialSet<P, Type, OrderSpec1>
unite( Poly1<P, Type, OrderSpec1> const& pset1,
       Poly2<P, Type, OrderSpec2> const& pset2 )
{
    typedef PolynomialSet<P, Type, OrderSpec1> res_type;
    typedef typename res_type::value_type value_type;

    if(pset1.order()!=pset2.order() || pset1.coeff().size2()!=pset2.coeff().size2())
        throw std::invalid_argument("polynomial union requires matching expansion orders and extents");

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

    ublas::matrix<value_type> m ( svdOfA.rank(), M.size2() );

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
    res_type result( pset1, res_type::polyset_type::toType( m ), true );
    if constexpr ( res_type::is_order_dynamic )
        result.setOrder( pset1.order() );
    return result;
}

}

#endif /* __operations_H */
