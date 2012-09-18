/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-02-28

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
   \file testlagrange.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-02-28
 */
#ifndef __TestLagrange_H
#define __TestLagrange_H 1

#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/im.hpp>

using namespace Feel;


/*!
  \class TestLagrange
  \brief Test Lagrange Polynomials

  Test some properties of the Lagrange polynomials

  @author Christophe Prud'homme
  @see
*/
template<typename FE>
class TestLagrange
{
public:


    /** @name Typedefs
     */
    //@{

    typedef typename FE::value_type value_type;
    typedef typename FE::matrix_type matrix_type;

    //@}


    /** @name Constructors, destructor
     */
    //@{

    TestLagrange( value_type eps = value_type( 3000 )*Feel::type_traits<value_type>::epsilon() )
        :
        M_eps( eps )
    {}
    TestLagrange( TestLagrange const & tl )
        :
        M_eps( tl.M_eps )
    {}
    ~TestLagrange() {}

    //@}

    /** @name Operator overloads
     */
    //@{
    void operator()() const
    {
        std::ostringstream os;
        os << "lagrange." << FE::nDim << "." << FE::nOrder;

        //BOOST_CHECK( fe.familyName() == os.str() );
        if ( fe.name() != os.str() )
        {
            std::cout << "[FAILURE] invalid name : " << fe.name() << " instead of " << os.str() << "\n";
        }

        checkIdentity();
        //BOOST_MESSAGE( "checkDiff()\n" );
        checkDiff();
    }


    //@}

protected:

    void
    checkDiff() const
    {
        BOOST_TEST_MESSAGE( "Checking Lagrange polynomials differentiation" );

        ublas::vector<value_type> error( FE::nDim );
        ublas::vector<matrix_type> der( fe.derivate( fe.points() ) );

        for ( int i = 0; i < FE::nDim; ++i )
            error[i] = ublas::norm_frobenius( der[i] - fe.derivate( i ).evaluate( fe.points() ) );

#if defined(  USE_TEST )
        BOOST_TEST_MESSAGE( "[checkDiff] Lagrange " << " dim : " << FE::nDim << " order : " << FE::nOrder << " error: " << error << " eps: " << M_eps  );
        BOOST_CHECK( ublas::norm_inf( error ) < M_eps );
#endif
    }

    void
    checkIdentity() const
    {
        BOOST_TEST_MESSAGE( "Checking Lagrange polynomials identity" );

        matrix_type eval_at_pts( FE::polyset_type::toMatrix( fe.evaluate( fe.points() ) ) );
        BOOST_TEST_MESSAGE( "eval_at_pts = " << eval_at_pts  );
        value_type error = ublas::norm_frobenius( eval_at_pts-ublas::identity_matrix<value_type>( eval_at_pts.size1() ) );


#if defined(  USE_TEST )
        BOOST_TEST_MESSAGE( "[checkIdentity] Lagrange "
                            << " dim : " << FE::nDim
                            << " order : " << FE::nOrder
                            << " epsilon = " << M_eps
                            << " error = " << error );
        BOOST_CHECK( error < M_eps );
#endif
    }
private:
    FE fe;
    value_type M_eps;
};
#endif /* __TestLagrange_H */
