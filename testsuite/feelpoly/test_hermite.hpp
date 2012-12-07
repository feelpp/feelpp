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
   \file testhermite.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-03-04
 */
#ifndef __TestHermite_H
#define __TestHermite_H 1

#include <feel/feelpoly/hermite.hpp>
#include <feel/feelpoly/im.hpp>

using namespace Feel;


/*!
  \class TestHermite
  \brief Test Hermite Polynomials

  Test some properties of the Hermite polynomials

  @author Christophe Prud'homme
  @see
*/
template<typename FE>
class TestHermite
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

    TestHermite( value_type eps = value_type( 3000 )*Feel::type_traits<value_type>::epsilon() )
        :
        M_eps( eps )
    {}
    TestHermite( TestHermite const & tl )
        :
        M_eps( tl.M_eps )
    {}
    ~TestHermite() {}

    //@}

    /** @name Operator overloads
     */
    //@{
    void operator()() const
    {
        std::ostringstream os;
        os << "hermite." << FE::nDim << "." << FE::nOrder;

        //BOOST_CHECK( fe.familyName() == os.str() );
        if ( fe.name() != os.str() )
        {
            std::cout << "[FAILURE] invalid name : " << fe.name() << " instead of " << os.str() << "\n";
        }

        BOOST_MESSAGE( "checkIdentity()\n" );
        checkIdentity();
        BOOST_MESSAGE( "checkDiff()\n" );
        checkDiff();
    }


    //@}

protected:

    void
    checkDiff() const
    {
        const int nDim = fe.nDim;

        if ( nDim == 2 ) return;

        //BOOST_TEST_MESSAGE( "Checking Hermite polynomials differentiation" );

        ublas::vector<value_type> error( FE::nDim );
        ublas::vector<matrix_type> der( fe.derivate( fe.points() ) );
        glas::clean( der[0],1e-10 );
        BOOST_TEST_MESSAGE( "der[0]=" << der[0] << "\n" );

        matrix_type result( der[0].size1(), der[0].size2() );
        result = ublas::zero_matrix<value_type>( result.size1(), result.size2() );

        for ( int d = 0; d < nDim; ++d )
        {
            ublas::subrange( result, ( d+1 )*( nDim+1 ), ( d+2 )*( nDim+1 ), ( d )*( nDim+1 ), ( d+1 )*( nDim+1 ) ) = ublas::identity_matrix<value_type>( nDim+1 );
            ublas::subrange( result, ( d+1 )*( nDim+1 ), ( d+2 )*( nDim+1 ), ( d+1 )*( nDim+1 ), ( d+2 )*( nDim+1 ) ) = ublas::identity_matrix<value_type>( nDim+1 );
        }

        BOOST_TEST_MESSAGE( "result = " << result << "\n" );

        for ( int i = 0; i < FE::nDim; ++i )
            error[i] = ublas::norm_frobenius( der[i] - result );

        //error[i] = ublas::norm_frobenius( der[i] - fe.derivate( i ).evaluate( fe.points() ) );

#if defined(  USE_TEST )
        BOOST_TEST_MESSAGE( "[checkDiff] Hermite " << " dim : " << FE::nDim << " order : " << FE::nOrder << " error: " << error << " eps: " << M_eps  );

        for ( int i = 0; i < FE::nDim; ++i )
            BOOST_CHECK(  error[i] < M_eps );

#endif
    }

    void
    checkIdentity() const
    {
        BOOST_MESSAGE( "Checking Hermite polynomials identity" );
        const int nDim = fe.nDim;
        BOOST_TEST_MESSAGE( "points=" << fe.points() );
        BOOST_TEST_MESSAGE( "Hemite at points=" << fe.evaluate( fe.points() ) );
        matrix_type eval_at_pts( FE::polyset_type::toMatrix( fe.evaluate( fe.points() ) ) );
        matrix_type result( eval_at_pts.size1(), eval_at_pts.size2() );
        result = ublas::zero_matrix<value_type>( result.size1(), result.size2() );
        ublas::subrange( result, 0, nDim+1, 0, nDim+1 ) = ublas::identity_matrix<value_type>( nDim+1 );

        for ( int d = 0; d < nDim; ++d )
            ublas::subrange( result, 0, nDim+1, ( d+1 )*( nDim+1 ), ( d+2 )*( nDim+1 ) ) = ublas::identity_matrix<value_type>( nDim+1 );

        if ( nDim == 2 )
            result( fe.numPoints-1, fe.numPoints-1 ) = 1;

        BOOST_TEST_MESSAGE( "result = " << result << "\n" );
        BOOST_TEST_MESSAGE( "Hermite eval_at_pts = " << eval_at_pts  );
        value_type error = ublas::norm_frobenius( eval_at_pts-result );

#if defined(  USE_TEST )
        BOOST_TEST_MESSAGE( "[checkIdentity] Hermite "
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
#endif /* __TestHermite_H */

