/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-04-13

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
   \file test_dubiner.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-04-13
 */
#include <iostream>

#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/dubinerba.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/operations.hpp>

static const int order = 7;
using namespace Feel;


template<typename PolySet>
void
writeEvalPoly( std::string const& name, PolySet const& q, int npts = 100 )
{
    typedef typename  PolySet::points_type points_type;
    points_type I( 1, npts );
    ublas::row( I, 0 ) = glas::linspace( -1.0, 1.0, npts );
    ublas::matrix<double> res( q.evaluate( I ) );
    //std::cout << "I = " << I << "\n";
    //std::cout << "res = " << res << "\n";
    std::ofstream ofs( name.c_str() );

    for ( int i = 0; i < npts; ++i )
    {
        ofs << I( 0, i ) << " ";

        for ( int j = 0; j < res.size1(); ++j )
            ofs << res( j, i ) << " ";

        ofs << "\n";
    }

}
template<typename Mat>
void
write( std::string const& name, Mat const& mass )
{
    std::ofstream mfs( name.c_str() );

    for ( int i = 0; i < mass.size1(); ++i )
    {
        for ( int j = 0; j < mass.size2(); ++j )
        {
            mfs << mass( i, j ) << " ";
        }

        mfs << "\n";
    }

    mfs << "\n";
}
template<typename PolySet1, typename PolySet2>
void
writeP( std::string const& name, PolySet1 const& p, PolySet2 const& q  )
{
    typedef typename  PolySet1::points_type points_type;
    int npts = order+1;
    points_type I( 1, npts );
    ublas::row( I, 0 ) = glas::linspace( -1.0, 1.0, npts );
    ublas::matrix<double> mp( ublas::trans( p.evaluate( I ) ) );
    ublas::matrix<double> mq( ublas::trans( q.evaluate( I ) ) );

    LU<ublas::matrix<double> > lu( mq );
    ublas::matrix<double> C = lu.solve( mp );
    ublas::matrix<double> P = ublas::trans( C );
    glas::clean( P, 1e-14 );
    std::cout << "P = " << P << "\n";

    write( name, P );
}
template<typename PolySet>
void
check( PolySet const& b )
{
    typedef typename  PolySet::points_type points_type;
    points_type p( 1,2 );
    p( 0,0 ) = -1;
    p( 0,1 )=1;

    IM_PK<1,2*order> im;

    std::cout << "p (-1) and p(1) = " << b.evaluate( p ) << "\n";
    std::cout << "mean(p) = " << ublas::prod( im.weights(),ublas::trans( b.evaluate( im.points() ) ) ) << "\n";

    for ( int i = 0; i < order+1; ++i )
        for ( int j = 0; j < order+1; ++j )
            std::cout << "int p_" << i << " p_" << j << "  = "
                      << ublas::prod( im.weights(),
                                      ublas::trans( ublas::element_prod( b.polynomial( i ).evaluate( im.points() ),
                                                    b.polynomial( j ).evaluate( im.points() ) ) ) ) << "\n";

}
int main()
{

    typedef BoundaryAdaptedPolynomialSet<1, order, Scalar, double, Simplex> bad_type;
    typedef bad_type::points_type points_type;

    bad_type b;
    points_type p( 1,2 );
    p( 0,0 ) = -1;
    p( 0,1 )=1;

    check( b );

    IM_PK<1,2*order> im;
    ublas::matrix<double> mass1( order+1, order+1 );

    for ( int i = 0; i < order+1; ++i )
        for ( int j = 0; j < order+1; ++j )
            mass1( i, j ) = ublas::prod( im.weights(),
                                         ublas::trans( ublas::element_prod( b.polynomial( i ).evaluate( im.points() ),
                                                 b.polynomial( j ).evaluate( im.points() ) ) ) )( 0 );


    write( "mass.dat", mass1 );

    writeEvalPoly( "evalbad.dat", b );

    typedef DubinerBoundaryAdapted<1,order,Scalar,double> dba_type;
    dba_type q;
    std::cout << "p (-1) and p(1) = " << q.evaluate( p ) << "\n";


    writeEvalPoly( "evaldba.dat", q );

    ublas::matrix<double> mass( ublas::prod( q.coeff(), ublas::trans( q.coeff() ) ) );
    write( "massnew.dat", mass );
    ublas::matrix<double> stiff( ublas::prod( dx( q ).coeff(), ublas::trans( dx( q ).coeff() ) ) );
    write( "stiffnew.dat", stiff );


    check ( q );

    writeP( "Pdba.dat", b, q );
    typedef OrthonormalPolynomialSet<1,order,Scalar,double> o_type;
    o_type o;
    writeP( "Po.dat", b, o );
}
