/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2013-10-22

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file test_slepc.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-10-22
*/

#define BOOST_TEST_MODULE test_slepc
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <feel/feel.hpp>
#include <feel/feelalg/solvereigen.hpp>

/** use Feel namespace */
using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description slepc( "test slepc options" );
    slepc.add_options()
        ("solver-type" , Feel::po::value<int>()->default_value( 5 ), "solver eigen type " )
        ( "nev" ,  Feel::po::value<int>()->default_value( 1 ), "nev " )
        ( "ncv" ,  Feel::po::value<int>()->default_value( 3 ), "ncv " )
        ( "tol" ,  Feel::po::value<double>()->default_value( 1e-10 ), "tolerance " )
        ( "maxiter" ,  Feel::po::value<int>()->default_value( 10000 ), "maxiter " )
        ;
    return  slepc.add( Feel::feel_options() ) ;
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_slepc" ,
                     "test_slepc" ,
                     "0.2",
                     "nD(n=1,2,3) test context of functionspace",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}



template<int Dim, int Order>
void
testSlepc()
{
    auto mesh=unitHypercube<Dim>();
    auto Xh = Pch<Order>( mesh );
    LOG(INFO)<<"nDof : "<<Xh->nLocalDof();
#if 0
    auto backend = backend_type::build( BACKEND_PETSC );
#else
    auto backend = backend_type::build( );
#endif

    auto u = Xh->element();
    auto v = Xh->element();

    auto A = backend->newMatrix( _test=Xh, _trial=Xh );
    auto B = backend->newMatrix( _test=Xh, _trial=Xh );
    form2( _test=Xh, _trial=Xh, _matrix=A ) = integrate( _range= elements( mesh ), _expr= idt( u )*id(  v  ) );
    form2( _test=Xh, _trial=Xh, _matrix=B ) = integrate( _range= elements( mesh ), _expr= idt( u )*id(  v  ) );
    A->close();
    B->close();
    SolverEigen<double>::eigenmodes_type modes;

    modes=
        eigs( _matrixA=A,
              _matrixB=B,
              _solver=( EigenSolverType )ioption(_name="solver-type"),
              _spectrum=SMALLEST_REAL,
              _transform=SINVERT,
              _ncv=ioption(_name="ncv"),
              _nev=ioption(_name="nev"),
              _tolerance=doption(_name="tol"),
              _maxit=ioption(_name="maxiter")
              );

    double eigen_value = modes.begin()->second.template get<0>();
    BOOST_CHECK_SMALL( (eigen_value-1), 1e-12 );

    auto eigen_vector = modes.begin()->second.template get<2>();
    auto Aw =  backend->newVector( Xh );
    auto Bw =  backend->newVector( Xh );
    A->multVector( eigen_vector, Aw );
    B->multVector( eigen_vector, Bw );
    //we should have Aw = eigen_value Bw
    Bw->scale( eigen_value );
    double energyAwAw = A->energy( Aw , Aw );
    double energyAwBw = A->energy( Aw , Bw );
    double energyBwBw = A->energy( Bw , Bw );
    BOOST_CHECK_SMALL( math::abs(energyAwAw-energyAwBw) , 1e-12 );
    BOOST_CHECK_SMALL( math::abs(energyAwAw-energyBwBw) , 1e-12 );
}


/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( slepc )

BOOST_AUTO_TEST_CASE( test_1 )
{
    testSlepc<2,1>();
}

BOOST_AUTO_TEST_SUITE_END()


