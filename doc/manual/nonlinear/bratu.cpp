/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-01-09

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2013 Feel++ Consortium

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
   \file bratu.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-02-04
 */
#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description bratuoptions( "Bratu problem options" );
    bratuoptions.add_options()
    ( "lambda", Feel::po::value<double>()->default_value( 1 ), "exp() coefficient value for the Bratu problem" )

    ( "penalbc", Feel::po::value<double>()->default_value( 30 ), "penalisation parameter for the weak boundary conditions" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )

    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return bratuoptions.add( Feel::feel_options() );
}

/**
 * Bratu Problem
 *
 * solve \f$ -\Delta u + \lambda \exp(u) = 0, \quad u_\Gamma = 0\f$ on \f$\Omega\f$
 */
int
main( int argc, char** argv )
{

    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="bratu",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org"));
    auto mesh = unitSquare();
    auto Vh = Pch<3>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    double penalbc = option(_name="penalbc").as<double>();
    double lambda = option(_name="lambda").as<double>();

    auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            if (!J) J = backend()->newMatrix( Vh, Vh );
            auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
            a = integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) );
            a += integrate( elements( mesh ), lambda*( exp( idv( u ) ) )*idt( u )*id( v ) );
            a += integrate( boundaryfaces( mesh ),
                            ( - trans( id( v ) )*( gradt( u )*N() )
                              - trans( idt( u ) )*( grad( v )*N() )
                              + penalbc*trans( idt( u ) )*id( v )/hFace() ) );
        };
    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            auto u = Vh->element();
            u = *X;
            auto r = form1( _test=Vh, _vector=R );
            r = integrate( elements( mesh ), gradv( u )*trans( grad( v ) ) );
            r +=  integrate( elements( mesh ),  lambda*exp( idv( u ) )*id( v ) );
            r +=  integrate( boundaryfaces( mesh ),
                             ( - trans( id( v ) )*( gradv( u )*N() )
                               - trans( idv( u ) )*( grad( v )*N() )
                               + penalbc*trans( idv( u ) )*id( v )/hFace() ) );
        };
    u.zero();
    backend()->nlSolver()->residual = Residual;
    backend()->nlSolver()->jacobian = Jacobian;
    backend()->nlSolve( _solution=u );

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
}





