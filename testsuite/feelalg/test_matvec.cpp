/*
 * SPDX-FileContributor: Vincent Huber <vincent.huber@cemosis.fr>
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2008-2009 Université Joseph Fourier (Grenoble I)
 * SPDX-FileCopyrightText: 2013-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/**
   \file myblas.cpp
   \author Vincent HUBER <vincent.huber@cemosis.fr>
   \date 2013-12-06
 */

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelvf/vf.hpp>
#include <cmath>
using namespace Feel;

int main( int argc, char* argv[] )
{
    // initialize feel++ environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about( _name="myblas",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org") );

    auto mesh = unitSquare();

    auto m = mat<2,2>(cst(1.),cst(0.),cst(0.),cst(1.));
    auto v1 = vec(cst(1),cst(1));
    auto v2 = mat<2,1>(cst(1),cst(1));
    auto v3 = trans(v2); 
    auto v4 = trans(v1); 
    auto mv1 = m*v1; //OK
    auto mv2 = m*v2; //OK
    auto mv3 = v3*m; // row vector times matrix is the valid orientation
    auto mv4 = v4*m; // row vector times matrix is the valid orientation
    
    double int_1 = integrate( _range = elements( mesh ), _expr= trans(v1)*mv1 ).evaluate()(0,0);
    double int_2 = integrate( _range = elements( mesh ), _expr= trans(v2)*mv2 ).evaluate()(0,0);
    double int_3 = integrate( _range = elements( mesh ), _expr= mv3*trans(v3) ).evaluate()(0,0);
    double int_4 = integrate( _range = elements( mesh ), _expr= mv4*trans(v4) ).evaluate()(0,0);

    auto isClose = []( double value ) { return std::abs( value - 2.0 ) < 1e-12; };
    if ( !isClose( int_1 ) || !isClose( int_2 ) || !isClose( int_3 ) || !isClose( int_4 ) )
    {
        std::cerr << "unexpected row/column vector contraction values: "
                  << int_1 << '\t'
                  << int_2 << '\t'
                  << int_3 << '\t'
                  << int_4 << '\n';
        return 1;
    }
    return 0;
} 
