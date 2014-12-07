/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Guillaume Dollé <guillaume.dolle@math.unistra.fr>

  Date: 2013-02-11

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)
                2013 Université de Strasbourg

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
   \file myfunctionspace.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>,
           Guillaume Dollé <guillaume.dolle@math.unistra.fr>
   \date 2010-07-15

   The myfunctionspace application compute integrals over a domain
   \see the \ref ComputingIntegrals section in the tutorial
   @author Christophe Prud'homme
*/
//! [all]
#include <feel/feel.hpp>
using namespace Feel;

int main( int argc, char** argv )
{

    //Initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="myfunctionspace",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" )  );

    //! [mesh]
    // create the mesh
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    //! [mesh]

    //! [space]
    // function space \( X_h \) using order 2 Lagrange basis functions
    auto Xh = Pch<2>( mesh );
    //! [space]

    //! [expression]
    auto g = expr<4>( soption(_name="functions.g"));
    auto gradg = grad<3>(g);
    //! [expression]

    //! [interpolant]
    // elements of \( u,w \in X_h \)
    auto u = Xh->element( "u" );
    auto w = Xh->element( "w" );
    
		// build the interpolant of u
    u.on( _range=elements( mesh ), _expr=g );
    
		// build the interpolant of the interpolation error
    w.on( _range=elements( mesh ), _expr=idv( u )-g );

    // compute L2 norms \(||\cdot||_{L^2}\)
    double L2g = normL2( elements( mesh ), g );
    double H1g = normL2( elements( mesh ), _expr=g,_grad_expr=gradg );
    double L2uerror = normL2( elements( mesh ), ( idv( u )-g ) );
    double H1uerror = normH1( elements( mesh ), _expr=( idv( u )-g ),
                              _grad_expr=( gradv( u )-gradg ) );
    if ( Environment::isMasterRank() )
    {
        std::cout << "||u-g||_0 = " << L2uerror/L2g << std::endl;
        std::cout << "||u-g||_1 = " << H1uerror/H1g << std::endl;
    }
    //! [interpolant]

    //! [export]
    // export for post-processing
    auto e = exporter( _mesh=mesh );
   
	 	// save interpolant
    e->add( "g", u );
   
	 	// save interpolant of interpolation error
    e->add( "u-g", w );

    e->save();
    //! [export]

}
//! [all]
