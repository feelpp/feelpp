/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

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
   \file generecrash.cpp
   \author Vincent HUBER <vincent.huber@cemosis.fr>
   \date 2013-08-22
 */
#include <feel/feel.hpp>
using namespace Feel;

int
main( int argc, char** argv )
{
    // Initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about( _name="test_on_inside" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    // create the mesh (specify the dimension of geometric entity)
    typedef Mesh<Simplex<2>> mesh_type;
    typedef FunctionSpace<mesh_type,bases<Lagrange<2,Scalar>>,double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    auto mesh = loadMesh(_mesh = new mesh_type );
    space_ptrtype Xh;
    Xh = space_type::New(mesh);
    auto u = Xh->element( "u" );
    auto v = Xh->element( "v" );

    auto l = form1(_test=Xh);
    auto a = form2(_test=Xh, _trial=Xh);

    l += integrate(_range=elements( mesh ), _expr=cst(1.));
    a += integrate(_range=elements( mesh ), _expr=idt(u)*id(v));
    a += on(_range=markedfaces(mesh,"toto"),_expr=cst(1),_rhs=l,_element=u);
    double _a_ = integrate( _range = markedfaces(mesh,"toto"),_expr = cst(1.0) ).evaluate()( 0,0 );
    if(Feel::detail::Environment::rank() == 0) std::cout << "_a_ = " << _a_ << std::endl;

    a += on(_range=boundaryfaces(mesh),_expr=cst(0.1),_rhs=l,_element=u);
    a += on(_range=markedfaces(mesh,"toto"),_expr=cst(0.1),_rhs=l,_element=u);
    
}


