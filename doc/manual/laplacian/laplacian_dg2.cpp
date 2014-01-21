/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-03-14

  Copyright (C) 2013 Université de Strasbourg

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
   \file laplacian_dg.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-03-14
 */
#include <feel/feel.hpp>

namespace Feel
{
/**
   \page LaplacianDG Laplacian using Discontinous Galerkin
   \author Feel++ Consortium

   <br>
   <br>

   `feelpp_doc_laplacian_dg` solves for the Laplacian in a square using a modal
   basis (Dubiner) using a DG formulation.

   \section LaplacianDG_Implementation
   the implementation is available in \ref doc/manual/laplacian/laplacian_dg.cpp
   \snippet laplacian_dg.cpp marker1
 */
}
int main(int argc, char**argv )
{
    /// [marker1]
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="laplacian_dg2",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    //auto mesh = unitSquare();
    auto mesh = loadMesh( _mesh=new Mesh<Hypercube<2>> );
    //auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    //auto mesh = unitCube();

    //auto Vh = Odh<1>( mesh );
    auto Vh = Pdh<1>( mesh, true );
    auto Xh = Pch<1>( mesh );

    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v));

    auto c = form2( _trial=Vh, _test=Vh,
                    _pattern=size_type(Pattern::EXTENDED) );
    //int ninternalfaces = nelements(internalfaces(mesh));
    c =integrate( internalfaces( mesh ),
                  + trans( jumpt( cst(1.0)/4 ) )*jump( cst( 1.0 )/4 ) / measFace() );
    // \int_Fint [ mean(u) ] \cdot [ mean(v) ] = \int_Fint [ 1/4 ] \cdot [ 1/4 ] / | F |
    // \int_Fint 1 / |F| = 1 ! = \sum_{F \in Fint} \int_F 1/|F|  = #{F\in Fint}
    if ( Environment::numberOfProcessors() == 1 )
        c.matrixPtr()->printMatlab( "c.m" );
    auto a = form2( _trial=Vh, _test=Vh,
                    _pattern=size_type(Pattern::EXTENDED) );
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a +=integrate( internalfaces( mesh ),
                   + trans( jumpt( cst(1.0)/4 ) )*jump( cst( 1.0 )/4 ) / (measFace()) );
    a += on( _range=boundaryfaces(mesh), _element=u, _rhs=l, _expr=cst(0.));

    a.solve(_rhs=l,_solution=u);

    auto p = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2 );
    auto uc = p->project( idv(u) );

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->add( "uc", uc );
    e->save();

    return 0;
}
