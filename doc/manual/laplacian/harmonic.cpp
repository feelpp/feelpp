/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-03-13

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
   \file harmonic.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-03-13
 */
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _directory=".",
                     _about=about(_name="harmonic",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<2,2> >,
                                _desc=domain( _name="kovaznay",
                                              _usenames=false,
                                              _shape="hypercube",
                                              _h=option(_name="mesh2d.hsize").as<double>(),
                                              _xmin=-0.5, _xmax=1,
                                              _ymin=-0.5, _ymax=1.5 ) );
    auto Vh = Pchv<2>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=trace(gradt(u)*trans(grad(v))) );
    a+=on(_range=markedfaces(mesh,1), _rhs=l, _element=u,
          _expr=zero<2,1>() );
    a+=on(_range=markedfaces(mesh,3), _rhs=l, _element=u,
          _expr=zero<2,1>() );
    a+=on(_range=markedfaces(mesh,4), _rhs=l, _element=u,
          _expr=zero<2,1>() );
    a+=on(_range=markedfaces(mesh,2), _rhs=l, _element=u,
          _expr=vec(cst(0.),0.08*(Px()+0.5)*(Px()-1)*(Px()*Px()-1)));
    a.solve(_rhs=l,_solution=u);


    auto m1 = lagrangeP1(_space=Vh)->mesh();
    auto XhVisu = Pchv<1>(m1);

    auto opIVisu = opInterpolation(_domainSpace=Vh,
                                   _imageSpace=XhVisu,
                                   _type=InterpolationNonConforme(false,true,false) );
    auto uVisu = opIVisu->operator()(u);
    auto e = exporter( _mesh=m1, _name="initial" );
    e->step(0)->setMesh( m1 );
    e->step(0)->add( "u", uVisu );
    e->save();

    meshMove( m1, uVisu );

    auto e1 = exporter( _mesh=m1, _name="moved" );
    e1->step(0)->setMesh( m1  );
    e1->step(0)->add( "u", uVisu );
    e1->save();

#if 0
	auto e = exporter( _mesh=m1, _name="initial" );
    e->step(0)->setMesh( m1 );
    e->step(0)->add( "u", u );
    e->save();

    meshMove( mesh, u );

    auto Vh2 = Pchv<2>( mesh );
    auto m2 = lagrangeP1(_space=Vh2)->mesh();
    auto e1 = exporter( _mesh=m2, _name="moved" );
    e1->step(0)->setMesh( m2  );
    e1->step(0)->add( "u", u );
    e1->save();
#endif


    return 0;
}


