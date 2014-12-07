/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Tue Feb 25 12:13:15 2014

   Copyright (C) 2014 Feel++ Consortium

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
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    //! [marker1]
    using namespace Feel;
	po::options_description stokesoptions( "Stokes options" );
	stokesoptions.add_options()
		( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
		;
	Environment env( _argc=argc, _argv=argv,
                     _desc=stokesoptions,
                     _about=about(_name="qs_stokes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    double meshSize = option(_name="gmsh.hsize").as<double>();
    GeoTool::Rectangle R( meshSize,"myRectangle",GeoTool::Node(0,0),GeoTool::Node(5,1));
    R.setMarker(_type="line",_name="inlet",_marker4=true);
    R.setMarker(_type="line",_name="outlet",_marker2=true);
    R.setMarker(_type="line",_name="wall",_marker1=true,_marker3=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);

    auto mesh = R.createMesh(_mesh=new Mesh<Simplex<2>>,_name="qs_stokes");

    auto g = expr<2,1>( soption(_name="functions.g") );
    auto Vh = THch<1>( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto v = V.element<0>(g,"poiseuille");
    auto p = U.element<1>();
    auto q = V.element<1>();

    auto deft = gradt( u );
    auto def = grad( v );
    double mu = doption(_name="mu");

    auto l = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate( _range=elements( mesh ), _expr=mu*inner( deft,def ) );
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) + divt( u )*id( q ) );

    a+=on(_range=markedfaces(mesh,"wall"), _rhs=l, _element=u,
          _expr=zero<2,1>() ) ;
    a+=on(_range=markedfaces(mesh,"inlet"), _rhs=l, _element=u,
          _expr=g );

    a.solve(_rhs=l,_solution=U);

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->add( "v", v );
    e->add( "p", p );
    e->save();
    //! [marker1]

    return 0;
}
