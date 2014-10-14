/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-04-17

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
#include <feel/feelalg/solvereigen.hpp>
int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	po::options_description stokesoptions( "Stokes options" );
	stokesoptions.add_options()
		( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
		;
	Environment env( _argc=argc, _argv=argv,
                     _desc=stokesoptions,
                     _about=about(_name="ge_stokes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));


    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);

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
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) - divt( u )*id( q ) );
    a += integrate( _range=elements( mesh ), _expr=1e-6*idt(p)*id(q) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=zero<2,1>() ) ;

    auto b = form2( _trial=Vh, _test=Vh);
    b = integrate( _range=elements( mesh ), _expr=trans(idt(u))*id(v)+0.*idt(p)*id(p) );

    auto modes= veigs( _formA=a, _formB=b );


    auto e = exporter( _mesh=mesh );
    int i = 0;
    for( auto const& mode: modes )
    {
        LOG_IF(WARNING, (i==0)&&( math::abs( mode.first - 52.34 ) > 1e-1 ) )
            << "Invalid stokes first eigen value " << mode.first << " should be " << 52.34;
        u = mode.second.element<0>();
        p = mode.second.element<1>();
        e->add( ( boost::format( "mode-u-%1%" ) % i ).str(), u );
        e->add( ( boost::format( "mode-p-%1%" ) % i ).str(), p );
        auto normL2Div = normL2( _range=elements(mesh), _expr=divv(u) );
        if ( Environment::isMasterRank() )
        {
            std::cout << "Lambda_" << i << " = " <<  mode.first << " Divergence = " << normL2Div << "\n";
        }
        i++;
    }
    e->save();

    return 0;
}
