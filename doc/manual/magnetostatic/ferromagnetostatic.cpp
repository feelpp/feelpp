/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-28

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
#include <feel/feeldiscr/ned1h.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
    po::options_description magnetooptions( "Ferro-Magnetostatic options" );
    magnetooptions.add_options()
        ( "weakdirichlet", po::value<bool>()->default_value( false ), "false if strong Dirichlet, true if weak Dirichlet treatment" )
        ;

    /// [marker1]
    Environment env( _argc=argc, _argv=argv,
                     _desc=magnetooptions,
                     _about=about(_name="ferro-magnetostatic",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<FEELPP_DIM>> );

    auto Xh = Pchv<1>( mesh );
    auto Nh = Ned1h<0>( mesh );

    auto a = form2( _trial=Nh, _test=Nh );
    auto c = doption("parameters.c");
    auto f = expr<FEELPP_DIM,1>(soption("functions.f"), "f");
    auto e = expr<FEELPP_DIM,1>(soption("functions.e"), "e");
    //std::cout << "curl(curl(E)) + E = " << curl(curl(e)) << "\n";
    auto u = Nh->element();
    auto v = Nh->element(e);
    auto w = Xh->element(e);
    auto z = Xh->element();
    double penaldir=30;

    for( auto i : { 1, 2 ,3, 4, 5 } )
    {
#if FEELPP_DIM == 2
        a = integrate(_range=elements(mesh), _expr=c*trans(idt(u))*id(v)+curlxt(u)*curlx(v));

        if ( boption("weakdirichlet") ) //weak Dirichlet
            a += integrate(boundaryfaces(mesh), -curlxt(u)*(cross(N(),id(u)) )
                           - curlx(u)*(cross(N(),idt(u)) )
                           + penaldir*trans(cross(idt(u),N()))*cross(id(u),N())/hFace() );
#else //Dim = 3
        a = integrate(_range=elements(mesh), _expr=c*trans(idt(u))*id(v)+trans(curlt(u))*curl(v));

        if ( boption("weakdirichlet") ) //weak Dirichlet
            a += integrate(boundaryfaces(mesh), -trans(curlt(u))*(cross(N(),id(u)) )
                           - trans(curl(u))*(cross(N(),idt(u)) )
                           + penaldir*trans(cross(idt(u),N()))*cross(id(u),N())/hFace() );
#endif

        auto l = form1( _test=Nh );
        l = integrate( _range=elements(mesh), _expr=trans(f)*id(v));
        if ( boption("weakdirichlet") ) //weak Dirichlet
        {
#if FEELPP_DIM == 2
            l += integrate(boundaryfaces(mesh), - curlx(u)*(cross(N(),e) )
                           + penaldir*trans(cross(e,N()))*cross(id(u),N())/hFace() );
#else //Dim = 3
            l += integrate(boundaryfaces(mesh), - trans(curl(u))*(cross(N(),e) )
                           + penaldir*trans(cross(e,N()))*cross(id(u),N())/hFace() );
#endif
        }
        else //strong Dirichlet
        {
            a += on(_range=boundaryfaces(mesh), _rhs=l, _element=u,_expr=e );
        }

        a.solve(_rhs=l,_solution=u);
    }

    auto b = form2( _test=Xh, _trial=Xh );
    b = integrate(_range=elements(mesh), _expr=trans(idt(w))*id(w));
    auto ll = form1( _test=Xh );
    ll = integrate( _range=elements(mesh), _expr=trans(idv(u))*id(w));
    b.solve( _solution=z, _rhs=ll, _rebuild=true );

    auto ex = exporter( _mesh=mesh );
    ex->add( "u", u );
    ex->add( "Ihe", v );
    ex->add( "e", w );
    ex->add( "ul2", z );
    ex->save();
}
