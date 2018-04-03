//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Vincent Chabannes <vincent.chabannes@feelpp.org>
//! @date 26 July 2017
//! @copyright 2017 Feel++ Consortium
//!

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description laplacianoptions( "Laplacian options" );
	laplacianoptions.add_options()
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
        ( "weakdir", po::value<bool>()->default_value( true ), "use weak Dirichlet condition" )
        ( "penaldir", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
		;

	Environment env( _argc=argc, _argv=argv,
                     _desc=laplacianoptions,
                     _about=about(_name="laplacian_meshsupport",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    typedef Mesh<Simplex<2,1>> mesh_type;
    auto mesh = loadMesh(_mesh=new mesh_type );
    typedef FunctionSpace< mesh_type,bases<Lagrange<2,Scalar> > > space_type;
    //auto Vh = Pch<2>( mesh );
    auto myrange = markedelements(mesh,"mat2");
    auto Vh= space_type::New( _mesh=mesh,_range=myrange );
    auto u = Vh->element("u");
    auto v = Vh->element();
    auto g = expr( soption(_name="functions.g"), "g" );

    auto l = form1( _test=Vh );
    l = integrate(_range=myrange,
                  _expr=id(v));

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=myrange,
                  _expr=gradt(u)*trans(grad(v)) );

    auto rangebf = markedfaces(mesh,{"line-h2","line-h3","line-vl2","line-vr2"});
    if ( !boption(_name="weakdir") )
    {
        a+=on(_range=rangebf, _rhs=l, _element=u, _expr=g );
    }
    else
    {
        double penaldir = doption(_name="penaldir");
        a += integrate( _range=rangebf,
                        _expr=  ( -( gradt( u )*N() )*id( v )
                                  -( grad( v )*N() )*idt( u )
                                  +penaldir*id( v )*idt( u )/hFace() ) );

        l += integrate( _range=rangebf,
                        _expr=g*( -grad( v )*N()
                                  + penaldir*id( v )/hFace() ) );

    }

    a.solve(_rhs=l,_solution=u);

    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "u", u );
    e->save();

    return 0;
}

