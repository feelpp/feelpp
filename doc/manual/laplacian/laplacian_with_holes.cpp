// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/exporter.hpp>


int main(int argc, char**argv )
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _directory=".",
                     _about=about(_name="qs_laplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    typedef Mesh<Simplex<2> > mesh_type;

    GeoTool::Rectangle R1( 0.05,"R1",GeoTool::Node( 0,0 ),GeoTool::Node( 1,1 ) );
    GeoTool::Circle C1( 0.05,"C1",GeoTool::Node( 0.5,0.5 ),GeoTool::Node( 0.75,0.75 ) );

    auto R1mesh = R1.createMesh(_mesh=new mesh_type,_name="R1" );
    auto C1mesh = C1.createMesh(_mesh=new mesh_type,_name="C1" );
    auto R1mC1mesh = ( R1-C1 ).createMesh(_mesh=new mesh_type,_name="R1-C1" );

    // auto mesh=R1mesh;
	// auto mesh=C1mesh;
    auto mesh=R1mC1mesh;

    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v));


    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=constant(0.) );
    a.solve(_rhs=l,_solution=u);

	cout<<"equazione risolta"<<endl;

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();

    return 0;
}

