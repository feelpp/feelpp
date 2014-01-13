// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="qs_laplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    //# endmarker1 #

    //# marker2 #
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    //# endmarker2 #

    //# marker3 #
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v));

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=constant(0.) );
    if ( Environment::numberOfProcessors() == 1 )
    {
        a.matrixPtr()->printMatlab("A.m");
        l.vectorPtr()->printMatlab("b.m");
    }
    a.solve(_rhs=l,_solution=u);
    //# endmarker3 #
    if ( Environment::numberOfProcessors() == 1 )
    {
        Vh->dof()->pointIdToDofRelation( "feelpp2msh.m" );
        u.printMatlab( "ug.m", true );
        u.printMatlab( "uf.m", false );
    }
    //# marker4 #
    auto e = exporter( _mesh=mesh );

    e->add( "u", u );

    e->save();
    return 0;
    //# endmarker4 #
}
