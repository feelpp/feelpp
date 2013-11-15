// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="ls_laplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto syms = symbols<3>();
    auto g = option(_name="functions.g").as<std::string>();
    auto laplacian_g = laplacian( g, syms  );

    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=-expr( laplacian_g, syms )*id(v));

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=expr( g, syms ) );
    a.solve(_rhs=l,_solution=u);

    LOG(INFO) << "L2 error norm : " << normL2( _range=elements(mesh), _expr=idv(u)-expr( g, syms ) );
}
