// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
	po::options_description laplacianoptions( "Laplacian options" );
	laplacianoptions.add_options()
		( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
		;
	Environment env( _argc=argc, _argv=argv,
                     _desc=laplacianoptions.add(feel_options()),
                     _about=about(_name="lapsubmesh",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    //# endmarker1 #

    //# marker2 #
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto submesh = createSubmesh( mesh, elements(mesh) );
    auto Vh = Pch<2>( mesh );
    auto Wh = Pch<2>( submesh );
    auto u = Vh->element();
    auto v = Wh->element();
    auto g = expr( option(_name="functions.g").as<std::string>() );
    auto l = form1( _test=Wh );
    l = integrate(_range=elements(submesh),
                  _expr=id(v));
    l+=integrate(_range=boundaryfaces(submesh),
                 _expr=g*(grad(v)*N()+id(v)/hFace()));

    auto a = form2( _trial=Vh, _test=Wh);
    a = integrate(_range=elements(submesh),
                  _expr=gradt(u)*trans(grad(v)) );

    a+= integrate( _range=boundaryfaces(submesh),
                   _expr=-gradt(u)*N()*id(v)+grad(v)*N()*idt(u)+idt(u)*id(v)/hFace());
    a.solve(_rhs=l,_solution=u);
    //# endmarker2 #


    //# marker4 #
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
    return 0;
    //# endmarker4 #
}
