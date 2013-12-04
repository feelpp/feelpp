// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>
#include <feel/feelvf/on.hpp>

int main(int argc, char**argv )
{
    //# marker1 #
    using namespace Feel;
    po::options_description desc("laplacian_hypercube");
    desc.add(feel_options()).add_options()( "mu",po::value<double>() -> default_value(1.),"mu" );
	Environment env( _argc=argc, _argv=argv,
                     _desc=desc,
                     _about=about(_name="laplacian_hypercube",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    //# endmarker1 #

    //# marker2 #
    auto mesh = loadMesh(_mesh=new Mesh<Hypercube<2>>);
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    //# endmarker2 #

    auto g = option(_name="functions.g").as<std::string>();
    auto vars = Symbols{"x","y"};
    //# marker3 #
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=-expr<6>( GiNaC::laplacian(g,vars),vars)*id(v));
                  //_expr=id(v));

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=expr( g, vars ) );
    if ( Environment::numberOfProcessors() == 1 )
        a.matrixPtr()->printMatlab("A.m");
    a.solve(_rhs=l,_solution=u);
    //# endmarker3 #

    auto l2 = normL2( _range=elements(mesh), _expr=idv(u)-expr<6>( g, vars ) );
    auto semih1 = normL2( _range=elements(mesh), _expr=gradv(u)-expr<1,2,6>( GiNaC::grad(g,vars), vars ) );
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "L2 : " << l2 << "\n";
        std::cout << "H1 : " << math::sqrt(l2*l2+semih1*semih1) << "\n";
    }
    on2(_range=boundaryfaces(mesh),_element=u,_expr=2*idv(u)).apply();
    //# marker4 #
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    v = vf::project( _space=Vh, _range=elements(mesh), _expr= expr( g, vars ) );
    e->add( "uexact", v );
    e->save();
    return 0;
    //# endmarker4 #
}
