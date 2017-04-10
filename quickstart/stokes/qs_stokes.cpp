/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description laplacianoptions( "Stokes options" );
	laplacianoptions.add_options()
        ( "mu", po::value<double>()->default_value( 1.0 ), "viscosity" )
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
		;

	Environment env( _argc=argc, _argv=argv,
                   _desc=laplacianoptions,
                   _about=about(_name="qs_stokes",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    tic();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    toc("loadMesh");

    tic();
    auto Vh = THch<1>( mesh );
    auto U = Vh->element("U");
    auto u = U.element<0>();
    auto p = U.element<1>();
    auto v = U.element<0>();
    auto q = U.element<1>();
    auto mu = doption(_name="mu");
    auto f = expr<FEELPP_DIM,1>( soption(_name="functions.f"), "f" );
    auto g = expr<FEELPP_DIM,1>( soption(_name="functions.g"), "g" );
    toc("Vh");

    tic();
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=inner(f,id(v)));
    toc("l");

    tic();
    auto a = form2( _trial=Vh, _test=Vh);
    auto Id = eye<FEELPP_DIM,FEELPP_DIM>();
    auto deft = sym(gradt(u));
    auto sigmat = -idt(p)*Id + 2*mu*deft;
    a = integrate(_range=elements(mesh),
                  _expr=inner( sigmat, grad(v) ) );
    a += integrate(_range=elements(mesh),
                   _expr=id(q)*divt(u) );
    a+=on(_range=markedfaces(mesh,"inlet"), _rhs=l, _element=u, _expr=g );
    a+=on(_range=markedfaces(mesh,{"wall","letters"}), _rhs=l, _element=u, _expr=zero<FEELPP_DIM,1>() );
    toc("a");

    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
        a.solve(_rhs=l,_solution=U);
    toc("a.solve");

    tic();
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "u", u );
    e->add( "p", p );
    e->save();
    toc("Exporter");
    return 0;
}
