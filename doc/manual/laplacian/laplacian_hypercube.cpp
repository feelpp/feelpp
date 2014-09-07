// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>


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
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    //# endmarker2 #

    auto g = option(_name="functions.g").as<std::string>();
    auto vars = Symbols{"x","y"};
    auto lap = GiNaC::laplacian(g,vars);
    double intg = integrate( _range=elements(mesh), _expr=expr<6>( g, vars ) ).evaluate()(0,0);
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "lap : " << lap <<"\n";
        std::cout << "integral : " << intg <<  "\n";
    }
    //# marker3 #
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=-expr<6>( lap,vars)*id(v));
    //_expr=id(v));


    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)));
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=expr( g, vars ) );
    if ( Environment::numberOfProcessors() == 1 )
        a.matrixPtr()->printMatlab("A.m");

    auto b = form2( _trial=Vh, _test=Vh);
    b = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    if ( Environment::numberOfProcessors() == 1 )
        b.matrixPtr()->printMatlab("B.m");
    auto c = form2( _trial=Vh, _test=Vh);
    c = integrate(_range=elements(mesh),
                  _expr=idt(u)*id(v));
    auto d = form2( _trial=Vh, _test=Vh);
    d = integrate(_range=elements(mesh),
                  _expr=idt(u)*id(v)+gradt(u)*trans(grad(v)));

    auto m = form1( _test=Vh );
    m = integrate(_range=elements(mesh), _expr=expr( g,vars)*id(v) );
    c.solve(_rhs=m,_solution=v );
    auto l2proj = normL2( _range=elements(mesh), _expr=idv(v)-expr<6>( g, vars ) );
    d.solve(_rhs=m,_solution=v );
    auto h1proj = normL2( _range=elements(mesh), _expr=idv(v)-expr<6>( g, vars ) );

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "l2 proj : " << l2proj <<"\n";
        std::cout << "h1 proj : " << h1proj <<"\n";
    }

    a.solve(_rhs=l,_solution=u);
    //# endmarker3 #

    v = project( _space=Vh, _range=elements(mesh), _expr=expr(g,vars));
    auto interpl2 = normL2( _range=elements(mesh), _expr=idv(v)-expr<6>( g, vars ));
    double bvv = b(v,v);
    double buu = b(u,u);
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "b(v,v) : " << bvv <<"\n";
        std::cout << "b(u,u) : " << buu <<"\n";
    }

    auto l2 = normL2( _range=elements(mesh), _expr=idv(u)-expr<6>( g, vars ));
    auto semih1 = normL2( _range=elements(mesh), _expr=gradv(u)-expr<1,2,6>( GiNaC::grad(g,vars), vars ));
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "interpL2 : " << interpl2 << "\n";
        std::cout << "L2 : " << l2 << "\n";
        std::cout << "H1 : " << math::sqrt(l2*l2+semih1*semih1) << "\n";
    }
    //on2(_range=boundaryfaces(mesh),_element=u,_expr=2*idv(u)).apply();
    //# marker4 #
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    v = vf::project( _space=Vh, _range=elements(mesh), _expr= expr( g, vars ) );
    e->add( "uexact", v );
    if ( Vh->nLocalDof() < 30 )
    {
        // save the basis function to display them
        for( auto d : Vh->dof()->localDof() )
        {
            v.zero();
            v[d.second.index()] =1;
            e->add( (boost::format("basis-%1%") % d.second.index()).str(), v );
        }
    }
    e->save();
    return 0;
    //# endmarker4 #
}
