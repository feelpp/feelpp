#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    // tag::env[]
    using namespace Feel;
    using Feel::cout;
	po::options_description geimcaseoptions( "GEIM CASE options" );
	geimcaseoptions.add_options()
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
        ( "alpha", po::value<double>()->default_value( 1.0 ), "alpha" )
        ( "beta", po::value<double>()->default_value( 1.0 ), "beta" )
        ( "gamma", po::value<double>()->default_value( 1.0 ), "gamma" )
		;

	Environment env( _argc=argc, _argv=argv,
                   _desc=geimcaseoptions,
                   _about=about(_name="qs_geim_case",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));
    // end::env[]

    // tag::mesh_space[]
    tic();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2,1>>);
    toc("loadMesh");

    tic();
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element("u");
    auto v = Vh->element("v");
    auto g = expr(soption(_name="functions.g"), "g" );
    auto alp = doption(_name="alpha");
    auto bet = doption(_name="beta");
    auto gam = doption(_name="gamma");
    toc("Vh");
    // end::mesh_space[]
    // tag::forms[]
    tic();
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),_expr=id(v));
    l += integrate(_range=markedelements(mesh,"Omega1"), _expr=(alp*sin(Px())+bet*cos(gam*Py()*pi))*id(v),_quad=4);
    toc("l");

    tic();
    auto a = form2( _trial=Vh, _test=Vh);
    tic();
    a = integrate(_range=elements(mesh),
                  _expr=inner(gradt(u),grad(v)),_quad=_Q<>(3) );
    toc("a");
    a+=on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );

    // end::forms[]

    // tag::solve[]
    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
        a.solve(_rhs=l,_solution=u);
    toc("a.solve");
    // end::solve[]

    // tag::export[]
    tic();
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "uh", u );
    e->save();
    toc("Exporter");
    // end::export[]

    return 0;

}
// end::global[]
