#include <feel/feel.hpp>

using namespace Feel;

inline po::options_description makeOptions()
{
	po::options_description laplacianoptions( "Laplacian options" );
	laplacianoptions.add_options()
		( "nu", po::value<double>()->default_value( 1.0 ), "viscosity" )
		;
	return laplacianoptions.add( Feel::feel_options() );
}

int main(int argc, char**argv )
{
	using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
			_desc=makeOptions(),
			_about=about(_name="qs_laplacian_onefeel",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));

	double nu = option(_name="nu").as<double>();
	
	auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
	auto Vh = Pch<2>( mesh );
	auto u = Vh->element();
	auto v = Vh->element();

	auto l = form1( _test=Vh );
	l = integrate(_range=elements(mesh),
			_expr=id(v));

	auto a = form2( _trial=Vh, _test=Vh);
	a = integrate(_range=elements(mesh),
			_expr=cst(nu)*gradt(u)*trans(grad(v)) );
	a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
			_expr=expr( option(_name="functions.g").as<std::string>(), symbols<2>() ) );
	a.solve(_rhs=l,_solution=u);

	auto e = exporter( _mesh=mesh );
	e->add( "u", u );
	e->save();
	return 0;
}
