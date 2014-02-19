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
                     _about=about(_name="qs_stokes",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    double meshSize = option(_name="gmsh.hsize").as<double>();
    GeoTool::Rectangle R( meshSize,"myRectangle",GeoTool::Node(0,0),GeoTool::Node(5,1));
    R.setMarker(_type="line",_name="inlet",_marker4=true);
    R.setMarker(_type="line",_name="outlet",_marker2=true);
    R.setMarker(_type="line",_name="wall",_marker1=true,_marker3=true);
    R.setMarker(_type="surface",_name="Omega",_markerAll=true);

    auto mesh = R.createMesh(_mesh=new Mesh<Simplex<2>>,_name="qs_stokes");

    auto Vh = THch<1>( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();

    auto deft = gradt( u );
    auto def = grad( v );
    double mu = option(_name="mu").as<double>();

    auto l = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh);
    a = integrate( _range=elements( mesh ), _expr=mu*inner( deft,def ) );
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) + divt( u )*id( q ) );

    a+=on(_range=markedfaces(mesh,"wall"), _rhs=l, _element=u,
          _expr=0*one() );
    a+=on(_range=markedfaces(mesh,"inlet"), _rhs=l, _element=u,
          _expr=-expr( option(_name="functions.g").as<std::string>() )*N() );

    a.solve(_rhs=l,_solution=U);

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->add( "p", p );
    e->save();

    return 0;
}

