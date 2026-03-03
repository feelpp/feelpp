// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

template<int Dim>
Feel::Symbols makeSymbols()
{
    if constexpr ( Dim == 2 )
        return Feel::Symbols{"x","y"};
    else
        return Feel::Symbols{"x","y","z"};
}

template<int Dim, int Order>
void runLaplacianHypercube()
{
    using namespace Feel;

    auto mesh = loadMesh( _mesh=new Mesh<Hypercube<Dim>> );
    auto Vh = Pch<Order,double,PointSetEquiSpaced>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto g = option( _name="functions.g" ).as<std::string>();
    auto vars = makeSymbols<Dim>();
    auto lap = GiNaC::laplacian( g, vars );
    auto gradg = GiNaC::grad( g, vars );

    auto cacheSeed = g + "_d" + std::to_string( Dim ) + "_p" + std::to_string( Order );
    auto cacheId = std::to_string( std::hash<std::string>{}( cacheSeed ) );
    auto gexpr = expr<8>( g, vars, "gexpr_" + cacheId );
    auto lapexpr = expr<1,1,8>( lap, vars, "lapg", "lapexpr_" + cacheId );
    auto gradexpr = expr<1,Dim,8>( gradg, vars, "gradg", "gradexpr_" + cacheId );
    auto gradexprvisu = expr<Dim,1,8>( gradg, vars, "gradg_visu", "gradvisu_" + cacheId );

    auto l = form1( _test=Vh );
    l = integrate( _range=elements(mesh), _expr=-lapexpr*id(v), _quad=_Q<8>() );

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements(mesh), _expr=gradt(u)*trans(grad(v)), _quad=_Q<8>() );
    a += on( _range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=gexpr, _type="elimination" );
    a.solve( _rhs=l, _solution=u );

    auto l2 = normL2( _range=elements(mesh), _expr=idv(u)-gexpr, _quad=_Q<8>() );
    auto semih1 = normL2( _range=elements(mesh), _expr=gradv(u)-gradexpr, _quad=_Q<8>() );
    auto h1 = math::sqrt( l2*l2 + semih1*semih1 );
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "dim/order             : " << Dim << "D/Q" << Order << "\n";
        std::cout << "Exact solution g      : " << g << "\n";
        std::cout << "lap(g)                : " << lap << "\n";
        std::cout << "grad(g)               : " << gradg << "\n";
        std::cout << "L2 error  ||u-g||_L2  = " << l2 << "\n";
        std::cout << "H1 error  ||u-g||_H1  = " << h1 << "\n";
    }

    auto e = exporter( _mesh=mesh );
    auto lapg = Vh->element();
    lapg.on( _range=elements(mesh), _expr=lapexpr );
    auto Vhv = Pchv<Order,PointSetEquiSpaced>( mesh );
    auto gradgvisu = Vhv->element();
    gradgvisu.on( _range=elements(mesh), _expr=gradexprvisu );
    e->add( "u", u );
    v = vf::project( _space=Vh, _range=elements(mesh), _expr=gexpr );
    e->add( "uexact", v );
    e->add( "lapg", lapg );
    e->add( "gradg", gradgvisu );
    e->save();
}

int main(int argc, char**argv )
{
    using namespace Feel;
    po::options_description desc("laplacian_hypercube");
    desc.add( feel_options() ).add_options()
        ( "case.dim", po::value<int>()->default_value( 2 ), "dimension: 2 or 3" )
        ( "case.order", po::value<int>()->default_value( 1 ), "polynomial order: 1 (Q1) or 2 (Q2)" );
    Environment env( _argc=argc, _argv=argv,
                     _desc=desc,
                     _about=about(_name="laplacian_hypercube",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    int dim = option( _name="case.dim" ).as<int>();
    int order = option( _name="case.order" ).as<int>();

    if ( dim == 2 && order == 1 )
        runLaplacianHypercube<2,1>();
    else if ( dim == 2 && order == 2 )
        runLaplacianHypercube<2,2>();
    else if ( dim == 3 && order == 1 )
        runLaplacianHypercube<3,1>();
    else if ( dim == 3 && order == 2 )
        runLaplacianHypercube<3,2>();
    else
        CHECK( false ) << "Unsupported case.dim/case.order = " << dim << "/" << order << " (expected dim in {2,3}, order in {1,2})";

    return 0;
}
