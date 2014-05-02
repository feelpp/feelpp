/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel++ library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date     : Tue Feb 25 12:13:15 2014

   Copyright (C) 2014 Feel++ Consortium

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#include <feel/feel.hpp>

using namespace Feel;

class Oseen  : public Simget
{
public:
    Oseen( Application* a ) : app( a ) {}
    void run();
    void setOutputs( std::vector<std::string> const& o ) { M_outputs = o; }
    Application* app;
    std::vector<std::string> M_outputs;
};

void Oseen::run()
{
    boost::timer ti;
    ti.restart();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM>>);

    if ( Environment::isMasterRank() )
    {
        std::cout << "mesh loaded time:  " << ti.elapsed() << "s\n";
    }
    ti.restart();
    auto Vh = THch<1>( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();

    double mu = doption( "mu" );
    double rho = doption( "rho" );
    //auto mu = option(_name="parameters.mu").as<double>();
    //auto rho = option(_name="parameters.rho").as<double>();

    if ( Environment::isMasterRank() )
    {
        std::cout << "Total number of dof : " << Vh->nDof() << "\n";
        std::cout << "Total number of dof(local) : " << Vh->nLocalDof() << "\n";
        std::cout << "Velocity number of dof : " << Vh->functionSpace<0>()->nDof() << "\n";
        std::cout << "Velocity number of dof(local) : " << Vh->functionSpace<0>()->nLocalDof() << "\n";
        std::cout << "Pressure number of dof : " << Vh->functionSpace<1>()->nDof() << "\n";
        std::cout << "Pressure number of dof(local) : " << Vh->functionSpace<1>()->nLocalDof() << "\n";

        std::cout << "function space time:  " << ti.elapsed() << "s\n";
    }
    ti.restart();
    auto deft = sym(gradt( u ));
    auto def = sym(grad( v ));

    auto g = expr<FEELPP_DIM,1>( soption(_name="functions.g"), "g" );

    auto intUz = integrate(_range=markedfaces(mesh,"inlet"), _expr=g ).evaluate()(0,0) ;
    auto aireIn = integrate(_range=markedfaces(mesh,"inlet"),_expr=cst(1.)).evaluate()(0,0);
    auto meanU = intUz/aireIn;
    auto reynolds = mean(_range=markedfaces(mesh,"inlet"), _expr=rho*trans(g)*N()*h()/mu).norm();
    auto flow = integrate(_range=markedfaces(mesh,"inlet"), _expr=inner(g,N())).evaluate()(0,0) ;
    if ( Environment::isMasterRank() )
    {
        std::cout<<"  Integrale U = "<< intUz << "\n";
        std::cout<<"   Area Inlet = "<< aireIn << "\n";
        std::cout<<"       Mean U = "<< meanU << "\n";
        std::cout<<"         Flow = "<< flow << "\n";
        std::cout<<"     Reynolds = "<< reynolds <<"\n";
        std::cout << " time:  " << ti.elapsed() << "s\n";
    }
    ti.restart();

    auto mybdf = bdf( _space=Vh, _name="mybdf" );

    auto ft = form1( _test=Vh );

    auto a = form2( _trial=Vh, _test=Vh), at = form2( _trial=Vh, _test=Vh);
    a = integrate( _range=elements( mesh ), _expr=2*mu*inner( deft,def ) + mybdf->polyDerivCoefficient(0)*trans(rho*idt(u))*id(u) );
    a +=integrate( _range=elements( mesh ), _expr=-div( v )*idt( p ) - divt( u )*id( q ) );

    if ( Environment::isMasterRank() )
    {
        std::cout << "assembly a time:  " << ti.elapsed() << "s\n";
    }
    auto e = exporter( _mesh=mesh );
    if (Environment::worldComm().isMasterRank())
        std::cout<<"Nb of elements in throat surface " <<nelements(markedfaces(mesh,"throat"))<<" \n";

    for ( mybdf->start();  mybdf->isFinished() == false; mybdf->next(U) )
    {
        this->setMeshSize( mybdf->time() );
        ti.restart();

        auto bdf_poly = mybdf->polyDeriv();
        auto rhsu =  bdf_poly.element<0>();
        auto extrap = mybdf->poly();
        auto extrapu = extrap.element<0>();
        ft = integrate( _range=elements(mesh), _expr=(rho*trans(idv(rhsu))*id(u) ) );

        M_stats.put( "t.assembly.rhs", ti.elapsed() );

        ti.restart();

        at = a;
        at += integrate( _range=elements( mesh ), _expr= trans(rho*gradt(u)*idv(extrapu))*id(v) );

        M_stats.put( "t.assembly.lhs", ti.elapsed() );

        ti.restart();
        at+=on(_range=markedfaces(mesh,"wall"), _rhs=ft, _element=u,
               _expr=zero<FEELPP_DIM>() );
        at+=on(_range=markedfaces(mesh,"inlet"), _rhs=ft, _element=u,
               _expr=g );

        M_stats.put( "t.assembly.on", ti.elapsed() );

        ti.restart();

        auto r = at.solve(_rhs=ft,_solution=U);

        M_stats.put( "t.solve.total", ti.elapsed() );
        M_stats.put( "d.solve.bool.converged",r.isConverged() );
        M_stats.put( "d.solve.int.nit",r.nIterations() );
        M_stats.put( "d.solve.double.residual",r.residual() );

        ti.restart();

        for( auto marker : M_outputs )
        {
            ti.restart();

            auto intUz= integrate(_range=markedfaces(mesh,marker),_expr=idv(u)).evaluate()(0,0);
            auto area= integrate(_range = markedfaces(mesh,marker),_expr=cst(1.)).evaluate()(0,0);
            auto meanU = mean(_range=markedfaces(mesh,marker), _expr=idv(u)).norm();
            auto reynolds = mean(_range=markedfaces(mesh,marker), _expr=rho*idv(u)*h()/mu).norm();
            auto flowrate = integrate(_range=markedfaces(mesh,marker), _expr=inner(idv(u),N())).evaluate()(0,0) ;

            std::string key = (boost::format("d.%1%")%marker).str();
            std::string key2 = (boost::format("t.integrate.%1%")%marker).str();
            M_stats.put( key2, ti.elapsed() );
            M_stats.put( key+".double.area", area );
            M_stats.put( key+".double.intUz", intUz );
            M_stats.put( key+".double.meanU", meanU );
            M_stats.put( key+".double.reynolds", reynolds );
            M_stats.put( key+".double.flowrate", flowrate );
        }

        ti.restart();
        e->step(mybdf->time())->add( "u", u );
        e->step(mybdf->time())->add( "p", p );
        e->save();
        M_stats.put( "t.export.total", ti.elapsed() );

        M_stats.put( "h", mybdf->time() );
        M_stats.put( "level", mybdf->iteration() );
        app->storeStats( this->name(), M_stats );
        if ( Environment::isMasterRank() )
            app->printStats( std::cout, Application::ALL );
        M_stats.clear();
    }
}

int main( int argc, char** argv )
{
    po::options_description nsoseenoptions( "Navier-Stokes Oseen options" );
    nsoseenoptions.add_options()
        ( "mu", Feel::po::value<double>()->default_value( 1. ), "Dynamic viscosity" )
        ( "rho", Feel::po::value<double>()->default_value( 1000. ), "Fluid density" )
        ( "outputs", Feel::po::value<std::string>(), "list of face markers (space separated) on which some statistics are computed" )
        ;

    boost::mpi::timer ti;
	Environment env( _argc=argc, _argv=argv,
                     _desc=nsoseenoptions,
                     _about=about(_name="ns_oseen",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    if ( Environment::isMasterRank() )
    {
        std::cout << "Environment ok time:  " << ti.elapsed() << "s\n";
    }
    std::ofstream out;
    if ( env.isMasterRank() )
        out.open( (boost::format("res-%1%.dat") % env.numberOfProcessors() ).str().c_str() );

    Application benchmark;
    auto oseen = new Oseen(  &benchmark );
    benchmark.add( oseen );
    std::vector<std::string> keys = { "t.assembly",
                                      "t.assembly",
                                      "t.solve",
                                      "d.solve",
                                      "t.export"};
    std::vector<std::string> SplitVec;
    std::string outputs = soption( "outputs" );
    algorithm::split( SplitVec, outputs, algorithm::is_any_of(" ") );
    oseen->setOutputs( SplitVec );
    if ( SplitVec.size() )
        keys.push_back( "t.integrate" );
    for( auto marker :  SplitVec)
    {
        keys.push_back(  (boost::format("d.%1%")%marker).str() );
    }
    benchmark.setStats( keys );
    benchmark.run();
    benchmark.printStats( std::cout );
    benchmark.printStats( out );

}
