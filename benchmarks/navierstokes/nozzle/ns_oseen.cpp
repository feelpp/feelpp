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

int main(int argc, char**argv )
{
    boost::timer ti;
    using namespace Feel;
    //# marker1 #
    /* po::options_description nsoseenoptions( "Navier-Stokes Oseen options" );
    nsoseenoptions.add_options()
        ( "mu", Feel::po::value<double>()->default_value( 1. ), "Dynamic viscosity" )
        ( "rho", Feel::po::value<double>()->default_value( 1000. ), "Fluid density" );
        //return nsoseenoptions.add(Feel::feel_options()).add(Feel::backend_options( "oseen" ) );*/

	Environment env( _argc=argc, _argv=argv,
                     //_desc=nsoseenoptions,
                     _about=about(_name="ns_oseen",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    if ( Environment::isMasterRank() )
    {
        std::cout << "Environment ok time:  " << ti.elapsed() << "s\n";
    }
    ti.restart();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
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

    double mu = 0.0035;
    double rho = 1056;
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

    auto g = expr( option(_name="functions.g").as<std::string>(), "g" );

    auto intUz = integrate(_range=markedfaces(mesh,"inlet"), _expr=g ).evaluate()(0,0) ;
    auto aireIn = integrate(_range=markedfaces(mesh,"inlet"),_expr=cst(1.)).evaluate()(0,0);
    auto meanU = intUz/aireIn;
    auto flow = integrate(_range=markedfaces(mesh,"inlet"), _expr=inner(g*N(),N())).evaluate()(0,0) ;
    auto reynolds = meanU*0.012/(mu/rho);
    if ( Environment::isMasterRank() )
    {
        std::cout<<"  Integrale U = "<< intUz << "\n";
        std::cout<<"   Aire Inlet = "<< aireIn << "\n";
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
        ti.restart();
        if (Environment::worldComm().isMasterRank())
        {
            std::cout << "============================================================\n";
            std::cout << "Time : " << mybdf->time() << "s\n";
        }
        auto bdf_poly = mybdf->polyDeriv();
        auto rhsu =  bdf_poly.element<0>();
        auto extrap = mybdf->poly();
        auto extrapu = extrap.element<0>();
        ft = integrate( _range=elements(mesh), _expr=(rho*trans(idv(rhsu))*id(u) ) );
        if (Environment::worldComm().isMasterRank())
        {
            std::cout << "assembly rhs time:  " << ti.elapsed() << "s\n";
        }
        ti.restart();

        at = a;
        at += integrate( _range=elements( mesh ), _expr= trans(rho*gradt(u)*idv(extrapu))*id(v) );
        if (Environment::worldComm().isMasterRank())
        {
            std::cout << "assembly convect time:  " << ti.elapsed() << "s\n";
        }
        ti.restart();
        at+=on(_range=markedfaces(mesh,"wall"), _rhs=ft, _element=u,
               _expr=0*one() );
        at+=on(_range=markedfaces(mesh,"inlet"), _rhs=ft, _element=u,
               _expr=-g*N() );
        if (Environment::worldComm().isMasterRank())
        {
            std::cout << "assembly on time:  " << ti.elapsed() << "s\n";
        }
        ti.restart();
        at.solve(_rhs=ft,_solution=U);
        if (Environment::worldComm().isMasterRank())
        {
            std::cout << "solve time:  " << ti.elapsed() << "s\n";
        }
        ti.restart();
        auto intUzt= integrate(_range=markedfaces(mesh,"throat"),_expr=idv(u)).evaluate()(0,0);
        auto aireInThroat= integrate(_range = markedfaces(mesh,"throat"),_expr=cst(1.)).evaluate()(0,0);
        auto meanUThroat = mean(_range=markedfaces(mesh,"throat"), _expr=idv(u));
        auto reynoldsThroat = meanUThroat*0.004/(mu/rho);
        auto flowT = integrate(_range=markedfaces(mesh,"throat"), _expr=inner(idv(u),N())).evaluate()(0,0) ;

        if (Environment::worldComm().isMasterRank())
        {
            std::cout<<"Integrale U throat = "<< intUzt << "\n";
            std::cout<<"   Airea of the throat's inlet = "<< aireInThroat << "\n";
            std::cout<<"   Mean U at throat's inlet  = "<< meanUThroat(2,0) <<"\n";
            std::cout<<"   Flow at the throat's inlet= "<< flowT<<"\n";
            std::cout<<"   Reynolds at the throat inlet = "<<reynoldsThroat(2,0)<< "\n";
            std::cout << "        ==========================  \n";
        }

        auto intUzOut= integrate(_range=markedfaces(mesh,"outlet"),_expr=idv(u)).evaluate()(0,0);
        auto aireOut= integrate(_range = markedfaces(mesh,"outlet"),_expr=cst(1.)).evaluate()(0,0);
        auto meanUOut= mean(_range=markedfaces(mesh,"outlet"), _expr=idv(u));
        auto reynoldsOut = meanUOut*0.012/(mu/rho);
        auto flowOut = integrate(_range=markedfaces(mesh,"outlet"), _expr=inner(idv(u),N())).evaluate()(0,0) ;

        if (Environment::worldComm().isMasterRank())
        {
            std::cout<<"Integrale U outlet = "<< intUzOut;
            std::cout<<"   Airea of the outlet = "<< aireOut;
            std::cout<<"   Mean U at the outlet  = "<< meanUOut(2,0);
            std::cout<<"   Flow at the outlet = "<< flowOut<<"\n";
            std::cout<<"   Reynolds at the outlet = "<<reynoldsOut(2,0);
            std::cout << "postproc time:  " << ti.elapsed() << "s\n";
        }
        ti.restart();
        e->step(mybdf->time())->add( "u", u );
        e->step(mybdf->time())->add( "p", p );
        e->save();
        if (Environment::worldComm().isMasterRank())
        {
            std::cout << "export time:  " << ti.elapsed() << "s\n";
        }

    }
    return 0;
}
