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
    //# marker1 #
    using namespace Feel;
	//po::options_description qsnsoptions( "Navier-Stokes Oseen options" );
	//qsnsoptions.add_options()( "mu", po::value<double>()->default_value( 1.0 ), "coeff" );
	Environment env( _argc=argc, _argv=argv,
                     //_desc=qsnsoptions,
                     _about=about(_name="ns_oseen",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));


    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);

    auto Vh = THch<1>( mesh );
    auto U = Vh->element();
    auto V = Vh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();

    auto deft = sym(gradt( u ));
    auto def = sym(grad( v ));
    double mu = 0.0035;
    double rho = 1056;
    //double nu = option(_name="nu").as<double>();

    mpi::timer ti;
    auto g = expr( option(_name="functions.g").as<std::string>(), "g" );

    auto intUz = integrate(_range=markedfaces(mesh,"inlet"), _expr=g ).evaluate()(0,0) ;
    auto aireIn = integrate(_range=markedfaces(mesh,"inlet"),_expr=cst(1.)).evaluate()(0,0);
    auto meanU = intUz/aireIn;
    auto flow = integrate(_range=markedfaces(mesh,"inlet"), _expr=inner(g*N(),N())).evaluate()(0,0) ;
    auto reynolds = meanU*0.012/(0.0035/1056);
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
        std::cout<<"Nb elexments in throat surface" <<nelements(markedfaces(mesh,"throat"))<<"\n";

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
        auto meanUThroat=intUzt/aireInThroat;
        auto reynoldsThroat = meanUThroat*0.004/(mu*rho);
        auto flowT = integrate(_range=markedfaces(mesh,"throat"), _expr=inner(idv(u),N())).evaluate()(0,0) ;

        if (Environment::worldComm().isMasterRank())
        {
            std::cout<<"Integrale U throat = "<< intUzt;
            std::cout<<"   Airea of the throat's inlet = "<< aireInThroat;
            std::cout<<"   Mean U at throat's inlet  = "<< meanUThroat;
            std::cout<<"   Reynolds at the throat inlet = "<<reynoldsThroat;
            std::cout<<"   Flow at the throat's inlet= "<< flowT<<"\n";
            std::cout << "postproc time:  " << ti.elapsed() << "s\n";
        }
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
