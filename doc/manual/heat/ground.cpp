/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2011-07-25

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
/**
   \file ground.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-07-25
 */

#include <feel/feel.hpp>

//# marker1 #
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description groundoptions("ground options");
    groundoptions.add_options()
        // mesh parameters
        ("D", Feel::po::value<double>()->default_value( 1 ),"depth")
        ("W", Feel::po::value<double>()->default_value( 0.5 ),"width")
        ("soil.0.k", Feel::po::value<double>()->default_value( 0.2 ),"conductivity")
        ("soil.0.c", Feel::po::value<double>()->default_value( 1 ),"capacity")
        ("soil.0.rho", Feel::po::value<double>()->default_value( 1 ),"density")
        ("soil.1.k", Feel::po::value<double>()->default_value( 0.2 ),"conductivity")
        ("soil.1.c", Feel::po::value<double>()->default_value( 1 ),"capacity")
        ("soil.1.rho", Feel::po::value<double>()->default_value( 1 ),"density")
        ("T0", Feel::po::value<std::string>()->default_value( "TR + TA*sin(omega*t):x:y:z:t:TR:TA:omega" ),"T0")
        ("TA", Feel::po::value<double>()->default_value( 1 ),"TA")
        ("TR", Feel::po::value<double>()->default_value( 0 ),"TR")
        ("omega", Feel::po::value<double>()->default_value( 2*M_PI ),"w")

        ("gamma", Feel::po::value<double>()->default_value( 20 ),"gamma");


    return groundoptions;
}
//# endmarker1 #

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "ground" ,
                           "ground" ,
                           "0.1",
                           "nD(n=1,2,3) heating and cooling of the ground due to surface temperature",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2014 Feel++ Consortium");
    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "");
    return about;
}

int main(int argc, char** argv )
{
    using namespace Feel;

    const int Dim = FEELPP_DIM;
    const int Order = FEELPP_ORDER;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    Environment::changeRepository( boost::format( "%1%/%2%D/P%3%/h_%4%" )
                                   % Environment::about().appName() % Dim % Order
                                   % option(_name="gmsh.hsize").as<double>() );

    auto W = option(_name="W").as<double>();
    auto D = option(_name="D").as<double>();
    auto k0 = option(_name="soil.0.k").as<double>();
    auto c0 = option(_name="soil.0.c").as<double>();
    auto rho0 = option(_name="soil.0.rho").as<double>();
    auto k1 = option(_name="soil.1.k").as<double>();
    auto c1 = option(_name="soil.1.c").as<double>();
    auto rho1 = option(_name="soil.1.rho").as<double>();
    auto TR = option(_name="TR").as<double>();
    auto TA = option(_name="TA").as<double>();
    auto frequency = option(_name="omega").as <double>();
    auto gamma = option(_name="gamma").as <double>();

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "meshSize = " << option(_name="gmsh.hsize").as<double>() << "\n"
                  << "D = "<<D<<"\n"
                  << "W = " << W << "\n"
                  << "TR = " << TR << "\n"
                  << "TA = " << TA << "\n"
                  << "rho0 = " << rho0 << "\n"
                  << "k0 = " << k0 << "\n"
                  << "c0 = " << c0 << "\n"
                  << "rho1 = " << rho1 << "\n"
                  << "k1 = " << k1 << "\n"
                  << "c1 = " << c1 << "\n"
                  << "frequency = " << frequency << "\n";

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<Dim>> );
    auto Xh = Pch<Order>( mesh );
    auto ts = bdf( _space=Xh, _name="Temperature" );
    auto T = Xh->element();
    auto v = Xh->element();


    auto a = form2( Xh, Xh );
    a = integrate( _range= markedelements(mesh,"soil.0"), _expr= k0*gradt(T)*trans(grad(v)) );
    a += integrate( _range= markedelements(mesh,"soil.1"), _expr= k1*gradt(T)*trans(grad(v)) );
    a += integrate( _range= markedfaces(mesh, "ground"),
                    _expr=k0*(-gradt(T)*N()*id(v)-grad(T)*N()*idt(v)+gamma*idt(T)*id(v)/hFace() ) );
    a += integrate( _range=markedelements(mesh, "soil.0"), _expr=rho0*c0*idt(T)*id(v)*ts->polyDerivCoefficient( 0 ) )
        + integrate( _range=markedelements(mesh, "soil.1"), _expr=rho1*c1*idt(T)*id(v)*ts->polyDerivCoefficient( 0 ) );

    T = project( _space=Xh, _expr=cst(TR) );
    auto T0=expr(option(_name="T0").as<std::string>());
    ts->initialize(T);

    auto e = exporter( _mesh=mesh );
    for ( ts->start(); ts->isFinished()==false; ts->next(T) )
    {
        if ( Environment::worldComm().isMasterRank() )
        {
            std::cout << "============================================================\n";
            std::cout << "t =" << ts->time() << "s " << " = " << ts->time()/60/60 << "h\n";
        }
        T0.setParameterValues( {{"t",ts->time()}} );
        // update right hand side with time dependent terms
        auto tsrhs = ts->polyDeriv();
        auto l = form1( _test=Xh);
        l = integrate( _range= markedfaces(mesh, "ground"), _expr=T0*k0*(-grad(v)*N()+gamma*id(v)/hFace() ) )+
            integrate( _range=markedelements(mesh, "soil.0"), _expr=rho0*c0*idv(tsrhs)*id(v)) +
            integrate( _range=markedelements(mesh, "soil.1"), _expr=rho1*c1*idv(tsrhs)*id(v) );

        a.solve( _rhs=l, _solution=T );

        e->step(ts->time())->setMesh( mesh );
        e->step(ts->time())->add( "Temperature", T );
        e->save();

     }

    if ( Environment::worldComm().isMasterRank() )
        std::cout << "Resolution ended, export done \n";

}
