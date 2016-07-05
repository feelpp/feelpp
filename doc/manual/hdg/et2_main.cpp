/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 26 Feb 2016

 Copyright (C) 2016 Feel++ Consortium

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
#include <et2.hpp>

using namespace Feel;
inline
po::options_description
makeOptions()
{
    po::options_description testhdivoptions( "Electro-Thermal Model options" );
    testhdivoptions.add_options()
        ( "sigma", po::value<std::string>()->default_value( "1" ), "electric conductivity" )
        ( "k", po::value<std::string>()->default_value( "1" ), "thermal conductivity" )
        ( "h", po::value<double>()->default_value( 1 ), "heat transfert coefficient" )
        ( "Tw", po::value<double>()->default_value( 1 ), "water temperature" )
        ( "V", po::value<double>()->default_value( 1 ), "electric potential" )
        ( "I_target", po::value<double>()->default_value( 1 ), "target current" )
        ( "tau_constant", po::value<double>()->default_value( 1.0 ), "stabilization constant for hybrid methods" )
        ( "tau_order", po::value<int>()->default_value( 0 ), "order of the stabilization function on the selected edges"  ) // -1, 0, 1 ==> h^-1, h^0, h^1
        ( "picard.itol", po::value<double>()->default_value( 1e-4 ), "tolerance" )
        ( "picard.itmax", po::value<int>()->default_value( 10 ), "iterations max" )
        ("Kp", po::value<double>()->default_value(100.), "PID proportional coefficient")
        ("Ki", po::value<double>()->default_value(0.), "PID integral coefficient")
        ("Kd", po::value<double>()->default_value(0.), "PID derivative coefficient")
        ;
    return testhdivoptions;
}

inline
AboutData
makeAbout()
{
    AboutData about( "electro-thermal" ,
                     "electro-thermal" ,
                     "0.1",
                     "Electro-Thermal Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Christophe Trophime", "developer", "christophe.trophime@lncmi.cnrs.fr", "" );
    about.addAuthor( "Daniele Prada", "developer", "", "" );
    return about;

}

int main(int argc, char *argv[])
{
    using namespace Feel;
    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=makeOptions() );
    FeelModels::ElectroThermal<2,1> ET;
    ET.run();
    return 0;
}
