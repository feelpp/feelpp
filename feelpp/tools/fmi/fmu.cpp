/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Jean-Baptiste Wahl
            Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 01 Dec 2019

 Copyright (C) 2018-2019 Feel++ Consortium

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
#include <iostream>
#include <string>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfmi/fmi4cpp.hpp>

using namespace fmi4cpp;

const double stop = 0.01;
const double stepSize = 1E-3;

int main( int argc, char* argv[] )
{
    using namespace Feel;
    Environment env( _argc = argc, _argv = argv );

    const std::string fmuPath = "../resources/fmus/2.0/cs/20sim/4.6.4.8004/"
                                "ControlledTemperature/ControlledTemperature.fmu";
    fmi2::fmu fmu( fmuPath );
    auto cs_fmu = fmu.as_cs_fmu();
    auto md = cs_fmu->get_model_description();

    auto var = md->model_variables->getByValueReference( 47 ).as_real();
    std::cout << "Name=" << var.name() << ", start=" << var.start().value_or( 0 ) << std::endl;

    auto slave1 = cs_fmu->new_instance();
    auto slave2 = cs_fmu->new_instance();

    std::cout << "model_identifier=" << slave1->get_model_description()->model_identifier << std::endl;

    slave1->setup_experiment();
    slave1->enter_initialization_mode();
    slave1->exit_initialization_mode();

    slave2->setup_experiment();
    slave2->enter_initialization_mode();
    slave2->exit_initialization_mode();

    std::vector<fmi2Real> ref( 2 );
    std::vector<fmi2ValueReference> vr = { md->get_variable_by_name( "Temperature_Reference" ).value_reference,
                                           md->get_variable_by_name( "Temperature_Room" ).value_reference };

    double t;
    while ( ( t = slave1->get_simulation_time() ) <= stop )
    {

        if ( !slave1->step( stepSize ) )
        {
            break;
        }
        if ( !slave1->read_real( vr, ref ) )
        {
            break;
        }
        std::cout << "t=" << t << ", Temperature_Reference=" << ref[0] << ", Temperature_Room=" << ref[1] << std::endl;
    }

    std::cout << "FMU '" << fmu.model_name() << "' terminated with success: " << ( slave1->terminate() == 1 ? "true" : "false" ) << std::endl;
    std::cout << "FMU '" << fmu.model_name() << "' terminated with success: " << ( slave2->terminate() == 1 ? "true" : "false" ) << std::endl;

    return 0;
}