/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-02-06

  Copyright (C) 2010-2011 Université Joseph Fourier (Grenoble I)

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
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

#include <opusdata.hpp>
#include <opusmodelbase.hpp>
#include <opusmodelfactory.hpp>

#include <eadsmfemapp.hpp>

namespace Feel
{
void
EadsMFemApp::init()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    if ( this->vm()["steady"].as<bool>() )
    {
        std::cout << "[EadsMFemApp::init] steady case" << std::endl;
        this->changeRepository( boost::format( "%1%/P%2%P%3%P%4%/D_%5%/e_%6%/h_%7%/stab_%8%/steady" )
                                % this->about().appName()
                                % this->vm()["order-u"].as<int>() % this->vm()["order-p"].as<int>() % this->vm()["order-temp"].as<int>()
                                % this->vm()["fluid.flow-rate"].as<double>()
                                % this->vm()["air.e"].as<double>()
                                % this->vm()["hsize"].as<double>()
                                % this->vm()["stab"].as<bool>()
                              );
    }

    else
    {
        std::cout << "[EadsMFemApp::init] unsteady case" << std::endl;
        this->changeRepository( boost::format( "%1%/P%2%P%3%P%4%/D_%5%/h_%6%/stab_%7%/to_%8%_dt_%9%" )
                                % this->about().appName()
                                % this->vm()["order-u"].as<int>() % this->vm()["order-p"].as<int>()  % this->vm()["order-temp"].as<int>()
                                % this->vm()["fluid.flow-rate"].as<double>()
                                % this->vm()["hsize"].as<double>()
                                % this->vm()["stab"].as<bool>()
                                % this->vm()["bdf.order"].as<int>()
                                % this->vm()["bdf.time-step"].as<double>()
                              );
    }

    std::cout << "[EadsMFemApp::init] changing repository done" << std::endl;
    std::cout << "[EadsMFemApp::init] building model" << std::endl;
    opus = OpusModelFactory::New( this->vm() );
    std::cout << "[EadsMFemApp] building model done\n";
}
void
EadsMFemApp::run()
{
    std::cout << "[EadsMFemApp::run] running model\n";
    std::vector<double> X( 7 ),Y( 2 );
    X[0]=this->vm()["ic1.k"].as<double>(); // kic
    X[1]=this->vm()["fluid.flow-rate"].as<double>(); // fluid flow rate
    X[2]=this->vm()["ic1.Q"].as<double>(); // Q
    X[3]=100;
    X[4]=this->vm()["air.e"].as<double>(); // ea
    X[5]=this->vm()["hsize"].as<double>(); // h
    X[6]=this->vm()["order-temp"].as<int>(); // polynomial order
    opus->run( X.data(), X.size(), Y.data(), Y.size() );
    //M_opus->run();
    std::cout << "[EadsMFemApp::run] running model done\n";
}
void
EadsMFemApp::run( const double * X, unsigned long N,
                  double * Y, unsigned long P )
{

    opus->run( X, N, Y, P );
}

}
