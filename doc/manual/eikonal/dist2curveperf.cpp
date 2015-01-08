/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file dist2curveperf.cpp
   \author Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
   \date 2014-01-21
 */

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feelpde/reinit_fms.hpp>
#include <fstream>
#include <sstream>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

int main( int argc, char** argv )
{
    const int dim = 3;
    const int order = 1;

    boost::timer chrono;

    Feel::Environment env( argc, argv );

    auto mesh = unitHypercube<dim>();
    auto Xh = Pch<order>(mesh);

    // some ellipse parameters:
    const double x0=0.5;
    const double y0=0.5;
    const double z0 = dim==2 ? 0 : 0.5;
    const double r0 = 0.2;
    auto X0 = Px() - x0;
    auto Y0 = Py() - y0;
    auto Z0 = Pz() - z0;

    auto initShape = vf::project(Xh, elements(mesh),
                                 sqrt( X0*X0 + Y0*Y0 + Z0*Z0 ) - r0 );

    chrono.restart();
    auto fm = fms( Xh );
    const double timeInitFms = chrono.elapsed();

    chrono.restart();
    auto phi = fm->march( initShape );
    const double timeMarch = chrono.elapsed();

    const double eL2 = std::sqrt( integrate(elements(mesh),
                                            (idv( phi ) - idv(initShape)) * (idv( phi ) - idv(initShape)) ).evaluate()(0,0) );

    const double eH1 = std::sqrt( integrate(elements(mesh),
                                            (idv( phi ) - idv(initShape)) * (idv( phi ) - idv(initShape))
                                            + (gradv( phi ) - gradv( initShape )) * trans( gradv( phi ) - gradv( initShape ) ) ).evaluate()(0,0) );

    const double eSH1 = std::sqrt( integrate(elements(mesh),
                                             (gradv( phi ) - gradv( initShape )) * trans( gradv( phi ) - gradv( initShape ) ) ).evaluate()(0,0) );

    const double area = integrate(elements(mesh), idv(phi) < 0 ).evaluate()(0,0);


    std::ostringstream resFileName;
    resFileName<< "res_order_" <<order
               << "_h_"
               << option("gmsh.hsize").as<double>();

    std::fstream resFile( resFileName.str(), std::fstream::out );
    resFile << std::left<<std::setw(15) << "#h"
            << std::left<<std::setw(15) << "np"
            << std::left<<std::setw(15) << "ndof"
            << std::left<<std::setw(15) << "el2"
            << std::left<<std::setw(15) << "esh1"
            << std::left<<std::setw(15) << "eh1"
            << std::left<<std::setw(15) << "area"
            << std::left<<std::setw(15) << "tinitfms"
            << std::left<<std::setw(15) << "tmarch"
            << std::endl

            << std::left<<std::setw(15) << option("gmsh.hsize").as<double>()
            << std::left<<std::setw(15) << Environment::worldComm().size()
            << std::left<<std::setw(15) << Xh->dof()->nDof()
            << std::left<<std::setw(15) << eL2
            << std::left<<std::setw(15) << eSH1
            << std::left<<std::setw(15) << eH1
            << std::left<<std::setw(15) << area
            << std::left<<std::setw(15) << timeInitFms
            << std::left<<std::setw(15) << timeMarch
            << std::endl;


    if ( option("exporter.export").as<bool>() )
      {
        auto Xhv = Pchv<order>(mesh);
        auto gradphi = vf::project(Xhv, elements(mesh), trans(gradv(phi)) );
        auto gradinitshape = vf::project(Xhv, elements(mesh), trans(gradv(initShape)) );

        auto diffphi2 = vf::project(Xh, elements(mesh), (idv( phi ) - idv(initShape)) * (idv( phi ) - idv(initShape)) );

        auto exp = exporter(_mesh=mesh, _name="dist2curveperf");
        exp->step(0)->add("initShape", initShape);
        exp->step(0)->add("phi", phi);
        exp->step(0)->add("diffphi2", diffphi2);
        exp->step(0)->add("gradphi", gradphi);
        exp->step(0)->add("gradinitshape", gradinitshape);
        exp->save();
      }

    return 0;
}
