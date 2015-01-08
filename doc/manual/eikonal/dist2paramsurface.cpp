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
   \file dist2paramcurve.cpp
   \author Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
   \date 2014-02-01
 */

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feelpde/reinit_fms.hpp>
#include <feel/feelpde/disttocurve.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;


int main( int argc, char** argv )
{
    const int dim = 3;
    Feel::Environment env( argc, argv );

    auto mesh = unitHypercube<dim>();
    auto Xh1 = Pch<1>(mesh);
    auto Xh0 = Pdh<0>(mesh);

    // distance to curve
    auto disttocurve = distToCurve( Xh0, Xh1 );

    // fast marching
    auto fm = fms( Xh1 );

    // ------------- ellipse --------------
#if 0
    const double a_ell=0.2;
    const double b_ell=0.3;
    const double c_ell=0.4;

    auto x_ell = [&](double t1, double t2 ) -> double { return 0.5 + a_ell * cos(t1)*sin(2*t2+3.14159); };
    auto y_ell = [&](double t1, double t2 ) -> double { return 0.5 + b_ell * sin(t1)*sin(2*t2+3.14159); };
    auto z_ell = [&](double t1, double t2 ) -> double { return 0.5 + c_ell * cos(t2); };

    auto ellipse = disttocurve->fromParametrizedCurve(x_ell, y_ell, z_ell,
                                                      0,
                                                      6.29,
                                                      option("gmsh.hsize").as<double>() / 10.,
                                                      -3.14159/2,
                                                      3.14159/2,
                                                      option("gmsh.hsize").as<double>() / 16.,
                                                      true,option("gmsh.hsize").as<double>() / 2,true,"ellipses");

    *ellipse = fm->march( ellipse, true );

#endif

    //sickle_cell

    const double H=0.8;
    const double R=0.75;
    const double r=0.1;


    const double t2max=asin(H/(2*R));
    const double l=option("gmsh.hsize").as<double>()*2.;
    const double dt2=l/R;


    auto x_scell = [&](double t1, double t2 ) -> double {
        if ((-t2max<=t2)&&(t2<-t2max+dt2))
            return  (r * sin((t2+t2max)*3.14159/(2*t2max)) * cos(t1)) + 2*R * (cos(-t2max+dt2) - cos(t2max)) +0.5 ;
        else if ((-t2max+dt2<=t2)&&(t2<t2max-dt2))
            return  (r * sin((t2+t2max)*3.14159/(2*t2max)) * cos(t1)) + 2*R * (cos(t2) - cos(t2max)) +0.5 ;
        else if((t2max-dt2<t2)&&(t2<=t2max))
            return (r * sin((t2+t2max)*3.14159/(2*t2max)) * cos(t1)) + 2*R * (cos(t2max-dt2) - cos(t2max)) +0.5 ;
    };

    auto y_scell = [&](double t1, double t2 ) -> double { return  (r * sin((t2+t2max)*3.14159/(2*t2max)) * sin(t1)) + 0.5; };

    auto z_scell = [&](double t1, double t2 ) -> double {
        if ((-t2max<=t2)&&(t2<-t2max+dt2))
            return  (R * sin(-t2max+dt2)) + 0.5;
        else if ((-t2max+dt2<=t2)&&(t2<t2max-dt2))
            return  (R * sin(t2)) + 0.5;
        else if((t2max-dt2<t2)&&(t2<=t2max))
            return  (R * sin(t2max-dt2)) + 0.5;
    };


    auto sickle_cell = disttocurve->fromParametrizedCurve(x_scell, y_scell, z_scell,
                                                          0,
                                                          6.29,
                                                          option("gmsh.hsize").as<double>() / 8.,
                                                          -t2max,
                                                          t2max,
                                                          option("gmsh.hsize").as<double>() / 16.,
                                                          true,option("gmsh.hsize").as<double>() / 3,true,"cell");

    *sickle_cell = fm->march( sickle_cell, true );
  // ------------------------------------

    auto mark2=vf::project(Xh0,marked2elements(mesh,1),cst(1));

    auto exp = exporter(_mesh=mesh, _name="dist2paramcurve");
    //    exp->step(0)->add("ellipse", *ellipse);
    exp->step(0)->add("mark2",mark2);
    exp->step(0)->add("sickle_cell",*sickle_cell);

    exp->save();


}
