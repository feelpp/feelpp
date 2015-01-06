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
#include <feel/feelpde/curveparametrizations.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;


int main( int argc, char** argv )
{
  const int dim = 2;

  Feel::Environment env( argc, argv );

  auto mesh = unitHypercube<dim>();
  auto Xh1 = Pch<1>(mesh);
  auto Xh0 = Pdh<0>(mesh);

  // distance to curve
  auto disttocurve = distToCurve( Xh0, Xh1 );

  // fast marching
  auto fm = fms( Xh1 );


  // ------------- ellipse --------------
  const double a_ell=0.2;
  const double b_ell=0.3;
  auto x_ell = [&](double t ) -> double { return 0.5 + a_ell * cos(t); };
  auto y_ell = [&](double t ) -> double { return 0.5 + b_ell * sin(t); };

  auto ellipse = disttocurve->fromParametrizedCurve(x_ell, y_ell,
                                                    0,    /*tStart*/
                                                    6.29, /*tEnd*/
                                                    option("gmsh.hsize").as<double>() / 2. /*dt*/ );
  *ellipse = fm->march( ellipse, true );
  std::cout<<"ellipse done"<<std::endl;
  // ------------------------------------




  // ------------- pre defined ellipse --------------
  auto ellipseParam = CurveParametrization::ellipse( 0.25, 0.28, 0.5, 0.5 );
  auto predefEllipse = disttocurve->fromParametrizedCurve(ellipseParam,
                                                          option("gmsh.hsize").as<double>() / 2. );
  *predefEllipse = fm->march( predefEllipse, true );
  std::cout<<"predef ellipse done"<<std::endl;
  // ------------------------------------------------




  // ------------- epitrochoid --------------
  const double a_epi=0.1; // 1 / nb_branch
  const double b_epi=0.8;
  auto x_epi = [&](double t) -> double { return ((1+a_epi) * cos(a_epi*t) - a_epi*b_epi * cos( (1+a_epi) * t )) / 4 + 0.5; };
  auto y_epi = [&](double t) -> double { return ((1+a_epi) * sin(a_epi*t) - a_epi*b_epi * sin( (1+a_epi) * t )) / 4 + 0.5; };

  auto epitro = disttocurve->fromParametrizedCurve( x_epi, y_epi,
                                                    0, 100, option("gmsh.hsize").as<double>()/2. );
  *epitro = fm->march( epitro, true );
  std::cout<<"epitrochoid done"<<std::endl;
  // ----------------------------------------




  // ------------- Sickle cell --------------
  const double R1=0.5;
  const double R2=0.21;
  const double H=0.2; // height of the cell
  const double l=option("gmsh.hsize").as<double>(); // length of the cut at the peak (top and bottom)

  // position of the sickle cell (the center of the first circle)
  const double dx = 0.;
  const double dy = 0.5;

  const double dt1=l/R1;
  const double dt2=l/R2;

  const double t1=asin(H/R1);
  const double t2=asin(H/R2);

  const double X=R1*cos(t1)-R2*cos(t2);

  auto x_sc = [&](double t) -> double {
      if ((-t1+dt1 <= t) && (t <= t1-dt1))
          return R1*cos(t)+dx;

      else if ((t1 - dt1 < t) && (t <= t1 + dt2))
          return R1*cos(t1-dt1) + dx +  (t-(t1-dt1))/(t1+dt2-(t1-dt1)) *(R2*cos(t2-dt2)+X+dx-(R1*cos(t1-dt1)+dx ));

      else if ((t1+dt2 < t) && (t <= t1+2*t2-dt2))
          return R2*cos(-(t-t1-t2))+X+dx;

      else if ((t1+2*t2-dt2 < t) && (t <= t1+2*t2+dt1))
          return R2*cos(-(t2-dt2))+X+dx + (t-(t1+2*t2-dt2))/(t1+2*t2+dt1-(t1+2*t2-dt2))*(R1*cos(-t1+dt1)+dx-(R2*cos(-(t2-dt2))+X+dx ));

      else
          return 0.5;
  };

  auto y_sc = [&](double t) -> double {
      if ((-t1+dt1 <= t) && (t <= t1-dt1))
          return R1*sin(t)+dy;

      else if ((t1 - dt1 < t) && (t <= t1 + dt2))
          return R1*sin(t1-dt1)+dy +  (t-(t1-dt1))/(t1+dt2-(t1-dt1))*(R2*sin(t2-dt2)+dy-(R1*sin(t1-dt1)+dy ));

      else if ((t1+dt2 < t) && (t <= t1+2*t2-dt2))
          return R2*sin(-(t-t1-t2))+dy;

      else if ((t1+2*t2-dt2 < t) && (t <= t1+2*t2+dt1))
          return R2*sin(-(t2-dt2)) + dy + (t-(t1+2*t2-dt2))/(t1+2*t2+dt1-(t1+2*t2-dt2))*(R1*sin(-t1+dt1)+dy-(R2*sin(-(t2-dt2))+dy));

      else
          return 0.5;
  };


  auto sickle_cell = disttocurve->fromParametrizedCurve(/* x(t), y(t) */
                                                        x_sc, y_sc,
                                                        /* tStart, tEnd, dt */
                                                       -t1+dt1, t1+2*t2+dt1, option("gmsh.hsize").as<double>()/8.,
                                                        /* broaden points for detection, broadening thickness, export points*/
                                                        true, option("gmsh.hsize").as<double>() / 10., true );
  *sickle_cell = fm->march( sickle_cell, true );
  std::cout<<"sc done"<<std::endl;
  // ----------------------------------------




  // ------------- pre-defined Sickle cell  --------------
  auto scpredef = CurveParametrization::crescent(R1, R2, H, l, 0.5, 0.5);

  auto siclkeCellPredefied = disttocurve->fromParametrizedCurve( get<0>(scpredef),
                                                                 get<1>(scpredef),
                                                                 get<2>(scpredef),
                                                                 get<3>(scpredef),
                                                                 option("gmsh.hsize").as<double>()/8.,
                                                                 true, option("gmsh.hsize").as<double>() / 10., true, "sc_pred" );

  *siclkeCellPredefied = fm->march( siclkeCellPredefied, true );
  std::cout<<"predef sc done"<<std::endl;
  // ----------------------------------------




  auto exp = exporter(_mesh=mesh, _name="dist2paramcurve");
  exp->step(0)->add("ellipse", *ellipse);
  exp->step(0)->add("predefEllipse", *predefEllipse);
  exp->step(0)->add("epitro", *epitro);
  exp->step(0)->add("sickle_cell",*sickle_cell);
  exp->step(0)->add("siclkeCellPredefied",*siclkeCellPredefied);
  exp->save();

}
