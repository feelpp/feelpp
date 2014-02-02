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

using namespace Feel;
using namespace Feel::vf;

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
  // ------------------------------------




  // ------------- epitrochoid --------------
  const double a_epi=0.1; // 1 / nb_branch
  const double b_epi=0.8;
  auto x_epi = [&](double t) -> double { return ((1+a_epi) * cos(a_epi*t) - a_epi*b_epi * cos( (1+a_epi) * t )) / 10 + 0.5; };
  auto y_epi = [&](double t) -> double { return ((1+a_epi) * sin(a_epi*t) - a_epi*b_epi * sin( (1+a_epi) * t )) / 10 + 0.5; };

  auto epitro = disttocurve->fromParametrizedCurve( x_epi, y_epi,
                                                    0, 100, option("gmsh.hsize").as<double>()/2. );
  *epitro = fm->march( epitro, true );
  // ----------------------------------------



  // // ------------- heart --------------------
  // auto x_heart = [&](double t) -> double { return cos(t) / 2.; };
  // auto y_heart = [&](double t) -> double { return (sin(t) + std::sqrt( std::abs( cos(t) ) )) / 2.; } ;

  // auto heart = disttocurve->fromParametrizedCurve( x_heart, y_heart,
  //                                                  0, 100, option("gmsh.hsize").as<double>()/5. );
  // *heart = fm->march( heart, true );
  // // ----------------------------------------



  auto exp = exporter(_mesh=mesh, _name="dist2paramcurve");
  exp->step(0)->add("ellipse", *ellipse);
  exp->step(0)->add("epitro", *epitro);
  //  exp->step(0)->add("heart", *heart);
  exp->save();



}
