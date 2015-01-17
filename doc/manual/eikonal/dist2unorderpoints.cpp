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
    /*
      This test works only if the macro DISTANCE_FROM_UNORDERED_POINTS
      in "feel/feelpde/disttocurve.hpp" is defined
    */

  const int dim = 2;

  Feel::Environment env( argc, argv );

  auto mesh = unitHypercube<dim>();
  auto Xh1 = Pch<1>(mesh);
  auto Xh0 = Pdh<0>(mesh);

  // distance to curve
  auto disttocurve = distToCurve( Xh0, Xh1 );

  // fast marching
  auto fm = fms( Xh1 );

  // ------------- pre defined ellipse --------------
  auto ellipseParam = CurveParametrization::ellipse( 0.25, 0.28, 0.5, 0.5 );
  auto predefEllipse = disttocurve->fromParametrizedCurveDisordered(ellipseParam,
                                                                    option("gmsh.hsize").as<double>() / 16.,
                                                                    true,
                                                                    option("gmsh.hsize").as<double>() / 2.,
                                                                    true
                                                                    );
  std::cout<<"predef ellipse on boundary done"<<std::endl;

  *predefEllipse = fm->march( predefEllipse, true );

  std::cout<<"predef ellipse done"<<std::endl;
  // ------------------------------------------------

  auto mark2 = vf::project(Xh0, marked2elements(mesh, 1), cst(1) );

  auto exp = exporter(_mesh=mesh, _name="dist2paramcurve");
  exp->step(0)->add("predefEllipse", *predefEllipse);
  // exp->step(0)->add("Ellipse_resigned", *Ellipse_resigned);
  exp->step(0)->add("mark2", mark2);
  exp->save();

}
