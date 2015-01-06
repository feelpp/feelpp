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
   \file dist2curve.cpp
   \author Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
   \date 2014-01-21
 */

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelpde/reinit_fms.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

int main( int argc, char** argv )
{
    const int dim = 2;
    const int order = 3;
    const bool isperiodic = true;

    typedef Mesh< Simplex<dim> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef mpl::if_c< isperiodic,
                       Periodic<>,
                       NoPeriodicity >::type periodicity_type;

    typedef bases<Lagrange<order, Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type, Feel::Periodicity <Feel::NoPeriodicity> > space_type;

    Feel::Environment env( argc, argv );

    // mesh and high order space
    auto mesh = loadMesh( _mesh=new mesh_type, _filename=option("gmsh.filename").as<std::string>(),
                          _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER);

    // defines the translation
    node_type transl( dim );
    transl[0]=1.;

    // defines the periodicity
    periodicity_type theperiodicity(4,2, transl);

    auto Xh = space_type::New(_mesh=mesh );

    // operator lagrange P1 and associated space
    auto opLagP1 = lagrangeP1( Xh );

     auto XhP1 = Pch<1>( opLagP1->mesh() );

    // interpolation operator
    auto opIntHoToP1 = opInterpolation(_domainSpace = Xh, _imageSpace = XhP1, _type=InterpolationNonConforme(false) );
    auto opIntP1ToHo = opInterpolation(_domainSpace = XhP1, _imageSpace = Xh, _type=InterpolationNonConforme(false) );

    // some ellipse parameters:
    const double x0=0.6;
    const double y0=0.5;
    const double z0 = dim==2 ? 0 : 0.5;
    const double aAxis = 0.1;
    const double bAxis = 0.3;

    auto X0 = Px() - x0;
    auto Y0 = Py() - y0;
    auto Z0 = Pz() - z0;

    // ellipse function (not exactly a distance function)
    auto ellipseShape = vf::project(Xh, elements(mesh),
                                    vf::sqrt( (X0/aAxis) * (X0/aAxis)
                                              + (Y0/bAxis) * (Y0/bAxis)
                                              + (Z0/bAxis) * (Z0/bAxis) ) - 1 );

    // interface local projection method on high order space
    auto ilpEllipse = vf::project(Xh, elements(mesh),
                                  idv(ellipseShape)
                                  / sqrt( inner( gradv(ellipseShape), gradv(ellipseShape) ) ) );

    auto ilpEllipseP1 = XhP1->element();
    opIntHoToP1->apply( ilpEllipse, ilpEllipseP1);

    // XhP1 is not periodic, the periodicity information is given directly to the fm
    auto fm = fms( XhP1, theperiodicity );
    auto phiP1 = fm->march( ilpEllipseP1 );

    auto phiHo = Xh->element();
    opIntP1ToHo->apply( phiP1, phiHo );

    auto expHo = exporter(_mesh=mesh, _name="dist2curveHo");
    expHo->step(0)->add("phiHo", phiHo);
    expHo->step(0)->add("ilpEllipse", ilpEllipse);
    expHo->step(0)->add("ellipseShape",ellipseShape);
    expHo->save();

    // auto expP1 = exporter(_mesh=opLagP1->mesh(), _name="dist2curveP1");
    // expP1->step(0)->add("phiP1", phiP1);
    // expP1->step(0)->add("ilpEllipseP1", ilpEllipseP1);
    // expP1->save();

    return 0;
}
