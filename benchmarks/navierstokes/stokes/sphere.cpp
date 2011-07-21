/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2011-07-19

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
   \file sphere.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-07-19
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

template<typename MeshType>
void myexport( std::string const& name, boost::shared_ptr<MeshType> mesh )
{
    using namespace Feel;
    using namespace Feel::vf;

    typedef Exporter<MeshType> exporter_type;
    auto exporter = exporter_type::New( "gmsh", name );
    exporter->step(0)->setMesh( mesh );
    exporter->save();
}

int main(int argc, char** argv)
{
    double hsize = 1;
    double R = 1;
    double U = 1;
    double mu = 1;
    bool straighten = true;
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("geomap", po::value<int>()->default_value( 0 ), "straighten mesh (0=ho,1=p1,2=opt)")
        ("straighten", po::value<bool>(&straighten)->default_value( true ), "straighten mesh")
        ("hsize", po::value<double>(&hsize)->default_value( 1 ), "h size")
        ("R", po::value<double>(&R)->default_value( 1 ), "sphere radius")
        ("U", po::value<double>(&U)->default_value( 1 ), "characteristic fluid velocity")
        ("mu", po::value<double>(&U)->default_value( 1 ), "fluid viscosity")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    const int Order=2;
    using namespace Feel;
    using namespace Feel::vf;
    Feel::Environment env(argc, argv );
    typedef Mesh<Simplex<3,Order> > mesh_type;

    GeoTool::Cube R1(hsize,"R1",GeoTool::Node(-3,-3,-3),GeoTool::Node(3,3,3));
    R1.setMarker(_type="surface",_name="Cube",_markerAll=true);
    R1.setMarker(_type="volume",_name="Omega",_markerAll=true);

    GeoTool::Sphere C1(hsize,"C1",GeoTool::Node(0,0,0),GeoTool::Node(R,0,0));
    C1.setMarker(_type="surface",_name="Sphere",_markerAll=true);

    auto R1mC1mesh = (R1-C1).createMesh<mesh_type>( "R1-C1" );

    myexport<mesh_type>("R1-C1",R1mC1mesh);

    auto mesh = R1mC1mesh;
    if ( straighten )
        mesh = straightenMesh( _mesh=R1mC1mesh );


    typedef FunctionSpace<mesh_type,bases<Lagrange<Order,Vectorial> > > Vh_t;
    typedef FunctionSpace<mesh_type,bases<Lagrange<Order,Scalar> > > Ph_t;
    auto Vh = Vh_t::New( mesh );
    auto Ph = Ph_t::New( mesh );

    auto r = sqrt(Px()*Px()+Py()*Py()+Pz()*Pz());
    auto rxy = sqrt(Px()*Px()+Py()*Py());
    auto sintheta = sin(acos(Pz()/r));
    auto costheta = Pz()/r;

    auto er=vec(Px()/r,Py()/r,Pz()/r);
    auto etheta=vec(Px()*Pz()/(r*rxy),Py()*Pz()/(r*rxy),-rxy/r);
    auto ephi=vec(-Py()/rxy,Px()/rxy,cst(0.));

    auto phi = vf::project( _space=Ph, _expr=-R*2*U*0.5*pow((r/R)*sintheta,2)*(1-1.5*R/r+0.5*pow(R/r,3) ));
    auto w = vf::project( _space=Vh, _expr=1.5*(U/R)*pow(R/r,2)*sintheta*ephi);
    auto ur = -2*U*(1-1.5*(R/r)+.5*pow(R/r,3))*costheta/R;
    auto utheta = (2*U*r*sintheta*sintheta*(1-1.5*(R/r)+.5*pow(R/r,3)) +
                   U*(r*sintheta)*(r*sintheta)*(1.5*(R/(r*r))-1.5*pow(R,3)/(r*r*r*r)))/(R*r*sintheta);
    auto u = vf::project( _space=Vh, _expr=ur*er+utheta*etheta);
    auto p = vf::project( _space=Ph, _expr=1.5*(mu*U/R)*(R/r)*(R/r)*costheta);
    auto st = vf::project( _space=Ph, _expr=sintheta);
    auto uur = vf::project( _space=Ph, _expr=ur);
    auto uut = vf::project( _space=Ph, _expr=utheta);

    auto sigmaN=-idv(p)*N()+.5*mu*(gradv(u)+trans(gradv(u)))*N();
    auto volume = integrate( _range=elements(mesh), _expr=cst(1.0),_quad=_Q<6>(), _geomap=(GeomapStrategyType)vm["geomap"].as<int>()  ).evaluate();
    auto surface = integrate( _range=markedfaces(mesh,"Sphere"), _expr=cst(1.0),_quad=_Q<6>() ).evaluate();
    auto force = integrate( _range=markedfaces(mesh,"Sphere"), _expr=sigmaN,_quad=_Q<6>() ).evaluate();
    Eigen::Vector3d force_exact;
    force_exact << 6*M_PI*mu*U*R, 0, 0;

    std::cout << "volume = " << volume << "\n";
    double ve = 6*6*6-4*M_PI*R*R*R/3;
    std::cout << "volume exact = " << ve << "\n";
    std::cout << "error volume = " << math::abs(volume(0,0)-ve) << "\n";

    std::cout << "surface = " << surface << "\n";
    std::cout << "surface exact = " << 4*M_PI*R*R << "\n";
    std::cout << "error surface = " << math::abs(surface(0,0)-4*M_PI*R*R) << "\n";

    std::cout << "force = " << force << "\n";
    std::cout << "force_exacte = " << force_exact << "\n";
    std::cout << "error=" << (force-force_exact).norm() << "\n";

    typedef Exporter<mesh_type,Order> exporter_type;
    auto exporter = exporter_type::New( "gmsh", "sphere" );
    exporter->step(0)->setMesh( mesh );
    exporter->step(0)->add( "phi", phi );
    exporter->step(0)->add( "w", w );
    exporter->step(0)->add( "u", u );
    exporter->step(0)->add( "sin(theta)", st );
    exporter->step(0)->add( "ur", uur );
    exporter->step(0)->add( "utheta", uut );
    exporter->step(0)->add( "ux", u.comp(X) );
    exporter->step(0)->add( "uy", u.comp(Y) );
    exporter->step(0)->add( "uy", u.comp(Z) );
    exporter->step(0)->add( "p", p );
    exporter->save();

}


