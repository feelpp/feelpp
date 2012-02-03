/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-03-04

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file convection_other.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-03-04
 */
#include <boost/lexical_cast.hpp>

#include "convection.hpp"

// Gmsh geometry/mesh generator
#include <feel/feelfilters/gmsh.hpp>

// gmsh importer
#include <feel/feelfilters/gmsh.hpp>





// ****** CONSTRUCTEURS ****** //
template <int Order_s, int Order_p, int Order_t>
Convection<Order_s,Order_p,Order_t>::Convection(int argc,
                                                char** argv,
                                                AboutData const& ad,
                                                po::options_description const& od)
	:
    super(argc,argv,ad,od),
    M_backend( backend_type::build( this->vm() ) ),
	exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{
    Parameter h(_name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.01:0.2:0.3" );
    this->
        addParameter( Parameter(_name="order-u",_type=DISC_ATTR,_latex="N_u",_values=boost::lexical_cast<std::string>( Order_s ).c_str() ) )
        .addParameter( Parameter(_name="order-p",_type=DISC_ATTR,_latex="N_p",_values=boost::lexical_cast<std::string>( Order_p ).c_str() ) )
        .addParameter( Parameter(_name="order-t",_type=DISC_ATTR,_latex="N_T",_values=boost::lexical_cast<std::string>( Order_t ).c_str() ) )
        .addParameter( Parameter(_name="gr",_type=CONT_ATTR,_latex="\\mbox{Gr}",_values="1:100:10000000") )
        .addParameter( Parameter(_name="pr",_type=CONT_ATTR,_latex="\\mbox{Pr}",_values="0.01:10:10000000") )
        .addParameter( h );

    Log() << "parameter added\n";
    std::vector<Parameter> depend;
    std::vector<std::string> funcs;
    std::vector<Parameter> depend2;
    depend2.push_back(h);
    std::vector<std::string> funcs2;
    funcs2.push_back("h**2");
    this->
          addOutput( Output(_name="AverageT",_latex="T",_dependencies=depend,_funcs=funcs) )
         .addOutput( Output(_name="FlowRate",_latex="D",_dependencies=depend,_funcs=funcs) )
         .addOutput( Output(_name="norm_L2",_latex="\\left\\| . \\right\\|_{L^2}",_dependencies=depend2,_funcs=funcs2) );

    Log() << "output added\n";
}

template <int Order_s, int Order_p, int Order_t>
typename Convection<Order_s, Order_p, Order_t>::mesh_ptrtype
Convection<Order_s, Order_p, Order_t>::createMesh()
{
    double h = this->vm()["hsize"].template as<double>();
    double l = this->vm()["length"].template as<double>();

    timers["mesh"].first.restart();
    mesh_ptrtype mesh( new mesh_type );
    //mesh->setRenumber( false );

    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = " << 2 << ";\n"
         << "a=" << 0 << ";\n"
         << "b=" << l << ";\n"
         << "c=" << 0 << ";\n"
         << "d=" << 1 << ";\n"
         << "h=" << h << ";\n"
         << "Point(1) = {a,c,0.0,h};\n"
         << "Point(2) = {b,c,0.0,h};\n"
         << "Point(3) = {b,d,0.0,h};\n"
         << "Point(4) = {a,d,0.0,h};\n"
         << "Point(5) = {b/2,d,0.0,h};\n"
         << "Point(6) = {b/2,d/2,0.0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,5};\n"
         << "Line(4) = {5,4};\n"
         << "Line(5) = {4,1};\n"
         << "Line(6) = {5,6};\n"
         << "Line Loop(7) = {1,2,3,4,5,6};\n"
        //<< "Line Loop(7) = {1,2,3,4,5};\n"
         << "Plane Surface(8) = {7};\n"
         << "Physical Line(\"Tinsulated\") = {1,3,4};\n"
         << "Physical Line(\"Tfixed\") = {5};\n"
         << "Physical Line(\"Tflux\") = {2};\n"
         << "Physical Line(\"Fflux\") = {6};\n"
         << "Physical Surface(\"domain\") = {8};\n";
    Gmsh gmsh;
    std::string fname = gmsh.generate( "domain", ostr.str()  );

    ImporterGmsh<mesh_type> import( fname );
    import.setVersion( "2.0" );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    Log() << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;

}



template <int Order_s, int Order_p, int Order_t>
void
Convection<Order_s,Order_p,Order_t> ::solve( sparse_matrix_ptrtype & J ,
                                             element_type & u,
                                             vector_ptrtype& F )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;

    this->updateResidual( U, F );
    this->updateJacobian( U, J );

    M_backend->nlSolve( J, U, F, 1e-10, 10 );
    u = *U;


};

template <int Order_s, int Order_p, int Order_t>
void Convection<Order_s,Order_p,Order_t> ::exportResults( boost::format fmt, element_type& U, double t)
{
    exporter->addPath( fmt );
    exporter->step(t)->setMesh( U.functionSpace()->mesh() );
    exporter->step(t)->add( "u", U.template element<0>() );
    exporter->step(t)->add( "p", U.template element<1>() );
    exporter->step(t)->add( "T", U.template element<2>() );
    exporter->save();
};

// instantiation
template class Convection<2,1,2>;
