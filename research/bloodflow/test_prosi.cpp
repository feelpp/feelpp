/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-02-16

  Copyright (C) 2006 EPFL

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
   \file main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-02-16
 */
#include <map>
#include <string>
#include <sstream>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/bdf.hpp>

// petsc may not be available
#if defined(FEELPP_HAS_PETSC_H)

#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/solverlinearpetsc.hpp>

#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/importergambit.hpp>
#include <feel/feelfilters/exporterensight.hpp>

#include <feel/feelvf/vf.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>

#define dot_gradt_grad( u, v ) (dxt(u)*dx(v) + dyt(u)*dy(v) + dyt(u)*dy(v))
//#define dot_idv_gradt( u, v ) (idv(u,0)*dxt(v) + idv(u,1)*dyt(v) + idv(u,2)*dzt(v))
//#define dot_idv_gradt( u, v ) (idv(u)*dxt(v) + idv(u)*dyt(v) + idv(u)*dzt(v))
#define dot_idv_gradt( u, v ) 0.0

namespace Feel
{
namespace ublas = boost::numeric::ublas;

const uint16_type DIM = 3;
using namespace vf;

Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "twodomainsmt" ,
                           "twodomainsmt" ,
                           "0.1",
                           "Two domain mass transport",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2006 EPFL" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Martin Prosi", "developer", "martin.prosi@mate.polimi.it", "" );
    return about;

}

const double u_filt = 1.78E-6;
const double U_0 = 15.7;
const double L_0 = 0.67;
const double Klag = 0.1164;
const double porosity = 0.15;

const double Pr = 2.0e-8; // permeability
const double kappa = 0.6176;
const double s = 2.134e-3;
const double D_L = 2.87e-7;
const double D_w = 8.0e-9;
const double k_w = -3.197e-4;

// Media
const uint16_type LUMEN = 1;
const uint16_type WALL = 2;

// BC ids
const uint16_type ADVENDITIA = 3;
const uint16_type ENDOTHELIUM = 4;
const uint16_type INTERNAL_LUMEN = 5;
const uint16_type CUTTING_PLANES = 6;
const uint16_type INFLOW_LUMEN = 7;
const uint16_type OUTFLOW_LUMEN = 8;

class TwoDomainsMTApp
    :
public Application
{
    typedef Application super;
public:
    typedef double value_type;



    static const uint16_type feOrder = 1;



    /*matrix*/
    typedef MatrixPetsc<value_type> sparse_matrix_type;
    typedef VectorPetsc<value_type> vector_type;
    //typedef MatrixGmm<value_type, gmm::row_major> sparse_matrix_type;
    //typedef VectorUblas<value_type> vector_type;


    typedef Mesh<GeoEntity<Simplex<DIM, 1> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /* spaces */
    // concentration
    typedef FunctionSpace<mesh_type, fem::Lagrange<DIM, 1, Scalar, Continuous> > concentration_space_type;
    typedef boost::shared_ptr<concentration_space_type> concentration_space_ptrtype;
    // pressure
    typedef concentration_space_type pressure_space_type;
    typedef boost::shared_ptr<pressure_space_type> pressure_space_ptrtype;

    // velocity
    typedef FunctionSpace<mesh_type, fem::Lagrange<DIM, 1, Vectorial, Continuous> > velocity_lumen_space_type;
    typedef boost::shared_ptr<velocity_lumen_space_type> velocity_lumen_space_ptrtype;

    //typedef FunctionSpace<mesh_type, fem::RaviartThomas<DIM, 2> > velocity_wall_space_type;
    typedef FunctionSpace<mesh_type, fem::Lagrange<DIM, 1, Vectorial, Continuous> > velocity_wall_space_type;
    typedef boost::shared_ptr<velocity_wall_space_type> velocity_wall_space_ptrtype;

    /* time discretizations */
    typedef Bdf<concentration_space_type> bdf_concentration_type;

    /*quadrature*/
    typedef IM_PK<DIM, 2, value_type> im1_type;
    typedef IM_PK<DIM, 4, value_type> im2_type;

    TwoDomainsMTApp( int argc, char** argv, AboutData const& ad )
        :
        super( argc, argv, ad ),
        M_mesh_c( new mesh_type ),
        M_mesh_u_lumen( new mesh_type ),
        M_mesh_u_wall( new mesh_type )
    {}
    TwoDomainsMTApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_mesh_c( new mesh_type ),
        M_mesh_u_lumen( new mesh_type ),
        M_mesh_u_wall( new mesh_type )
    {}

    void setFunctionSpaces();
    void solveSystem();
    void concentrationInit( concentration_space_type::element_type&,
                            concentration_space_type::element_type& );
    void concentrationStep( value_type dt,
                            bdf_concentration_type& time_discr_lumen,
                            bdf_concentration_type& time_discr_wall,
                            concentration_space_type::element_type& c_l,
                            velocity_lumen_space_type::element_type& u_l,
                            concentration_space_type::element_type& c_w,
                            velocity_wall_space_type::element_type& u_w );

    void run()
    {
        boost::timer t, t1;

        setFunctionSpaces();

        VLOG(1) << "time elapsed in TwoDomainsMTApp::setFunctionSpaces() : " << t1.elapsed() << "\n";
        t1.restart();

        solveSystem();
        VLOG(1) << "time elapsed in TwoDomainsMTApp::solveSystem() : " << t1.elapsed() << "\n";

        VLOG(1) << "time elapsed in TwoDomainsMTApp::run() : " << t.elapsed() << "\n";
    }

private:
    mesh_ptrtype M_mesh_c;
    mesh_ptrtype M_mesh_u_lumen;
    mesh_ptrtype M_mesh_u_wall;

    concentration_space_ptrtype Ch_lumen, Ch_wall;
    velocity_lumen_space_ptrtype Uh_lumen;
    velocity_wall_space_ptrtype Uh_wall;
    pressure_space_ptrtype Ph_lumen, Ph_wall;

    std::map<std::string, concentration_space_type::element_type> c_map;

};

void
TwoDomainsMTApp::setFunctionSpaces()
{
    po::variables_map vm = this->vm();
    std::ostringstream ostr, meshstr;
    meshstr << "wall_lumen_" <<  this->vm()["nel"].as<int>() << ".neu";

    ImporterGambit<mesh_type> import( meshstr.str() );
    M_mesh_c->accept( import );

    // lumen concentration
    Ch_lumen = concentration_space_type::New( M_mesh_c );
    VLOG(1) << "creating concentration lumen space done\n"
            << " - dimension = " << Ch_lumen->nDof() << "\n";

    // wall concentration
    Ch_wall = concentration_space_type::New( M_mesh_c );
    VLOG(1) << "creating concentration wall space done\n"
            << " - dimension = " << Ch_wall->nDof() << "\n";

    // lumen fluid
    //import.setFilename( u_lumen );
    //M_mesh_u_lumen->accept( import );
    //Uh_lumen = velocity_lumen_space_type::New( M_mesh_u_lumen );
    M_mesh_c->createSubmesh( *M_mesh_u_lumen,
                              M_mesh_c->beginElementWithMarker( LUMEN  ),
                              M_mesh_c->endElementWithMarker( LUMEN  ) );
    Uh_lumen = velocity_lumen_space_type::New( M_mesh_u_lumen );
    Ph_lumen = pressure_space_type::New( M_mesh_u_lumen );
    VLOG(1) << "creating velocity/pressurelumen space done\n"
            << " - dimension Uh_lumen = " << Uh_lumen->nDof() << "\n"
            << " - dimension Ph_lumen = " << Ph_lumen->nDof() << "\n";

    // lumen concentration
    //import.setFilename( u_wall );
    //M_mesh_u_wall->accept( import );
    //Uh_wall = velocity_wall_space_type::New( M_mesh_u_wall );
    M_mesh_c->createSubmesh( *M_mesh_u_wall,
                              M_mesh_c->beginElementWithMarker( WALL ),
                              M_mesh_c->endElementWithMarker( WALL  ) );
    Uh_wall = velocity_wall_space_type::New( M_mesh_u_wall );
    Ph_wall = pressure_space_type::New( M_mesh_u_wall );
    VLOG(1) << "creating velocity/pressure wall space done\n"
            << " - dimension Uh_wall = " << Uh_wall->nDof() << "\n"
            << " - dimension Ph_wall = " << Ph_wall->nDof() << "\n";

}

void
TwoDomainsMTApp::concentrationInit( concentration_space_type::element_type& c_l,
                                    concentration_space_type::element_type& c_w )
{
    c_w = vf::project( Ch_wall, markedelements( *M_mesh_c, WALL ), constant( 0.01 ) );
    c_l = vf::project( Ch_lumen, markedelements( *M_mesh_c, LUMEN ), constant( 1.00 ) );

}
void
TwoDomainsMTApp::concentrationStep( value_type dt,
                                    bdf_concentration_type& time_discr_lumen,
                                    bdf_concentration_type& time_discr_wall,
                                    concentration_space_type::element_type& c_l,
                                    velocity_lumen_space_type::element_type& u_l,
                                    concentration_space_type::element_type& c_w,
                                    velocity_wall_space_type::element_type& u_w )
{
    concentration_space_type::element_type c_w_prev( Ch_wall ), c_w_v( Ch_wall );
    concentration_space_type::element_type c_l_prev( Ch_lumen ), c_l_v( Ch_lumen );


    concentration_space_type::element_type alpha_l( Ch_lumen ), beta_l( Ch_lumen );
    concentration_space_type::element_type alpha_w( Ch_wall ), beta_w( Ch_wall );



    SolverLinearPetsc<double> solver;


    double lhs_der_time = time_discr_lumen.derivateCoefficient( 0 );

    const int N_jacobi = 4;

    // Jacobi iteration
    for ( int iter = 1; iter < N_jacobi; ++iter )
    {
        alpha_l = ( ublas::scalar_vector<value_type>( alpha_l.size(), u_filt ) - c_map["P_l"] ) / D_L;
        beta_l = ublas::element_prod( c_map["P_l"], c_w ) / ( D_L*porosity );

        c_l_prev = time_discr_lumen.derivate();

        //
        // Lumen
        //

        vector_type F_l;
        form( Ch_lumen, F_l ) = ( integrate( elements( *M_mesh_c ), im1_type(),
                                             idv( c_l_prev ) * id( c_l_v ) ) +
                                  integrate( markedfaces( *M_mesh_c, ENDOTHELIUM ), im1_type(),
                                             dt * D_L*idv( beta_l ) * id( c_l_v ) )
                                );
        F_l.close();


        sparse_matrix_type M_l;
        form( Ch_lumen, Ch_lumen, M_l ) = ( integrate( elements( *M_mesh_c ), im1_type(),
                                            lhs_der_time*idt( c_l )*id( c_l_v )+
                                            dt * D_L*( dot_gradt_grad( c_l, c_l_v ) ) +
                                            dt * dot_idv_gradt( u_l, c_l )*id( c_l_v ) )  +

                                            integrate( markedfaces( *M_mesh_c, ENDOTHELIUM ), im1_type(),
                                                    dt * D_L*idv( alpha_l ) * idt( c_l ) * id( c_l_v ) ) +

                                            on( markedfaces( *M_mesh_c, INFLOW_LUMEN ), c_l, F_l, constant( 1.0 ) ) +
                                            on( markedfaces( *M_mesh_c, INTERNAL_LUMEN ), c_l, F_l, constant( 1.0 ) ) );
        M_l.close();

        //timer.restart();
        vector_type C_l( Ch_lumen->dof()->nDof(), Ch_lumen->dof()->nLocalDof() );
        solver.solve( M_l, C_l, F_l, 1e-10, 1000 );

        //copy C_l in c_l
        for ( size_type i = 0; i < c_l.size(); ++i )
        {
            c_l( i ) = C_l( i );
        }

        //VLOG(1) << "[timer] solver time : " << timer.elapsed() << "\n";

        alpha_w = - c_map["P_w"] / ( porosity * D_w ) - ublas::scalar_vector<value_type>( c_w.size(), Klag*u_filt/D_w );
        beta_w = ublas::element_prod( c_map["P_w"], c_l ) / D_w;

        //
        // Wall
        //
        c_w_prev = time_discr_wall.derivate();
        vector_type F_w;
        form( Ch_wall, F_w ) = ( integrate( elements( *M_mesh_c ), im1_type(),
                                            idv( c_w_prev ) * id( c_w_v ) ) +
                                 integrate( markedfaces( *M_mesh_c, ENDOTHELIUM ), im1_type(),
                                            dt * D_w*idv( beta_w ) * id( c_w_v ) )
                               );
        F_w.close();


        sparse_matrix_type M_w;
        form( Ch_wall, Ch_wall, M_w ) = ( integrate( elements( *M_mesh_c ), im1_type(),
                                          ( lhs_der_time+dt*k_w )*idt( c_w )*id( c_w_v )+
                                          dt * D_w*dot_gradt_grad( c_w, c_w_v ) +
                                          dt* Klag * dot_idv_gradt( u_w, c_w )*id( c_w_v ) ) +

                                          integrate( markedfaces( *M_mesh_c, ENDOTHELIUM ), im1_type(),
                                                  dt * D_w*idv( alpha_w ) * idt( c_w ) * id( c_w_v ) ) +

                                          on( markedfaces( *M_mesh_c, ADVENDITIA ), c_w, F_w, constant( 0.01 ) ) );
        M_w.close();

        //timer.restart();
        vector_type C_w( Ch_wall->dof()->nDof(), Ch_wall->dof()->nLocalDof() );
        solver.solve( M_w, C_w, F_w, 1e-10, 1000 );

        //copy C_w in c_w
        for ( size_type i = 0; i < c_w.size(); ++i )
        {
            c_w( i ) = C_w( i );
        }
    }

    time_discr_lumen.shiftRight( c_l );
    time_discr_wall.shiftRight( c_w );

}
void
TwoDomainsMTApp::solveSystem()
{
    concentration_space_type::element_type c_w( Ch_wall );
    concentration_space_type::element_type c_l( Ch_lumen );

    velocity_lumen_space_type::element_type u_l( Uh_lumen );
    velocity_wall_space_type::element_type u_w( Uh_wall );

    concentrationInit( c_l, c_w );

    // time discretizations setup
    bdf_concentration_type time_discr_lumen( Ch_lumen, BDF_ORDER_ONE );
    bdf_concentration_type time_discr_wall( Ch_wall, BDF_ORDER_ONE );
    time_discr_lumen.shiftRight( c_l );
    time_discr_wall.shiftRight( c_w );


    value_type T = 1;
    value_type dt = 0.1;

    //
    // save the results
    //
    typedef Exporter<mesh_type> export_type;
    typedef Exporter<mesh_type>::timeset_type timeset_type;

    ExporterEnsight<mesh_type> wall_lumen( "wall_lumen" );
    export_type::timeset_ptrtype ts_wl( new timeset_type( "wall_lumen" ) );
    ts_wl->setTimeIncrement( dt );
    wall_lumen.addTimeSet( ts_wl );

    ExporterEnsight<mesh_type> wall( "wall" );
    export_type::timeset_ptrtype ts_w( new timeset_type( "wall" ) );
    ts_w->setTimeIncrement( dt );
    wall.addTimeSet( ts_w );

    ExporterEnsight<mesh_type> lumen( "lumen" );
    export_type::timeset_ptrtype ts_l( new timeset_type( "lumen" ) );
    ts_l->setTimeIncrement( dt );
    lumen.addTimeSet( ts_l );

    velocity_lumen_space_type::P1Lagrange Uh_lumen_p1( Uh_lumen );
    velocity_wall_space_type::P1Lagrange Uh_wall_p1( Uh_wall );

    for ( double t = dt; t < T; t += dt )
    {
        // velocity in the Lumen (could have come from an NS code )
        u_l = vf::project( Uh_lumen,
                           elements( *M_mesh_u_lumen ),
                           oneX()* ( 2.0/L_0 )*u_filt*( 2.0-( 4.0/( L_0*L_0 ) )*( Px()*Px()+Py()*Py() ) )*Px() +
                           oneY()* ( 2.0/L_0 )*u_filt*( 2.0-( 4.0/( L_0*L_0 ) )*( Px()^( 2 )+Py()^( 2 ) )*Py() ) +
                           oneZ()* 2.0*U_0*( 1.0-( 4.0/( L_0*L_0 ) )*( Px()^( 2 )+Py()^( 2 ) ) )*( 1.0-( 4.0*u_filt/( U_0*L_0 ) )*Pz() ) );

        // velocity in the Arterial wall (could have come from an Darcy code )
        u_w = vf::project( Uh_wall,
                           elements( *M_mesh_u_wall ),
                           oneX()* ( Klag/porosity )*u_filt*( L_0/2.0 )*Px()/( Px()*Px()+Py()*Py() ) +
                           oneY()* ( Klag/porosity )*u_filt*( L_0/2.0 )*Py()/( Px()*Px()+Py()*Py() ) );


        // update P (in the future may want to update D, k, f .... whatever)
        c_map["P_w"] = vf::project( Ch_wall, elements( *M_mesh_c ), constant( Pr ) );
        c_map["P_l"] = vf::project( Ch_lumen, elements( *M_mesh_c ), constant( Pr ) );

        //concentrationStep( dt, time_discr_lumen, time_discr_wall, c_l, u_l, c_w, u_w );


        velocity_lumen_space_type::P1Lagrange::p1_type::element_type u_l_p1 = Uh_lumen_p1( u_l );
        velocity_wall_space_type::P1Lagrange::p1_type::element_type u_w_p1 = Uh_wall_p1( u_w );

        timeset_type::step_ptrtype ts_step_wl = ts_wl->step( t );
        ts_step_wl->setMesh( M_mesh_c );
        ts_step_wl->addNodalScalar( "c_wall_lumen", c_l.size(), c_l.begin(), c_l.end() );
        //ts_step_wl->addNodalScalar( "c_w", c_l.size(), c_l.begin(), c_l.end() );

        timeset_type::step_ptrtype ts_step_l = ts_l->step( t );
        ts_step_l->setMesh( Uh_lumen_p1.mesh() );
        ts_step_l->addNodalVector( "u_lumen", u_l_p1.size(), u_l_p1.begin(), u_l_p1.end() );
        //ts_step_l->addNodalVector( "u_l", u.size(), u.begin(), u.end() );

        timeset_type::step_ptrtype ts_step_w = ts_w->step( t );
        ts_step_w->setMesh( Uh_wall_p1.mesh() );
        ts_step_w->addNodalVector( "u_wall", u_w_p1.size(), u_w_p1.begin(), u_w_p1.end() );
        //ts_step_w->addNodalVector( "u_l", u.size(), u.begin(), u.end() );

        wall_lumen.save();
        wall.save();
        lumen.save();
    }
}
} // Feel

int main( int argc,  char** argv )
{
    using namespace Feel;

    Feel::po::options_description test( "twodomainsmt options" );
    test.add_options()
    ( "dt", Feel::po::value<double>()->default_value( 0.1 ), "time step value" )
    ( "T", Feel::po::value<double>()->default_value( 1.0 ), "final Time" )
    ( "nel", Feel::po::value<int>()->default_value( 9473 ), "number of elements" )
    ;

    TwoDomainsMTApp app( argc, argv, makeAbout(), test );
    VLOG(1) << "N process: " << Application::nProcess() << "\n"
            << "Id : " << Application::processId() << "\n";

    app.run();
}
#else
int main( int argc,  char** argv )
{
    // do nothing here
}
#endif // FEELPP_HAS_PETSC_H
