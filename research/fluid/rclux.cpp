/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-03-17

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file rclux.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-03-17
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feeldiscr/advreact.hpp>

#include <feel/feelvf/vf.hpp>



std::pair<std::string, std::string> createMesh( double hsize );

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description rcluxoptions( "Rclux options" );
    rcluxoptions.add_options()
    // mesh parameters
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )

    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "peps", Feel::po::value<double>()->default_value( 0. ), "epsilon for pressure term" )
    ( "divtol", Feel::po::value<double>()->default_value( 1.e10 ), "divergence tolerance" )

    ( "dt", Feel::po::value<double>()->default_value( 0.05 ), "time step value" )
    ( "ft", Feel::po::value<double>()->default_value( 1 ), "Final time value" )

    ( "umax", Feel::po::value<double>()->default_value( 0.1 ), "amplitude of the poiseuille profile (m/s)" )
    ( "c0", Feel::po::value<double>()->default_value( 100 ), "concentration at inflow (mol/m^3)" )

    // default values: water
    ( "viscosity", Feel::po::value<double>()->default_value( 0.001 ), "viscosity value ()" )
    ( "density", Feel::po::value<double>()->default_value( 1000 ), "density value ()" )

    //
    ( "transport", "add transport equation" )

    // export in matlab
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return rcluxoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "rclux" ,
                           "rclux" ,
                           "0.1",
                           "2D",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

namespace Feel
{
/**
 * \class Rclux
 * \brief
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class Rclux
    :
public Application
{
    typedef Application super;
public:

#define Entity Simplex

    // -- TYPEDEFS --
    static const uint16_type Dim = 2;
    static const uint16_type Order = 2;
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    /*matrix*/
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<2, 1> entity_type;
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<fem::Lagrange<Dim, 0, Scalar, Discontinuous> > > p0_space_type;
    typedef p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<fem::Lagrange<Dim, Order, Scalar, Continuous, double, Entity> > concentration_basis_type;
    typedef fusion::vector<fem::Lagrange<Dim, Order, Vectorial, Continuous, double, Entity> > velocity_basis_type;
    typedef fusion::vector<fem::Lagrange<Dim, Order-1, Scalar, Continuous, double, Entity> > pressure_basis_type;

    /*concentration*/
    typedef FunctionSpace<mesh_type, concentration_basis_type, value_type> concentration_space_type;
    typedef boost::shared_ptr<concentration_space_type> concentration_space_ptrtype;
    typedef concentration_space_type::element_type concentration_type;

    /*velocity*/
    typedef FunctionSpace<mesh_type, velocity_basis_type, value_type> velocity_space_type;
    typedef boost::shared_ptr<velocity_space_type> velocity_space_ptrtype;
    typedef velocity_space_type::element_type velocity_type;

    /*pressure*/
    typedef FunctionSpace<mesh_type, pressure_basis_type, value_type> pressure_space_type;
    typedef boost::shared_ptr<pressure_space_type> pressure_space_ptrtype;
    typedef pressure_space_type::element_type pressure_type;


    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    //typedef IM<Dim, imOrder, value_type, Entity> im_c_type;
    typedef IM<Dim, 3*Order-1, value_type, Simplex> im_u_type;
    typedef IM<Dim, 2*Order,   value_type, Simplex> im_p_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef Exporter<mesh_type>::timeset_type timeset_type;

    typedef AdvReact<concentration_space_type,imOrder,Entity> advreact_type;
    typedef boost::shared_ptr<advreact_type> advreact_ptrtype;

    Rclux( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backendu( backend_type::build( this->vm() ) ),
        M_backendp( backend_type::build( this->vm() ) ),
        M_backendc( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        bcCoeff( this->vm()["bccoeff"].as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "rclux" ) )
    {
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "rclux" );

        this->changeRepository( boost::format( "%1%/%2%/%3%/" )
                                % this->about().appName()
                                % entity_type::name()
                                % this->vm()["hsize"].as<double>()
                              );
        using namespace Feel::vf;


        /*
         * First we create the mesh
         */
        mesh = createMesh( meshSize );

        Ch = concentration_space_type::New( mesh );
        Vh = velocity_space_type::New( mesh );
        Ph = pressure_space_type::New( mesh );

        M_advreact = advreact_ptrtype( new advreact_type( Ch, M_backendc ) );




    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve the system
     */
    void solve( backend_ptrtype& backend, sparse_matrix_ptrtype& D, vector_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, concentration_type& c, velocity_type& v, pressure_type& p );

private:

    backend_ptrtype M_backendu;
    backend_ptrtype M_backendp;
    backend_ptrtype M_backendc;

    double meshSize;
    double bcCoeff;

    mesh_ptrtype mesh;
    concentration_space_ptrtype Ch;
    velocity_space_ptrtype Vh;
    pressure_space_ptrtype Ph;

    sparse_matrix_ptrtype Dcst;
    vector_ptrtype Fcst;

    boost::shared_ptr<export_type> exporter;
    export_type::timeset_ptrtype timeSet;


    advreact_ptrtype M_advreact;
};


Rclux::mesh_ptrtype
Rclux::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    Gmsh gmsh;
    std::string mesh_name, mesh_desc;
    boost::tie( mesh_name, mesh_desc ) = ::createMesh( meshSize );
    std::string fname = gmsh.generate( mesh_name, mesh_desc );

    ImporterGmsh<mesh_type> import( fname );
    import.setVersion( "2.0" );
    mesh->accept( import );
    return mesh;
} // Rclux::createMesh


void
Rclux::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    using namespace Feel::vf;

    concentration_type c( Ch, "c" );
    concentration_type d( Ch, "d" );
    concentration_type cn( Ch, "cn" );

    velocity_type u( Vh, "u" );
    velocity_type v( Vh, "v" );
    velocity_type un( Vh, "un" );

    pressure_type p( Ph, "p" );
    pressure_type q( Ph, "q" );
    pressure_type pn( Ph, "pn" );


    double r1 = 2e-3;

    /*
     * Construction of the right hand side
     *
     * \f$ f = \int_\Omega g * v \f$ where \f$ g \f$ is a vector
     * directed in the \f$ y \f$ direction.
     */
    vector_ptrtype F( M_backendu->newVector( Vh ) );
    vector_ptrtype fp( M_backendp->newVector( Ph ) );

    // --- Construction of regularized laplacian on pressure space Ap
    sparse_matrix_ptrtype Ap( M_backendp->newMatrix( Ph, Ph ) );
    double peps = this->vm()["peps"].as<double>();

    // Note: peps = 0 leads to a singular matrix, because there it contains
    //       a null space of dimension one (the constant pressures)
    //       However, this is not a problem for the Krylov solver, because
    //       the Krylov subspace built up will not contain the constant
    //       unless you start with it as an initial vector. The incomplete
    //       factorization preconditioner doesn't have a problem either, as
    //       it would encounter a zero pivot only at the last step, where the
    //       factorization would in fact be complete. So peps = 0 is safe.
    form2( Ph, Ph, Ap, _init=true ) = integrate( elements( mesh ), im_p_type(),
                                      gradt( p )*trans( grad( q ) )+ peps*idt( p )*id( p )
                                               )
                                      +
                                      integrate( markedfaces( mesh,mesh->markerName( "outflow" ) ), im_p_type(),
                                              - trans( gradt( p )*N() )*id( p )
                                              - trans( grad( p )*N() )*idt( p )
                                              + bcCoeff/hFace() * trans( idt( p ) ) * id( p ) );

    Ap->close();

    if ( this->vm().count( "export-matlab" ) )
        Ap->printMatlab( "Ap.m" );

    LOG(INFO) << "pressure matrix assembled\n";
    // --- Construction of velocity mass matrix Mu
#if VELOCITY_UPDATE
    sparse_matrix_ptrtype Mu( M_backendu->newMatrix( Vh, Vh ) );
    form2( Vh, Vh, Mu, _init=true ) = integrate( elements( mesh ), im_u_type(),
                                      trans( idt( v ) )*id( u )
                                               );
    Mu->close();
#endif

    double viscosity = this->vm()["viscosity"].as<double>();
    double density = this->vm()["density"].as<double>();
    double dt = this->vm()["dt"].as<double>();
    double time  = dt;

    // --- Construction of velocity diffusion-reaction matrix Lu
    sparse_matrix_ptrtype Lu( M_backendu->newMatrix( Vh, Vh ) );
    form2( Vh, Vh, Lu, _init=true ) =
        integrate( elements( mesh ), im_u_type(),
                   viscosity*dt*trace( gradt( u )*trans( grad( u ) ) )+
                   density*trans( idt( u ) )*id( u ) )+
        integrate( markedfaces( mesh,mesh->markerName( "wall" ) ), im_u_type(),
                   viscosity*dt*( - trans( gradt( u )*N() )*id( u )
                                  - trans( grad( u )*N() )*idt( u )
                                  + bcCoeff/hFace() * trans( idt( u ) ) * id( u ) )
                 )+
        integrate( markedfaces( mesh,mesh->markerName( "inflow" ) ), im_u_type(),
                   viscosity*dt*( - trans( gradt( u )*N() )*id( u )
                                  - trans( grad( u )*N() )*idt( u )
                                  + bcCoeff/hFace() * trans( idt( u ) ) * id( u ) )
                 );
    Lu->close();

    if ( this->vm().count( "export-matlab" ) )
        Lu->printMatlab( "Lu.m" );

    LOG(INFO) << "velocity diffusion reaction matrix assembled\n";

    // --- Time loop
    for ( int iter = 0;
            time < this->vm()["ft"].as<double>();
            ++iter, time += dt )
    {
        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "Time: " << time << "s" << " dt=" << dt << " ft=" << this->vm()["ft"].as<double>() << "\n";
        double divTol = this->vm()["divtol"].as<double>();
        double divError = 2*divTol + 1;
        int subiter = 0;

        while ( divError > divTol )
        {
            ++subiter;

            // --- update rhs for solve for intermediate ux and uy
            form1( Vh, F, _init=true ) = integrate( elements( mesh ), im_u_type(),
                                                    ( density*trans( idv( un ) )-dt*gradv( pn ) )*id( u )
                                                  );

            // --- Construction of convection operator on velocity
            //     space
            sparse_matrix_ptrtype A( M_backendu->newMatrix( Vh, Vh ) );
            form2( Vh, Vh, A, _init=true ) =
                integrate( elements( mesh ), im_u_type(),
                           density*dt*trans( gradt( u )*idv( un ) )*id( u ) );
            A->close();


            // --- Construction of convection-diffusion-reaction
            A->addMatrix( 1.0, *Lu );

            double umax = this->vm()["umax"].as<double>();
            AUTO( inflow, vec( constant( 0. ), umax*( 1.0-( Px()^( 2 ) )/( r1*r1 ) ) ) );
            form1( Vh, F ) += integrate( markedfaces( mesh,mesh->markerName( "inflow" ) ), im_u_type(),
                                         dt*viscosity*trans( inflow )* ( bcCoeff/hFace()*id( u )- grad( u )*N() ) );

            if ( this->vm().count( "export-matlab" ) )
            {
                A->printMatlab( "A.m" );
                F->printMatlab( "F.m" );
            }

            /*
             * Solution phase
             */
            // --- solve for u
            this->solve( M_backendu, A, un, F );

            LOG(INFO) << "solve for cdr matrix done\n";

            // rhs for pressure increment
            form1( Ph, fp, _init=true ) = integrate( elements( mesh ), im_p_type(), -divv( un )*id( p )/dt );

            if ( this->vm().count( "export-matlab" ) )
            {
                fp->printMatlab( "fp.m" );
            }

            // --- solve for pressure increment phi
            this->solve( M_backendp, Ap, p, fp );
            LOG(INFO) << "solve for pressure matrix done\n";

            divError = math::sqrt( integrate( elements( mesh ), im_p_type(),
                                              divv( un )^2 ).evaluate()( 0,0 ) );
            LOG(INFO) << "[Splitting] ||div u||_2 = " << divError << "\n";

            // --- update pressure
            pn += p;

#if VELOCITY_UPDATE
            // --- update velocity
            form1( Vh, F ) = integrate( elements( mesh ), im_u_type(),
                                        ( trans( idv( un ) )-dt*gradv( phi ) )*id( U )
                                      );

            this->solve( M_backendu, Mu, un, F );

            divError = std::sqrt( integrate( elements( mesh ), im_p_type(),
                                             divv( un )^2 ).evaluate() );
            LOG(INFO) << "[Splitting] ||div u||_2 = " << divError << "\n";
#endif

        } // inner loop

        LOG(INFO) << "inner divergence loop done\n";
        u = un;

        if ( this->vm().count( "transport" ) )
        {
            double c0 = this->vm()["c0"].as<double>();
            //vec(constant(0.), (r1-Px())*(r1+Px())/(r1*r1) ) );
            //AUTO( g, constant(1.0)*(chi(time < T/2.0) ) );
            M_advreact->update( /* sigma */ constant( 1. )/dt, //+k*I*chi( emarker()==mesh->markerName( "" ) ),
                                            /* beta */  idv( u ),
                                            /* f */     idv( cn )/dt,
                                            /* g */     constant( c0 ) );

            M_advreact->solve();

            cn = M_advreact->phi();
        }

        this->exportResults( time, cn, un, pn );
    }


} // Rclux::run

void
Rclux::solve( backend_ptrtype& backend, sparse_matrix_ptrtype& D, vector_type& u, vector_ptrtype& F  )
{
    boost::timer ti;
    LOG(INFO) << "Solving for matrix size " << D->size1() << " u.size: " << u.size() << " F.size: " << F->size() << "\n";
    vector_ptrtype U( backend->newVector( u.map() ) );
    backend->solve( D, D, U, F );
    u = *U;
    LOG(INFO) << "Solving for matrix size " << D->size1() << " u.size: " << u.size() << " F.size: " << F->size() << " done in " << ti.elapsed() << "s.\n";
} // Rclux::solve


void
Rclux::exportResults( double time,
                      concentration_type& c,
                      velocity_type& v,
                      pressure_type& p                      )
{
    boost::timer ti;
    LOG(INFO) << "exporting results\n";
    timeset_type::step_ptrtype timeStep = timeSet->step( time );
    timeStep->setMesh( c.functionSpace()->mesh() );
    timeStep->add( "pid",
                   regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( c.functionSpace()->mesh() ) ) ) );
    timeStep->add( "concentration", c );
    timeStep->add( "velocity", v );
    timeStep->add( "pressure", p );
    exporter->save();
    LOG(INFO) << "exporting results done in " << ti.elapsed() << "s.\n";
} // Rclux::export
} // Feel


int
main( int argc, char** argv )
{
    using namespace Feel;

    /* define and run application */
    Rclux rclux( argc, argv, makeAbout(), makeOptions() );

    rclux.run();
}




