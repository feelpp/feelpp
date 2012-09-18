/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2007-01-18

  Copyright (C) 2007 EPFL

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file levelset_trilinos.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2007-01-18
 */

#ifndef _LEVELSET_TRILINOS_HPP_
#define _LEVELSET_TRILINOS_HPP_

#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelfilters/exporterensight.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/advreact.hpp>
#include <feel/feelalg/backendtrilinos.hpp>
#include <feel/feelalg/backend_adaptive_reuse_pc.hpp>

#include "reinit_fms.hpp"
#include "reinit_ilp.hpp"
#include "indicator.hpp"

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description levelsetoptions( "LevelSet options" );
    levelsetoptions.add_options()
    ( "dt", Feel::po::value<double>()->default_value( 0.1 ),
      "time step value" )
    ( "ft", Feel::po::value<double>()->default_value( 1 ),
      "Final time value" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ),
      "first h value to start convergence" )
    ( "export", "export results(ensight, data file(1D)" )
    ( "stabcoeff", Feel::po::value<double>()->default_value( 0.1 ),
      "interior penalty stabilization coefficient" )
    ;

    Feel::po::options_description solveroptions( "algebraic solver options" );
    solveroptions.add_options()
    ( "tolerance", Feel::po::value<double>()->default_value( 2.e-10 ),
      "solver tolerance" )
    ( "verbose", Feel::po::value<int>()->default_value( 0 ),
      "(=0,1,2) print solver iterations" )
    ( "maxiter", Feel::po::value<int>()->default_value( 1000 ),
      "set maximum number of iterations" )
    ( "fillin", Feel::po::value<int>()->default_value( 2 ),
      "fill-in for incomplete factorizations" )
    ( "threshold", Feel::po::value<double>()->default_value( 1.e-3 ),
      "threshold for incomplete factorizations" )
    ;
    return levelsetoptions.add( solveroptions );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "levelset" ,
                           "levelset" ,
                           "0.1",
                           "2D and 3D Level Set Test Problem",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2006 EPFL" );

    about.addAuthor( "Christoph Winkelmann", "developer",
                     "christoph.winkelmann@epfl.ch", "" );
    return about;

}


namespace Feel
{
template<int Dim>
class LevelSet
    :
public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type pOrder = 1;

    typedef double value_type;

    /* entity */
#define ENTITY Simplex

    /* mesh */
    typedef Mesh<GeoEntity<ENTITY<Dim, 1, Dim> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    /* bases */
    typedef fusion::vector<fem::Lagrange<Dim, pOrder,
            Scalar, Continuous,
            double, ENTITY> >
            basis_p_type;
    typedef fusion::vector<fem::Lagrange<Dim, 0,
            Scalar, Discontinuous,
            double, ENTITY> >
            basis_i_type;

    /* spaces */
    typedef FunctionSpace<mesh_type, basis_p_type, value_type> space_p_type;
    typedef FunctionSpace<mesh_type, basis_i_type, value_type> space_i_type;
    typedef boost::shared_ptr<space_p_type> space_p_ptrtype;
    typedef boost::shared_ptr<space_i_type> space_i_ptrtype;
    typedef typename space_p_type::element_type element_p_type;
    typedef typename space_i_type::element_type element_i_type;

    /* quadrature for postprocessing */
    static const uint16_type imOrder = 20; // 2*pOrder
    typedef IM<Dim, imOrder, value_type, ENTITY> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef typename Exporter<mesh_type>::timeset_type timeset_type;

    LevelSet( int argc, char** argv, AboutData const& ad )
        :
        super( argc, argv, ad ),
        M_meshSize( this->vm()["hsize"].template as<double>() ),
        M_exporter( new ExporterEnsight<mesh_type>( "levelset" ) ),
        M_timeSet( new timeset_type( "levelset" ) ),
        M_timers(),
        M_im()
    {
        VLOG(1) << "[LevelSet] hsize = " << M_meshSize << "\n";
        VLOG(1) << "[LevelSet] export = "
                << this->vm().count( "export" ) << "\n";

        M_timeSet->setTimeIncrement( this->vm()["dt"].template as<double>() );
        M_exporter->addTimeSet( M_timeSet );
    }

    LevelSet( int argc,
              char** argv,
              AboutData const& ad,
              po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_meshSize( this->vm()["hsize"].template as<double>() ),
        M_exporter( new ExporterEnsight<mesh_type>( "levelset" ) ),
        M_timeSet( new timeset_type( "levelset" ) ),
        M_timers(),
        M_im()
    {
        VLOG(1) << "[LevelSet] hsize = " << M_meshSize << "\n";
        VLOG(1) << "[LevelSet] export = "
                << this->vm().count( "export" ) << "\n";

        M_timeSet->setTimeIncrement( this->vm()["dt"].template as<double>() );
        M_exporter->addTimeSet( M_timeSet );
    }

    LevelSet( LevelSet const& tc )
        :
        super( tc ),
        M_meshSize( tc.M_meshSize ),
        M_exporter( new ExporterEnsight<mesh_type>( "levelset" ) ),
        M_timeSet( new timeset_type( "levelset" ) ),
        M_timers( tc.M_timers ),
        M_im()
    {
        VLOG(1) << "[LevelSet] hsize = " << M_meshSize << "\n";
        VLOG(1) << "[LevelSet] export = "
                << this->vm().count( "export" ) << "\n";

        M_timeSet->setTimeIncrement( this->vm()["dt"].template as<double>() );
        M_exporter->addTimeSet( M_timeSet );
    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptr_type createMesh( double meshSize );

    /**
     * alias for run()
     */
    void operator()()
    {
        run();
    }

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time,
                        element_p_type& phi,
                        element_p_type& psi,
                        element_i_type& kappa );

private:

    double M_meshSize;

    boost::shared_ptr<export_type> M_exporter;
    typename export_type::timeset_ptrtype M_timeSet;

    std::map<std::string,std::pair<boost::timer,double> > M_timers;

    im_type M_im;
}; // LevelSet

template<int Dim>
typename LevelSet<Dim>::mesh_ptr_type
LevelSet<Dim>::createMesh( double meshSize )
{
    M_timers["mesh"].first.restart();
    mesh_ptr_type mesh( new mesh_type );
    std::string fname;

    if ( Application::processId() == 0 )
    {
        GmshHypercubeDomain<Dim,1,ENTITY> td;
        td.setCharacteristicLength( meshSize );
        fname = td.generate( ENTITY<Dim,1,Dim>::name() );
    }

    Application::Broadcast( fname, 0 );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    if ( Application::nProcess() != 1 )
        mesh->partition();

    M_timers["mesh"].second = M_timers["mesh"].first.elapsed();
    VLOG(1) << "[LevelSet] createMesh(): "
            << M_timers["mesh"].second << "\n";

    return mesh;

} // LevelSet::createMesh


template<int Dim>
void
LevelSet<Dim>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    this->changeRepository( boost::format( "%1%/h_%2%" )
                            % this->about().appName()
                            % M_meshSize );
    this->setLogs();

    using namespace Feel::vf;

    /*
     * First we create the mesh
     */
    mesh_ptr_type mesh = createMesh( M_meshSize );

    /*
     * The function spaces and some associate elements are then defined
     */
    M_timers["init"].first.restart();
    space_p_ptrtype space_p = space_p_type::New( mesh );
    space_i_ptrtype space_i = space_i_type::New( mesh );
    //space_i->dof()->showMe();
    element_p_type phi( space_p, "phi" );
    element_p_type psi( space_p, "psi" );
    element_i_type kappa( space_i, "kappa" );

    typedef BackendTrilinos backendS_type;
    typedef boost::shared_ptr<backendS_type> backendS_ptrtype;

    typedef BackendAdaptiveReusePC<backendS_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    backend_ptrtype backend( new backend_type );
    backend->set_noisy( this->vm()["verbose"].template as<int>() );
    backend->set_maxiter( this->vm()["maxiter"].template as<int>() );
    backend->set_tol( this->vm()["tolerance"].template as<double>() );
    backend->set_fillin( this->vm()["fillin"].template as<int>() );
    backend->set_threshold( this->vm()["threshold"].template as<double>() );

    backendS_ptrtype backendSymmP1( new backendS_type );
    backendSymmP1->set_noisy( this->vm()["verbose"].template as<int>() );
    backendSymmP1->set_maxiter( this->vm()["maxiter"].template as<int>() );
    backendSymmP1->set_tol( this->vm()["tolerance"].template as<double>() );
    backendSymmP1->set_fillin( this->vm()["fillin"].template as<int>() );
    backendSymmP1->set_threshold( this->vm()["threshold"].template as<double>() );
    backendSymmP1->set_symmetric( true );

    backendS_ptrtype backendSymmP0( new backendS_type );
    backendSymmP0->set_noisy( this->vm()["verbose"].template as<int>() );
    backendSymmP0->set_maxiter( this->vm()["maxiter"].template as<int>() );
    backendSymmP0->set_tol( this->vm()["tolerance"].template as<double>() );
    backendSymmP0->set_fillin( this->vm()["fillin"].template as<int>() );
    backendSymmP0->set_threshold( this->vm()["threshold"].template as<double>() );
    backendSymmP0->set_symmetric( true );

    typedef __typeof__( elements( *mesh ) ) IteratorRange;
    ReinitializerFMS<space_p_type, IteratorRange>
    reinitializerFMS( space_p, elements( *mesh ) );
    ReinitializerILP<space_p_type, backendS_type, ENTITY>
    reinitializerILP( space_p, backendSymmP1 );
    Indicator<space_p_type, backendS_type, ENTITY>
    indicator( space_p, space_i, backendSymmP0 );

    // -- initial condition
    value_type x0 = 0.5;
    value_type y0 = 0.5; // 0.2
    value_type radius = 0.25; // 0.1
    value_type slotWidth = 0.1;
    value_type slotDepth = 0.0;
    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );

    // bubble
    //     psi = project( space_p, elements(*mesh),
    //                    radius - sqrt( pow( Px()-x0, 2.0 ) + pow( Py()-y0, 2.0 ) )
    //                    );
    //     psi = project( space_p, elements(*mesh),
    //                    pow(radius,2.0) - pow(Px()-x0,2.0) - pow(Py()-y0,2.0)
    //                    );

    // Zalesak slotted disk
    psi = project( space_p, elements( *mesh ),
                   max( ( min( slotWidth+( Px()-x0 ),
                               ( min( slotWidth-( Px()-x0 ),
                                      slotDepth-( Py()-y0 ) ) ) ) ),
                        sqrt( pow( Px()-x0,2.0 )+pow( Py()-y0,2.0 ) )-radius ) );
    value_type mass0 = ( pi-std::asin( slotWidth/radius ) )*std::pow( radius,2.0 )
                       - slotWidth *
                       ( std::sqrt( std::pow( radius,2.0 )-std::pow( slotWidth,2.0 ) ) +
                         2.0*slotDepth );

    // circle
    //     psi = project( space_p, elements(*mesh),
    //                    2*(sqrt(pow(Px()-x0,2.0)+pow(Py()-y0,2.0))-radius) );
    //     psi = project( space_p, elements(*mesh),
    //                    pow( pow(Px()-x0,2.0) + pow(Py()-y0,2.0), 2.0*Px()+0.5 ) -
    //                    pow( pow(radius ,2.0)                   , 2.0*Px()+0.5 ) );
    //     value_type mass0 = pi*std::pow(radius,2.0);

    // straight line
    //     psi = project( space_p, elements(*mesh), (Px()-x0)*exp(10*Py()) );
    //     psi = project( space_p, elements(*mesh), 2*(Px()-x0) );
    //     func phi1 = Px()-x0;
    //     value_type mass0 = x0;

    // square
    // psi = project( space_p, elements(*mesh),
    //                2*max( (max(-radius-(Px()-x0),-radius+(Px()-x0))),
    //                          (max(-radius-(Py()-y0),-radius+(Py()-y0))) ) );
    //func phi1 = max(max(-radius-(Px()-x0),-radius+(Px()-x0)),
    //                max(-radius-(Py()-y0),-radius+(Py()-y0)));
    //value_type mass0 = 4*pow(radius,2.0);

    // two circles intersecting
    //     psi = project( space_p, elements(*mesh),
    //                    min(sqrt( pow(Px()-x0+radius/2,2.0) + pow(Py()-y0,2.0) )
    //                           -radius,
    //                           sqrt( pow(Px()-x0-radius/2,2.0) + pow(Py()-y0,2.0) )
    //                           -radius) );
    //     func phi1 = min(sqrt( pow(Px()-x0+radius/2,2.0) + pow(Py()-y0,2.0) )
    //                        -radius,
    //                        sqrt( pow(Px()-x0-radius/2,2.0) + pow(Py()-y0.2.0) )
    //                        -radius);
    // value_type mass0 = pow(radius,2.0)*(4*pi/3+sqrt(3.)/2);

    double mass = integrate( elements( *mesh ), M_im,
                             chi( idv( psi )<0 )
                           ).evaluate()( 0,0 );
    VLOG(1) << "[LevelSet] mass before reinit = " << mass << "\n";
    VLOG(1) << "[LevelSet]   rel. mass error  = " << mass/mass0-1.0 << "\n";

    indicator.update( psi );
    kappa = indicator.indicatorGamma();
    phi = reinitializerILP( psi, kappa );
    phi = reinitializerFMS( phi );

    mass = integrate( elements( *mesh ), M_im,
                      chi( idv( phi )<0 )
                    ).evaluate()( 0,0 );
    VLOG(1) << "[LevelSet] mass after  reinit = " << mass << "\n";
    VLOG(1) << "[LevelSet]   rel. mass error  = " << mass/mass0-1.0 << "\n";

    mass = integrate( elements( *mesh ), M_im,
                      chi( idv( phi )*idv( psi )<0 )
                    ).evaluate()( 0,0 );
    VLOG(1) << "[LevelSet] sign change error  = " << mass << "\n";
    VLOG(1) << "[LevelSet] sign chg. error/h  = " << mass/M_meshSize << "\n";

    mass = std::sqrt( integrate( elements( *mesh ), M_im,
                                 pow( sqrt( trans( gradv( phi ) )*gradv( phi ) ) - 1.0,
                                      2.0 )
                               ).evaluate()( 0,0 ) );
    VLOG(1) << "[LevelSet] distance from dist = " << mass << "\n";

    M_timers["init"].second = M_timers["init"].first.elapsed();

    AdvReact<space_p_type, imOrder, backend_type, ENTITY>
    advreact( space_p, backend );
    advreact.set_stabcoeff( this->vm()["stabcoeff"].template as<double>() );

    double time        = 0;
    this->exportResults( time, phi, psi, kappa );
    double dt          = this->vm()["dt"].template as<double>();
    time               = dt;

    AUTO( betax, 2*pi*( 0.5-Py() ) );
    AUTO( betay, 2*pi*( Px()-0.5 ) );
    //     AUTO( betax,  sin(2.0*pi*Py())*sin(2.0*pi*Px()) );
    //     AUTO( betay, -sin(3.0*pi*Px())*sin(1.0*pi*Py()) );

    AUTO( beta, betax*oneX()+betay*oneY() );

    //     element_p_type vx = project( space_p, elements(*mesh), betax );
    //     element_p_type vy = project( space_p, elements(*mesh), betay );
    //     AUTO( beta, idv(vx)*oneX()+idv(vy)*oneY() );

    const value_type theta = 0.5;

    // --- Time loop
    for ( int iter = 0;
            time-dt/2 < this->vm()["ft"].template as<double>();
            ++iter, time += dt )
    {
        std::cout << "[LevelSet] update advreact\n" << std::flush;
        M_timers["update"].first.restart();
        advreact.update( /* sigma = */ 1.0/dt,
                                       /* beta  = */ theta*beta,
                                       /* f     = */ idv( phi )/dt
                                       /*         */ - ( 1-theta )*( trans( beta )*gradv( phi ) ),
                                       /* g     = */ idv( phi ),
                                       /* updtJ = */ !backend->reusePC() );
        M_timers["update"].second += M_timers["update"].first.elapsed();
        VLOG(1) << "[LevelSet] assembly time: "
                << M_timers["update"].first.elapsed() << "\n";

        std::cout << "[LevelSet] solve  advreact\n" << std::flush;
        M_timers["solve"].first.restart();
        advreact.solve();
        M_timers["solve"].second += M_timers["solve"].first.elapsed();
        VLOG(1) << "[LevelSet] solving  time: "
                << M_timers["solve"].first.elapsed() << "\n";

        psi = advreact.phi();

        mass = integrate( elements( *mesh ), M_im,
                          chi( idv( psi )<0 )
                        ).evaluate()( 0,0 );
        VLOG(1) << "[LevelSet] mass before reinit = " << mass << "\n";

        std::cout << "[LevelSet] reinitialize\n" << std::flush;
        M_timers["reinit"].first.restart();
        indicator.update( psi );
        kappa = indicator.indicatorGamma();
        phi = reinitializerILP( psi, kappa );
        phi = reinitializerFMS( phi );
        M_timers["reinit"].second += M_timers["reinit"].first.elapsed();

        VLOG(1) << "[LevelSet] t = " << time << "\n";
        std::cout << "[LevelSet] t = " << time << "\n";
        mass = integrate( elements( *mesh ), M_im,
                          chi( idv( phi )<0 )
                        ).evaluate()( 0,0 );
        VLOG(1) << "[LevelSet] mass after  reinit = " << mass << "\n";

        mass = integrate( elements( *mesh ), M_im,
                          sqrt( trans( gradv( phi ) )*gradv( phi ) )
                        ).evaluate()( 0,0 );
        VLOG(1) << "[LevelSet] mean gradient magnitude = " << mass << "\n";

        this->exportResults( time, phi, psi, kappa );

    } // time loop

    VLOG(1) << "[Levelset] total timings:\n";
    VLOG(1) << "[Levelset]   init:     " << M_timers["init"].second << "\n";
    VLOG(1) << "[Levelset]   mesh:     " << M_timers["mesh"].second << "\n";
    VLOG(1) << "[Levelset]   assembly: " << M_timers["update"].second << "\n";
    VLOG(1) << "[Levelset]   solving:  " << M_timers["solve"].second << "\n";
    VLOG(1) << "[Levelset]   reinit.:  " << M_timers["reinit"].second << "\n";

} // LevelSet::run

template<int Dim>
void
LevelSet<Dim>::exportResults( double time,
                              element_p_type& phi,
                              element_p_type& psi,
                              element_i_type& kappa )
{
    M_timers["export"].first.restart();

    // -- EXPORT --
    if ( this->vm().count( "export" ) )
    {
        typename timeset_type::step_ptrtype
        timeStep = M_timeSet->step( time );
        timeStep->setMesh( phi.functionSpace()->mesh() );
        timeStep->add( "phi", phi );
        timeStep->add( "psi", psi );
        timeStep->add( "kappa", kappa );
        M_exporter->save();
    } // export

    M_timers["export"].second += M_timers["export"].first.elapsed();
    VLOG(1) << "[LevelSet] exporting time: "
            << M_timers["export"].second << "\n";
} // LevelSet::export

} // Feel

#endif /* _LEVELSET_TRILINOS_HPP_ */
