/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-06

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
   \file twofluid.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-06
 */
#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>


#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporterensight.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/oseen.hpp>
#include <feel/feeldiscr/advreact.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description twofluidoptions( "TwoFluid options" );
    twofluidoptions.add_options()
    ( "dt", Feel::po::value<double>()->default_value( 0.1 ),
      "time step value" )
    ( "ft", Feel::po::value<double>()->default_value( 1.0 ),
      "Final time value" )
    ( "mu+", Feel::po::value<double>()->default_value( 1.0 ),
      "viscosity of fluid +" )
    ( "mu-", Feel::po::value<double>()->default_value( 1.0 ),
      "viscosity of fluid -" )
    ( "rho+", Feel::po::value<double>()->default_value( 1.0 ),
      "density of fluid +" )
    ( "rho-", Feel::po::value<double>()->default_value( 2.0 ),
      "density of fluid -" )
    ( "g", Feel::po::value<double>()->default_value( 9.81 ),
      "gravitation" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ),
      "first h value to start convergence" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ),
      "coefficient for weak Dirichlet conditions" )
    ( "fixpointtol", Feel::po::value<double>()->default_value( 1.e6 ),
      "Convergence tolerance for fixed point sub-iterations" )
    ( "divcomp", Feel::po::value<double>()->default_value( 0.5 ),
      "divergence compensation coefficient in advection (0,0.5,1)" )
    ( "stabcoeff", Feel::po::value<double>()->default_value( 0.5 ),
      "streamline diffusion stabilization coefficient" )
    ( "export", "export results(ensight, data file(1D)" )
    ;

    return twofluidoptions;
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "twofluid" ,
                           "twofluid" ,
                           "0.1",
                           "two fluid problem",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2006 EPFL" );

    about.addAuthor( "Christoph Winkelmann", "developer",
                     "christoph.winkelmann@epfl.ch", "" );
    return about;

}


namespace Feel
{
template<int Dim>
class TwoFluid
    :
public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type uOrder = 2;
    static const uint16_type pOrder = 1;

    typedef double value_type;

    /* mesh */
    typedef Mesh<GeoEntity<Simplex<Dim, 1> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    /* bases */
    typedef fusion::vector<fem::Lagrange<Dim, uOrder,
            Vectorial, Continuous, double> >
            basis_U_type;
    typedef fusion::vector<fem::Lagrange<Dim, uOrder,
            Scalar, Continuous, double> >
            basis_u_type;
    typedef fusion::vector<fem::Lagrange<Dim, pOrder,
            Scalar, Continuous, double> >
            basis_p_type;

    typedef fusion::vector<basis_U_type, basis_p_type> basis_type;

    /* spaces */
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_U_type;
    typedef boost::shared_ptr<space_U_type> space_U_ptrtype;
    typedef typename space_U_type::element_type element_U_type;

    /* quadrature for postprocessing */
    typedef IM<Dim, 2*uOrder, value_type, Simplex> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef typename Exporter<mesh_type>::timeset_type timeset_type;

    TwoFluid( int argc, char** argv, AboutData const& ad )
        :
        super( argc, argv, ad ),
        M_meshSize( this->vm()["hsize"].template as<double>() ),
        M_bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        M_rhoM( 0.5 * ( this->vm()["rho+"].template as<double>() +
                        this->vm()["rho-"].template as<double>() ) ),
        M_rhoD( 0.5 * ( this->vm()["rho+"].template as<double>() -
                        this->vm()["rho-"].template as<double>() ) ),
        M_muM ( 0.5 * ( this->vm()["mu+"].template as<double>() +
                        this->vm()["mu-"].template as<double>() ) ),
        M_muD ( 0.5 * ( this->vm()["mu+"].template as<double>() -
                        this->vm()["mu-"].template as<double>() ) ),
        M_g ( this->vm()["g"].template as<double>() ),
        M_exporter( new ExporterEnsight<mesh_type>( "twofluid" ) ),
        M_timeSet( new timeset_type( "twofluid" ) ),
        M_timers(),
        M_stats()
    {
        VLOG(1) << "[TwoFluid] hsize = " << M_meshSize << "\n";
        VLOG(1) << "[TwoFluid] bccoeff = " << M_bcCoeff << "\n";
        VLOG(1) << "[TwoFluid] rho- = " << M_rhoM-M_rhoD << "\n";
        VLOG(1) << "[TwoFluid] rho+ = " << M_rhoM+M_rhoD << "\n";
        VLOG(1) << "[TwoFluid] mu- = " << M_muM-M_muD << "\n";
        VLOG(1) << "[TwoFluid] mu+ = " << M_muM+M_muD << "\n";
        VLOG(1) << "[TwoFluid] export = "
                << this->vm().count( "export" ) << "\n";

        M_timeSet->setTimeIncrement( this->vm()["dt"].template as<double>() );
        M_exporter->addTimeSet( M_timeSet );
    }

    TwoFluid( int argc,
              char** argv,
              AboutData const& ad,
              po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_meshSize( this->vm()["hsize"].template as<double>() ),
        M_bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        M_rhoM( 0.5 * ( this->vm()["rho+"].template as<double>() +
                        this->vm()["rho-"].template as<double>() ) ),
        M_rhoD( 0.5 * ( this->vm()["rho+"].template as<double>() -
                        this->vm()["rho-"].template as<double>() ) ),
        M_muM ( 0.5 * ( this->vm()["mu+"].template as<double>() +
                        this->vm()["mu-"].template as<double>() ) ),
        M_muD ( 0.5 * ( this->vm()["mu+"].template as<double>() -
                        this->vm()["mu-"].template as<double>() ) ),
        M_g ( this->vm()["g"].template as<double>() ),
        M_exporter( new ExporterEnsight<mesh_type>( "twofluid" ) ),
        M_timeSet( new timeset_type( "twofluid" ) ),
        M_timers(),
        M_stats()
    {
        VLOG(1) << "[TwoFluid] hsize = " << M_meshSize << "\n";
        VLOG(1) << "[TwoFluid] bccoeff = " << M_bcCoeff << "\n";
        VLOG(1) << "[TwoFluid] rho- = " << M_rhoM-M_rhoD << "\n";
        VLOG(1) << "[TwoFluid] rho+ = " << M_rhoM+M_rhoD << "\n";
        VLOG(1) << "[TwoFluid] mu- = " << M_muM-M_muD << "\n";
        VLOG(1) << "[TwoFluid] mu+ = " << M_muM+M_muD << "\n";
        VLOG(1) << "[TwoFluid] export = "
                << this->vm().count( "export" ) << "\n";

        M_timeSet->setTimeIncrement( this->vm()["dt"].template as<double>() );
        M_exporter->addTimeSet( M_timeSet );
    }

    TwoFluid( TwoFluid const& tc )
        :
        super( tc ),
        M_meshSize( tc.M_meshSize ),
        M_bcCoeff( tc.M_bcCoeff ),
        M_rhoM( tc.M_rhoM ),
        M_rhoD( tc.M_rhoD ),
        M_muM( tc.M_muM ),
        M_muD( tc.M_muD ),
        M_g ( tc.g ),
        M_exporter( new ExporterEnsight<mesh_type>( "twofluid" ) ),
        M_timeSet( new timeset_type( "twofluid" ) ),
        M_timers( tc.M_timers ),
        M_stats( tc.M_stats )
    {
        VLOG(1) << "[TwoFluid] hsize = " << M_meshSize << "\n";
        VLOG(1) << "[TwoFluid] bccoeff = " << M_bcCoeff << "\n";
        VLOG(1) << "[TwoFluid] rho- = " << M_rhoM-M_rhoD << "\n";
        VLOG(1) << "[TwoFluid] rho+ = " << M_rhoM+M_rhoD << "\n";
        VLOG(1) << "[TwoFluid] mu- = " << M_muM-M_muD << "\n";
        VLOG(1) << "[TwoFluid] mu+ = " << M_muM+M_muD << "\n";
        VLOG(1) << "[TwoFluid] export = "
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
                        element_U_type& U,
                        element_p_type& phi );

private:

    double M_meshSize;
    double M_bcCoeff;
    double M_rhoM;
    double M_rhoD;
    double M_muM;
    double M_muD;
    double M_g;

    boost::shared_ptr<export_type> M_exporter;
    typename export_type::timeset_ptrtype M_timeSet;

    std::map<std::string,std::pair<boost::timer,double> > M_timers;
    std::map<std::string,double> M_stats;
}; // TwoFluid

template<int Dim>
typename TwoFluid<Dim>::mesh_ptr_type
TwoFluid<Dim>::createMesh( double meshSize )
{
    M_timers["mesh"].first.restart();
    mesh_ptr_type mesh( new mesh_type );

    Gmsh __gmsh;
    __gmsh.setOrder( GMSH_ORDER_ONE );
    std::ostringstream ostr;
    std::ostringstream nameStr;
    std::string fname;

    switch ( Dim )
    {
    case 2:
        fname = __gmsh.generateSquare( "twofluid2d", meshSize );
        break;

    case 3:
        fname = __gmsh.generateCube( "twofluid3d", meshSize );
        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << Dim;
        throw std::logic_error( os.str() );
    }

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    M_timers["mesh"].second = M_timers["mesh"].first.elapsed();
    VLOG(1) << "[timer] createMesh(): " << M_timers["mesh"].second << "\n";
    return mesh;
} // TwoFluid::createMesh


template<int Dim>
void
TwoFluid<Dim>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    this->changeRepository( boost::format( "%1%/h_%2%/Re_%3%" )
                            % this->about().appName()
                            % M_meshSize
                            % ( 1./M_muM ) );
    this->setLogs();

    using namespace Feel::vf;

    /*
     * First we create the mesh
     */
    mesh_ptr_type mesh = createMesh( M_meshSize );
    M_stats["nelt"] = mesh->elements().size();

    /*
     * The function spaces and some associate elements are then defined
     */
    M_timers["init"].first.restart();
    space_U_ptrtype Uh = space_U_type::New( mesh );
    space_u_ptrtype Xh = space_u_type::New( mesh );
    space_p_ptrtype Yh = space_p_type::New( mesh );
    //Xh->dof()->showMe();
    element_U_type U( Uh, "U" );
    element_u_type ux( Xh, "ux" );
    element_u_type uy( Xh, "uy" );
    element_u_type uxn( Xh, "uxn" );
    element_u_type uyn( Xh, "uyn" );
    element_p_type p( Yh, "p" );
    element_p_type phi( Yh, "phi" );
    element_p_type phin( Yh, "phin" );
    M_timers["init"].second = M_timers["init"].first.elapsed();
    M_stats["ndof"] = Dim * Xh->nDof() + 2*Yh->nDof();

    // -- initial condition : bubble
    value_type xb = 0.5;
    value_type yb = 0.2;
    value_type rb = 0.1;
    phi = project( Yh, elements( *mesh ),
                   rb - sqrt( ( ( Px()-xb )^2.0 ) + ( ( Py()-yb )^2.0 ) )
                 );
    phin = phi;
    double time = 0;
    this->exportResults( time, U, ux, uy, p, phi );

    const uint16_type imOrder = 20; // 20, to have many points, costly
    Oseen<space_u_type, space_p_type, imOrder> oseen( Xh, Yh );
    AdvReact<space_p_type, imOrder> advReact( Yh );

    oseen.set_noisy( this->vm()["verbose"].template as<int>() );
    oseen.set_maxiter( this->vm()["maxiter"].template as<int>() );
    oseen.set_fillin( this->vm()["fillin"].template as<int>() );
    oseen.set_threshold( this->vm()["threshold"].template as<double>() );
    oseen.set_bccoeff( M_bcCoeff );
    oseen.set_tol( this->vm()["tolerance"].template as<double>() );

    advReact.set_noisy( this->vm()["verbose"].template as<int>() );
    advReact.set_maxiter( this->vm()["maxiter"].template as<int>() );
    advReact.set_fillin( this->vm()["fillin"].template as<int>() );
    advReact.set_threshold( this->vm()["threshold"].template as<double>() );
    advReact.set_tol( this->vm()["tolerance"].template as<double>() );
    advReact.set_stabcoeff( this->vm()["stabcoeff"].template as<double>() );

    double dt          = this->vm()["dt"].template as<double>();
    time               = dt;
    double fixpointTol = this->vm()["fixpointtol"].template as<double>();
    double divcomp     = this->vm()["divcomp"].template as<double>();

    // --- Time loop
    for ( int iter = 0;
            time-dt/2 < this->vm()["ft"].template as<double>();
            ++iter, time += dt )
    {
        double fixpointErr = 2*fixpointTol+1.0;
        uint16_type subiter;

        for ( subiter = 0; fixpointErr>fixpointTol; ++subiter )
        {
            /* rho = M_rhoM + M_rhoD*sign(phin) */
            std::cout << "[TwoFluid] update oseen\n" << std::flush;
            oseen.update( /* sigma = */ ( M_rhoM +
                                          M_rhoD*sign( idv( phin ) )
                                        ) / dt,
                                        /* nu    = */ M_muM + M_muD*sign( idv( phin ) ),
                                        /* betax = */ ( M_rhoM +
                                                M_rhoD*sign( idv( phin ) )
                                                      ) * idv( uxn ),
                                        /* betay = */ ( M_rhoM +
                                                M_rhoD*sign( idv( phin ) )
                                                      ) * idv( uyn ),
                                        /* fx    = */ ( M_rhoM +
                                                M_rhoD*sign( idv( phin ) )
                                                      ) * ( idv( ux )/dt ),
                                        /* fy    = */ ( M_rhoM +
                                                M_rhoD*sign( idv( phin ) )
                                                      ) * ( idv( uy )/dt - M_g ),
                                        /* gx    = */ 0.0,
                                        /* gy    = */ 0.0 );
            std::cout << "[TwoFluid] solve  oseen\n" << std::flush;
            oseen.solve();

            const element_u_type& uxnn = oseen.velocityX();
            const element_u_type& uynn = oseen.velocityY();

            std::cout << "[TwoFluid] update advreact\n" << std::flush;
            //                     advReact.update( /* sigma = */ 1.0/dt + 0.5*divcomp*
            //                                      /*         */ ( dxv(uxnn) + dyv(uynn) ),
            //                                      /* betax = */ 0.5*idv(uxnn),
            //                                      /* betay = */ 0.5*idv(uynn),
            //                                      /* f     = */ idv(phi)/dt - 0.5*divcomp*
            //                                      /*         */ ( dxv(ux) + dyv(uy) )
            //                                      /*         */ -0.5*(idv(ux)*dxv(phi)+
            //                                                          idv(uy)*dyv(phi)),
            //                                      /* g     = */ 0.0 // no inflow anyway!
            //                                      );
            advReact.update( /* sigma = */ 1.0/dt + divcomp*
                                           /*         */ ( dxv( uxnn ) + dyv( uynn ) ),
                                           /* betax = */ idv( uxnn ),
                                           /* betay = */ idv( uynn ),
                                           /* f     = */ idv( phi )/dt,
                                           /* g     = */ 0.0 // no inflow anyway!
                           );
            std::cout << "[TwoFluid] solve  advreact\n" << std::flush;
            advReact.solve();
            std::cout << "[TwoFluid] done\n" << std::flush;

            const element_p_type& phinn = advReact.phi();

            fixpointErr =
                std::sqrt( integrate( elements( *mesh ), im_type(),
                                      ( ( idv( uxn )-idv( uxnn ) )^2 )
                                      +( ( idv( uyn )-idv( uynn ) )^2 )
                                      +( ( idv( phin )-idv( phinn ) )^2 )
                                    ).evaluate()( 0,0 ) );

            VLOG(1) << "[TwoFluid] fixpoint iteration " << subiter
                    << "\n";
            VLOG(1) << "[TwoFluid] fixpoint error = " << fixpointErr
                    << "\n";

            uxn = uxnn;
            uyn = uynn;
            phin = phinn;

        }

        /*
         * Post processing phase
         */
        ux = uxn;
        uy = uyn;
        p = oseen.pressure();
        phi = phin;

        VLOG(1) << "[TwoFluid] t = " << time << "\n";
        VLOG(1) << "[TwoFluid] #subiter = " << subiter << "\n";

        double divError = std::sqrt( integrate( elements( *mesh ),
                                                im_type(),
                                                ( dxv( ux )+dyv( uy ) )^2
                                              ).evaluate()( 0,0 ) );
        VLOG(1) << "[TwoFluid] ||div u||_2 = " << divError << "\n";

        double mass = integrate( elements( *mesh ), im_type(),
                                 chi( idv( phi )>0 )
                               ).evaluate()( 0,0 );
        VLOG(1) << "[TwoFluid] mass = " << mass << "\n";

        U.comp( X ) = ux;
        U.comp( Y ) = uy;

        this->exportResults( time, U, phi );

    } // time loop

    VLOG(1) << "[timer] run():     init: " << M_timers["init"].second << "\n";
    VLOG(1) << "[timer] run(): assembly: " << M_timers["assembly"].second
            << "\n";

} // TwoFluid::run

template<int Dim>
void
TwoFluid<Dim>::exportResults( double time,
                              element_U_type& U,
                              element_p_type& phi )
{
    M_timers["export"].first.restart();

    // -- EXPORT --
    if ( this->vm().count( "export" ) )
    {
        typename timeset_type::step_ptrtype
        timeStep = M_timeSet->step( time );
        timeStep->setMesh( ux.functionSpace()->mesh() );
        timeStep->add( "U", U );
        timeStep->add( "phi", phi );
        M_exporter->save();
    } // export

    M_timers["export"].second = M_timers["export"].first.elapsed();
    VLOG(1) << "[timer] exportResults(): " << M_timers["export"].second << "\n";
} // TwoFluid::export
} // Feel
