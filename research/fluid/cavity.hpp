/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-06-15

  Copyright (C) 2006 EPFL
  Copyright (C) 2006,2007,2008 Université Joseph Fourier (Grenoble I)

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
   \file cavity.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-06-15
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelcore/timermap.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/oseen.hpp>

#include <feel/feelvf/vf.hpp>

// tag of the lid
const int LID = 1;
std::string createCavity( double );


/*
  makeOptions: serves to determine which options are available.
           (note: also for inline help  :-)
*/
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description cavityoptions( "Cavity options" );
    cavityoptions.add_options()
    ( "dt", Feel::po::value<double>()->default_value( 0.1 ), "time step value" )
    ( "ft", Feel::po::value<double>()->default_value( 10 ), "final time value" )

    ( "nu", Feel::po::value<double>()->default_value( 0.01 ), "viscosity value" )
    ( "rho", Feel::po::value<double>()->default_value( 1.0 ), "fluid density value" )

    ( "fixpointtol", Feel::po::value<double>()->default_value( 1.e6 ),
      "Convergence tolerance for fixed point sub-iterations" )
    ( "fixpointmaxiter", Feel::po::value<int>()->default_value( 10 ),
      "Maximum number of fixed point sub-iterations" )

    ( "peps", Feel::po::value<double>()->default_value( 0 ), "epsilong for pressure term" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "export", Feel::po::value<int>()->default_value( 1 ), "export results(ensight, data file(1D)" )
    ;
    return cavityoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "cavity" ,
                           "cavity" ,
                           "0.1",
                           "2D and 3D Driven Cavity",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2007,2008 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


namespace Feel
{
template<int Dim>
class Cavity : public Application
{
    typedef Application super;
public:

    // constants
    static const uint16_type imOrder = 5; // interpolation order
    static const uint16_type uOrder = 2;  // velocity field polynomial order
    static const uint16_type pOrder = 1;  // pressure  polynomial order

    // -- TYPEDEFS --
    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Simplex<Dim, 1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    /* Omega is "dim" dimensional.
       We want a lagrange vectorial (in "dim" dimensions") basis of continuous functions
       for the velocity and
       a lagrange scalar basis of continuous functions for the pressure
    */
    typedef Lagrange<uOrder, Vectorial> basis_u_type;
    typedef Lagrange<pOrder, Scalar> basis_p_type;
    typedef bases<basis_u_type,basis_p_type> basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    // element_type is the type of an element of $V_h$
    typedef typename space_type::element_type element_type;
    // $V_h^u$: velocity function
    typedef typename element_type::template sub_element<0>::type element_0_type;
    // $V_h^p$: pressure function
    typedef typename element_type::template sub_element<1>::type element_1_type;
    /*quadrature*/
    typedef IM<Dim, imOrder, value_type, Simplex> im_type;
    // this is an example of template metaprogramming
    template<int O>
    struct MyIm
    {
        typedef IM<Dim, O, value_type, Simplex> type;
    };


    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef typename Exporter<mesh_type>::timeset_type timeset_type;

    Cavity( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backendNS( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        M_nu( this->vm()["nu"].template as<double>() ),
        M_rho( this->vm()["rho"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "cavity" ) ),
        timers()
    {
        this->printInfo();

        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "cavity" );
    }

    void printInfo() const
    {
        LOG(INFO) << "[Cavity] hsize = " << meshSize << "\n";
        LOG(INFO) << "[Cavity] M_nu = " << M_nu << "\n";
        LOG(INFO) << "[Cavity] M_rho = " << M_rho << "\n";
        LOG(INFO) << "[Cavity] export = " << this->vm().count( "export" ) << "\n";
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
     * export results to ensight format (enabled by  --export cmd line options)
     * here u is both the velocity and the pressure
     */
    void exportResults( double time, element_type& u );

private:

    backend_ptrtype M_backendNS;
    double meshSize;
    double M_nu;
    double M_rho;

    boost::shared_ptr<export_type> exporter;
    typename export_type::timeset_ptrtype timeSet;

    // list of execution times for different tasks (e.g., mesh, solver,...)
    std::map<std::string,std::pair<boost::timer,double> > timers;
}; // Cavity

template<int Dim>
typename Cavity<Dim>::mesh_ptrtype
Cavity<Dim>::createMesh( double meshSize )
{
    timers["mesh"].first.restart();
    mesh_ptrtype mesh( new mesh_type );

    Gmsh gmsh;
    ImporterGmsh<mesh_type> import( gmsh.generate( entity_type::name().c_str(), createCavity( meshSize ) ) );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    LOG(INFO) << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // Cavity::createMesh


template<int Dim>
void
Cavity<Dim>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;
    this->changeRepository( boost::format( "%1%/%2%/P%3%P%4%/h_%5%/nu_%6%" )
                            % this->about().appName()
                            % entity_type::name()
                            % 2 % 1
                            % this->vm()["hsize"].template as<double>()
                            % this->vm()["nu"].template as<double>()
                          );
    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );
    /*
     * The function space and some associate elements are then defined
     */
    space_ptrtype Xh = space_ptrtype( space_type::New( mesh ) );
    element_type Un( Xh, "Un" );


    /*
     * a quadrature rule for numerical integration
     */
    im_type im;

    // Note: peps = 0 leads to a singular matrix, because there it contains
    //       a null space of dimension one (the constant pressures)
    //       However, this is not a problem for the Krylov solver, because
    //       the Krylov subspace built up will not contain the constant
    //       unless you start with it as an initial vector. The incomplete
    //       factorization preconditioner doesn't have a problem either, as
    //       it would encounter a zero pivot only at the last step, where the
    //       factorization would in fact be complete. So peps = 0 is safe.
    const double peps = this->vm()["peps"].template as<double>();
    const double dt = this->vm()["dt"].template as<double>();
    // Temporal horizon, i.e., final time
    const double ft = this->vm()["ft"].template as<double>();
    double fixpointTol = this->vm()["fixpointtol"].template as<double>();
    double fixpointMaxiter = this->vm()["fixpointmaxiter"].template as<int>();
    boost::timer ti;

    std::set<flag_type> dirichletFlags;
    std::set<flag_type> neumannFlags;

    for ( flag_type flag=1; flag<=2; ++flag )
    {
        dirichletFlags.insert( flag );
    }

    Oseen<space_type, imOrder, Simplex> oseen( Xh, M_backendNS, dirichletFlags, neumannFlags, this->vm() );


    TimerMap timers;
    timers["timer overall time"].restart();

    for ( double t = dt; t <= ft; t += dt )
    {
        timers["timer per iteration"].restart();

        LOG(INFO) << "--------------------------------------------------------------------------------\n";
        LOG(INFO) << "T = " << t << "s\n";
        LOG(INFO) << "ft = " << ft << "s, dt = " << dt << ", M_nu = " << M_nu << "\n";

        timers["timer oseen update per iteration"].reset();
        timers["timer oseen update"].restart();

        // store the last time step approximation
        //Un = oseen.solution();

        oseen.update( /* itRan = */ elements( *mesh ),
                                    /* sigma = */ constant( 1.0/dt ),
                                    /* nuInc = */ M_nu,
                                    /* nuAbs = */ M_nu,
                                    ///* beta  = */ idv( oseen.velocity() ),
                                    /* beta  = */ constant( 0 )*oneX(),
                                    /* f     = */ idv( oseen.velocity() )/dt,
                                    /* c     = */ constant( peps ),
                                    /* g     = */ chi( ( Py() > 0.99 ) && ( Px() > 0.01 && Px() < 0.99 ) )*oneX(),
                                    /* noSlip= */ 1.0,
                                    /* updtJ = */ true  );
        timers["timer oseen update per iteration"].accumulate();
        timers["timer oseen update"].accumulate();

        timers["timer solve per iteration"].reset();
        timers["timer solve"].restart();
        oseen.solve();
        timers["timer solve per iteration"].accumulate();
        timers["timer solve"].accumulate();

        LOG(INFO) << "[Cavity] t = " << t << "\n";
        //LOG(INFO) << "[Cavity] #subiter = " << subiter << "\n";

        this->exportResults( t, oseen.solution() );

        timers["timer overall time"].accumulate();

        timers.report( std::string( "cavity" ) );
    }

    timers.report( std::string( "cavity" ) );

} // Cavity::run

template<int Dim>
void
Cavity<Dim>::exportResults( double time, element_type& u )
{
    double dt = this->vm()["dt"].template as<double>();

    if ( int( time/dt ) % this->vm()["export"].template as<int>() == 0 )
    {
        typename timeset_type::step_ptrtype timeStep = timeSet->step( time );
        timeStep->setMesh( u.functionSpace()->mesh() );
        timeStep->add( "Velocity", u.template element<0>() );
        timeStep->add( "Presssure", u.template element<1>() );
        exporter->save();

    }
} // Cavity::export
} // Feel

/*
  Reminder:
  idt, gradt, ... is used to refer to a trial function (or vector field)
  eg.: idt(u) represents _any_ function of type u
  id, grad, ...   is used to refer to a test function (or vector field)
  idv, gradv, ... is the _value_ of the argument
  eg.: idv(u) is the value of the instance u

  trans is the transpose of the argument

  in this example, you have that
  idt(u) \equiv idt(v)
  id(u)  \equiv id(v)
  BUT idv(u) \neq idv(v) !!!
*/
