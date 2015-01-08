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
#include <feel/feel.hpp>
/** error management **/
// #include <research/hifimagnet/applications/Tools/error.hpp>
#include "error.hpp"

/** use Feel namespace */
using namespace Feel;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
  po::options_description
makeOptions()
{
  po::options_description laplacian_parabolic_eqoptions( "Laplacian_parabolic options" );
  laplacian_parabolic_eqoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
    ( "geofile", Feel::po::value<std::string>()->default_value( "" ), "name of the geofile input")
    ( "geo_depends", Feel::po::value<std::string>()->default_value(""), "list of dependants file")
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
    ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
    ( "penaldir", Feel::po::value<double>()->default_value( 10 ),
      "penalisation parameter for the weak boundary Dirichlet formulation" )
    ( "initial_u", po::value<std::string>()->default_value( "0.0"), "Value or function to initialize T (instationnary mode)")
    ;
  return laplacian_parabolic_eqoptions.add( Feel::feel_options() );
}


/**
 * \class Laplacian_parabolic
 *
 * Laplacian equation with instationnary term solver using continuous approximation spaces
 * solve \f$ \dfrac{\partial u}{\partial t} - nu*\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem = 2
 */
template<int Dim, int Order = 1>
class Laplacian_parabolic
:
  public Simget
{
  typedef Simget super;
  public:

  //! numerical type is double
  typedef double value_type;

  //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 1
  typedef Simplex<Dim> convex_type;
  //! mesh type
  typedef Mesh<convex_type> mesh_type;
  //! mesh shared_ptr<> type
  typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

  //! the basis type of our approximation space
  typedef bases<Lagrange<Order,Scalar> > basis_type;
  //! the approximation function space type
  typedef FunctionSpace<mesh_type, basis_type> space_type;
  //! the approximation function space type (shared_ptr<> type)
  typedef boost::shared_ptr<space_type> space_ptrtype;
  //! an element type of the approximation function space
  typedef typename space_type::element_type element_type;

  //! the exporter factory type
  typedef Exporter<mesh_type> export_type;
  //! the exporter factory (shared_ptr<> type)
  typedef boost::shared_ptr<export_type> export_ptrtype;

  /* BDF discretization */
  typedef Bdf<space_type>  bdf_type;
  typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

  /* ErrorBase */
  typedef ErrorBase<Dim, Order> error_type;
  typedef boost::shared_ptr<error_type> error_ptrtype;

  /**
   * Constructor
   */
  Laplacian_parabolic()
    :
      super(),
      meshSize( doption("hsize") ),
      shape( soption("shape") )
  {
  }

  void run();

  private:

  //! mesh characteristic size
  double meshSize;

  //! shape of the domain
  std::string shape;

  mesh_ptrtype mesh;
}; // Laplacian_parabolic


/**
 *	\brief Defines the edp to be solved in the instationnary case : \f$ dfrac{\partial u}{\partial t} - nu*\Delta u = f \f$
 */
struct transient_edp {
  ex operator()(ex u, std::vector<symbol> vars, std::vector<symbol> p) const
  {
    return (diff(u, p[0], 1) - laplacian(u, vars) );
  };
};

/**
 *	\brief Defines the edp to be solved in the stationnary case : \f$ - nu*\Delta u = f \f$
 */
struct steady_edp {
  ex operator()(ex u, std::vector<symbol> vars) const
  {
    return (-laplacian(u, vars) );
  };
};

/**
 *	\brief Function to compute the equation and find the unknown
 */
template<int Dim, int Order>
  void
Laplacian_parabolic<Dim,Order>::run()
{
  LOG(INFO) << "------------------------------------------------------------\n";
  LOG(INFO) << "Execute Laplacian_parabolic<" << Dim << ">\n";

  std::cout << "------------------------------------------------------------\n";
  std::cout << "Execute Laplacian_parabolic<" << Dim << ">\n";

  Environment::changeRepository( boost::format( "doc/manual/laplacian/%1%/D%2%/P%3%/h_%4%/" )
      % this->about().appName()
      % Dim
      % Order
      % meshSize );

  /**
   * Loading variables from cfg file
   */
  bool weak_dirichlet = ioption("weakdir");
  value_type penaldir = doption("penaldir");
  std::string geofile = soption("geofile");

  // loading exact and rhs
  std::string exact  = soption("error.exact");
  std::string rhs  = soption("error.rhs");
  std::string params = soption("error.params");
  value_type nu = doption("nu");

  // loading time loop variables
  bool steady = boption("bdf.steady");
  double t_final = doption("bdf.time-final");
  double t0 = this->vm()["bdf.time-initial"].template as <double>();
  double dt = this->vm()["bdf.time-step"].template as <double>();

  // loadgin initial temperature expression
  std::string initial_u = soption("initial_u");


  ///**
  // * Creation of a new mesh depending on the information of the geofile
  // */
  //// access to the geofile
  //std::string access_geofile = ( boost::format( "%1%.geo" ) % geofile ).str();

  //gmsh_ptrtype desc_geo;
  //desc_geo=geo(_filename = access_geofile, _dim = Dim, _order = 1, _h = meshSize);

  //// creation of a new mesh depending on the information of the geofile (available in desc_geo)
  //mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
  //    _desc=desc_geo,
  //    _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES );
  ///** \code */
  mesh_ptrtype mesh = loadMesh(_mesh=new mesh_type,_filename=soption("geofile")); 
  /** \endcode */

#if 0
  //auto eit = mesh->beginElementWithProcessId( this->comm().rank() );
  //for( int ne = 0; ne < 10; ++ne )
  mesh->eraseElement( mesh->beginElementWithProcessId( this->comm().rank() ) );
  mesh->eraseElement( mesh->beginElementWithProcessId( this->comm().rank() ) );
  mesh->eraseElement( mesh->beginElementWithProcessId( this->comm().rank() ) );
  mesh->eraseElement( boost::prior(mesh->endElementWithProcessId( this->comm().rank() ) ) );
  std::cout << "n elements after: "  << mesh->numElements() << "\n";
#endif


  /**
   * The function space and some associated elements (functions) are then defined
   */
  /** \code */
  space_ptrtype Xh = space_type::New( mesh );

  // print some information (number of local/global dof in logfile)
  Xh->printInfo();

  element_type u( Xh, "u" );
  element_type v( Xh, "v" );
  element_type Rhs( Xh, "rhs" );
  element_type gproj(Xh, "exact_proj");
  /** \endcode */


  error_ptrtype cvg(new error_type(this->vm(), ""));

  if( !exact.empty() )
  {
    //! [marker11]
    if ( !params.empty() )
      cvg->setParams ( params );
    //! [marker11]
    LOG(INFO) << "Loading function : " << exact << std::endl;
    std::cout << "Loading function : " << exact << std::endl;
    //! [marker12]
    cvg->setSolution(exact, params);
    //! [marker12]
    cvg->print();
  }

  /**
   * Add extra parameters ( t for example )
   */
  /** \code */
  auto parameters = cvg->getParams(); //symbols<Dim>();
  /** \endcode */
  if( parameters.size())
    std::cout << "WARNING -- " << parameters.size() << " parameters have been defined\n";

  auto vars = cvg->getVars(); //symbols<Dim>();

  /// [marker13]
  // if rhs should be computed
  if(cvg->computedrhs())
  {
    // if it is not defined
    if( rhs.empty() )
    {
      // compute it from the exact function
      FEELPP_ASSERT( exact.size() ).error( "no exact function defined, no rhs function defined. Aborting...");

      // with the edp of the parabolic equation
      if ( !steady )
      {
        typedef typename boost::function<ex (ex u, std::vector<symbol> vars, std::vector<symbol> p)> t_edp_type;
        t_edp_type t_edp = transient_edp();
        cvg->setRhs(&t_edp);
      }
      // or with the normal laplacian equation
      else
      {
        typedef typename boost::function<ex (ex u, std::vector<symbol> vars)> edp_type;
        edp_type edp =  steady_edp();
        cvg->setRhs(&edp);
      }
      LOG(INFO) << "computed rhs is : " << cvg->getRhs() << "\n";

    }
    // else it is defined, so ok
    else
    {
      cvg->setRhs();
      LOG(INFO) << "rhs is : " << cvg->getRhs() << "\n";
    }
    std::cout << "rhs is : " << cvg->getRhs() << "\n";
  }
  /// [marker13]

  /**
   * Initializing u, g and f from initial temperature expression
   */

  // initialising T (to 303K if nothing is specified)
  ex initial_u_expr;

  LOG(INFO) << "Loading initial function : " << initial_u << std::endl;
  /** \code */
  initial_u_expr = parse(initial_u, vars, parameters);
  /** \endcode */
  LOG(INFO) << "initial function is : " << initial_u_expr << std::endl;

  if ( !steady && parameters.size() )
  {
    u = cvg->exact_project(Xh, t0);
    gproj = cvg->exact_project(Xh, t0);
  }
  else
    gproj = project(Xh, elements(mesh), expr(initial_u_expr,vars) );

  if(cvg->computedrhs())
    if ( parameters.size() )
      Rhs = cvg->rhs_project(Xh, t0);
    else
      Rhs = cvg->rhs_project(Xh);
  else
    Rhs = vf::project( Xh, elements(mesh), cst(0.0) );



  /**
   * BDF implementation
   */
  // set geometry exporting static
  auto exp = exporter(_mesh = mesh, _geo = "static");

  /** \code */
  /// [marker1]
  // create the BDF structure
  bdf_ptrtype M_bdf = bdf( _space=Xh, _name="u" );

  // start the time at time-initial
  M_bdf->start();
  // create the finite difference polynome of the unknown
  M_bdf->initialize(u);
  /// [marker1]
  /** \endcode */

  // print some information
  if( !steady )
  {
    std::cout<< " \n ****** Time loop ****** \n"
      << "Initial time : "<< t0 << "\n"
      << "Final time : "<< t_final << "\n"
      << "Time step : "<< dt << "\n"
      << "Initial condition : "<< initial_u << "\n"
      << "BDF order : " << M_bdf->timeOrder() << "\n" << std::endl;
  }
  else
    std::cout<< "\n****** Steady state ******"<< std::endl;


  std::cout << "t \t L2error \t H1error" << std::endl;

  double L2Time_error = 0.0;

  // time not depending terms

  /**
   * create the matrix that will hold the algebraic representation
   * of the left hand side (only stationnary terms)
   */
  //! [marker3]
  /** \code */
  auto D = backend()->newMatrix( _test=Xh, _trial=Xh  );
  /** \endcode */

  //! assemble \(\int_\Omega \nu \nabla u \cdot \nabla v\)
  /** \code */
  auto a = form2( _test=Xh, _trial=Xh, _matrix=D );
  a = integrate( _range=elements( mesh ), _expr=nu*gradt( u )*trans( grad( v ) ) );
  /** \endcode */
  //! [marker3]

  if ( weak_dirichlet )
  {
    /** weak dirichlet conditions treatment for the boundaries marked 1 and 3
     * -# assemble \(\int_{\partial \Omega} -\nabla u \cdot \mathbf{n} v\)
     * -# assemble \(\int_{\partial \Omega} -\nabla v \cdot \mathbf{n} u\)
     * -# assemble \(\int_{\partial \Omega} \frac{\gamma}{h} u v\)
     */
    /** \code */
    //! [marker10]
    a += integrate( _range=markedfaces( mesh,"Dirichlet" ),
        _expr= nu * ( -( gradt( u )*vf::N() )*id( v )
          -( grad( v )*vf::N() )*idt( u )
          +penaldir*id( v )*idt( u )/hFace() ) );
    //! [marker10]
    /** \endcode */
  }

  //! assemble \(int_\Omega u^{n+1} v\)
  /** \code */
  /// [marker8]
  if( !steady )
    a += integrate( _range=elements( mesh ), _expr = cst(M_bdf->polyDerivCoefficient( 0 ))*idt(u)*id(v) );
  /// [marker8]

  // time depending term

  do{ // temporal loop start

    if(cvg->computedrhs())
      if ( parameters.size() )
        Rhs = cvg->rhs_project(Xh, t0);
      else
        Rhs = cvg->rhs_project(Xh);
    else
      Rhs = vf::project( Xh, elements(mesh), cst(0.0) );

    if ( !steady && parameters.size() && !cvg->getExactSolution().empty() )
      gproj = cvg->exact_project(Xh, M_bdf->time() );
    else
      gproj = cvg->exact_project(Xh);


    /**
     * Construction of the right hand side. F is the vector that holds
     * the algebraic representation of the right habd side of the
     * problem
     */
    /** \code */
    //# marker2 #
    auto F = backend()->newVector( Xh );
    auto l = form1( _test=Xh, _vector=F );
    l = integrate( _range=elements( mesh ), _expr=idv(Rhs)*id( v ) );
    //# endmarker2 #
    if ( weak_dirichlet )
    {
      //# marker41 #
      l += integrate( _range=markedfaces( mesh,"Dirichlet" ),
          _expr=nu*idv(gproj)*( -grad( v )*vf::N()
            + penaldir*id( v )/hFace() ) );
      //# endmarker41 #
    }

    /** \endcode */

    //! add temporal term to the lhs and the rhs
    /** \code */
    /// [marker9]
    if( !steady )
    {
      auto bdf_poly = M_bdf->polyDeriv();
      l += integrate( _range = elements(mesh), _expr = idv( bdf_poly )*id(v));
    }
    /// [marker9]
    /** \endcode */


    /**
     * add time depending terms for the left hand side
     */
    /** \code */
    auto Dt = backend()->newMatrix( _test=Xh, _trial=Xh );
    auto at = form2( _test=Xh, _trial=Xh, _matrix=Dt );
    /** \endcode */
    if( !weak_dirichlet )
    {
      /** strong(algebraic) dirichlet conditions treatment for the boundaries marked 1 and 3
       * -# first close the matrix (the matrix must be closed first before any manipulation )
       * -# modify the matrix by cancelling out the rows and columns of D that are associated with the Dirichlet dof
       */
      /** \code */
      //# marker5 #
      at = on( _range=markedfaces( mesh, "Dirichlet" ),
          _element=u, _rhs=F, _expr=idv(gproj) );
      //# endmarker5 #
      /** \endcode */

    }

    /** \code */
    at += a;
    /** \endcode */

    //! solve the system
    /** \code */
    //# marker6 #
    backend( _rebuild=true  )->solve( _matrix=Dt, _solution=u, _rhs=F );
    //# endmarker6 #
    /** \endcode */

    /**
     * Computing L2 and H1 error
     */
    /** \code */
    //# marker7 #

    if ( !cvg->getExactSolution().empty() )
    {
      ex solution;
      if( !steady && parameters.size() )
        solution  = cvg->getSolution(M_bdf->time() );
      else
        solution = cvg->getSolution();

      auto g = expr(solution,vars);
      auto gradg = expr<1,Dim,2>(grad(solution,vars), vars );

      /// [marker14]
      double L2error = normL2( _range=elements( mesh ),_expr=( idv( u )-idv(gproj) ) );
      double H1seminorm = math::sqrt( integrate( elements(mesh), (gradv(u) - gradg)*trans(gradv(u) - gradg) ).evaluate()(0,0) );
      double H1error = math::sqrt( L2error*L2error + H1seminorm*H1seminorm);
      /// [marker14]

      LOG(INFO) << "||error||_L2=" << L2error << "\n";
      LOG(INFO) << "||error||_H1=" << H1error << "\n";
      std::cout <<  M_bdf->time() <<"\t" << L2error << "\t" << H1error << "\n";

      L2Time_error = L2Time_error + L2error*L2error;
    }

    //# endmarker7 #
    /** \endcode */

    //! save the results
    /** \code */
    //! exporting exact and approached solution at t
    if ( exp->doExport() )
    {
      LOG(INFO) << "exportResults starts at " << M_bdf->time() << "\n";

      exp->step( M_bdf->time() )->add( "solution", u );
      exp->step( M_bdf->time() )->add( "exact", gproj );
      exp->step( M_bdf->time() )->add( "rhs", Rhs );

      exp->save();
      LOG(INFO) << "exportResults done\n";
    }
    /** \endcode */

    /// [marker15]
    M_bdf->shiftRight(u);
    /// [marker15]

    /// [marker16]
    M_bdf->next();
    /// [marker16]
  }
  while( M_bdf->isFinished() == false );


  if (!steady )
  {
    L2Time_error = math::sqrt(M_bdf->timeStep()*L2Time_error);

    std::cout << "\n***************************\n";
    std::cout << "BDF order : " << M_bdf->timeOrder() << "\t";
    std::cout << "Time Step : "<< M_bdf->timeStep() << "\t";
    std::cout << "Hsize : " << meshSize << "\t";
    std::cout << "L2 Time error = " << L2Time_error << std::endl;
  }
  /* end of bdf implementation */

} // Laplacian_parabolic::run

/**
 * main function: entry point of the program
 */
  int
main( int argc, char** argv )
{
  /**
   * Initialize Feel++ Environment
   */
  Environment env( _argc=argc, _argv=argv,
      _desc=makeOptions(),
      _about=about(_name="parabolic_eq",
        _author="Christophe Prud'homme",
        _email="christophe.prudhomme@feelpp.org") );

  Application app;

  //if ( app.nProcess() == 1 )
  //    app.add( new Laplacian_parabolic_<1>() );
  app.add( new Laplacian_parabolic<2>() );
  //app.add( new Laplacian_parabolic<3>() );

  app.run();

}
