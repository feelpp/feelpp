/* -*- mode: c++; coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-07

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

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
   \file residualestimator.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-07-15
 */
/** include predefined life command line options */
#include <life/options.hpp>

/** include linear algebra backend */
#include <life/lifealg/backend.hpp>

/** include function space class */
#include <life/lifediscr/functionspace.hpp>

/** include helper function to define \f$P_0\f$ functions associated with regions  */
#include <life/lifediscr/region.hpp>

/** include integration methods */
#include <life/lifepoly/im.hpp>

/** include gmsh mesh importer */
#include <life/lifefilters/gmsh.hpp>

/** include exporter factory class */
#include <life/lifefilters/exporter.hpp>

/** include  polynomialset header */
#include <life/lifepoly/polynomialset.hpp>

/** include  the header for the variational formulation language (vf) aka FEEL++ */
#include <life/lifevf/vf.hpp>

/** use Life namespace */
using namespace Life;
using namespace Life::vf;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Life::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description residualestimatoroptions("ResidualEstimator options");
    residualestimatoroptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.2), "mesh size")
        ("dim", po::value<int>()->default_value( 0 ), "dimension of the geometry( 0: all three, 1, 2 or 3")
        ("shape", Life::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        ("weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ("penaldir", Life::po::value<double>()->default_value( 20 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ("alpha", Life::po::value<double>()->default_value( 3 ), "Regularity coefficient for function f")
        ;
    return residualestimatoroptions.add( Life::life_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Life::Application subclass.
 *
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
{
    AboutData about( "residualestimator" ,
                     "residualestimator" ,
                     "0.2",
                     "nD(n=1,2,3) Residual Estimator on Laplacian equation on simplices or simplex products",
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


/**
 * \class ResidualEstimator
 *
 * Laplacian Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
template<int Dim>
class ResidualEstimator
    :
    public Simget
{
    typedef Simget super;
public:

    //! Polynomial order \f$P_1\f$
    static const uint16_type Order = 1;

    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    //! sparse matrix type associated with backend
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    //! sparse matrix type associated with backend (shared_ptr<> type)
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    //! vector type associated with backend
    typedef typename backend_type::vector_type vector_type;
    //! vector type associated with backend (shared_ptr<> type)
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 1
    typedef Simplex<Dim> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! function space that holds piecewise constant (\f$P_0\f$) functions (e.g. to store material properties or partitioning
    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar,Discontinuous> > > p0_space_type;
    typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
    //! an element type of the \f$P_0\f$ discontinuous function space
    typedef typename p0_space_type::element_type p0_element_type;


    //! the basis type of our approximation space
    typedef bases<Lagrange<1,Scalar> > p1_basis_type;
    typedef FunctionSpace<mesh_type, p1_basis_type> p1_space_type;
    typedef typename p1_space_type::element_type p1_element_type;


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

    /**
     * Constructor
     */
    ResidualEstimator( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        M_backendP1( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( "gmsh", this->about().appName() ) ),
        dim( this->vm()["dim"].template as<int>() ),
        shape( this->vm()["shape"].template as<std::string>() )
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:

    /**
     * solve the system \f$D u = F\f$
     *
     * \param in D sparse matrix
     * \param inout u solution of the system
     * \param in F vector representing the right hand side of the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );
    void solveP1( sparse_matrix_ptrtype& D, p1_element_type& u, vector_ptrtype& F );

private:

    //! linear algebra backend
    backend_ptrtype M_backend;
    backend_ptrtype M_backendP1;

    //! mesh characteristic size
    double meshSize;

    int dim;

    //! shape of the domain
    std::string shape;

    export_ptrtype exporter;

    //! Piecewise constant functions space
    p0_space_ptrtype P0h;

}; // ResidualEstimator

template<int Dim> const uint16_type ResidualEstimator<Dim>::Order;

template<int Dim>
void
ResidualEstimator<Dim>::run()
{
    if ( dim && dim != Dim ) return ;
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute ResidualEstimator<" << Dim << ">\n";
    std::vector<double> X( 2 );
    X[0] = meshSize;
    if ( shape == "hypercube" )
        X[1] = 1;
    else // default is simplex
        X[1] = 0;
    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<int Dim>
void
ResidualEstimator<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long n )
{
    if ( X[1] == 0 ) shape = "simplex";
    if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "doc/tutorial/%1%/%2%-%3%/P%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=(boost::format( "%1%-%2%" ) % shape % Dim).str() ,
                                                      _usenames=true,
                                                      _shape=shape,
                                                      _dim=Dim,
                                                      _h=X[0] ) );

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    P0h = p0_space_type::New( mesh );
    std::cout<<"P0h->nLocalDof() = "<<P0h->nLocalDof()<<std::endl;

    space_ptrtype Xh = space_type::New( mesh );
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    element_type gproj( Xh, "v" );
    /** \endcode */

    /** define \f$g\f$ the expression of the exact solution and
     * \f$f\f$ the expression of the right hand side such that \f$g\f$
     * is the exact solution
     */
    /** \code */
    //# marker1 #
    value_type alpha = this->vm()["alpha"].template as<double>();
    value_type pi = M_PI;
    //! deduce from expression the type of g (thanks to keyword 'auto')
    //auto g = (1-Px()*Px())*(1-Py()*Py())*(1-Pz()*Pz())*pow(trans(vf::P())*vf::P(),(alpha/2.0));
    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());
    gproj = vf::project( Xh, elements(mesh), g );

    //! deduce from expression the type of f (thanks to keyword 'auto')
    auto f = pi*pi*Dim*g;
    //# endmarker1 #
    /** \endcode */

    bool weakdir = this->vm()["weakdir"].template as<int>();
    value_type penaldir = this->vm()["penaldir"].template as<double>();
    value_type hsize = this->vm()["hsize"].template as<double>();

    using namespace Life::vf;

    /**
     * Construction of the right hand side. F is the vector that holds
     * the algebraic representation of the right habd side of the
     * problem
     */
    /** \code */
    //# marker2 #
    vector_ptrtype F( M_backend->newVector( Xh ) );
    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements(mesh), f*id(v) )+
        integrate( markedfaces( mesh, mesh->markerName("Neumann") ), gradv(gproj)*vf::N()*id(v) );
    //# endmarker2 #
    if ( this->comm().size() != 1 || weakdir )
        {
            //# marker41 #
            form1( _test=Xh, _vector=F ) +=
                integrate( markedfaces(mesh,mesh->markerName("Dirichlet")), g*(-grad(v)*vf::N()+penaldir*id(v)/hFace()) );
            //# endmarker41 #
        }
    F->close();

    /** \endcode */

    /**
     * create the matrix that will hold the algebraic representation
     * of the left hand side
     */
    //# marker3 #
    /** \code */
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );
    /** \endcode */

    //! assemble $\int_\Omega \nu \nabla u \cdot \nabla v$
    /** \code */
    form2( Xh, Xh, D, _init=true ) =
        integrate( elements(mesh), gradt(u)*trans(grad(v)) );
    /** \endcode */
    //# endmarker3 #

    if ( this->comm().size() != 1 || weakdir )
        {
            /** weak dirichlet conditions treatment for the boundaries marked 1 and 3
             * -# assemble \f$\int_{\partial \Omega} -\nabla u \cdot \mathbf{n} v\f$
             * -# assemble \f$\int_{\partial \Omega} -\nabla v \cdot \mathbf{n} u\f$
             * -# assemble \f$\int_{\partial \Omega} \frac{\gamma}{h} u v\f$
             */
            /** \code */
            //# marker10 #
            form2( Xh, Xh, D ) +=
                integrate( markedfaces(mesh,mesh->markerName("Dirichlet")),
                           -(gradt(u)*vf::N())*id(v)
                           -(grad(v)*vf::N())*idt(u)
                           +penaldir*id(v)*idt(u)/hFace());
            D->close();
            //# endmarker10 #
            /** \endcode */
        }
    else
        {
            /** strong(algebraic) dirichlet conditions treatment for the boundaries marked 1 and 3
             * -# first close the matrix (the matrix must be closed first before any manipulation )
             * -# modify the matrix by cancelling out the rows and columns of D that are associated with the Dirichlet dof
             */
            /** \code */
            //# marker5 #
            D->close();
            form2( Xh, Xh, D ) +=
                on( markedfaces(mesh, mesh->markerName("Dirichlet")), u, F, g );
            //# endmarker5 #
            /** \endcode */

        }
    /** \endcode */


    //! solve the system
    /** \code */
    this->solve( D, u, F );
    /** \endcode */

    //! compute the \f$L_2$ norm of the error
    /** \code */
    //# marker7 #
    double L2error2 =integrate(elements(mesh),(idv(u)-g)*(idv(u)-g) ).evaluate()(0,0);
    double L2error =   math::sqrt( L2error2 );
    double semiH1error2 = integrate(elements(mesh),(gradv(u)-gradv(gproj))*(gradv(u)-gradv(gproj))).evaluate()(0,0);
    double H1error = math::sqrt(L2error2+semiH1error2);


    /*******************residual estimator**********************/

    auto estimatorP0 =   integrate(elements(mesh), vf::h()*vf::h()*trace(hessv(u))*trace(hessv(u)) ).broken(P0h);

     auto estimatorP0_internalfaces = integrate(internalfaces(mesh),
                                              vf::hFace()* (jumpv(gradv(u))) * (jumpv(gradv(u))) ).broken( P0h );

     auto estimatorP0_Neumann = integrate( markedfaces(mesh,mesh->markerName("Neumann")),
                                          vf::hFace()* (gradv(u)*vf::N()-idv(gproj)) * (gradv(u)*vf::N()-idv(gproj)) ).broken( P0h );





     //     auto estimatorP0_Neumann = integrate( markedfaces(mesh,mesh->markerName("Neumann")),
     //           print(vf::hFace()* (gradv(u)*vf::N()-idv(gproj)) * (gradv(u)*vf::N()-idv(gproj)),"neuman:") ).broken( P0h );



    //double number_elem=mesh->numElements();
    int number_elem=P0h->nLocalDof();
    /*    LIFE_ASSERT( P0h->nLocalDof() == number_elem )( P0h->nLocalDof() )( number_elem ).warn( "invalid p0 space" );
    LIFE_ASSERT( estimatorP0.size() == number_elem )( estimatorP0.size() )( number_elem ).warn( "invalid p0 function" );
    LIFE_ASSERT( estimatorP0_internalfaces.size() == number_elem )( estimatorP0_internalfaces.size() )( number_elem ).warn( "invalid p0 internal function" );
    LIFE_ASSERT( estimatorP0_Neumann.size() == number_elem )( estimatorP0_Neumann.size() )( number_elem ).warn( "invalid p0 neumann function" );
    */


    if(Dim==2){
      for(int i=0;i<number_elem;i++){
         if(estimatorP0(i)<0) std::cout<<"estimatorP0("<<i<<") = "<<estimatorP0(i)<<std::endl;
	 if(estimatorP0_internalfaces(i)<0) std::cout<<"faces("<<i<<") = "<<estimatorP0_internalfaces(i)<<std::endl;
	 if(estimatorP0_Neumann(i)<0) std::cout<<"Neumann("<<i<<") = "<<estimatorP0_Neumann(i)<<std::endl;
      }
    }

     estimatorP0+=estimatorP0_internalfaces;
    estimatorP0+=estimatorP0_Neumann;
    estimatorP0.printMatlab("estimator_based_on_residuals.m");
    double estimator =  math::sqrt(estimatorP0.sum()) ;



    std::cout<<"Number of elements : "<<number_elem<<std::endl;
    std::cout<<"real error in L2 norm : "<<L2error<<std::endl;
    std::cout<<"estimated error based on residuals for L2 norm : " << estimator*hsize <<std::endl ;
    std::cout<<"real error in H1 norm : "<<H1error<<std::endl;
    std::cout<<"estimated error based on residuals for H1 norm : " << estimator <<std::endl ;



    std::ofstream file ("FICHIER.txt",std::ios_base::app) ;
    file << "hsize = "<<hsize<<" errorL2 = "<<L2error<<" estimatedL2 = "<<estimator*hsize<<" errorH1 = "<<H1error<<" estimatedH1 = "<<estimator
	 <<" N = "<<number_elem<<" et DIM = "<<Dim<<"\n";

    /*
    //Now we will project the estimatorP0 on a P1 space
    auto P1h = p1_space_type::New( mesh, MESH_CHECK );
    auto M1 = M_backendP1->newMatrix( P1h, P1h );
    auto v1 = P1h->element("v1");
    auto estimatorP1 = P1h->element("estimatorP1");
    form2( P1h, P1h, M1, _init=true ) = integrate( elements(mesh), idt(estimatorP1)*id(v1));
    auto F1 = M_backendP1->newVector( P1h );
    form1( P1h, F1, _init=true ) = integrate( elements(mesh), vf::sqrt(idv(estimatorP0))*id(v1));
    this->solveP1( M1, estimatorP1, F1 );
    */






    //estimatorP1 is now the representation at nodes of the error estimator
    /*******************************/


    //Log() << "||error||_L2=" << L2error << "\n";
    //# endmarker7 #
    /** \endcode */

    //! save the results
    /** \code */
    //! project the exact solution
    element_type exact_solution( Xh, "exact_solution" );
    exact_solution = vf::project( Xh, elements(mesh), g );

    auto u_minus_exact = u;
    u_minus_exact -= exact_solution;
    for(int i=0;i<u_minus_exact.size();i++){
      u_minus_exact(i)=math::abs(u_minus_exact(i));
    }

    auto meas = P0h->element();
    meas = vf::project( P0h, elements(mesh), vf::meas() );
    meas.printMatlab( "meas.m" );
    /*
    export_ptrtype exporter( export_type::New( this->vm(),
                                               (boost::format( "%1%-%2%-%3%" )
                                                % this->about().appName()
                                                % shape
						% Dim).str() ) );*/
    if ( exporter->doExport() )
    {
        Log() << "exportResults starts\n";

        exporter->step(0)->setMesh( mesh );
        exporter->step(0)->add( "unknown", u );
        exporter->step(0)->add( "exact solution", exact_solution);
        //exporter->step(0)->add( "estimated error",estimatorP1);
        //exporter->step(0)->add( "u - exact solution", u_minus_exact) ;
        //exporter->step(0)->add( "measure", meas) ;

        exporter->save();
        Log() << "exportResults done\n";
    }
    /** \endcode */
} // ResidualEstimator::run

//# marker6 #
template<int Dim>
void
ResidualEstimator<Dim>::solve( sparse_matrix_ptrtype& D,
                       element_type& u,
                       vector_ptrtype& F )
{
    //! solve the system, first create a vector U of the same size as
    //! u, then call solve,
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    //! call solve, the second D is the matrix which will be used to
    //! create the preconditionner
    M_backend->solve( D, D, U, F );
    //M_backend->solve( _matrix=D, _solution=U, _rhs=F, _rtolerance=1e-8 );
    //! copy U in u
    u = *U;
} // ResidualEstimator::solve
//# endmarker6 #

template<int Dim>
void
ResidualEstimator<Dim>::solveP1( sparse_matrix_ptrtype& D,
                                 p1_element_type& u,
                            vector_ptrtype& F )
{
  //! solve the system, first create a vector U of the same size as
  //! u, then call solve,
  vector_ptrtype U = M_backendP1->newVector( u.functionSpace() );
  //! call solve, the second D is the matrix which will be used to
  //! create the preconditionner
  M_backendP1->solve( D , D, U, F );
  //M_backendP1->solve( _matrix=D, _solution=U, _rhs=F, _rtolerance=1e-8 );
  //! copy U in u
  u = *U;
} // ElectroHeat::solve


/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    /**
     * create an application
     */
    /** \code */
    Application app( argc, argv, makeAbout(), makeOptions() );
    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }
    /** \endcode */

    /**
     * register the simgets
     */
    /** \code */
    app.add( new ResidualEstimator<1>( app.vm(), app.about() ) );
    app.add( new ResidualEstimator<2>( app.vm(), app.about() ) );
    app.add( new ResidualEstimator<3>( app.vm(), app.about() ) );
    /** \endcode */

    /**
     * run the application
     */
    /** \code */
    app.run();
    /** \endcode */
}





