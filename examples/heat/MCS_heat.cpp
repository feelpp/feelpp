/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Julien Calandreau et Anna Effeindzourou M2 MCS
       Date: 2011-11-16

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
   \file equation_chaleur.cpp
   \author Julien Calandreau et Anna Effeindzourou M2 MCS
   \date 2011-11-16
 */
/** include predefined feel command line options */
#include <feel/options.hpp>

/** include linear algebra backend */
#include <feel/feelalg/backend.hpp>

/** include function space class */
#include <feel/feeldiscr/functionspace.hpp>

/** include helper function to define \f$P_0\f$ functions associated with regions  */
#include <feel/feeldiscr/region.hpp>

/** include integration methods */
#include <feel/feelpoly/im.hpp>

/** include gmsh mesh importer */
#include <feel/feelfilters/gmsh.hpp>

/** include exporter factory class */
#include <feel/feelfilters/exporter.hpp>

/** include  polynomialset header */
#include <feel/feelpoly/polynomialset.hpp>

/** include  the header for the variational formulation language (vf) aka FEEL++ */
#include <feel/feelvf/vf.hpp>

#include <feel/feelfilters/geotool.hpp>


/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

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
    po::options_description laplacianoptions("Laplacian options");
    laplacianoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ("shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)")
        ("nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient")
        ("weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ("penaldir", Feel::po::value<double>()->default_value( 10 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ;
    return laplacianoptions.add( Feel::feel_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
{
    AboutData about( "laplacian" ,
                     "laplacian" ,
                     "0.2",
                     "nD(n=1,2,3) Laplacian on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


/**
 * \class Laplacian
 *
 * Laplacian Solver using continuous approximation spaces
 * solve \f$ \frac{\partial T}{\partial t}-\Delta T = f\f$ on \f$\I \times Omega\f$ and \f$T= h1\f$ on \f$\Gamma1\f$
 * and \f$ \frac{\partial T}{\partial n}=h2 \f$  on \f$\Gamma2\f$ and \f$ \frac{\partial T}{\partial n}=h3\f$ \f$\Gamma3\f$
 * and \f$T(0,x)=h4 on \f$\ Omega\f$
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
template<int Dim>
class BlocHeat
    :
    public Simget
{
    typedef Simget super;
public:

    //! Polynomial order \f$P_2\f$
    static const uint16_type Order = 2;

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
    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar, Discontinuous> > > p0_space_type;
    //! an element type of the \f$P_0\f$ discontinuous function space
    typedef typename p0_space_type::element_type p0_element_type;

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
    BlocHeat( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! shape of the domain
    std::string shape;
}; // BlocHeat

template<int Dim> const uint16_type BlocHeat<Dim>::Order;

template<int Dim>
void
BlocHeat<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute BlocHeat<" << Dim << ">\n";
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
BlocHeat<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    //if ( X[1] == 0 ) shape = "simplex";
    //if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "MCS/heat/%1%/%2%-%3%/P%4%/h_%5%/" ) //chemin où seront générés les résultats
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc = geo(_filename="MCS_heat.geo") );
    //auto mesh = GeoTool::createMeshFromGeoFile<mesh_type>("MCS_heat.geo", "nameExport",meshSize);

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    space_ptrtype Xh = space_type::New( mesh );
    element_type T( Xh, "T" );
    element_type v( Xh, "v" );


    /** \endcode */

    /** define \f$g\f$ the expression of the exact solution and
     * \f$f\f$ the expression of the right hand side such that \f$g\f$
     * is the exact solution
     */
    /** \code */
    //# marker1 #

    //! deduce from expression the type of f (thanks to keyword 'auto')
    auto f = 0;
    //# endmarker1 #
    /** \endcode */

    bool weakdir = this->vm()["weakdir"].template as<int>();
    value_type penaldir = this->vm()["penaldir"].template as<double>();
    value_type nu = this->vm()["nu"].template as<double>();
    using namespace Feel::vf;

    /**
     * Construction of the right hand side. F is the vector that holds
     * the algebraic representation of the right hand side of the
     * problem
     */
    /** \code */
    //# marker2 #
    auto F = M_backend->newVector( Xh );

    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( markedfaces(mesh,11) ,0.1*id(v) );


    //# endmarker2 #
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
        integrate( elements(mesh), nu*gradt(T)*(trans(grad(v)) ))
        + integrate( elements(mesh), idt(T)*id(v) ) ;
    /** \endcode */
    //# endmarker3 #

    /** strong(algebraic) dirichlet conditions treatment for the boundaries marked 1 and 3
     * -# first close the matrix (the matrix must be closed first before any manipulation )
     * -# modify the matrix by cancelling out the rows and columns of D that are associated with the Dirichlet dof
     */
    /** \code */
    //# marker5 #
    D->close();
    form2( Xh, Xh, D ) +=
        on( markedfaces(mesh,10), T, F, cst(1.) );
    //# endmarker5 #
    /** \endcode */

    /** \endcode */

    //! solve the system
    /** \code */
	//# marker6 #
    backend_type::build()->solve( _matrix=D, _solution=T, _rhs=F );
	//# endmarker6 #
    /** \endcode */

    export_ptrtype exporter( export_type::New( this->vm(),
                                               (boost::format( "%1%-%2%-%3%" )
                                                % this->about().appName()
                                                % shape
                                                % Dim).str() ) );
    if ( exporter->doExport() )
    {
        Log() << "exportResults starts\n";

        exporter->step(0)->setMesh( mesh );

        exporter->step(0)->add( "T", T );

        exporter->save();
        Log() << "exportResults done\n";
    }
    /** \endcode */
} // BlocHeat::run

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    Environment env( argc, argv );
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
    app.add( new BlocHeat<2>( app.vm(), app.about() ) );
    /** \endcode */

    /**
     * run the application
     */
    /** \code */
    app.run();
    /** \endcode */
}





