/* -*- mode: c++ coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-07

  Copyright (C) 2008-2009 Universite Joseph Fourier (Grenoble I)

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
   \file laplacian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-07
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

/** include gmsh generator for tensorized domains */
#include <life/lifefilters/gmshtensorizeddomain.hpp>

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
    po::options_description laplacianoptions("Laplacian options");
    laplacianoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ("nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient")
        ("weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ("penaldir", Life::po::value<double>()->default_value( 10 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ;
    return laplacianoptions.add( Life::life_options() );
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
    AboutData about( "laplacian" ,
                     "laplacian" ,
                     "0.2",
                     "nD(n=1,2,3) Laplacian on simplices or simplex products",
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


/**
 * \class Laplacian
 *
 * Laplacian Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
template<int Dim>
class Laplacian
    :
    public Application
{
    typedef Application super;
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
    typedef FunctionSpace<mesh_type, bases<Lagrange<0,Scalar> >, Discontinuous > p0_space_type;
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
    Laplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
    {
    }

    /**
     * create the mesh using mesh size \c meshSize
     *
     * \param meshSize the mesh characteristic size
     * \return a mesh of an hyper cube of dimension Dim
     */
    mesh_ptrtype createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve the system \f$D u = F\f$
     *
     * \param in D sparse matrix
     * \param inout u solution of the system
     * \param in F vector representing the right hand side of the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     *
     * \param u the function to save in ensight format.
     * \param e the exact solution projected in the finite element space to save in ensight format.
     */
    void exportResults( element_type& u, element_type& e );

private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! exporter factory
    export_ptrtype exporter;
}; // Laplacian

template<int Dim> const uint16_type Laplacian<Dim>::Order;

template<int Dim>
typename Laplacian<Dim>::mesh_ptrtype
Laplacian<Dim>::createMesh( double meshSize )
{
    /** instantiate a new mesh */
    /** \code */
    mesh_ptrtype mesh( new mesh_type );
    /** \endcode */

    //! generate a tensorized domain (hyper cube of dimension Dim)
    /** \code */
    GmshTensorizedDomain<convex_type::nDim,
                         convex_type::nOrder,
                         convex_type::nDim,
                         Simplex> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( -1, 1 ) );
    td.setY( std::make_pair( -1, 1 ) );
    std::string fname = td.generate( convex_type::name().c_str() );
    /** \endcode */

    //! importer the mesh generated by td
    /** \code */
    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    /** \endcode */

    return mesh;
} // Laplacian::createMesh


template<int Dim>
void
Laplacian<Dim>::run()
{
    /**
     * print help if --help is passed to the command line
     */
    /** \code */
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }
    /** \endcode */

    /**
     * we change to the directory where the results and logs will be
     * stored
     */
    /** \code */
    this->changeRepository( boost::format( "doc/tutorial/%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % convex_type::name()
                            % Order
                            % this->vm()["hsize"].template as<double>()
                            );
    /** \endcode */
    /**
     * First we create the mesh
     */
    /** \code */
    mesh_ptrtype mesh = createMesh( meshSize );
    /** \endcode */

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    space_ptrtype Xh = space_type::New( mesh );
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    /** \endcode */

    /** define \f$g\f$ the expression of the exact solution and
     * \f$f\f$ the expression of the right hand side such that \f$g\f$
     * is the exact solution
     */
    /** \code */
    //# marker1 #
    value_type pi = M_PI;
    //! deduce from expression the type of g (thanks to keyword 'auto')
    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());
    //! deduce from expression the type of f (thanks to keyword 'auto')
    auto f = pi*pi*Dim*g;
    //# endmarker1 #
    /** \endcode */

    bool weakdir = this->vm()["weakdir"].template as<int>();
    value_type penaldir = this->vm()["penaldir"].template as<double>();


    /**
     * Construction of the right hand side. F is the vector that holds
     * the algebraic representation of the right habd side of the
     * problem
     */
    /** \code */
    //# marker2 #
    vector_ptrtype F( M_backend->newVector( Xh ) );

    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements(mesh), f*id(v) );
    //# endmarker2 #
    if ( Application::nProcess () != 1 || weakdir )
        {
            //# marker41 #
            form1( _test=Xh, _vector=F ) +=
                integrate( markedfaces(mesh,1), g*(-grad(v)*N()+penaldir*id(v)/hFace()) ) +
                integrate( markedfaces(mesh,3), g*(-grad(v)*N()+penaldir*id(v)/hFace()) );
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

    value_type nu = this->vm()["nu"].template as<double>();

    //! assemble $\int_\Omega \nu \nabla u \cdot \nabla v$
    /** \code */
    form2( Xh, Xh, D, _init=true ) =
        integrate( elements(mesh), nu*gradt(u)*trans(grad(v)) );
    /** \endcode */
    //# endmarker3 #

    if ( Application::nProcess () != 1 || weakdir )
        {
            /** weak dirichlet conditions treatment for the boundaries marked 1 and 3
             * -# assemble \f$\int_{\partial \Omega} -\nabla u \cdot \mathbf{n} v\f$
             * -# assemble \f$\int_{\partial \Omega} -\nabla v \cdot \mathbf{n} u\f$
             * -# assemble \f$\int_{\partial \Omega} \frac{\gamma}{h} u v\f$
             */
            /** \code */
            //# marker10 #
            form2( Xh, Xh, D ) +=
                integrate( markedfaces(mesh,1),

                           -(gradt(u)*N())*id(v)

                           -(grad(v)*N())*idt(u)
                           +penaldir*id(v)*idt(u)/hFace()) +
                integrate( markedfaces(mesh,3),
                           -(gradt(u)*N())*id(v)
                           -(grad(v)*N())*idt(u)
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
                on( markedfaces(mesh, 1), u, F, g )+
                on( markedfaces(mesh, 3), u, F, g );
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
    double L2error2 =integrate(elements(mesh),
                               (idv(u)-g)*(idv(u)-g) ).evaluate()(0,0);
    double L2error =   math::sqrt( L2error2 );


    Log() << "||error||_L2=" << L2error << "\n";
    //# endmarker7 #
    /** \endcode */

    //! save the results
    /** \code */
    //! project the exact solution
    element_type e( Xh, "e" );
    e = vf::project( Xh, elements(mesh), g );
    //std::cout << "e=" << e << "\n";
    this->exportResults( u, e );
    /** \endcode */
} // Laplacian::run

//# marker6 #
template<int Dim>
void
Laplacian<Dim>::solve( sparse_matrix_ptrtype& D,
                       element_type& u,
                       vector_ptrtype& F )
{
    //! solve the system, first create a vector U of the same size as
    //! u, then call solve,
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    //! call solve, the second D is the matrix which will be used to
    //! create the preconditionner
    M_backend->solve( D, D, U, F );
    //! copy U in u
    u = *U;
} // Laplacian::solve
//# endmarker6 #

template<int Dim>
void
Laplacian<Dim>::exportResults( element_type& U, element_type& E )
{
    if ( exporter->doExport() )
    {
        Log() << "exportResults starts\n";

        exporter->step(0)->setMesh( U.functionSpace()->mesh() );

        exporter->step(0)->add( "pid",regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );
        exporter->step(0)->add( "u", U );
        exporter->step(0)->add( "g", E );

        exporter->save();
    }
} // Laplacian::export



/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    /**
     * intantiate a Laplacian<Dim> class with Dim=2 (e.g. geometric dimension is 2)
     */
    /** \code */
    Laplacian<2> laplacian( argc, argv, makeAbout(), makeOptions() );
    /** \encode */

    /**
     * run the application
     */
    /** \code */
    laplacian.run();
    /** \endcode */
}





