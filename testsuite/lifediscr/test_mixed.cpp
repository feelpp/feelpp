/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-07-07

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_mixed.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-07-07
 */
#include <life/options.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>


#include <life/lifediscr/region.hpp>

#include <life/lifepoly/im.hpp>
#include <life/lifefilters/importergmsh.hpp>

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
    typedef bases<Lagrange<Order,Vectorial> > v_basis_type;
    typedef bases<Lagrange<Order,Scalar> > p_basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, v_basis_type> v_space_type;
    typedef FunctionSpace<mesh_type, p_basis_type> p_space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<v_space_type> v_space_ptrtype;
    typedef boost::shared_ptr<p_space_type> p_space_ptrtype;
    //! an element type of the approximation function space
    typedef typename v_space_type::element_type v_element_type;
    typedef typename p_space_type::element_type p_element_type;

    /**
     * Constructor
     */
    Laplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() )
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

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

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
    v_space_ptrtype Xh = v_space_type::New( mesh );
    p_space_ptrtype Yh = p_space_type::New( mesh );
    v_element_type u( Xh, "u" );
    p_element_type p( Yh, "p" );
    /** \endcode */


    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Yh ) );
    form2( _trial=Xh, _test=Yh, _matrix=D, _init=true)= integrate( elements(mesh), _Q<2>(), divt(u)*id(p) );

    D->close();
    D->printMatlab( "D.m" );
} // Laplacian::run



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






