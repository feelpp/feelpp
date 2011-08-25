/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2011-08-24

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file mortar.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-08-24
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

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
        ("shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
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
                     "nD(n=1,2,3) Laplacian using mortar",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier");

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
    Laplacian( po::variables_map const& vm, AboutData const& about )
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
}; // Laplacian

template<int Dim> const uint16_type Laplacian<Dim>::Order;

template<int Dim>
void
Laplacian<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute Laplacian<" << Dim << ">\n";
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
Laplacian<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
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

    mesh_ptrtype mesh1 = createGMSHMesh( _mesh=new mesh_type,
                                         _desc=domain( _name=(boost::format( "%1%-%2%" ) % shape % Dim).str() ,
                                                       _usenames=true,
                                                       _shape=shape,
                                                       _dim=Dim,
                                                       _h=X[0],
                                                       _xmin=-1,_xmax=0 ) );
    mesh_ptrtype mesh2 = createGMSHMesh( _mesh=new mesh_type,
                                         _desc=domain( _name=(boost::format( "%1%-%2%" ) % shape % Dim).str() ,
                                                       _usenames=true,
                                                       _shape=shape,
                                                       _dim=Dim,
                                                       _h=X[0],
                                                       _xmin=0,_xmax=1) );

    /**
     * The function space and some associated elements(functions) are then defined
     */
    space_ptrtype Xh1 = space_type::New( mesh1 );
    auto u1 = Xh1->element();
    auto v1 = Xh1->element();

    lagmult_space_ptrtype Lh1 = lagmult_space_type::New( mesh1->trace( markedfaces(mesh,"gamma") ));
    auto mu = Lh1->element();
    auto nu = Lh1->element();

    space_ptrtype Xh2 = space_type::New( mesh2 );
    auto u2 = Xh2->element();
    auto v2 = Xh2->element();

    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());
    auto f = pi*pi*Dim*g;

    bool weakdir = this->vm()["weakdir"].template as<int>();
    value_type penaldir = this->vm()["penaldir"].template as<double>();
    value_type nu = this->vm()["nu"].template as<double>();

    auto F1 = M_backend->newVector( Xh1 );
    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements(mesh), f*id(v) )+
        integrate( markedfaces( mesh, "Neumann" ),
				  nu*gradv(gproj)*vf::N()*id(v) );
    form1( _test=Xh, _vector=F ) +=
        integrate( markedfaces(mesh,"Dirichlet"),
                   g*(-grad(v)*vf::N()+penaldir*id(v)/hFace()) );

    auto D = M_backend->newMatrix( Xh, Xh );


    form2( Xh1, Xh1, D1, _init=true ) =
        integrate( elements(mesh), nu*gradt(u1)*trans(grad(v1)) );
    form2( _trial=Xh1, _test=Lh, B1, _init=true ) +=
        integrate( markedfaces(mesh,"gamma"), idt(u1)*id(nu) );

    form2( Xh1, Xh1, D1 ) +=
        integrate( markedfaces(mesh,"Dirichlet"),
                   -(gradt(u1)*vf::N())*id(v1)
                   -(grad(v1)*vf::N())*idt(u1)
                   +penaldir*id(v1)*idt(u1)/hFace());

    form2( Xh2, Xh2, D2, _init=true ) =
        integrate( elements(mesh), nu*gradt(u2)*trans(grad(v2)) );

    form2( _trial=Xh2, _test=Lh, B2, _init=true ) +=
        integrate( markedfaces(mesh,"gamma"), -idt(u2)*id(nu) );


    form2( Xh2, Xh2, D2 ) +=
        integrate( markedfaces(mesh,"Dirichlet"),
                   -(gradt(u2)*vf::N())*id(v2)
                   -(grad(v2)*vf::N())*idt(u2)
                   +penaldir*id(v2)*idt(u2)/hFace());
    D->close();

    backend_type::build()->solve( _matrix=D, _solution=u, _rhs=F );

    double L2error2 =integrate(elements(mesh),
                               (idv(u)-g)*(idv(u)-g) ).evaluate()(0,0);
    double L2error =   math::sqrt( L2error2 );


    Log() << "||error||_L2=" << L2error << "\n";

    element_type e( Xh, "e" );
    e = vf::project( Xh, elements(mesh), g );

    export_ptrtype exporter( export_type::New( this->vm(),
                                               (boost::format( "%1%-%2%-%3%" )
                                                % this->about().appName()
                                                % shape
                                                % Dim).str() ) );
    if ( exporter->doExport() )
    {
        Log() << "exportResults starts\n";

        exporter->step(0)->setMesh( mesh );

        exporter->step(0)->add( "u", u );
        exporter->step(0)->add( "g", e );

        exporter->save();
        Log() << "exportResults done\n";
    }
    /** \endcode */
} // Laplacian::run

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
    app.add( new Laplacian<1>( app.vm(), app.about() ) );
    app.add( new Laplacian<2>( app.vm(), app.about() ) );
    app.add( new Laplacian<3>( app.vm(), app.about() ) );
    /** \endcode */

    /**
     * run the application
     */
    /** \code */
    app.run();
    /** \endcode */
}






