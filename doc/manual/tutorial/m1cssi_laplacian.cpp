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
/**
   \file laplacian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-07-15
 */

#include <feel/feel.hpp>


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
    po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
        ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
        ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
        ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ( "penaldir", Feel::po::value<double>()->default_value( 10 ),
          "penalisation parameter for the weak boundary Dirichlet formulation" )
        ( "exact1D", po::value<std::string>()->default_value( "sin(2*Pi*x)" ), "exact 1D solution" )
        ( "exact2D", po::value<std::string>()->default_value( "sin(2*Pi*x)*cos(2*Pi*y)" ), "exact 2D solution" )
        ( "exact3D", po::value<std::string>()->default_value( "sin(2*Pi*x)*cos(2*Pi*y)*cos(2*Pi*z)" ), "exact 3D solution" )
        ( "rhs1D", po::value<std::string>()->default_value( "" ), "right hand side 1D" )
        ( "rhs2D", po::value<std::string>()->default_value( "" ), "right hand side 2D" )
        ( "rhs3D", po::value<std::string>()->default_value( "" ), "right hand side 3D" )
        ;
    return laplacianoptions.add( Feel::feel_options() );
}


/**
 * \class Laplacian
 *
 * Laplacian Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
template<int Dim, int Order = 1>
class Laplacian
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

    /**
     * Constructor
     */
    Laplacian()
        :
        super(),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
        {
        }

    void run();

private:

    //! mesh characteristic size
    double meshSize;

    //! shape of the domain
    std::string shape;

    mesh_ptrtype mesh;
}; // Laplacian

template<int Dim, int Order>
void
Laplacian<Dim,Order>::run()
{
    LOG(INFO) << "------------------------------------------------------------\n";
    LOG(INFO) << "Execute Laplacian<" << Dim << ">\n";

    Environment::changeRepository( boost::format( "doc/manual/tutorial/%1%/%2%-%3%/P%4%/h_%5%/" )
                                   % this->about().appName()
                                   % shape
                                   % Dim
                                   % Order
                                   % meshSize );

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                                      _usenames=true,
                                                      _shape=shape,
                                                      _h=meshSize,
                                                      _xmin=-1,
                                                      _ymin=-1 ) );


    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    space_ptrtype Xh = space_type::New( mesh );

    // print some information (number of local/global dof in logfile)
    Xh->printInfo();

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    /** \endcode */

    bool weak_dirichlet = this->vm()["weakdir"].template as<int>();
    value_type penaldir = this->vm()["penaldir"].template as<double>();
    value_type nu = this->vm()["nu"].template as<double>();

    auto g = sin(2*pi*Px());
    auto gradg= 2*pi*cos(2*pi*Px());
    auto f = 4*pi*pi*sin(2*pi*Px());
    /** \code */

    auto F = backend()->newVector( Xh );
    auto l = form1( _test=Xh, _vector=F );
    l = integrate( _range=elements( mesh ), _expr=f*id( v ) );

    auto D = backend()->newMatrix( _test=Xh, _trial=Xh  );
    //! assemble $\int_\Omega \nu \nabla u \cdot \nabla v$
    auto a = form2( _test=Xh, _trial=Xh, _matrix=D );
    a = integrate( _range=elements( mesh ), _expr=nu*gradt( u )*trans( grad( v ) ) );

    // apply Dirichlet boundary conditions
    a += on( _range=boundaryfaces(mesh),_element=u, _rhs=F, _expr=g );

    // solve the system
    backend( _rebuild=true  )->solve( _matrix=D, _solution=u, _rhs=F );
    //! compute the \f$L_2$ norm of the error
    double L2error =normL2( _range=elements( mesh ),_expr=( idv( u )-g ) );

    LOG(INFO) << "||error||_L2=" << L2error << "\n";


    //! save the results
    //! project the exact solution
    element_type e( Xh, "e" );
    e = project( _space=Xh, _range=elements( mesh ), _expr=g );

    export_ptrtype exporter( export_type::New() );
    if ( exporter->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( mesh );
        exporter->step( 0 )->add( "solution", u );
        exporter->step( 0 )->add( "exact", e );

        exporter->save();
        LOG(INFO) << "exportResults done\n";
    }

    /** \endcode */
} // Laplacian::run

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
                     _about=about(_name="laplacian",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Application app;


    app.add( new Laplacian<1>() );
    app.run();

}





