/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2012 Feel++ Consortium

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
/** include predefined feel command line options */
#include <feel/feel.hpp>

#include <boost/tokenizer.hpp>

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
    po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
    ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
    ( "penaldir", Feel::po::value<double>()->default_value( 10 ),
      "penalisation parameter for the weak boundary Dirichlet formulation" )
        ( "filename", po::value<std::string>()->default_value( ""  ), "mesh filename to load (generate mesh if empty)" )
        ( "Neumann", po::value<std::string>()->default_value( ""  ), "Neumann boundaries list" )
        ( "Dirichlet", po::value<std::string>()->default_value( ""  ), "Dirichlet boundaries list" )
    ;
    return laplacianoptions;
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
}; // Laplacian

template<int Dim> const uint16_type Laplacian<Dim>::Order;

template<int Dim>
void
Laplacian<Dim>::run()
{
    LOG(INFO) << "------------------------------------------------------------\n";
    LOG(INFO) << "Execute Laplacian<" << Dim << ">\n";

    Environment::changeRepository( boost::format( "applications/check/%1%/%2%-%3%/P%4%/h_%5%/" )
                                   % this->about().appName()
                                   % shape
                                   % Dim
                                   % Order
                                   % meshSize );

    mesh_ptrtype mesh;
    if ( this->vm().count( "filename" ) )
    {
        LOG(INFO) << "Loading " << this->vm()["filename"].template as<std::string>() << "..." << std::endl;
        mesh = loadGMSHMesh( _mesh=new mesh_type,
                             _filename=this->vm()["filename"].template as<std::string>(),
                             _rebuild_partitions=(Environment::worldComm().size() > 1),
                             _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    }
    else
        mesh = createGMSHMesh( _mesh=new mesh_type,
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
    element_type gproj( Xh, "v" );
    /** \endcode */

    /** define \f$g\f$ the expression of the exact solution and
     * \f$f\f$ the expression of the right hand side such that \f$g\f$
     * is the exact solution
     */
    /** \code */
    //# marker1 #
    value_type pi = 2*M_PI;
    //! deduce from expression the type of g (thanks to keyword 'auto')
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    gproj = vf::project( Xh, elements( mesh ), g );

    //! deduce from expression the type of f (thanks to keyword 'auto')
    auto f = pi*pi*Dim*g;
    //# endmarker1 #
    /** \endcode */

    bool weak_dirichlet = this->vm()["weakdir"].template as<int>();
    value_type penaldir = this->vm()["penaldir"].template as<double>();
    value_type nu = this->vm()["nu"].template as<double>();

    using namespace Feel::vf;

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",-;|");
    tokenizer neumann_tokens(this->vm()["Neumann"].template as<std::string>(), sep);
    std::for_each( neumann_tokens.begin(), neumann_tokens.end(),[]( std::string const& s ) { LOG(INFO) << "token: " << s << std::endl; } );
    tokenizer dirichlet_tokens(this->vm()["Dirichlet"].template as<std::string>(), sep);
    std::for_each( dirichlet_tokens.begin(), dirichlet_tokens.end(),[]( std::string const& s ) { LOG(INFO) << "token: " << s << std::endl; } );

    /**
     * Construction of the right hand side. F is the vector that holds
     * the algebraic representation of the right habd side of the
     * problem
     */
    /** \code */
    //# marker2 #
    auto F = backend( _vm=this->vm() )->newVector( Xh );
    auto rhs = form1( _test=Xh, _vector=F, _init=true );
    rhs = integrate( _range=elements( mesh ), _expr=f*id( v ) );

    std::for_each( neumann_tokens.begin(), neumann_tokens.end(),
                  [&]( std::string const& mark )
                  {
                      LOG(INFO) << "Integrating on Neumann boundary '" << mark << "'\n" ;
                      rhs += integrate( _range=markedfaces( mesh, mark ),
                                        _expr=nu*gradv( gproj )*vf::N()*id( v ) );
                  } );
    //# endmarker2 #
    if ( weak_dirichlet )
    {
        //# marker41 #
        std::for_each( dirichlet_tokens.begin(), dirichlet_tokens.end(),
                      [&]( std::string const& mark )
                      {
                          LOG(INFO) << "Integrating on Dirichlet boundary '" << mark << "'\n" ;
                          rhs += integrate( _range=markedfaces( mesh,mark ),
                                            _expr=g*( -grad( v )*vf::N()+penaldir*id( v )/hFace() ) );
                      } );
        //# endmarker41 #
    }

    /** \endcode */

    /**
     * create the matrix that will hold the algebraic representation
     * of the left hand side
     */
    //# marker3 #
    /** \code */
    auto D = backend()->newMatrix( _test=Xh, _trial=Xh  );
    /** \endcode */

    //! assemble $\int_\Omega \nu \nabla u \cdot \nabla v$
    /** \code */
    auto a = form2( _test=Xh, _trial=Xh, _matrix=D );
    a = integrate( _range=elements( mesh ), _expr=nu*gradt( u )*trans( grad( v ) ) );
    /** \endcode */
    //# endmarker3 #

    if ( weak_dirichlet )
    {
        /** weak dirichlet conditions treatment for the boundaries marked 1 and 3
         * -# assemble \f$\int_{\partial \Omega} -\nabla u \cdot \mathbf{n} v\f$
         * -# assemble \f$\int_{\partial \Omega} -\nabla v \cdot \mathbf{n} u\f$
         * -# assemble \f$\int_{\partial \Omega} \frac{\gamma}{h} u v\f$
         */
        /** \code */
        //# marker10 #
        std::for_each( dirichlet_tokens.begin(), dirichlet_tokens.end(),
                      [&]( std::string const& mark )
                      {
                          LOG(INFO) << "Integrating on Dirichlet boundary '" << mark << "'\n" ;
                          a += integrate( _range=markedfaces( mesh,mark ),
                                          _expr= ( -( gradt( u )*vf::N() )*id( v )
                                                   -( grad( v )*vf::N() )*idt( u )
                                                   +penaldir*id( v )*idt( u )/hFace() ) );
                      } );
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
        std::for_each( dirichlet_tokens.begin(), dirichlet_tokens.end(),
                       [&]( std::string const& mark )
                       {
                           LOG(INFO) << "Integrating on Dirichlet boundary '" << mark << "'\n" ;
                           a += on( _range=markedfaces( mesh, mark ),
                                    _element=u, _rhs=F, _expr=g );
                       } );
        //# endmarker5 #
        /** \endcode */

    }

    /** \endcode */

    //! solve the system
    /** \code */
    //# marker6 #
    backend( _rebuild=true )->solve( _matrix=D, _solution=u, _rhs=F );
    //# endmarker6 #
    /** \endcode */

    //! compute the \f$L_2$ norm of the error
    /** \code */
    //# marker7 #
    double L2error = normL2( _range=elements( mesh ), _expr=( idv( u )-g ) );

    LOG(INFO) << "||error||_L2=" << L2error << "\n";

    std::ofstream res(this->vm()["result-file"].template as<std::string>() );
    res << "L2="<< L2error << "\n";

    //# endmarker7 #
    /** \endcode */

    //! save the results
    /** \code */
    //! project the exact solution
    element_type e( Xh, "e" );
    e = vf::project( Xh, elements( mesh ), g );

    export_ptrtype exporter( export_type::New( this->vm(),
                             ( boost::format( "%1%-%2%-%3%" )
                               % this->about().appName()
                               % shape
                               % Dim ).str() ) );

    if ( exporter->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( mesh );
        exporter->step( 0 )->addRegions();
        exporter->step( 0 )->add( "exact_solution", e );
        exporter->step( 0 )->add( "solution", u );

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
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="check",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Application app;
    app.add( new Laplacian<3>() );
    app.run();

}





