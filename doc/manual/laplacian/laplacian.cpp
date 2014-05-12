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
    return laplacianoptions.add( backend_options( "Laplacian1D") ).add( backend_options( "Laplacian2D") ).add( backend_options( "Laplacian3D") );
}


/**
 * \class Laplacian
 *
 * Laplacian Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
template<int Dim, int Order = 2>
class Laplacian
    :
public Simget
{
public:

    void run();

}; // Laplacian

template<int Dim, int Order>
void
Laplacian<Dim,Order>::run()
{
    LOG(INFO) << "------------------------------------------------------------\n";
    LOG(INFO) << "Execute Laplacian<" << Dim << ">\n";

    Environment::changeRepository( boost::format( "doc/manual/laplacian/%1%/D%2%/P%3%/h_%4%/" )
                                   % this->about().appName()
                                   % Dim
                                   % Order
                                   % doption("gmsh.hsize") );

    auto mesh = loadMesh( new Mesh<Simplex<Dim>> );
    auto Xh = Pch<Order>( mesh );
    Xh->printInfo();

    auto u = Xh->element();
    auto v = Xh->element();
    auto gproj = Xh->element();

    bool weak_dirichlet = ioption( _name="weakdir" );
    double penaldir = doption(_name="penaldir");
    double nu = doption(_name="nu");
    auto exact = expr(soption(_name=(boost::format("exact%1%D")%Dim).str()));
    std::string rhs_str = soption(_name=(boost::format("rhs%1%D")%Dim).str());
    auto f = -nu*laplacian(exact);
    auto gradg = grad<Dim,2>(exact);
    if ( Environment::isMasterRank() )
    {
        std::cout << "Looking for u such that -Delta u = f on Omega\n"
                  << "  with \n"
                  << "    - u = g on Dirichlet and \n"
                  << "    - \\partial_n u = \\partial_n g on Neumann\n";
        //std::cout << "      f = " << f << std::endl;
        std::cout << "    rhs = " << rhs_str << std::endl;
        std::cout << "      g = " << exact << std::endl;
        std::cout << " grad g = " << gradg << std::endl;
    }

    // build Xh-interpolant of g
    gproj = project( _space=Xh, _range=elements( mesh ), _expr=exact );


    auto l = form1( _test=Xh );
    if ( rhs_str.empty() )
        l = integrate( _range=elements( mesh ), _expr=f*id( v ) );
    else
        l = integrate( _range=elements( mesh ), _expr=expr(rhs_str)*id( v ) );
    l += integrate( _range=markedfaces( mesh, "Neumann" ),
                    _expr=nu*gradg*N()*id( v ) );

    //# endmarker2 #
    if ( weak_dirichlet )
    {
        //# marker41 #
        l += integrate( _range=markedfaces( mesh,"Dirichlet" ),
                          _expr=nu*exact*( -grad( v )*N()
                                      + penaldir*id( v )/hFace() ) );
        //# endmarker41 #
    }

    /** \endcode */

    //! assemble $\int_\Omega \nu \nabla u \cdot \nabla v$
    /** \code */
    auto a = form2( _test=Xh, _trial=Xh );
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
        a += integrate( _range=markedfaces( mesh,"Dirichlet" ),
                        _expr= nu * ( -( gradt( u )*N() )*id( v )
                                      -( grad( v )*N() )*idt( u )
                                      +penaldir*id( v )*idt( u )/hFace() ) );
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
        a += on( _range=markedfaces( mesh, "Dirichlet" ),
                 _element=u, _rhs=l, _expr=exact );
        //# endmarker5 #
        /** \endcode */

    }


    a.solve( _rhs=l, _solution=u, _name=(boost::format("Laplacian%1%D")%Dim).str() );


    //! compute the \f$L_2$ norm of the error
    double L2error =normL2( _range=elements( mesh ),_expr=( idv( u )-exact ) );

    LOG(INFO) << "||error||_L2=" << L2error << "\n";
    if ( Environment::isMasterRank() )
    {
        std::cout << "||error||_L2=" << L2error << "\n";
    }

    //! save the results
    /** \code */
    //! project the exact solution
    v = project( _space=Xh, _range=elements( mesh ), _expr=exact );

    auto e = exporter( _mesh = mesh );
    e->add( "solution", u );
    e->add( "exact", v );
    e->save();

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


    if ( app.nProcess() == 1 )
        app.add( new Laplacian<1>() );
    app.add( new Laplacian<2>() );
    app.add( new Laplacian<3>() );

    app.run();

}
