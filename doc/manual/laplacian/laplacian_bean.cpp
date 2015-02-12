/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2015-02-11

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
   \file laplacian_bean.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2015-02-11
 */

#include <feel/feel.hpp>
#include <feel/feelfilters/straightenmesh_impl.hpp>


/** use Feel namespace */
using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
        ( "geofile", Feel::po::value<std::string>()->default_value( "bean" ), "name of the geofile input")
        ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
        ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ( "penaldir", Feel::po::value<double>()->default_value( 10 ),
          "penalisation parameter for the weak boundary Dirichlet formulation" )
        ( "exact", po::value<std::string>()->default_value( "sin(2*Pi*x)*cos(2*Pi*y)" ), "exact solution" )
        ( "rhs", po::value<std::string>()->default_value( "" ), "right hand side" )
        ;
    return laplacianoptions.add( backend_options("LaplacianG1") ).add( backend_options("LaplacianG2") );
}


/**
 * \class Laplacian
 *
 * Laplacian Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=2 or 3)
 */
template<int Dim, int POrder = 2, int GOrder = 1>
class Laplacian
    :
public Simget
{
public:

    void run();

}; // Laplacian

template<int Dim, int POrder, int GOrder>
void
Laplacian<Dim,POrder,GOrder>::run()
{
    LOG(INFO) << "------------------------------------------------------------\n";
    LOG(INFO) << "Execute Laplacian<" << Dim << ">\n";

    Environment::changeRepository( boost::format( "doc/manual/laplacian_bean/%1%/D%2%/P%3%/G%4%/h_%5%/" )
                                   % this->about().appName()
                                   % Dim
                                   % POrder
                                   % GOrder
                                   % doption("gmsh.hsize") );

    auto mesh = loadMesh( _mesh=new Mesh<Hypercube<Dim,GOrder,Dim>>, _filename=soption("geofile") );
    auto Xh = Pch<POrder>( mesh );
    Xh->printInfo();

    auto u = Xh->element();
    auto v = Xh->element();
    auto gproj = Xh->element();

    bool weak_dirichlet = ioption( _name="weakdir" );
    double penaldir = doption(_name="penaldir");
    double nu = doption(_name="nu");
    auto exact = expr(soption(_name="exact"));
    std::string rhs_str = soption(_name="rhs");
    auto f = -nu*laplacian(exact);
    auto gradg = grad<Dim,2>(exact);
    if ( Environment::isMasterRank() )
    {
        std::cout << "Looking for u such that -Delta u = f on Omega\n"
                  << "  with \n"
                  << "    - u = g on Dirichlet and \n";

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

    if ( weak_dirichlet )
    {
        l += integrate( _range=boundaryfaces( mesh ),
                        _expr=nu*exact*( -grad( v )*N()
                                         + penaldir*id( v )/hFace() ) );
    }

    auto a = form2( _test=Xh, _trial=Xh );
    a = integrate( _range=elements( mesh ), _expr=nu*gradt( u )*trans( grad( v ) ) );

    if ( weak_dirichlet )
    {
        a += integrate( _range=boundaryfaces( mesh ),
                        _expr= nu * ( -( gradt( u )*N() )*id( v )
                                      -( grad( v )*N() )*idt( u )
                                      +penaldir*id( v )*idt( u )/hFace() ) );
    }
    else
    {
        a += on( _range=boundaryfaces( mesh ),
                 _element=u, _rhs=l, _expr=exact );
    }


    a.solve( _rhs=l, _solution=u, _name=(boost::format("LaplacianG%1%")%GOrder).str() );


    //! compute the \f$L_2$ norm of the error
    double L2error =normL2( _range=elements( mesh ),_expr=( idv( u )-exact ) );

    LOG(INFO) << "||error||_L2=" << L2error << "\n";
    if ( Environment::isMasterRank() )
    {
        std::cout << "||error||_L2=" << L2error << "\n";
    }

    //! project the exact solution
    v = project( _space=Xh, _range=elements( mesh ), _expr=exact );

    auto e = exporter( _mesh = mesh );
    e->add( "solution", u );
    e->add( "exact", v );
    e->save();

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
                     _about=about(_name="laplacian_bean",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org") );

    Application app;
    app.add( new Laplacian<2,2,1>() );
    app.add( new Laplacian<2,2,2>() );
    app.run();

}
