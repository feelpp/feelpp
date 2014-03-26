/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake1@e.ujf-grenoble.fr>
       Date: 2011-08-20

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
   \file test_trace.cpp
   \author Abdoulaye Samake <abdoulaye.samake.@imag.fr>
   \date 2011-08-20
 */

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatortrace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

inline
po::options_description
makeOptions()
{
    po::options_description testoptions( "Test options" );
    testoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.05 ), "mesh size" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
    ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
    ( "penaldir", Feel::po::value<double>()->default_value( 10 ),
      "penalisation parameter for the weak boundary Dirichlet formulation" )
    ;
    return testoptions.add( Feel::feel_options() );
}

template<typename T, int Dim, int Order = 1>
struct imesh
{
    typedef Simplex<Dim, Order> convex_type;
    typedef Mesh<convex_type, T > type;
    typedef boost::shared_ptr<type> ptrtype;
};

template<int Dim, int Order>
class Test
    :
public Simget
{
    typedef Simget super;
public:

    typedef double double_type;

    typedef typename imesh<double_type,Dim,Order>::convex_type convex_type;
    typedef typename imesh<double_type,Dim,Order>::type mesh_type;
    typedef typename mesh_type::trace_mesh_type trace_mesh_type;
    typedef FunctionSpace<mesh_type, bases<Lagrange<Order, Scalar> >, double> space_type;
    typedef Backend<double_type> backend_type;
    typedef Exporter<trace_mesh_type,Order> trace_export_type;

    /**
     * Constructor
     */
    Test()
        :
        super(),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
    {}

    void run();

private:

    //! linear algebra backend
    backend_ptrtype M_backend;
    //! mesh characteristic size
    double meshSize;
    //! shape of the domain
    std::string shape;
}; // Test

template<int Dim, int Order>
void
Test<Dim,Order>::run()
{

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute Test<" << Dim << "," << Order << ">\n";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/P%4%G%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                              _usenames=true,
                                              _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                              _shape=shape,
                                              _order=Order,
                                              _dim=Dim,
                                              _h=meshSize ) );

    auto Xh = space_type::New(_mesh=mesh);
    auto u = Xh->element();
    auto v = Xh->element();
    auto gproj = Xh->element();

    double_type pi = M_PI;

    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    gproj = vf::project( _space=Xh, _range=elements( mesh ), _expr=g );

    auto f = pi*pi*Dim*g;

    auto trace_mesh = mesh->trace( boundaryfaces( mesh ) );
    auto Th = Xh->trace( boundaryfaces( mesh ) ) ;
    auto t = vf::project( _space=Th, _range=elements( Th->mesh() ), _expr=g );

    auto op_trace = operatorTrace(Xh);
    //auto g_trace = op_trace->trace( _range=boundaryfaces( mesh ), _expr=g );
    auto g_trace = op_trace->trace( _expr=g );

    double measure = integrate( _range=elements( Th->mesh() ), _expr=cst( 1.0 ), _quad=_Q<5>() ).evaluate()( 0,0 );
    std::cout << " -- measure  =" << measure << "\n";
    double e1 = integrate( _range=elements( Th->mesh() ), _expr=idv( t ), _quad=_Q<15>() ).evaluate()( 0,0 );
    std::cout << " -- |trace|  =" << e1 << "\n";
    double e2 = integrate( _range=elements( Th->mesh() ), _expr=g, _quad=_Q<15>() ).evaluate()( 0,0 );
    std::cout << " -- |g|_bdy  =" << e2 << "\n";
    double e3 = integrate( _range=boundaryfaces( mesh ), _expr=g, _quad=_Q<15>() ).evaluate()( 0,0 );
    std::cout << " -- |g|  =" << e3 << "\n";
    double error = integrate( _range=elements( Th->mesh() ), _expr=idv( t )-idv( g_trace ) ).evaluate()( 0,0 );
    std::cout << " -- |op_trace-trace|  =" << error << "\n";

    std::vector<std::string> bdynames = boost::assign::list_of( "Dirichlet" )( "Neumann" );
    BOOST_FOREACH( auto bdy, bdynames )
    {
        std::cout << "============================================================\n";
        std::cout << "Boundary "  << bdy << "\n";
        auto Lh = Xh->trace( markedfaces( mesh,bdy ) ) ;
        auto l = vf::project( _space=Lh, _range=elements( Lh->mesh() ), _expr=g );


        double measure = integrate( _range=elements( Lh->mesh() ), _expr=cst( 1.0 ), _quad=_Q<5>() ).evaluate()( 0,0 );
        std::cout << " -- measure(" << bdy << ")  =" << measure << "\n";
        double e1 = integrate( _range=elements( Lh->mesh() ), _expr=idv( l ), _quad=_Q<15>() ).evaluate()( 0,0 );
        std::cout << " -- |trace|(" << bdy << ")  =" << e1 << "\n";
        double e2 = integrate( _range=elements( Lh->mesh() ), _expr=g, _quad=_Q<15>() ).evaluate()( 0,0 );
        std::cout << " -- |g|_bdy(" << bdy << ")  =" << e2 << "\n";
        double e3 = integrate( _range=markedfaces( mesh,bdy ), _expr=g, _quad=_Q<15>() ).evaluate()( 0,0 );
        std::cout << " -- |g|(" << bdy << ")  =" << e3 << "\n";

        auto g2 = trans( vf::P() )*vf::P();
        auto l2 = vf::project( _space=Lh, _range=elements( Lh->mesh() ), _expr=g2 );
        double p1 = integrate( _range=elements( Lh->mesh() ), _expr=idv( l2 ), _quad=_Q<15>() ).evaluate()( 0,0 );
        std::cout << " -- proj(x^2+y^2)(" << bdy << ")  =" << p1 << "\n";
        double p2 = integrate( _range=elements( Lh->mesh() ),  _expr=g2, _quad=_Q<15>() ).evaluate()( 0,0 );
        std::cout << " -- x^2+y^2(" << bdy << ")  =" << p2 << "\n";
        double p3 = integrate( _range=markedfaces( Xh->mesh(), bdy ),  _expr=g2, _quad=_Q<15>() ).evaluate()( 0,0 );
        std::cout << " -- 2D(x^2+y^2)(" << bdy << ")  =" << p3 << "\n";

    }

    auto trace_exporter = trace_export_type::New( this->vm(),
                                                  ( boost::format( "trace-%1%-%2%-%3%" )
                                                    % this->about().appName()
                                                    % shape
                                                    % Dim ).str() );


    if ( trace_exporter->doExport() )
    {
        LOG(INFO) << "trace export starts\n";
        trace_exporter->step( 0 )->setMesh( trace_mesh );
        trace_exporter->step( 0 )->add( "trace_g", t );
        trace_exporter->save();
        LOG(INFO) << "trace export done\n";
    }

} // Test::run

int
main( int argc, char** argv )
{

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="test_trace",
                                  _author="Abdoulaye Samake",
                                  _email="abdoulaye.samake@imag.fr") );
    Application app;

    app.add( new Test<2,3>() );
    app.run();

}
