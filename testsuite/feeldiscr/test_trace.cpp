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
   \author Abdoulaye Samake <abdoulaye.samake1.prudhomme@e.ujf-grenoble.fr>
   \date 2011-08-20
 */

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatortrace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;
using namespace Feel::vf;

inline
po::options_description
makeOptions()
{
    po::options_description testoptions("Test options");
    testoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.05 ), "mesh size")
        ("shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        ("nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient")
        ("weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ("penaldir", Feel::po::value<double>()->default_value( 10 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ;
    return testoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_trace" ,
                     "test_trace" ,
                     "0.2",
                     "nD(n=1,2,3) Test Trace operations on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

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
    typedef typename imesh<double_type,Dim,Order>::ptrtype mesh_ptrtype;
    typedef typename mesh_type::trace_mesh_type trace_mesh_type;
    typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;


    typedef FunctionSpace<mesh_type, bases<Lagrange<Order, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef Backend<double_type> backend_type;

    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef Exporter<mesh_type,Order> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef Exporter<trace_mesh_type,Order> trace_export_type;
    typedef boost::shared_ptr<trace_export_type> trace_export_ptrtype;




    /**
     * Constructor
     */
    Test( po::variables_map const& vm, AboutData const& about )
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
}; // Test

template<int Dim,int Order>
void
Test<Dim,Order>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute Test<" << Dim << "," << Order << ">\n";
    std::vector<double> X( 2 );
    X[0] = meshSize;
    if ( shape == "hypercube" )
        X[1] = 1;
    else // default is simplex
        X[1] = 0;
    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<int Dim, int Order>
void
Test<Dim,Order>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) shape = "simplex";
    if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/P%4%G%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=(boost::format( "%1%-%2%" ) % shape % Dim).str() ,
                                                      _usenames=true,
                                                      _convex=(convex_type::is_hypercube)?"Hypercube":"Simplex",
                                                      _shape=shape,
                                                      _order=Order,
                                                      _dim=Dim,
                                                      _h=X[0] ) );

    space_ptrtype Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();
    auto gproj = Xh->element();

    double_type pi = M_PI;

    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());
    gproj = vf::project( Xh, elements(mesh), g );

    auto f = pi*pi*Dim*g;
#if 0
    bool weakdir = this->vm()["weakdir"].template as<int>();
    double_type penaldir = this->vm()["penaldir"].template as<double>();
    double_type nu = this->vm()["nu"].template as<double>();

    using namespace Feel::vf;

    auto F = M_backend->newVector( Xh );
    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( elements(mesh), f*id(v) )+
        integrate( markedfaces( mesh, "Neumann" ),
                   nu*gradv(gproj)*vf::N()*id(v) );

    if ( this->comm().size() != 1 || weakdir )
    {

        form1( _test=Xh, _vector=F ) +=
            integrate( markedfaces(mesh,"Dirichlet"),
                       g*(-grad(v)*vf::N()+penaldir*id(v)/hFace()) );

    }
    F->close();


    auto D = M_backend->newMatrix( Xh, Xh );

    form2( Xh, Xh, D, _init=true ) =
        integrate( elements(mesh), nu*gradt(u)*trans(grad(v)) );

    if ( this->comm().size() != 1 || weakdir )
        {

            form2( Xh, Xh, D ) +=
                integrate( markedfaces(mesh,"Dirichlet"),
                           -(gradt(u)*vf::N())*id(v)
                           -(grad(v)*vf::N())*idt(u)
                           +penaldir*id(v)*idt(u)/hFace());
            D->close();

        }
    else
        {

            D->close();
            form2( Xh, Xh, D ) +=
                on( markedfaces(mesh, "Dirichlet"), u, F, g );

        }

    backend_type::build()->solve( _matrix=D, _solution=u, _rhs=F );

    double L2error2 =integrate(elements(mesh),
                               (idv(u)-g)*(idv(u)-g) ).evaluate()(0,0);
    double L2error =   math::sqrt( L2error2 );


    Log() << "||error||_L2=" << L2error << "\n";

    auto e = Xh->element();
    e = vf::project( Xh, elements(mesh), g );


    export_ptrtype exporter( export_type::New( this->vm(),
                                               (boost::format( "%1%-%2%-%3%" )
                                                % this->about().appName()
                                                % shape
                                                % Dim).str() ) );

    if ( exporter->doExport() )
    {
        Log() << "export starts\n";

        exporter->step(0)->setMesh( mesh );

        exporter->step(0)->add( "u", u );
        exporter->step(0)->add( "g", e );

        exporter->save();
        Log() << "export done\n";
    }
#endif

    auto trace_mesh = mesh->trace( boundaryfaces(mesh) );
    auto Th = Xh->trace( boundaryfaces(mesh)) ;
    auto t = vf::project( Th, elements( Th->mesh() ), g );

    auto op_trace = operatorTrace( Th, M_backend );

    auto g_trace = op_trace->trace( _range=elements( Th->mesh() ), _expr=g);

    double error = integrate( elements( Th->mesh() ), idv(t)-idv(g_trace) ).evaluate()(0,0);

    std::cout << "error  =" << error << "\n";




    trace_export_ptrtype trace_exporter( trace_export_type::New( this->vm(),
                                                                 (boost::format( "trace-%1%-%2%-%3%" )
                                                                  % this->about().appName()
                                                                  % shape
                                                                  % Dim).str() ) );

    if ( trace_exporter->doExport() )
    {
        Log() << "trace export starts\n";

        trace_exporter->step(0)->setMesh( trace_mesh );
        trace_exporter->step(0)->add( "trace_g", t );
        trace_exporter->save();
        Log() << "trace export done\n";
    }

} // Test::run

int
main( int argc, char** argv )
{
    Environment env( argc, argv );

    Application app( argc, argv, makeAbout(), makeOptions() );
    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    // app.add( new Test<1>( app.vm(), app.about() ) );
    //app.add( new Test<2,1>( app.vm(), app.about() ) );
    //app.add( new Test<2,2>( app.vm(), app.about() ) );
    app.add( new Test<2,3>( app.vm(), app.about() ) );
    // app.add( new Test<3>( app.vm(), app.about() ) );

    app.run();

}





