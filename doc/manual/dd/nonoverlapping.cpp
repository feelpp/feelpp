/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Abdoulaye Samake <Abdoulaye.Samake@imag.fr>
   Date: 2012-01-01

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
   \file nonoverlap.cpp
   \author Abdoulaye Samake <Abdoulaye.Samake@imag.fr>
   \date 2012-01-01
*/

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/assign/std/vector.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/aitken.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>


namespace Feel
{
    gmsh_ptrtype nonOverlapGeometryLeft( int Dim, double hsize );
    gmsh_ptrtype nonOverlapGeometryRight( int Dim, double hsize );

using namespace Feel::vf;

inline
po::options_description
makeOptions()
{
    po::options_description relaxationoptions( "relaxation options" );
    relaxationoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.04 ), "mesh size" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
    ( "tol", Feel::po::value<double>()->default_value( 1e-06 ),  " tolerance " )
    ( "imax", Feel::po::value<double>()->default_value( 50 ), " maximum number of iteration" )
    ;
    return relaxationoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "nonoverlapping" ,
                     "nonoverlapping" ,
                     "0.2",
                     "2D, Nonoverlap on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Abdoulaye Samake", "developer", "Abdoulaye.Samake@imag.fr", "" );
    return about;

}

template<int Dim>
class ddmethod
    :
public Simget
{
    typedef Simget super;
public:

    static const uint16_type Order = 2;
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::vector_type vector_type;
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef  boost::numeric::ublas::vector<element_type> Vector_type ;

    /**
     * Constructor
     */
    ddmethod( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() ),
        M_firstExporter( export_type::New( this->vm(),
                                           ( boost::format( "%1%-%2%-%3%" )
                                             % this->about().appName()
                                             % Dim
                                             % int( 1 ) ).str() ) ),
        M_secondExporter( export_type::New( this->vm(),
                                            ( boost::format( "%1%-%2%-%3%" )
                                              % this->about().appName()
                                              % Dim
                                              % int( 2 ) ).str() ) ),
        timers()

    {}

    template<typename DirichletExpr,
             typename RhsExpr,
             typename InterfaceExpr>
    void localProblem( element_type& u,
                       std::vector<int> const& dirichletFlags, DirichletExpr gD,
                       RhsExpr f,
                       std::vector<int> const& interfaceFlags, InterfaceExpr w );

    void exportResults( element_type& u,element_type& v,double time );

    mesh_ptrtype createMesh(  element_type& u );

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:

    backend_ptrtype M_backend;
    double meshSize;
    std::string shape;
    mesh_ptrtype mesh1;
    mesh_ptrtype mesh2;
    export_ptrtype M_firstExporter;
    export_ptrtype M_secondExporter;

    std::map<std::string, std::pair<boost::timer, double> > timers;
    std::map<std::string,double> stats;

    std::vector<int> dirichletFlags1;
    std::vector<int> dirichletFlags2;

    std::vector<int> interfaceFlags1;
    std::vector<int> interfaceFlags2;

}; // nonoverlap

template<int Dim> const uint16_type ddmethod<Dim>::Order;

template<int Dim>
template<typename DirichletExpr,
         typename RhsExpr,
         typename InterfaceExpr>
void
ddmethod<Dim>::localProblem( element_type& u,
                             std::vector<int> const& dirichletFlags, DirichletExpr gD,
                             RhsExpr f,
                             std::vector<int> const& interfaceFlags, InterfaceExpr w )
{

    auto Xh=u.functionSpace();
    auto mesh=Xh->mesh();
    element_type v( Xh,"v" );

    auto B = M_backend->newVector( Xh );

    timers["assembly"].first.restart();

    form1( _test=Xh,_vector=B, _init=true ) =
        integrate( elements( mesh ), f*id( v ) );

    BOOST_FOREACH( int marker, interfaceFlags )
    {
        form1( _test=Xh,_vector=B ) +=
            integrate( markedfaces( mesh, marker ), w*id( v ) );
    }

    B->close();

    timers["assembly"].second = timers["assembly"].first.elapsed();
    timers["assembly_F"].second = timers["assembly"].first.elapsed();

    auto A = M_backend->newMatrix( Xh, Xh );

    timers["assembly"].first.restart();

    form2( _test=Xh, _trial=Xh, _matrix=A, _init=true ) =
        integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) );

    A->close();

    BOOST_FOREACH( int marker, dirichletFlags )
    {
        form2( Xh, Xh, A ) +=
            on( markedfaces( mesh, marker ) ,	u, B, gD );
    }

    timers["assembly"].second += timers["assembly"].first.elapsed();
    timers["assembly_D"].second = timers["assembly"].first.elapsed();

    timers["solver"].first.restart();

    backend_type::build()->solve( _matrix=A, _solution=u, _rhs=B );//, _reuse_prec=true );

    timers["solver"].second = timers["solver"].first.elapsed();

    Log() << "[timer] run():  assembly: " << timers["assembly"].second << "\n";
    Log() << "[timer] run():    o D : " << timers["assembly_D"].second << "\n";
    Log() << "[timer] run():    o F : " << timers["assembly_F"].second << "\n";
    Log() << "[timer] run():  solver: " << timers["solver"].second << "\n";
}

template<int Dim>
typename ddmethod<Dim>::mesh_ptrtype
ddmethod<Dim>::createMesh(  element_type& u )
{
    auto Xh=u.functionSpace();
    mesh_ptrtype mesh =Xh->mesh();
    return mesh;
}

template<int Dim>
void
ddmethod<Dim>::exportResults( element_type& u, element_type& v, double time )
{
    auto Xh1=u.functionSpace();
    auto Xh2=v.functionSpace();

    double pi = M_PI;
    using namespace vf;
    auto g = sin( pi*Px() )*cos( pi*Py() );

    auto proj1 = Xh1->element();
    auto proj2 = Xh2->element();

    proj1 = vf::project( Xh1, elements( createMesh( u ) ), g );
    proj2 = vf::project( Xh2, elements( createMesh( v ) ), g );

    Log() << "exportResults starts\n";
    timers["export"].first.restart();

    M_firstExporter->step( time )->setMesh( createMesh( u ) );
    M_firstExporter->step( time )->add( "solution", ( boost::format( "solution-%1%" ) % int( 1 ) ).str(), u );
    M_firstExporter->step( time )->add( "exact", ( boost::format( "exact-%1%" ) % int( 1 ) ).str(), proj1 );
    M_firstExporter->save();

    M_secondExporter->step( time )->setMesh( createMesh( v ) );
    M_secondExporter->step( time )->add( "solution",( boost::format( "solution-%1%" ) % int( 2 ) ).str(), v );
    M_secondExporter->step( time )->add( "exact",( boost::format( "exact-%1%" ) % int( 2 ) ).str(), proj2 );
    M_secondExporter->save();

    std::ofstream ofs( ( boost::format( "%1%.sos" ) % this->about().appName() ).str().c_str() );

    if ( ofs )
    {
        ofs << "FORMAT:\n"
            << "type: master_server gold\n"
            << "SERVERS\n"
            << "number of servers: " << int( 2 ) << "\n";

        for ( int j = 1; j <= 2; ++ j )
        {
            ofs << "#Server " << j << "\n";
            ofs << "machine id: " << mpi::environment::processor_name()  << "\n";
            ofs << "executable:\n";
            ofs << "data_path: .\n";
            ofs << "casefile: nonoverlapping-" << Dim << "-" << j << "-1_0.case\n";
        }
    }

    Log() << "exportResults done\n";
    timers["export"].second = timers["export"].first.elapsed();
    std::cout << "[timer] exportResults(): " << timers["export"].second << "\n";
} // ddmethod::export

template<int Dim>
void
ddmethod<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute ddmethod<" << Dim << ">\n";
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
ddmethod<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) shape = "simplex";

    if ( X[1] == 1 ) shape = "hypercube";

    value_type tol = this->vm()["tol"].template as<double>();
    value_type imax = this->vm()["imax"].template as<double>();

    Environment::changeRepository( boost::format( "doc/manual/dd/%1%/%2%-%3%/P%4%/h_%5%/" )
                                   % this->about().appName()
                                   %this->shape
                                   % Dim
                                   % Order
                                   %this->meshSize );


    mesh1 = createGMSHMesh( _mesh=new mesh_type,
                            _desc = nonOverlapGeometryLeft( Dim, this->meshSize ) );

    mesh2 = createGMSHMesh( _mesh=new mesh_type,
                            _desc = nonOverlapGeometryRight( Dim, this->meshSize ) );



    if ( Dim == 2 )
    {
        using namespace boost::assign;
        dirichletFlags1+= 1,2,4;
        dirichletFlags2+= 2,3,4;
        interfaceFlags1+= 3;
        interfaceFlags2+= 1;
    }

    else if ( Dim == 3 )
    {
        using namespace boost::assign;
        dirichletFlags1+= 1,2,3,4,5,6;
        dirichletFlags2+= 2,3,4,5,6;
        interfaceFlags1+= 7;
        interfaceFlags2+= 1;
    }


    auto Xh1 = space_type::New( mesh1 );
    auto Xh2 = space_type::New( mesh2 );
    element_type u1( Xh1, "u1" );
    element_type u2( Xh2, "u2" );
    auto uu = Xh1->element();
    auto uv = Xh2->element();

    auto lambda = Xh1->element();
    auto lambdaold = Xh1->element();
    lambda.zero();
    value_type pi = M_PI;
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    auto f = pi*pi*Dim*g;

    auto gradg = trans( pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() )*unitX()
                        -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() )*unitY()
                        -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() )*unitZ() );

    double L2erroru1 = 1.;
    double L2erroru2 = 1.;

    double H1erroru1 = 2.;
    double H1erroru2 = 2.;

    auto Ih12 = opInterpolation( _domainSpace=Xh1, _imageSpace=Xh2, _range=markedfaces( Xh2->mesh(),interfaceFlags2[0] ) );
    auto Ih21 = opInterpolation( _domainSpace=Xh2, _imageSpace=Xh1, _range=markedfaces( Xh1->mesh(),interfaceFlags1[0] ) );


    unsigned int cptExport = 0;

    std::string const fname1( "domain1.dat" );
    std::string const fname2( "domain2.dat" );
    std::ofstream history1( fname1.c_str() );
    std::ofstream history2( fname2.c_str(),std::ios::app );

    while ( ( ( L2erroru1 +L2erroru2 ) > tol ) && ( cptExport < imax ) )
    {
        ++cptExport;

        std::cout << "============================================================\n";
        std::cout << "iteration  : " << cptExport << "\n";
        std::cout << "L2erroru1  : " << L2erroru1  << "\n";
        std::cout << "L2erroru2  : " << L2erroru2  << "\n";
        std::cout << "H1erroru1  : " << H1erroru1  << "\n";
        std::cout << "H1erroru2  : " << H1erroru2  << "\n";

        localProblem( u1,
                      dirichletFlags1, /*dirichlet*/g,
                      /*rhs*/f,
                      interfaceFlags1,idv( lambda ) );

        Ih12->apply( lambda, uv );

        localProblem( u2,
                      dirichletFlags2, g,
                      f,
                      interfaceFlags2,-idv( uv ) );

        Ih21->apply( u2, uu );

        lambda.zero();
        lambda.add( 0.5,uu-u1 );
        lambda.add( 1.0,lambdaold );
        lambdaold = lambda;

        // compute L2error;
        double  L2error2u1 =integrate( elements( mesh1 ), ( idv( u1 )-g )*( idv( u1 )-g ) ).evaluate()( 0,0 );
        double L2error2u2 =integrate( elements( mesh2 ), ( idv( u2 )-g )*( idv( u2 )-g ) ).evaluate()( 0,0 );

        L2erroru1 = math::sqrt( L2error2u1 );
        L2erroru2 = math::sqrt( L2error2u2 );

        // compute H1error;
        double semi_H1error1 =integrate( elements( mesh1 ),
                                         ( gradv( u1 )-gradg )*trans( ( gradv( u1 )-gradg ) ) ).evaluate()( 0,0 );

        double semi_H1error2 =integrate( elements( mesh2 ),
                                         ( gradv( u2 )-gradg )*trans( ( gradv( u2 )-gradg ) ) ).evaluate()( 0,0 );


        H1erroru1 = math::sqrt( L2error2u1 + semi_H1error1 );
        H1erroru2 = math::sqrt( L2error2u2 + semi_H1error2 );

        if ( history1 )
        {
            history1.setf( std::ios::scientific );
            history1 << std::setw( 4 ) << cptExport << "  \t  ";
            history1<< std::setw( 8 ) << std::setprecision( 4 ) << L2erroru1 << "  \t  " << H1erroru1 <<"\n";
        }

        else
        {
            std::cerr << " convergence history filename " << fname1 << " could not be opened " << "\n";
        }

        if ( history2 )
        {
            history2.setf( std::ios::scientific );
            history2 << std::setw( 4 ) << cptExport << "  \t  ";
            history2<< std::setw( 8 ) << std::setprecision( 4 ) << L2erroru2 << "  \t  " << H1erroru2 <<"\n";
        }

        else
        {
            std::cerr << " convergence history filename " << fname2 << " could not be opened " << "\n";
        }

        this->exportResults( u1,u2, cptExport );

    }; // iteration loop

    std::cout << "-------------------------end iteration---------------\n";

    std::cout << "number of iteration  : " << cptExport << "\n";

    std::cout << "L2erroru1  : " << L2erroru1  << "\n";

    std::cout << "L2erroru2  : " << L2erroru2  << "\n";

    std::cout << "H1erroru1  : " << H1erroru1  << "\n";

    std::cout << "H1erroru2  : " << H1erroru2  << "\n";

} // nonoverlap::run
} // Feel

int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( argc, argv );
    Application app( argc, argv, makeAbout(), makeOptions() );

    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    ddmethod<2>  Relax( app.vm(), app.about() );

    Relax.run();
}





