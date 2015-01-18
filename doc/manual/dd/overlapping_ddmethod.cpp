/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Abdoulaye Samake <abdoulaye.samake@imag.fr>
   Date: 2011-11-14

   Copyright (C) 2008-2011 Universite Joseph Fourier (Grenoble I)

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
   \file ddmethod.cpp
   \author Abdoulaye Samake <abdoulaye.samake@imag.fr>
   \date 2011-11-14
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
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/geotool.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>

namespace Feel
{
gmsh_ptrtype ddmethodGeometryLeft( int RDim, double hsize );
gmsh_ptrtype ddmethodGeometryRight( int RDim, double hsize );
using namespace Feel::vf;
inline
po::options_description
makeOptions()
{
    po::options_description relaxationoptions( "relaxation options" );
    relaxationoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.04 ), "mesh size" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "maxIterations", po::value<double>()->default_value( 10 ), "maximal number of iterations" )
    ( "additive", po::value<int>()->default_value( 0 ), "use additive method" )
    ( "tolerance", Feel::po::value<double>()->default_value( 1e-08 ),  " tolerance " )
    ;
    return relaxationoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "overlapping_ddmethod" ,
                     "overlapping_ddmethod" ,
                     "0.2",
                     "nD(n=2,3) ddmethod on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Abdoulaye Samake", "developer", "abdoulaye.samake@imag.fr", "" );
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
    /**
     * Constructor
     */
    ddmethod()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        shape( soption("shape") ),
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

    double l2Error( element_type& u );
    double h1Error( element_type& u );
    void exportResults( element_type& u,element_type& v,double time );
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

}; // ddmethod

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
    auto v=Xh->element();
    auto B = M_backend->newVector( Xh );
    timers["assembly"].first.restart();
    form1( _test=Xh,_vector=B, _init=true ) =
        integrate( elements( mesh ), f*id( v ) );
    B->close();

    timers["assembly"].second = timers["assembly"].first.elapsed();
    timers["assembly_B"].second = timers["assembly"].first.elapsed();

    auto A = M_backend->newMatrix( Xh, Xh );
    timers["assembly"].first.restart();
    form2( _test=Xh, _trial=Xh, _matrix=A, _init=true ) =
        integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) );
    A->close();
    BOOST_FOREACH( int marker, dirichletFlags )
    {
        // std::cout << "apply strong dirichlet on   " << marker << std::endl;
        form2( Xh, Xh, _matrix=A ) +=
            on( markedfaces( mesh, marker ) ,	u, B, gD );
    }
    BOOST_FOREACH( int marker, interfaceFlags )
    {
        // std::cout << "apply interface condition on   " << marker << std::endl;
        form2( Xh, Xh, _matrix=A ) +=
            on( markedfaces( mesh, marker ) ,	u, B, w );
    }
    timers["assembly"].second += timers["assembly"].first.elapsed();
    timers["assembly_A"].second = timers["assembly"].first.elapsed();

    timers["solver"].first.restart();
    backend_type::build(soption("backend"))->solve( _matrix=A, _solution=u, _rhs=B );//, _reuse_prec=true );
    timers["solver"].second = timers["solver"].first.elapsed();

    LOG(INFO) << "[timer] run():  assembly: " << timers["assembly"].second << "\n";
    LOG(INFO) << "[timer] run():    o A : " << timers["assembly_A"].second << "\n";
    LOG(INFO) << "[timer] run():    o B : " << timers["assembly_B"].second << "\n";
    LOG(INFO) << "[timer] run():  solver: " << timers["solver"].second << "\n";
}

template<int Dim>
double
ddmethod<Dim>::l2Error( element_type& u )
{
    auto Xh=u.functionSpace();
    auto mesh=Xh->mesh();
    value_type pi = M_PI;
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    double L2error2 =integrate( elements( mesh ), ( idv( u )-g )*( idv( u )-g ) ).evaluate()( 0,0 );
    double error = math::sqrt( L2error2 );
    return error;
}

template<int Dim>
double
ddmethod<Dim>::h1Error( element_type& u )
{
    auto Xh=u.functionSpace();
    auto mesh=Xh->mesh();
    value_type pi = M_PI;
    auto gradg = trans( +pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() )*unitX()
                        -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() )*unitY()
                        -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() )*unitZ() );

    double semi_H1error =integrate( elements( mesh ),
                                    ( gradv( u )-gradg )*trans( ( gradv( u )-gradg ) ) ).evaluate()( 0,0 );
    double L2error2 = std::pow( l2Error( u ) , 2 );
    double error = math::sqrt( L2error2 + semi_H1error );
    return error;
}

template<int Dim>
void
ddmethod<Dim>::exportResults( element_type& u, element_type& v, double time )
{
    timers["export"].first.restart();

    M_firstExporter->step( time )->setMesh( u.functionSpace()->mesh() );
    M_firstExporter->step( time )->add( "solution", ( boost::format( "solution-%1%" ) % int( 1 ) ).str(), u );
    M_firstExporter->save();

    M_secondExporter->step( time )->setMesh( v.functionSpace()->mesh() );
    M_secondExporter->step( time )->add( "solution",( boost::format( "solution-%1%" ) % int( 2 ) ).str(), v );
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
            ofs << "casefile: ddmethod-" << Dim << "-" << j << "-1_0.case\n";
        }
    }

    std::cout << "exportResults done" << std::endl;
    timers["export"].second = timers["export"].first.elapsed();
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

    value_type tolerance = doption("tolerance");
    value_type maxIterations = doption("maxIterations");

    Environment::changeRepository( boost::format( "doc/manual/dd/%1%/%2%-%3%/P%4%/h_%5%/" )
                                   % this->about().appName()
                                   %this->shape
                                   % Dim
                                   % Order
                                   %this->meshSize );

    if ( Dim == 2 )
    {

        mesh1 = createGMSHMesh( _mesh=new mesh_type,
                                _desc = ddmethodGeometryLeft( Dim,this->meshSize ) );

        mesh2 = createGMSHMesh( _mesh=new mesh_type,
                                _desc = ddmethodGeometryRight( Dim,this->meshSize ) );

    }


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
        dirichletFlags1+= 2,3,4,5,6,7,8;
        dirichletFlags2+= 1,2,3,4,5;
        interfaceFlags1+= 1;
        interfaceFlags2+= 6;
    }

    auto Xh1 = space_type::New( mesh1 );
    auto Xh2 = space_type::New( mesh2 );
    auto u1 = Xh1->element();
    auto u2 = Xh2->element();
    auto u1old = Xh1->element();
    auto uv = Xh1->element();
    auto uu = Xh2->element();
    value_type pi = M_PI;
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    auto f = pi*pi*Dim*g;
    bool additive = ioption("additive");
    double L2erroru1 = 1.;
    double L2erroru2 = 1.;
    double H1erroru1 = 2.;
    double H1erroru2 = 2.;

    auto Ih12 = opInterpolation( _domainSpace=Xh1, _imageSpace=Xh2, _range=markedfaces( Xh2->mesh(),interfaceFlags2[0] ) );
    auto Ih21 = opInterpolation( _domainSpace=Xh2, _imageSpace=Xh1, _range=markedfaces( Xh1->mesh(),interfaceFlags1[0] ) );

    unsigned int cpt = 0;

    while ( ( L2erroru1 +L2erroru2 ) > tolerance && cpt <= maxIterations )
    {

        std::cout << "===============================\n";
        std::cout << "iteration  : " << cpt << std::endl;
        std::cout << "L2erroru1  : " << L2erroru1  << std::endl;
        std::cout << "L2erroru2  : " << L2erroru2  << std::endl;
        std::cout << "H1erroru1  : " << H1erroru1  << std::endl;
        std::cout << "H1erroru2  : " << H1erroru2  << std::endl;

        // u1old = u1;

        // if( additive )
        // {
        //     if(cpt==0) std::cout << " additive method" << std::endl;
        // }
        // else
        // {
        //     if(cpt==0) std::cout << " multiplicative method" << std::endl;
        // }
        localProblem( u1,
                      dirichletFlags1, g,
                      f,
                      interfaceFlags1,idv( uv ) );
        // if ( !additive )
        //    u1old = u1;
        Ih12->apply( u1, uu );

        localProblem( u2,
                      dirichletFlags2, g,
                      f,
                      interfaceFlags2,idv( uu ) );

        Ih21->apply( u2, uv );

        // compute L2error;
        L2erroru1 = l2Error( u1 );
        L2erroru2 = l2Error( u2 );
        // compute H1error;
        H1erroru1 = h1Error( u1 );
        H1erroru2 = h1Error( u2 );
        // export results
        this->exportResults( u1,u2, cpt );

        ++cpt;
    }; // iteration loop

    std::cout << "------end iteration--------\n";

    std::cout << "number of iteration  : " << cpt-1 << std::endl;

    std::cout << "L2erroru1  : " << L2erroru1  << std::endl;

    std::cout << "L2erroru2  : " << L2erroru2  << std::endl;

    std::cout << "H1erroru1  : " << H1erroru1  << std::endl;

    std::cout << "H1erroru2  : " << H1erroru2  << std::endl;

} // ddmethod::run
} // Feel
/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    Application app;

    ddmethod<2>  Relax;
    Relax.run();

}





