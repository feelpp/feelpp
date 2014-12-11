/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Abdoulaye Samake <abdoulaye.samake1@ujf-grenoble.fr>
   Date: 2011-04-14

   Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

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
   \file test_aitken.cpp
   \author Abdoulaye Samake <abdoulaye.samake1@ujf-grenoble.fr>
   \date 2011-05-13
*/

// Boost::Assign
#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>

/** include predefined feel command line options */
#include <feel/options.hpp>

/** include linear algebra backend */
#include <feel/feelalg/backend.hpp>

#include <feel/feelalg/aitken.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/geotool.hpp>

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
    po::options_description laplacianoptions( "Aitken testsuite options" );
    laplacianoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "relaxmethod", po::value<int>()->default_value( 0 ), "use DD (=0) or DN (=1) method" )
    ( "additive", po::value<int>()->default_value( 0 ), "use relax_aik additive method" )
    ( "x1max", Feel::po::value<double>()->default_value( 1.25 ),  " x1max for first subdomain" )
    ( "x2min", Feel::po::value<double>()->default_value( 1.25 ), " x2min for second subdomain" )
    ( "tol", Feel::po::value<double>()->default_value( 1e-06 ),  " tolerance " )
    ( "imax", Feel::po::value<double>()->default_value( 50 ), " maximum number of iteration" )
    ( "theta", Feel::po::value<double>()->default_value( 1. ), " relaxation parameter" )
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
    AboutData about( "test_aitken" ,
                     "test_aitken" ,
                     "0.2",
                     "nD(n=1,2,3) Aitken relaxation",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Abdoulaye Samake", "developer", "abdoulaye.samake@ujf-grenoble.fr", "" );
    return about;

}

enum DDMethod
{
    // Dirichlet-Dirichlet
    DD = 0,
    // Dirichlet-Neumann
    DN = 1
};

/**
 * \class TestAitken
 *
 * Relax_Aik Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
template<int Dim>
class TestAitken
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
    //! vector type associated with backend
    typedef typename backend_type::vector_type vector_type;

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
    TestAitken( )
        :
        super(  ),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() ),
        timers()
    {
    }

    template<typename DirichletExpr,
             typename RhsExpr,
             typename InterfaceExpr>
    void localProblem( element_type& u,
                       std::vector<int> const& dirichletFlags, DirichletExpr gD,
                       RhsExpr f,
                       std::vector<int> const& interfaceFlags, InterfaceExpr w,
                       DDMethod choice );
    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! shape of the domain
    std::string shape;

    std::map<std::string, std::pair<boost::timer, double> > timers;
    std::map<std::string,double> stats;

    // flags for dirichlet boundary conditions
    std::vector<int> dirichletFlags1;
    std::vector<int> dirichletFlags2;

    // flags for interface conditions
    std::vector<int> interfaceFlags1;
    std::vector<int> interfaceFlags2;

}; // TestAitken

template<int Dim> const uint16_type TestAitken<Dim>::Order;

template< int Dim>
template<typename DirichletExpr,
         typename RhsExpr,
         typename InterfaceExpr>
void
TestAitken<Dim>::localProblem( element_type& u,
                               std::vector<int> const& dirichletFlags, DirichletExpr gD,
                               RhsExpr f,
                               std::vector<int> const& interfaceFlags, InterfaceExpr w,
                               DDMethod choice )


{

    auto Xh=u.functionSpace();
    auto mesh=Xh->mesh();
    element_type v( Xh,"v" );

    auto B = M_backend->newVector( Xh );

    form1( _test=Xh,_vector=B, _init=true ) =
        integrate( elements( mesh ), f*id( v ) );

    if ( choice == DDMethod::DN )
    {
        BOOST_FOREACH( int marker, interfaceFlags )
        {
            form1( _test=Xh,_vector=B ) +=
                integrate( markedfaces( mesh, marker ), w*id( v ) );
        }
    }

    B->close();

    auto A = M_backend->newMatrix( Xh, Xh );

    timers["assembly"].first.restart();

    form2( _test=Xh, _trial=Xh, _matrix=A, _init=true ) =
        integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) );

    A->close();

    BOOST_FOREACH( int marker, dirichletFlags )
    {
        form2( _test=Xh, _trial=Xh, _matrix=A ) +=
            on( markedfaces( mesh, marker ) ,	u, B, gD );
    }

    if ( choice == DDMethod::DD )
    {
        BOOST_FOREACH( int marker, interfaceFlags )
        {
            form2( _test=Xh, _trial=Xh, _matrix=A ) +=
                on( markedfaces( mesh, marker ) ,	u, B, w );
        }
    }

    backend_type::build( soption( _name="backend" ) )->solve( _matrix=A, _solution=u, _rhs=B );

}

template<int Dim>
void
TestAitken<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute TestAitken<" << Dim << ">\n";
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
TestAitken<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) shape = "simplex";

    if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/P%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );


    uint32_type k0 = 0;
    uint32_type k1 = 1;

    value_type x1max = this->vm()["x1max"].template as<double>();
    value_type x2min = this->vm()["x2min"].template as<double>();
    value_type tol = this->vm()["tol"].template as<double>();
    value_type imax = this->vm()["imax"].template as<double>();
    value_type theta = this->vm()["theta"].template as<double>();


    mesh_ptrtype mesh1 = createGMSHMesh( _mesh=new mesh_type,
                                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                         _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" ) % shape % Dim % 1 % k0 ).str() ,
                                                 _addmidpoint=false,
                                                 _usenames=false,
                                                 _shape=shape,
                                                 _dim=Dim,
                                                 _h=X[0],
                                                 _xmin=0.,
                                                 _xmax=x1max,
                                                 _ymin=0.,
                                                 _ymax=2.,
                                                 _zmin=0.,
                                                 _zmax=2.
                                                     ) );

    mesh_ptrtype mesh2 = createGMSHMesh( _mesh=new mesh_type,
                                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                         _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" ) % shape % Dim % 1 % k1 ).str() ,
                                                 _addmidpoint=false,
                                                 _usenames=false,
                                                 _shape=shape,
                                                 _dim=Dim,
                                                 _h=X[0],
                                                 _xmin=x2min,
                                                 _xmax=2.,
                                                 _ymin=0.,
                                                 _ymax=2.,
                                                 _zmin=0.,
                                                 _zmax=2.
                                                     ) );



    if ( Dim == 1 )
    {
        using namespace boost::assign;
        dirichletFlags1 += 1;
        dirichletFlags2 += 3;
        interfaceFlags1 += 3;
        interfaceFlags2 += 1;
    }

    else if ( Dim == 2 )
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
        dirichletFlags1+= 6,15,19,23,28;
        dirichletFlags2+= 6,15,23,27,28;
        interfaceFlags1+= 27;
        interfaceFlags2+= 19;
    }

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    auto Xh1 = space_type::New( mesh1 );
    auto Xh2 = space_type::New( mesh2 );
    element_type u1( Xh1, "u1" );
    element_type u2( Xh2, "u2" );
    element_type u1old( Xh1, "u1old" );
    auto uu = Xh1->element();
    auto lambda = Xh2->element();
    auto residual = Xh2->element();

    AitkenType relaxmethod = ( AitkenType )this->vm()["relaxmethod"].template as<int>();
    auto aitkenRelax =  aitken( _space= Xh2,
                                _type=( relaxmethod == 0 ) ? "standard" : "method1",
                                _initial_theta=theta,
                                _tolerance=tol );
    aitkenRelax.initialize( residual, lambda );

    value_type pi = M_PI;

    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    auto f = pi*pi*Dim*g;

    auto gradg = trans( +pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() )*unitX()+
                        -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() )*unitY()+
                        -pi*sin( pi*Px() )*cos( pi*Py() )*sin( pi*Pz() )*unitZ() );


    bool additive = this->vm()["additive"].template as<int>();

    double err1 = 1.;
    double err2 = 1.;

    double error1 = 2.;
    double error2 = 2.;

    auto Ih12 = opInterpolation( _domainSpace=Xh1, _imageSpace=Xh2, _range=markedfaces( Xh2->mesh(), interfaceFlags2[0] ) );
    auto Ih21 = opInterpolation( _domainSpace=Xh2, _imageSpace=Xh1, _range=markedfaces( Xh1->mesh(), interfaceFlags1[0] ) );

    aitkenRelax.restart();

    while ( !aitkenRelax.isFinished() &&  aitkenRelax.nIterations() <= imax )
    {


        LOG(INFO) << "============================================================\n";
        LOG(INFO) << "iteration  : " << aitkenRelax.nIterations() << "\n";
        LOG(INFO) << "L2erroru1  : " << err1  << "\n";
        LOG(INFO) << "L2erroru2  : " << err2  << "\n";
        LOG(INFO) << "H1erroru1  : " << error1  << "\n";
        LOG(INFO) << "H1erroru2  : " << error2  << "\n";


        u1old = u1;


        if ( additive )
        {

            if ( aitkenRelax.nIterations()==1 ) LOG(INFO) << "test_aitken additive method" << "\n";
        }

        else
        {
            if ( aitkenRelax.nIterations()==1 ) LOG(INFO) << "test_aiken multiplicative method" << "\n";
        }


        Ih21->apply( u2, uu );

        localProblem( u1,
                      dirichletFlags1, /*dirichlet*/g,
                      /*rhs*/f,
                      /**/interfaceFlags1,idv( uu ),
                      DDMethod::DD );

        if ( !additive )
            u1old = u1;

        lambda = u2;

        localProblem( u2,
                      dirichletFlags2, /*dirichlet*/g,
                      /*rhs*/f,
                      /**/ interfaceFlags2,gradv( u1old )*vf::N(),
                      DDMethod::DN );


        residual = u2-lambda;

        u2 = aitkenRelax.apply( residual, u2  );

        aitkenRelax.printInfo();

        aitkenRelax.saveConvergenceHistory( std::string( "history.dat" ) );

        ++aitkenRelax;

        double L2error1 =integrate( elements( mesh1 ), ( idv( u1 )-g )*( idv( u1 )-g ) ).evaluate()( 0,0 );
        err1 = math::sqrt( L2error1 );

        double  L2error2 =integrate( elements( mesh2 ), ( idv( u2 )-g )*( idv( u2 )-g ) ).evaluate()( 0,0 );
        err2 = math::sqrt( L2error2 );

        double semi_H1error1 = 2.;
        double semi_H1error2 = 2.;


        semi_H1error1 =integrate( elements( mesh1 ), ( gradv( u1 )-gradg )*trans( ( gradv( u1 )-gradg ) ) ).evaluate()( 0,0 );
        semi_H1error2 =integrate( elements( mesh2 ), ( gradv( u2 )-gradg )*trans( ( gradv( u2 )-gradg ) ) ).evaluate()( 0,0 );

        error1 = math::sqrt( L2error1 + semi_H1error1 );
        error2 = math::sqrt( L2error2 + semi_H1error2 );

    }; // iteration loop

    export_ptrtype exporter( export_type::New( this->vm(),
                             ( boost::format( "%1%-%2%-%3%-%4%" )
                               % this->about().appName()
                               % shape
                               % k0
                               % Dim ).str() ) );

    export_ptrtype exporter1( export_type::New( this->vm(),
                              ( boost::format( "%1%-%2%-%3%-%4%" )
                                % this->about().appName()
                                % shape
                                % k1
                                % Dim ).str() ) );


    if ( exporter->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";
        exporter->step( 0 )->setMesh( mesh1 );
        exporter->step( 0 )->add( "u1", u1 );
        exporter->save();
        exporter1->step( 1 )->setMesh( mesh2 );
        exporter1->step( 1 )->add( "u2", u2 );
        exporter1->save();


        LOG(INFO) << "exportResults done\n";
    }

    /** \endcode */
} // TestAitken::run

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    Environment env( _argc=argc, _argv=argv, _desc=makeOptions(), _about=makeAbout() );
    /**
     * create an application
     */
    /** \code */
    Application app;

    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    /** \endcode */
    using namespace Feel::vf;

    /**
     * register the simgets
     */
    /** \code */
    app.add( new TestAitken<1>() );
    app.add( new TestAitken<2>() );
    // app.add( new TestAitken<3>( app.vm(), app.about() ) );
    /** \endcode */

    /**
     * run the application
     */
    /** \code */
    app.run();
    /** \endcode */
}





