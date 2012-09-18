//#include <boost/test/unit_test.hpp>
//using boost::unit_test::test_suite;


#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>

#define USE_GMM   1

#if USE_GMM
#include <gmm_iter_solvers.h>

//#include <feel/feelalg/matrixublas.hpp>
#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/matrixgmm.hpp>
#include <feel/feelcore/application.hpp>


#else
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
#include <feel/feelalg/solverlinearpetsc.hpp>
#endif

#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/lagrange.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporterensight.hpp>

#include <feel/feelvf/vf.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>


using boost::timer;

using namespace std;
using namespace Feel;

/** Streams **/

ofstream S_error( "Data/Multi/error_2D.dat" );
//ofstream S_cond("Data/Multi/cond_2D.dat");
ofstream S_timing( "Data/Multi/timing_2D.dat" );
ofstream S_dofs( "Data/Multi/dofs_2D.dat" );
ofstream S_iter( "Data/Multi/solver_iter_2D.dat" );

typedef int16_type int_type;

void timeout( double time );

#if 0
template<int_type N, typename T> T Poisson();



using boost::unit_test::test_suite;

template<int_type P>
void add_tests( test_suite* test )
{
    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Poisson<P,value_type> ) ) );

    add_tests<2*P>( test );
}

template<>
void add_tests<16>( test_suite* test )
{
    typedef double value_type;

    test->add( BOOST_TEST_CASE( ( Poisson<4,value_type> ) ) );
}


/* Main test function */

test_suite*
init_unit_test_suite( int /*argc*/, char** /*argv[]*/ )
{
    test_suite* test = BOOST_TEST_SUITE( "Spectral Accuracy (BA) test using Feel Language" );
#if defined( FEELPP_HAS_QD_REAL)
    unsigned int old_cw;
    fpu_fix_start( &old_cw );
#endif
    add_tests<16>( test );
    return test;
}
#endif

#define FUNC_EXPO 1


namespace Feel
{
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "pfem" ,
                           "pfem" ,
                           "0.1",
                           "Parallel Fem",
                           Feel::AboutData::License_LGPL,
                           "Copyright (c) 2005,2006 EPFL" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

template<int_type N, typename T = double>
class Poisson
    :
#if USE_GMM
public Application
#else
public Application
#endif
{
    typedef T value_type;
#if USE_GMM
    typedef Application super;
    /*matrix*/
    typedef MatrixGmm<value_type, gmm::row_major> sparse_matrix_type;
    typedef VectorUblas<value_type> vector_type;
#else
    typedef Application super;
    /*matrix*/
    typedef MatrixPetsc<value_type> sparse_matrix_type;
    typedef VectorPetsc<value_type> vector_type;
#endif

public:
    Poisson( int argc, char** argv, AboutData const& ad )
        :
        super( argc, argv, ad )
    {}
    Poisson( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od )
    {}

    std::pair<size_type, T>
    solve( sparse_matrix_type const& A,
           vector_type& X,
           vector_type const& F,
           double tol = 1e-15,
           int nits = 5000 )
    {
        po::variables_map vm = this->vm();;

#if USE_GMM
        double gmm_tol = tol;

        if ( vm.count( "gmm_tol" ) )
            gmm_tol = vm["gmm_tol"].as<double>();

        gmm::iteration projiter( gmm_tol );
        projiter.set_noisy( vm["gmm_noisy"].as<int>() );
        projiter.set_maxiter( nits );

        if ( vm["gmm_pc"].as<int>() == 0 )
        {
            gmm::identity_matrix P;
            gmm::cg( A.mat(), X, F, P, projiter );
        }

        else if ( vm["gmm_pc"].as<int>() == 1 )
        {
            gmm::diagonal_precond<typename sparse_matrix_type::matrix_type> P( A.mat() );
            gmm::cg( A.mat(), X, F, P, projiter );
        }

        else if ( vm["gmm_pc"].as<int>() == 2 )
        {
            int fill = vm["gmm_pc_ildltt_fill"].as<int>();
            double drop = vm["gmm_pc_ildltt_drop"].as<double>();

            gmm::ildltt_precond<typename sparse_matrix_type::matrix_type> P( A.mat(), fill, drop );
            gmm::cg( A.mat(), X, F, P, projiter );
        }

        if ( !projiter.converged() )
        {
            cout << "[gmm] Solver didn't converge" << endl;
        }

        else
        {
            cout << "[gmm] Solver converged in " << projiter.get_iteration() << " iterations ";
        }

        return std::make_pair( projiter.get_iteration(), tol );
#else
        SolverLinearPetsc<double> solver;
        size_type its;
        double resid;
        boost::tie( its,resid ) = solver.solve( A, X, F, tol, 1000 );

        std::vector<double> res;
        solver.getResidualHistory( res );

        for ( int i = 0; i < res.size(); ++i )
        {
            if ( i % 100 == 0 || i == res.size()-1 )
                cout << "iteration " << i << " residual " << res[i] << "\n";
        }

        return std::make_pair( its, resid );
#endif
    }
    void run()
    {
        using namespace std;
        cout << "N = " << N << std::endl;

        timer t_tot;
        timer t;

        Gmsh __gmsh;
        string fname;
        ostringstream ostr;

        po::variables_map vm = this->vm();
        double h= vm["hsize"].as<double>();

        if ( vm["domain"].as<int>()  == 1 )
        {
            ostr << "h=" << h << ";\n"
                 << "Point(1) = {-1, -1,0.0,h};\n"
                 << "Point(2) = { 1, -1,0.0,h};\n"
                 << "Point(3) = {-1,  1,0.0,h};\n"
                 << "Line(1) = {2,3};\n"
                 << "Line(2) = {3,1};\n"
                 << "Line(3) = {1,2};\n"
                 << "Line Loop(4) = {1,2,3};\n"
                 << "Plane Surface(5) = {4};\n"
                 << "Physical Surface(30) = {5};\n"
                 << "Physical Line(31) = {1};\n"
                 << "Physical Line(32) = {2};\n"
                 << "Physical Line(33) = {3};\n";
        }

        else if ( vm["domain"].as<int>()  == 2 )
        {
            ostr << "h=" << h << ";\n"
                 << "Point(1) = {-1, -1,0.0,h};\n"
                 << "Point(2) = { 1, -1,0.0,h};\n"
                 << "Point(3) = { 1,  1,0.0,h};\n"
                 << "Point(4) = {-1,  1,0.0,h};\n"
                 << "Line(1) = {1,2};\n"
                 << "Line(2) = {2,3};\n"
                 << "Line(3) = {3,4};\n"
                 << "Line(4) = {4,1};\n"
                 << "Line Loop(5) = {1,2,3,4};\n"
                 << "Plane Surface(6) = {5};\n"
                 << "Physical Surface(30) = {6};\n"
                 << "Physical Line(31) = {1};\n"
                 << "Physical Line(32) = {2};\n"
                 << "Physical Line(33) = {3};\n"
                 << "Physical Line(34) = {4};\n";
        }

        else
        {
            std::cerr << "Invalid domain\n";
            return;
        }

        t.restart();
        t_tot.restart();
        cout <<"Mesh generation ... ";
        fname = __gmsh.generate( "triangle", ostr.str() );

        /* Mesh */

        typedef Mesh<GeoEntity<Simplex<2, 1> > > mesh_type;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
        mesh_ptrtype mesh( new mesh_type );

        ImporterGmsh<mesh_type> import( fname );
        mesh->accept( import );
        mesh->partition();
        timeout( t.elapsed() );

        /** Approximation Space **/

        using namespace Feel::vf;

        cout << "Space and forms generation ...\n";
        t.restart();
        typedef BoundaryAdaptedPolynomialSet<2, N, Scalar, value_type, Simplex> basis_type;
        //typedef OrthonormalPolynomialSet<2, N, Scalar, value_type, Simplex> basis_type;
        //typedef OrthogonalPolynomialSet<2, N, Scalar, value_type, Simplex> basis_type;
        //typedef fem::Lagrange<2, N, Scalar, Continuous, value_type, Simplex> basis_type;

        typedef FunctionSpace<mesh_type, basis_type, value_type > space_type;
        boost::shared_ptr<space_type> Xh( new space_type( mesh ) );


        S_dofs << Xh.get()->nDof() << std::endl;

        IM_PK<2, 2*N, value_type, Gauss> im;

        typename space_type::element_type uM( Xh.get() );
        typename space_type::element_type vM( Xh.get() );

        sparse_matrix_type Mass;
        Mass.init( Xh->dof()->nDof(), Xh->dof()->nDof(),
                   Xh->dof()->nLocalDof(), Xh->dof()->nLocalDof() );
        BilinearForm<space_type, space_type, sparse_matrix_type> m( Xh, Xh, Mass );
        m =  integrate( elements( *mesh ), im, idt( uM )*id( vM ) );
        Mass.close();
        std::cout << "Mass matrix ... " << t.elapsed() << "s\n";
        t.restart();
        //Mass.printMatlab( "Mass.m" );

        typename space_type::element_type uG( Xh.get() );
        typename space_type::element_type vG( Xh.get() );

        sparse_matrix_type Stiff;
        Stiff.init( Xh->dof()->nDof(), Xh->dof()->nDof(),
                    Xh->dof()->nLocalDof(), Xh->dof()->nLocalDof() );
        BilinearForm<space_type, space_type, sparse_matrix_type> sti( Xh, Xh, Stiff );
        sti =  integrate( elements( *mesh ), im, dxt( uG )*dx( vG ) + dyt( uG )*dy( vG ) );

        Stiff.close();
        std::cout << "Stiff matrix ... " << t.elapsed() << "s\n";
        t.restart();
        //Stiff.printMatlab( "Stiff.m" );

        sparse_matrix_type Global;
        Global.init( Xh->dof()->nDof(), Xh->dof()->nDof(),
                     Xh->dof()->nLocalDof(), Xh->dof()->nLocalDof() );
        BilinearForm<space_type, space_type, sparse_matrix_type> g( Xh, Xh, Global );
        g =  integrate( elements( *mesh ), im, idt( uG )*id( vG ) + dxt( uG )*dx( vG ) + dyt( uG )*dy( vG ) );

        Global.close();
        std::cout << "Stiff+Mass matrix ... " << t.elapsed() << "s\n";
        t.restart();
        //Global.printMatlab( "Global.m" );


        // -- PROBLEM PARAMETERS --

#if FUNC_EXPO
#if 0 //sin(pi x)cos(pi y)
        value_type pi = 4.0 * math::atan( value_type( 1.0 ) );

        __typeof__( sin( pi*Px() )*cos( pi*Py() ) )
        exact_sol = sin( pi*Px() )*cos( pi*Py() );

        __typeof__( ( 2.0*pi*pi+1.0 )*sin( pi*Px() )*cos( pi*Py() ) )
        f = ( 2.0*pi*pi+1.0 )*sin( pi*Px() )*cos( pi*Py() );

        __typeof__( pi*cos( pi*Px() )*cos( pi*Py() ) )
        grad_x = pi*cos( pi*Px() )*cos( pi*Py() );

        __typeof__( -pi*sin( pi*Px() )*sin( pi*Py() ) )
        grad_y = -pi*sin( pi*Px() )*sin( pi*Py() );

#else // exp(-(x+y))
        __typeof__( exp( -Px()-Py() ) )
        exact_sol = exp( -Px()-Py() );

        __typeof__( -exp( -Px()-Py() ) )
        f = -exp( -Px()-Py() );

        __typeof__( -exp( -Px()-Py() ) )
        grad_x = -exp( -Px()-Py() );

        __typeof__( -exp( -Px()-Py() ) )
        grad_y = -exp( -Px()-Py() );
#endif
#endif

        //vector_type F( Xh->dof()->nDof(), Xh->dof()->nLocalDof() );
        vector_type F( Xh->dof()->nDof() );
        LinearForm<space_type, vector_type> l( Xh, F );

#if FUNC_EXPO
        l = integrate( elements( *mesh ), im, f*id( vG ) )
            + integrate( markedfaces( *mesh,31 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( vG ) )
            + integrate( markedfaces( *mesh,32 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( vG ) )
            + integrate( markedfaces( *mesh,33 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( vG ) );

        if ( vm["domain"].as<int>() == 2 )
            l += integrate( markedfaces( *mesh,34 ), im, ( Nx()*grad_x + Ny()*grad_y )*id( vG ) );

#else
        l = integrate( elements( *mesh ), im, id( vG ) );
#endif


        F.close();
        std::cout << "Rhs ... " << t.elapsed() << "s\n";
        t.restart();
        timeout( t.elapsed() );

        typename space_type::element_type vP( Xh.get() );
        //vector_type ProjF( Xh->dof()->nDof(), Xh->dof()->nLocalDof() );
        vector_type ProjF( Xh->dof()->nDof() );
        LinearForm<space_type, vector_type > pro( Xh, ProjF );
#if FUNC_EXPO
        pro = integrate( elements( *mesh ), im, exact_sol*id( vP ) );
#else
        pro = integrate( elements( *mesh ), im, id( vP ) );
#endif
        ProjF.close();

        typename space_type::element_type sol_proj( Xh.get() );

        //vector_type U( Xh->dof()->nDof(), Xh->dof()->nLocalDof() );
        vector_type U( Xh->dof()->nDof() );


        //cout << sol_proj.container() << std::endl;


        /** System Resolution **/

        cout << "System iterative resolution ...\n";

        typename space_type::element_type sol_ap( Xh.get() );

        //vector_type V( Xh->dof()->nDof(), Xh->dof()->nLocalDof() );
        vector_type V( Xh->dof()->nDof() );
        t.restart();
        size_type its;
        value_type resid;
        boost::tie( its,resid ) = this->solve( Global, V, F, 2e-15, 10000 );
        std::cout << "============================================================\n";
        //std::cout << "[solver type] s= " << solver2.solverType() << "\n";
        //std::cout << "[pc type]     p= " << solver2.preconditionerType() << "\n";
        std::cout << "[niter]   niter= " << its << "\n";
        std::cout << "[resid]   resid= " << resid << "\n";
        std::cout << "[nelts]   nelts= " << Xh->mesh()->numElements() << "\n";
        std::cout << "[ndofs]   ndofs= " << Xh->dof()->nDof() << "\n";
        std::cout << "[timer]       t= " << t.elapsed() << "\n";
        std::cout << "============================================================\n";

#if !USE_GMM
        cout << "Error Analysis ...\n";

        size_type s = U.lastLocalIndex() - U.firstLocalIndex();
        FEELPP_ASSERT( s != invalid_size_type_value ).error( "invalid size" );

        sol_proj.resize( Xh->dof()->nDof() );
        sol_ap.resize( Xh->dof()->nDof() );

        for ( size_type i = 0; i < s; ++i )
        {
            sol_proj( i ) = U( i );
            sol_ap( i ) = V( i );
        }

        VecAXPY( V.vec(), PetscScalar( -1 ), U.vec() );
        typename space_type::element_type error( Xh.get() );
        error = sol_proj - sol_ap;
        vector_type one( V.size() );
        typename mesh_type::element_iterator it = mesh->beginElement();
        typename mesh_type::element_iterator en = mesh->endElement();

        while ( it != en )
        {
            for ( int i = 0; i < 3; ++i )
                one.set( boost::get<0>( Xh->dof()->localToGlobal( it->id(), i ) ), 1 );

            ++it;
        }

        std::cout << "area(triangle/bf) = " << m.matrix().energy( one,one ) << "\n";
        PetscScalar res_1;
        VecDot( ProjF.vec(), one.vec(), &res_1 );

        std::cout << "area(triangle/lf) = " << res_1 << "\n";

        value_type ErrH1 = math::sqrt( g.matrix().energy( V,V ) );
        value_type ErrL2 = math::sqrt( m.matrix().energy( V,V ) );


        std::cout << "Error in H1-norm = " << ErrH1 << std::endl;
        std::cout << "Error in L2-norm = " << ErrL2 << std::endl;

        S_error << N << " " << ErrL2 << "  "  << ErrH1<< std::endl;
#endif // ! USE_GMM
        double time_tot( t_tot.elapsed() );

        cout << "Total resolution time = ";
        timeout( time_tot );

        S_timing << N << " " << time_tot << std::endl;

        if ( vm.count( "export" ) )
        {
            typename space_type::P1Lagrange p1lag( Xh );
            typename space_type::P1Lagrange::p1_type::element_type u_p1 = p1lag( sol_ap );
            typename space_type::P1Lagrange::p1_type::element_type ex_p1 = p1lag( sol_proj );
            //typename space_type::P1Lagrange::p1_type::element_type er_p1 = p1lag( error );

            //
            // save the results
            //
            {
                std::ostringstream ostr;
                ostr << "test_mult_petsc";
                typedef Exporter<mesh_type> export_type;
                typedef typename Exporter<mesh_type>::timeset_type timeset_type;

                boost::shared_ptr<export_type> __ensight;
                __ensight = boost::shared_ptr<export_type>( new ExporterEnsight<mesh_type>( ostr.str() ) );
                typename export_type::timeset_ptrtype __ts( new timeset_type( ostr.str(), true ) );
                __ts->setTimeIncrement( 1.0 );
                __ensight->addTimeSet( __ts );
                typename timeset_type::step_ptrtype __step = __ts->step( 1.0 );
                __step->setMesh( p1lag.mesh() );

                //__step->addElementScalar( "pfem_pid", pid.size(), pid.begin(), pid.end() );
                __step->addNodalScalar( "pfem_u_p1", u_p1.size(), u_p1.begin(), u_p1.end() );
                __step->addNodalScalar( "pfem_exact_p1", ex_p1.size(), ex_p1.begin(), ex_p1.end() );
                //__step->addNodalScalar( "pfem_error_p1", er_p1.size(), er_p1.begin(), er_p1.end() );

                __ensight->save();
            }
        }
    }


    void timeout( double time )
    {
        if ( time < 60 )
            cout << time << " seconds.\n";

        else if ( time < 3600 )
        {
            int_type min = ( int_type )floor( time )/60;
            cout << min << " minutes " << ( time-60*min ) << " seconds.\n";
        }

        else
        {
            int_type hour = ( int_type )floor( time )/3600;
            int_type min = ( int_type )floor( time-hour*3600 )/60;
            cout << hour << " hours " << min << " minutes " << ( time-3600*hour-60*min ) << " seconds.\n";
        }
    }
};
}

int main( int argc,  char** argv )
{
    using namespace Feel;

    Feel::po::options_description gmm_options( "gmm options" );
    gmm_options.add_options()
    ( "gmm_noisy", Feel::po::value<int>()->default_value( 0 ), "gmm noisy level" )
    ( "gmm_tol", Feel::po::value<double>()->default_value( 1e-12 ), "gmm tolerance" )
    ( "gmm_solver", Feel::po::value<int>()->default_value( 0 ), "gmm solver" )
    ( "gmm_pc", Feel::po::value<int>()->default_value( 0 ), "gmm preconditioner" )
    ( "gmm_pc_ildltt_fill", Feel::po::value<int>()->default_value( 2 ), "gmm ildltt fill" )
    ( "gmm_pc_ildltt_drop", Feel::po::value<double>()->default_value( 1e-3 ), "gmm ildltt drop" )
    ;
    Feel::po::options_description test( "pfem options" );
    test.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.2 ), "h value" )
    ( "beta", Feel::po::value<double>()->default_value( 0.0 ), "beta value in -Delta u + beta u = f" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "domain", Feel::po::value<int>()->default_value( 1 ), "type of domain (triangle=1, quadrangle=2)" )
    ( "N", Feel::po::value<int>()->default_value( 5 ), "Polynomial order (5,10,15,20)" )
    ( "testall", "run all test cases" )
    ( "export", "export results" )
    ;

    test.add( gmm_options );

    Poisson<5> app5( argc, argv, makeAbout(), test );
    Poisson<10> app10( argc, argv, makeAbout(), test );
    Poisson<15> app15( argc, argv, makeAbout(), test );
    Poisson<20> app20( argc, argv, makeAbout(), test );
    Poisson<25> app25( argc, argv, makeAbout(), test );

    VLOG(1) << "N process: " << Application::nProcess() << "\n"
            << "Id : " << Application::processId() << "\n";

    if ( app5.vm()["N"].as<int>()  == 5 )
    {
        app5.run();
    }

    if ( app5.vm()["N"].as<int>()  == 10 )
    {
        app10.run();
    }

    if ( app5.vm()["N"].as<int>()  == 15 )
    {
        app15.run();
    }

    if ( app5.vm()["N"].as<int>()  == 20 )
    {
        app20.run();
    }

    if ( app5.vm()["N"].as<int>()  == 25 )
    {
        app25.run();
    }

}
