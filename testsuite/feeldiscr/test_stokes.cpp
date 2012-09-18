/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-06-06

  Copyright (C) 2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_stokes.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-06-06
 */
//#define USE_BOOST_TEST 1
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
using boost::unit_test::test_suite;

#include <map>
#include <boost/lambda/bind.hpp>

#include <feel/feelconfig.h>

#include <feel/feelcore/application.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporterensight.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/matrixublas.hpp>
#include <feel/feelalg/matrixgmm.hpp>
#include <feel/feelalg/cg.hpp>
#include <feel/feelalg/preconditioner.hpp>

#include <feel/feelvf/vf.hpp>
#include <gmm_iter_solvers.h>

//
// global variables
//
int Nstep = 1;
const Feel::uint16_type gtOrder = 1;

//
// Global functions
//
template<typename VectorType>
void printVector( VectorType vector, std::string name,
                  std::ostream& os = std::cout )
{
    for ( int i=0; i<vector.size(); ++i )
    {
        os << name << "[" << i << "] = " << vector[i] << std::endl;
    }
}

template<typename MatrixType>
void printMatrix( MatrixType matrix, std::string name,
                  std::ostream& os = std::cout )
{
    os << name << " =" << std::endl;

    for ( int i=0; i<matrix.size1(); ++i )
    {
        for ( int j=0; j<matrix.size2(); ++j )
        {
            os << matrix( i,j ) << " ";
        }

        os << std::endl;
    }
}


inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description testoptions( "test convergence options" );
    testoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "nmesh", Feel::po::value<int>()->default_value( 3 ), "number of meshes to build " )
    ( "beta", Feel::po::value<double>()->default_value( 1.0 ), "beta value in -Delta u + beta u = f" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "Dirichlet condition type(0=elimination,1=penalisation, 2=weak" )
    ( "testall", "run all test cases" )
    ( "export", "export results(ensight, data file(1D)" )

    ;
    Feel::po::options_description solveroptions( "algebraic solver options" );
    solveroptions.add_options()
    ( "verbose", Feel::po::value<int>()->default_value( 0 ), "(=0,1,2) print solver iterations" )
    ( "maxiter", Feel::po::value<int>()->default_value( 1000 ), "set maximum number of iterations" )
    ;
    return testoptions.add( solveroptions );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_convergence" ,
                           "test_convergence" ,
                           "0.1",
                           "Convergence test in nD(n=1,2,3) for elliptic PDEs",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2005-2006 EPFL" );

    about.addAuthor( "Christoph Winkelmann", "developer", "christoph.winkelmann@epfl.ch", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


namespace Feel
{

template< typename value_type,
          typename Basis,
          int gtOrder,
          int imOrder>
class TestConvergence
    :
public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type dim = mpl::at_c<Basis,0>::type::nDim;
    static const uint16_type nDim = mpl::at_c<Basis,0>::type::nDim;
    static const uint16_type feOrder = mpl::at_c<Basis,0>::type::nOrder;

    /*matrix*/
    typedef MatrixGmm<value_type, gmm::row_major> sparse_matrix_type;
    /*mesh*/
    typedef Mesh<GeoEntity<Simplex<dim, 1> > > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptr_type;

    /*basis*/
    typedef Basis basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;

    /*quadrature*/
    typedef IM_PK<dim, imOrder, value_type> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef typename Exporter<mesh_type>::timeset_type timeset_type;

    TestConvergence( int argc, char** argv, AboutData const& ad )
        :
        super( argc, argv, ad ),
        initialMeshSize( this->vm()["hsize"].template as<double>() ),
        beta( this->vm()["beta"].template as<double>() ),
        bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        nMesh( this->vm()["nmesh"].template as<int>() ),
        exporter( new ExporterEnsight<mesh_type>( "test_convergence" ) ),
        timeSet( new timeset_type( "test_convergence" ) ),
        timers( nMesh )

    {
        VLOG(1) << "[TestConvergence] hsize = " << initialMeshSize << "\n";
        VLOG(1) << "[TestConvergence] beta = " << beta << "\n";
        VLOG(1) << "[TestConvergence] bccoeff = " << bcCoeff << "\n";
        VLOG(1) << "[TestConvergence] nmesh = " << nMesh << "\n";
        VLOG(1) << "[TestConvergence] bctype = " << this->vm()["bctype"].template as<int>() << "\n";
        VLOG(1) << "[TestConvergence] export = " << this->vm().count( "export" ) << "\n";

        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
    }

    TestConvergence( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        initialMeshSize( this->vm()["hsize"].template as<double>() ),
        beta( this->vm()["beta"].template as<double>() ),
        bcCoeff( this->vm()["bccoeff"].template as<double>() ),
        nMesh( this->vm()["nmesh"].template as<int>() ),
        exporter( new ExporterEnsight<mesh_type>( "test_convergence" ) ),
        timeSet( new timeset_type( "test_convergence" ) ),
        timers( nMesh )
    {
        VLOG(1) << "[TestConvergence] hsize = " << initialMeshSize << "\n";
        VLOG(1) << "[TestConvergence] beta = " << beta << "\n";
        VLOG(1) << "[TestConvergence] bccoeff = " << bcCoeff << "\n";
        VLOG(1) << "[TestConvergence] nmesh = " << nMesh << "\n";
        VLOG(1) << "[TestConvergence] bctype = " << this->vm()["bctype"].template as<int>() << "\n";
        VLOG(1) << "[TestConvergence] export = " << this->vm().count( "export" ) << "\n";

        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
    }

    TestConvergence( TestConvergence const& tc )
        :
        super( tc ),
        initialMeshSize( tc.initialMeshSize ),
        beta( tc.beta ),
        bcCoeff( tc.bcCoeff ),
        nMesh( tc.nMesh ),
        exporter( new ExporterEnsight<mesh_type>( "test_convergence" ) ),
        timeSet( new timeset_type( "test_convergence" ) ),
        timers( tc.timers )
    {
        VLOG(1) << "[TestConvergence] hsize = " << initialMeshSize << "\n";
        VLOG(1) << "[TestConvergence] beta = " << beta << "\n";
        VLOG(1) << "[TestConvergence] bccoeff = " << bcCoeff << "\n";
        VLOG(1) << "[TestConvergence] nmesh = " << nMesh << "\n";
        VLOG(1) << "[TestConvergence] bcweak = " << this->vm().count( "bcweak" ) << "\n";
        VLOG(1) << "[TestConvergence] export = " << this->vm().count( "export" ) << "\n";

        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptr_type createMesh( double meshSize );

    /**
     * solve the problem for the mesh size \c meshSize
     */
    value_type onceTest( int iMesh, double meshSize );

    /**
     * alias for run()
     */
    void operator()()
    {
        run();
    }

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve non-symmetric system
     */
    template<typename Mat, typename Vec1, typename Vec2>
    void solveNonSym( int iMesh, Mat const& D, Vec1& u, Vec2 const& F );

    /**
     * solve symmetric system
     */
    template<typename Mat, typename Vec1, typename Vec2>
    void solveSym( int iMesh, Mat const& D, Vec1& u, Vec2 const& F );

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( int iMesh, element_0_type const& u, element_1_type const& p );

private:

    double initialMeshSize;
    double beta;
    double bcCoeff;
    int nMesh;

    boost::shared_ptr<export_type> exporter;
    typename export_type::timeset_ptrtype timeSet;

    std::vector<std::map<std::string,std::pair<boost::timer,double> > > timers;
}; // TestConvergence

template< typename value_type,
          typename Basis,
          int gtOrder,
          int imOrder>
typename TestConvergence<value_type,Basis,gtOrder, imOrder>::mesh_ptr_type
TestConvergence<value_type,Basis,gtOrder, imOrder>::createMesh( double meshSize )
{
    mesh_ptr_type mesh( new mesh_type );

    Gmsh __gmsh;
    __gmsh.setOrder( GMSH_ORDER_ONE );
    std::ostringstream ostr;
    std::ostringstream nameStr;

    switch ( dim )
    {
    case 2:
        ostr << "h=" << meshSize << ";\n"
             << "Point(1) = {-1,-1,0,h};\n"
             << "Point(2) = {1,-1,0,h};\n"
             << "Point(3) = {1,1,0,h};\n"
             << "Point(4) = {-1,1,0,h};\n"
             << "Line(1) = {1,2};\n"
             << "Line(2) = {2,3};\n"
             << "Line(3) = {3,4};\n"
             << "Line(4) = {4,1};\n"
             << "Line Loop(4) = {1,2,3,4};\n"
             << "Plane Surface(5) = {4};\n"
             << "Physical Line(10) = {1,3};\n"
             << "Physical Line(20) = {2,4};\n"
             << "Physical Surface(6) = {5};\n";
        nameStr << "square." << meshSize;
        break;

    case 3:
        ostr << "h=" << meshSize << ";\n"
             << "Point(1) = {-1,-1,-1,h};\n"
             << "Point(2) = {-1, 1,-1,h};\n"
             << "Point(3) = { 1, 1,-1,h};\n"
             << "Point(4) = { 1,-1,-1,h};\n"
             << "Line(1) = {1,4};\n"
             << "Line(2) = {4,3};\n"
             << "Line(3) = {3,2};\n"
             << "Line(4) = {2,1};\n"
             << "Line Loop(5) = {3,4,1,2};\n"
             << "Plane Surface(6) = {5};\n"
             << "Extrude Surface {6, {0,0,2}};\n"
             << "Physical Surface(10) = {15,23,6,28};\n"
             << "Physical Surface(20) = {19,27};\n"
             << "Surface Loop(31) = {28,15,-6,19,23,27};\n"
             << "Volume(1) = {31};\n"
             << "Physical Volume(2) = {1};\n";
        nameStr << "cube." << meshSize;
        break;

    default:
        std::ostringstream os;
        os << "invalid dimension: " << dim;
        throw std::logic_error( os.str() );
    }

    std::string fname = __gmsh.generate( nameStr.str(), ostr.str() );
    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // TestConvergence::createMesh


template< typename value_type,
          typename Basis,
          int gtOrder,
          int imOrder>
void
TestConvergence<value_type,Basis,gtOrder, imOrder>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    std::vector<value_type> h( nMesh );
    std::vector<value_type> err( nMesh );
    std::vector<value_type> order( nMesh );
    std::vector<value_type> time( nMesh );

    value_type meshSize = initialMeshSize;

    for ( int iMesh = 0; iMesh<nMesh; ++iMesh, meshSize /= 2.0 )
    {
        h[iMesh] = meshSize;
        err[iMesh] = this->onceTest( iMesh, meshSize );

        if ( iMesh == 0 )
        {
            order[iMesh] = 0;
        }

        else
        {
            order[iMesh] = std::log( err[iMesh]/err[iMesh-1] ) /
                           std::log(  h[iMesh]/  h[iMesh-1] );
        }
    }

    int width = 20;
    std::cout << std::setw( width/2 ) << std::right << "h"
              << std::setw( width/2 ) << std::right << "P"
              << std::setw( width/1 ) << std::right << "rel. L2 error"
              << std::setw( width/1 ) << std::right << "order"
              << std::setw( width/1 ) << std::right << "ass. time(s)"
              << std::setw( width/1 ) << std::right << "sol. time(s)"
              << std::endl;

    for ( int iMesh = 0; iMesh<nMesh; ++iMesh )
    {
        std::cout << std::setw( width/2 ) << std::right << h[iMesh]
                  << std::setw( width/2 ) << std::right << feOrder
                  << std::setw( width/1 ) << std::right << err[iMesh]
                  << std::setw( width/1 ) << std::right << order[iMesh]
                  << std::setw( width/1 ) << std::right << timers[iMesh]["assembly"].second
                  << std::setw( width/1 ) << std::right << timers[iMesh]["solver"].second
                  << std::endl;
    }

    value_type orderTol = 0.09;
#if USE_BOOST_TEST
    BOOST_CHECK( order[nMesh-1] >= feOrder + 1 - orderTol );
#else

    if ( order[nMesh-1] < feOrder + 1 - orderTol )
        std::cout << "[FAILURE] wrong L2 convergence order for FE<"
                  << dim << "," << feOrder << ">\n"
                  << "[FAILURE] order should have been " << feOrder + 1 << "\n";

#endif
} // TestConvergence::run

template< typename value_type,
          typename Basis,
          int gtOrder,
          int imOrder>
value_type
TestConvergence<value_type,Basis,gtOrder, imOrder>::onceTest( int iMesh, double meshSize )
{
    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;

    value_type pi = 4.0 * math::atan( value_type( 1.0 ) );

    mesh_ptr_type mesh = createMesh( meshSize );


    // -- SPACE, FUNCTIONS AND QUADRATURE --
    timers[iMesh]["init"].first.restart();
    space_ptrtype Xh = space_type::New( mesh );
    //Xh->dof()->showMe();
    element_type U( Xh, "u" );
    element_type V( Xh, "v" );
    element_0_type u( U.template element<0>() );
    element_0_type v( U.template element<0>() );
    element_1_type p( U.template element<1>() );
    element_1_type q( U.template element<1>() );
    timers[iMesh]["init"].second = timers[iMesh]["init"].first.elapsed();

    IM_PK<dim, imOrder> im;

    // -- EXACT SOLUTION --
    __typeof__(     sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() ) )
    exact_sol = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );

    // -- RIGHT HAND SIDE --
    VectorUblas<value_type> F( u.size() );
    LinearForm<space_type, VectorUblas<value_type> > vf_F( Xh, F );

    timers[iMesh]["assembly"].first.restart();
    //vf_F =   integrate( elements(*mesh), im, 0 );
    timers[iMesh]["assembly"].second = timers[iMesh]["assembly"].first.elapsed();
    // -- MATRIX --
    sparse_matrix_type D;
    BilinearForm<space_type, space_type, sparse_matrix_type> vf_D( Xh, Xh, D );
    timers[iMesh]["assembly"].first.restart();
    vf_D =   integrate( elements( *mesh ), im,
                        //+ dot( gradt(u), grad(v) )
                        dxt( u )*dx( v )+ dyt( u )*dy( v )+dzt( u )*dz( v )
                        - div( v ) * idt( p )
                        - divt( u ) * id( q )
                        + 1e-6*idt( p )*id( q ) );

    if ( this->vm()[ "bctype" ].template as<int>() == 0 ||
            this->vm()[ "bctype" ].template as<int>() == 1 )
        vf_D +=
            on( markedfaces( *mesh,10 ), u, F, oneX(),
                on_strategy_type( this->vm()["bctype"].template as<int>() ) ) +
            on( markedfaces( *mesh,20 ), u, F, 0,
                on_strategy_type( this->vm()["bctype"].template as<int>() ) );

    D.close();
    timers[iMesh]["assembly"].second = timers[iMesh]["assembly"].first.elapsed();

    this->solveNonSym( iMesh, D, U, F );

    this->exportResults( iMesh,
                         U.template element<0>(),
                         U.template element<1>() );

    ++Nstep;

    //return errL2 / uExL2;
    return 1;

    //     std::cout << "h = " << meshSize << std::endl;
    //     std::cout << "L2 error = " << errL2 << std::endl;
    //     std::cout << "relative L2 error = " << errL2 / uExL2 << std::endl;

} // TestConvergence::onceTest

template< typename value_type,
          typename Basis,
          int gtOrder,
          int imOrder>
template<typename Mat, typename Vec1, typename Vec2>
void
TestConvergence<value_type,Basis,gtOrder, imOrder>::solveNonSym( int iMesh,
        Mat const& D,
        Vec1& u,
        Vec2 const& F )
{
    timers[iMesh]["solver"].first.restart();

    gmm::iteration iter( 2.0E-10 );
    iter.set_noisy( this->vm()["verbose"].template as<int>() );
    iter.set_maxiter( this->vm()["maxiter"].template as<int>() );
    // incomplete LU with k fill-in and threshold preconditioner.
    // Efficient but could be costly.
    gmm::ilut_precond<typename sparse_matrix_type::matrix_type> P( D.mat(), 2, 1e-3 );

    // Conjugate gradient
    gmm::bicgstab( D.mat(), u.container(), F, P, iter );

    if ( !iter.converged() )
    {
        std::cout << "stiffness solver didn't converge" << std::endl;
    }

    timers[iMesh]["solver"].second = timers[iMesh]["solver"].first.elapsed();
} // TestConvergence::solveNonSym
template< typename value_type,
          typename Basis,
          int gtOrder,
          int imOrder>
template<typename Mat, typename Vec1, typename Vec2>
void
TestConvergence<value_type,Basis,gtOrder, imOrder>::solveSym( int iMesh,
        Mat const& D,
        Vec1& u,
        Vec2 const& F )
{
    timers[iMesh]["solver"].first.restart();

    gmm::iteration iter( 2.0E-10 );
    iter.set_noisy( this->vm()["verbose"].template as<int>() );
    iter.set_maxiter( this->vm()["maxiter"].template as<int>() );
    // incomplete LU with k fill-in and threshold preconditioner.
    // Efficient but could be costly.
    gmm::ildltt_precond<typename sparse_matrix_type::matrix_type> P( D.mat(), 2, 1e-3 );

    // Conjugate gradient
    gmm::cg( D.mat(), u.container(), F, P, iter );

    if ( !iter.converged() )
    {
        std::cout << "mass solver didn't converge" << std::endl;
    }

    //std::cout << "iterations count = " << iter.get_iteration() << std::endl;

    timers[iMesh]["solver"].second = timers[iMesh]["solver"].first.elapsed();
}
template< typename value_type,
          typename Basis,
          int gtOrder,
          int imOrder>
void
TestConvergence<value_type,Basis,gtOrder, imOrder>::exportResults( int iMesh,
        element_0_type const& u,
        element_1_type const& p )

{
    timers[iMesh]["export"].first.restart();

    // -- EXPORT --
    if ( this->vm().count( "export" ) )
    {
#if 0
        // interpolate on the P1 associated space
        typename space_type::template P1Lagrange<0> p1lag( u.functionSpace() );
        typename space_type::template P1Lagrange<0>::p1_type::element_type u_p1 = p1lag( u );

        typename timeset_type::step_ptrtype timeStep = timeSet->step( float( Nstep ) );
        timeStep->setMesh( p1lag.mesh() );
        timeStep->addNodalVector( "u", u_p1.size(), u_p1.begin(), u_p1.end()  );
        //timeStep->addNodalScalar( "p", p.size(), p.begin(), p.end()  );
        exporter->save();
#endif
    } // export

    timers[iMesh]["export"].second = timers[iMesh]["export"].first.elapsed();
} // TestConvergence::export
} // Feel



#if USE_BOOST_TEST
using namespace Feel;
struct lagrange_test_suite : public test_suite
{
    template<int geoDim, int feOrder, int imOrder>
    struct tc
    {

        typedef TestConvergence<double,
                fusion::vector<fem::Lagrange<geoDim, feOrder, Vectorial, Continuous, double>,
                fem::Lagrange<geoDim, feOrder, Scalar, Continuous, double> >,
                gtOrder,
                imOrder> type;
    };
    template<typename AboutType, typename OptionsType>
    lagrange_test_suite( int argc, char** argv, AboutType const& about, OptionsType const& options )
        :
        test_suite( "lagrange basis convergence testsuite" )
    {
        // 2D
        add( BOOST_TEST_CASE( ( typename tc<2,1,2>::type( argc, argv, about, options ) ) ) );
        //add( BOOST_TEST_CASE( (typename tc<2,2,4>::type( argc, argv, about, options )) ) );
        //add( BOOST_TEST_CASE( (typename tc<2,5,10>::type( argc, argv, about, options )) ) );

    }
};

test_suite*
init_unit_test_suite( int argc, char* argv[] )
// main( int argc, char* argv[] )
{
    using namespace Feel;
    typedef double value_type;

    test_suite* test = BOOST_TEST_SUITE( "Elliptic finite element solver convergence test suite" );


    BOOST_MESSAGE( "Convergence test for scalar Lagrange basis functions" );
    test->add( new lagrange_test_suite( argc, argv, makeAbout(), makeOptions() ) );
    //BOOST_MESSAGE("Convergence test for vectorial Lagrange basis functions");
    //test->add( new lagrange_test_suite<Vectorial>( argc, argv, makeAbout(), makeOptions() ) );
    //typedef fem::Lagrange<dim, feOrder, Scalar, Continuous, value_type> basis_type;
    //typedef BoundaryAdaptedPolynomialSet<2, feOrder, Scalar, value_type, Simplex> basis_type;
    return test;
}
#else
int
main( int argc, char* argv[] )
{
    using namespace Feel;

#if 0

    try
    {
        TestConvergence<double, 2, 1, gtOrder, 2> app( argc, argv, makeAbout(), makeOptions() );
        app.run();
    }

    catch ( std::exception& e )
    {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }

    catch ( ... )
    {
        std::cerr << "Exception of unknown type!\n";
    }

#else
    typedef fusion::vector<fem::Lagrange<2, 2, Vectorial, Continuous, double>,
            fem::Lagrange<2, 1, Scalar, Continuous, double> > basis_type;
    TestConvergence<double, basis_type,gtOrder, 2> app( argc, argv, makeAbout(), makeOptions() );
    app.run();
#endif
    return 0;
}
#endif
