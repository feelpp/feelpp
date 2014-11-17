/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-07-07

  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_mixed.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-07-07
 */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE function space testsuite
#include <testsuite/testsuite.hpp>



#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>


#include <feel/feeldiscr/region.hpp>

#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */


namespace Feel
{
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
    ( "shape", po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (hypercube, simplex, ellipsoid)" )
    ( "nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient" )
    ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
    ( "penaldir", Feel::po::value<double>()->default_value( 10 ),
      "penalisation parameter for the weak boundary Dirichlet formulation" )
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
    AboutData about( "test_mixed" ,
                     "test_mixed" ,
                     "0.2",
                     "nD(n=1,2,3) Test_Mixed on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2010 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


/**
 * \class Test_Mixed
 *
 * TestMixed Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
template<int Dim, int Order>
class TestMixed
    :
public Application
{
    typedef Application super;
public:
    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    //! sparse matrix type associated with backend
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    //! sparse matrix type associated with backend (shared_ptr<> type)
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    //! vector type associated with backend
    typedef typename backend_type::vector_type vector_type;
    //! vector type associated with backend (shared_ptr<> type)
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

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
    typedef bases<Lagrange<Order,Vectorial> > v_basis_type;
    typedef bases<Lagrange<Order-1,Scalar> > p_basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, v_basis_type> v_space_type;
    typedef FunctionSpace<mesh_type, p_basis_type> p_space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<v_space_type> v_space_ptrtype;
    typedef boost::shared_ptr<p_space_type> p_space_ptrtype;
    //! an element type of the approximation function space
    typedef typename v_space_type::element_type v_element_type;
    typedef typename p_space_type::element_type p_element_type;

    /**
     * Constructor
     */
    TestMixed()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
    {
    }

    /**
     * run the convergence test
     */
    void run();

private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;
    //! shape of the domain
    std::string shape;

}; // TestMixed

template<int Dim, int Order>
void
TestMixed<Dim,Order>::run()
{
    /**
     * print help if --help is passed to the command line
     */
    /** \code */
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    /** \endcode */

    /**
     * we change to the directory where the results and logs will be
     * stored
     */
    /** \code */
    this->changeRepository( boost::format( "%1%/%2%/P%3%/%4%/h_%5%/" )
                            % this->about().appName()
                            % convex_type::name()
                            % Order
                            % shape
                            % meshSize
                          );
    /** \endcode */
    /**
     * First we create the mesh
     */
    /** \code */
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                                _usenames=true,
                                                _convex=( convex_type::is_hypercube )?"Hypercube":"Simplex",
                                                _shape=shape,
                                                _dim=Dim,
                                                _xmin=-1.,_ymin=-1.,_zmin=-1.,
                                                _h=meshSize ),
                                        _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
    /** \endcode */

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    v_space_ptrtype Xh = v_space_type::New( mesh );
    p_space_ptrtype Yh = p_space_type::New( mesh );
    v_element_type u( Xh, "u" );
    p_element_type p( Yh, "p" );
    /** \endcode */

    u = vf::project( Xh, elements( mesh ), P() );
    p = vf::project( Yh, elements( mesh ), cst( 1. ) );

    vector_ptrtype U = M_backend->newVector( Xh );
    *U = u;
    vector_ptrtype V = M_backend->newVector( Xh );
    vector_ptrtype Q = M_backend->newVector( Yh );
    vector_ptrtype P = M_backend->newVector( Yh );
    *P = p;
    U->close();
    P->close();
    //P->printMatlab( "P.m" );
    //U->printMatlab( "U.m" );

    double meas = integrate( elements( mesh ), cst( 1. ) ).evaluate()( 0, 0 );

    double res =  integrate( elements( mesh ), divv( u )*idv( p ) ).evaluate()( 0,0 ) ;
    BOOST_CHECK_CLOSE( res, Dim*meas, 5e-11 );
    BOOST_TEST_MESSAGE( "[(u,p)] res = " << res << " (must be equal to " << Dim*meas << ")\n" );

    res =  integrate( elements( mesh ), divv( u )*divv( u ) ).evaluate()( 0,0 ) ;
    BOOST_CHECK_CLOSE( res, Dim*Dim*meas, 5e-11 );
    BOOST_TEST_MESSAGE( "[(u,u)] res = " << res << " (must be equal to " << Dim*Dim*meas << ")\n" );

    res =  integrate( elements( mesh ), idv( p )*idv( p ) ).evaluate()( 0,0 ) ;
    BOOST_CHECK_CLOSE( res, meas, 5e-11 );
    BOOST_TEST_MESSAGE( "[(p,p)] res = " << res << " (must be equal to " << meas << ")\n" );

    {
        auto D = M_backend->newMatrix( Xh, Yh );
        BOOST_CHECK_EQUAL( D->size1(), Yh->nLocalDof() );
        BOOST_CHECK_EQUAL( D->size2(), Xh->nLocalDof() );
        form2( _trial=Xh, _test=Yh, _matrix=D, _init=true )= integrate( elements( mesh ), divt( u )*id( p ) );
        D->close();
        //D->printMatlab( "D.m" );
        D->multVector( U, Q );
        //Q->printMatlab( "Q.m" );
        res = inner_product( Q, P );
        BOOST_CHECK_CLOSE( res, Dim*meas, 5e-13 );
        BOOST_TEST_MESSAGE( "[D(p,u)] res = " << res << " (must be equal to " << Dim*meas << ")\n" );

        form2( _trial=Xh, _test=Yh, _matrix=D, _init=true )=
            integrate( elements( mesh ), -grad( p )*idt( u ) );

        D->close();
        //D->printMatlab( "Dg.m" );
        D->multVector( U, Q );
        res = inner_product( Q, P );
        BOOST_CHECK_SMALL( res, 1e-13 );
        BOOST_TEST_MESSAGE( "[Dg(p,u)] res = " << res << " (must be equal to " << 0 << ")\n" );

        form2( _trial=Xh, _test=Yh, _matrix=D, _init=true )=
            integrate( boundaryfaces( mesh ), trans( idt( u ) )*N()*id( p ) );

        D->close();
        //D->printMatlab( "Db.m" );
        D->multVector( U, Q );
        res = inner_product( Q, P );
        BOOST_CHECK_CLOSE( res, Dim*meas, 1e-12 );
        BOOST_TEST_MESSAGE( "[Db(p,u)] res = " << res << " (must be equal to " << Dim*meas << ")\n" );

        vector_ptrtype F( M_backend->newVector( Yh ) );
        BOOST_CHECK_EQUAL( F->localSize(), Yh->nLocalDofWithGhost() );
        form1( _test=Yh, _vector=F, _init=true )= integrate( elements( mesh ), divv( u )*id( p ) );

        F->close();
        //F->printMatlab( "Fu.m" );
        res = inner_product( F, P );
        BOOST_CHECK_CLOSE( res, Dim*meas, 5e-11 );
        BOOST_TEST_MESSAGE( "[Fu(p)] res = " << res << " (must be equal to " << Dim*meas << ")\n" );

    }
    {
        sparse_matrix_ptrtype D( M_backend->newMatrix( Yh, Xh ) );
        BOOST_CHECK_EQUAL( D->size1(), Xh->nLocalDof() );
        BOOST_CHECK_EQUAL( D->size2(), Yh->nLocalDof() );
        form2( _test=Xh, _trial=Yh, _matrix=D, _init=true )= integrate( elements( mesh ), div( u )*idt( p ) );

        D->close();
        //D->printMatlab( "Dt.m" );
        D->multVector( P, V );
        double res = inner_product( V, U );
        BOOST_CHECK_CLOSE( res, Dim*meas, 5e-11 );
        BOOST_TEST_MESSAGE( "[D(u,p)] res = " << res << " (must be equal to " << Dim*meas << ")\n" );


        vector_ptrtype F( M_backend->newVector( Xh ) );
        BOOST_CHECK_EQUAL( F->localSize(), Xh->nLocalDofWithGhost() );
        form1( _test=Xh, _vector=F, _init=true )= integrate( elements( mesh ), div( u )*idv( p ) );

        F->close();
        //F->printMatlab( "Fp.m" );
        res = inner_product( F, U );
        BOOST_CHECK_CLOSE( res, Dim*meas, 5e-11 );
        BOOST_TEST_MESSAGE( "[Fp(u)] res = " << res << " (must be equal to " << Dim*meas << ")\n" );


    }
    {
        sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );
        BOOST_CHECK_EQUAL( D->size1(), Xh->nLocalDof() );
        BOOST_CHECK_EQUAL( D->size2(), Xh->nLocalDof() );
        form2( _trial=Xh, _test=Xh, _matrix=D, _init=true )= integrate( elements( mesh ), divt( u )*div( u ) );
        D->close();
        //D->printMatlab( "divdiv.m" );
        D->multVector( U, V );
        double res = inner_product( V, U );
        BOOST_CHECK_CLOSE( res, Dim*Dim*meas, 5e-11 );
        BOOST_TEST_MESSAGE( "[(u,u)] res = " << res << " (must be equal to " << Dim*Dim*meas << ")\n" );

    }
    {
        sparse_matrix_ptrtype D( M_backend->newMatrix( Yh, Yh ) );
        BOOST_CHECK_EQUAL( D->mapRow().nLocalDofWithGhost(), Yh->nLocalDofWithGhost() );
        BOOST_CHECK_EQUAL( D->mapCol().nLocalDofWithGhost(), Yh->nLocalDofWithGhost() );

        form2( _trial=Yh, _test=Yh, _matrix=D, _init=true )= integrate( elements( mesh ), idt( p )*id( p ) );
        D->close();
        //D->printMatlab( "idid.m" );
        D->multVector( P, Q );
        double res = inner_product( Q, P );
        BOOST_CHECK_CLOSE( res, meas, 5e-11 );
        BOOST_TEST_MESSAGE( "[(p,p)] res = " <<  res << " (must be equal to " << meas << ")\n" );
    }
} // TestMixed::run
}
#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() );

BOOST_AUTO_TEST_SUITE( mixed )

BOOST_AUTO_TEST_CASE( test_mixed1_21 )
{
    BOOST_TEST_MESSAGE( "test_mixed1 (2D,Order 2,)" );
    Feel::TestMixed<2,2> t;
    t.run();
    BOOST_TEST_MESSAGE( "test_mixed1 (2D,Order 2,) done" );
}
#if 1
BOOST_AUTO_TEST_CASE( test_mixed1_23 )
{
    BOOST_TEST_MESSAGE( "test_mixed1 (2D,Order 3,)" );
    Feel::TestMixed<2,3> t;
    t.run();
    BOOST_TEST_MESSAGE( "test_mixed1 (2D,Order 3,) done" );
}

BOOST_AUTO_TEST_CASE( test_mixed2_31 )
{
    BOOST_TEST_MESSAGE( "test_mixed2 (3D,Order 2,)" );
    Feel::TestMixed<3,2> t;
    t.run();
    BOOST_TEST_MESSAGE( "test_mixed2 (3D,Order 2,) done" );
}

BOOST_AUTO_TEST_CASE( test_mixed2_33 )
{
    BOOST_TEST_MESSAGE( "test_mixed2 (3D,Order 3,)" );
    Feel::TestMixed<3,3> t;
    t.run();
    BOOST_TEST_MESSAGE( "test_mixed2 (3D,Order 3,) done" );
}
#endif
BOOST_AUTO_TEST_SUITE_END()


#if 0
int BOOST_TEST_CALL_DECL
main( int argc, char* argv[] )
{
    Feel::Environment env( argc, argv );
    Feel::Assert::setLog( "test_mixed.assert" );
    int ret = ::boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

    return ret;
}
#endif
#else

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    /**
     * intantiate a TestMixed<Dim> class with Dim=2 (e.g. geometric dimension is 2)
     */
    /** \code */
    Feel::TestMixed<2,2> testMixed( argc, argv, Feel::makeAbout(), Feel::makeOptions() );
    /** \encode */

    /**
     * run the application
     */
    /** \code */
    testMixed.run();
    /** \endcode */
}
#endif // USE_BOOST_TEST
