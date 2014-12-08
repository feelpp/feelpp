/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Author(s): Cecile Daversin  <cecile.daversin@lncmi.cnrs.fr>
       Date: 2011-12-07

  Copyright (C) 2011 UJF
  Copyright (C) 2011 CNRS

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
   \file test_hcurl.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2011-12-07
 */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE H_curl approximation
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>
#include <feel/feel.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelpoly/nedelec.hpp>
/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description testHcurlOptions( "test h_curl options" );
    testHcurlOptions.add_options()
        ( "meshes", po::value< std::vector<std::string> >(), "vector containing mesh names" )
        ;
    return testHcurlOptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_hcurl_oneelt" ,
                     "test_hcurl_oneelt" ,
                     "0.1",
                     "Test for h_curl space (ref elt)",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "cecile.daversin@lncmi.cnrs.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

class TestHCurlOneElt
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
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<2,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    //typedef bases<Nedelec<0> > basis_type;
    typedef bases<Nedelec<0,NedelecKind::NED1> > basis_type;
    typedef bases<Lagrange<1,Scalar> > lagrange_basis_s_type; //P1 scalar space
    typedef bases<Lagrange<1,Vectorial> > lagrange_basis_v_type; //P1 vectorial space

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef FunctionSpace<mesh_type, lagrange_basis_s_type> lagrange_space_s_type;
    typedef FunctionSpace<mesh_type, lagrange_basis_v_type> lagrange_space_v_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<lagrange_space_s_type> lagrange_space_s_ptrtype;
    typedef boost::shared_ptr<lagrange_space_v_type> lagrange_space_v_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    TestHCurlOneElt()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        std::cout << "[TestHCurlOneElt]\n";

        this->changeRepository( boost::format( "%1%" ) % this->about().appName() );
    }

    /**
     * run the application
     */
    void shape_functions( std::string one_element_mesh );
    void testProjector( std::string one_element_mesh );


private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! exporter factory
    export_ptrtype exporter;

}; //TestHCurlOneElt

void
TestHCurlOneElt::testProjector(std::string one_element_mesh )
{
    std::string mesh_name = one_element_mesh + ".msh";
    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type,
                                  _filename=mesh_name );

    space_ptrtype Nh = space_type::New( mesh );
    lagrange_space_v_ptrtype Yh_v = lagrange_space_v_type::New( mesh ); //lagrange vectorial space
    lagrange_space_s_ptrtype Yh_s = lagrange_space_s_type::New( mesh ); //lagrange scalar space

    auto E = Py()*unitX() + Px()*unitY();
    auto f = cst(0.); //curld[2d](E) = f

    // curl2D(u1,u2) = d(u2)/dx1 - d(u1)/dx2
    // 2D case : only curlx is initialized to curl2D(u1,u2) (curl is a vector with only one component initialized)

    // L2 projection (Lagrange)
    auto l2_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=L2 ); //l2 vectorial proj
    auto l2_lagS = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=L2 ); //l2 scalar proj
    auto E_pL2_lag = l2_lagV->project( _expr= E );
    auto error_pL2_lag = l2_lagS->project( _expr=curlxv(E_pL2_lag) - f );

    // H1 projection (Lagrange)
    auto h1_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=H1 ); //h1 vectorial proj
    auto E_pH1_lag = h1_lagV->project( _expr= E, _grad_expr=mat<2,2>(cst(0.),cst(1.),cst(1.),cst(0.)) );
    auto error_pH1_lag = l2_lagS->project( _expr=curlxv(E_pH1_lag) - f );

    // HCURL projection (Lagrange)
    auto hcurl_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=HCURL );
    auto E_pHCURL_lag = hcurl_lagV->project( _expr= E /*, _curl_expr=cst(0.)*/ );
    auto error_pHCURL_lag = l2_lagS->project( _expr=curlxv(E_pHCURL_lag) - f );

    // L2 projection (Nedelec)
    auto l2_ned = opProjection( _domainSpace=Nh, _imageSpace=Nh, _type=L2 );
    auto E_pL2_ned = l2_ned->project( _expr= E );
    auto error_pL2_ned = l2_lagS->project( _expr=curlxv(E_pL2_lag) - f );


    BOOST_TEST_MESSAGE("L2 projection [Lagrange]: error[curl(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pL2_lag, error_pL2_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pL2_lag, error_pL2_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("H1 projection [Lagrange]: error[curl(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pH1_lag, error_pH1_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pH1_lag, error_pH1_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("HCURL projection [Lagrange]: error[curl(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pHCURL_lag, error_pHCURL_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pHCURL_lag, error_pHCURL_lag ) ), 1e-13 );

    BOOST_TEST_MESSAGE("L2 projection [NED]: error[curl(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pL2_ned, error_pL2_ned ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pL2_ned, error_pL2_ned ) ), 1e-13 );

    std::string proj_name = "projection";
    export_ptrtype exporter_proj( export_type::New( this->vm(),
                                  ( boost::format( "%1%-%2%-%3%" )
                                    % this->about().appName()
                                    % ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str()
                                    % proj_name ).str() ) );

    exporter_proj->step( 0 )->setMesh( mesh );
    exporter_proj->step( 0 )->add( "proj_L2_E[Lagrange]", E_pL2_lag );
    exporter_proj->step( 0 )->add( "proj_H1_E[Lagrange]", E_pH1_lag );
    exporter_proj->step( 0 )->add( "proj_HCurl_E[Lagrange]", E_pHCURL_lag );
    exporter_proj->step( 0 )->add( "proj_L2_E[NED]", E_pL2_ned );
    exporter_proj->save();


}
void
TestHCurlOneElt::shape_functions( std::string one_element_mesh )
{
    auto mesh_name = one_element_mesh +".msh"; //create the mesh and load it
    fs::path mesh_path( mesh_name );

    mesh_ptrtype oneelement_mesh = loadMesh( _mesh=new mesh_type,
                                             _filename=mesh_name);

    auto refine_level = std::floor(1 - math::log( 0.1 ));
    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type,
                                      _filename=mesh_name,
                                      _refine=( int )refine_level);

    space_ptrtype Xh = space_type::New( oneelement_mesh );

    std::cout << "Family = " << Xh->basis()->familyName() << "\n"
              << "Dim    = " << Xh->basis()->nDim << "\n"
              << "Order  = " << Xh->basis()->nOrder << "\n"
              << "NDof   = " << Xh->nLocalDof() << "\n";

    element_type U_ref( Xh, "U" );
    element_type V_ref( Xh, "V" );

    // To store the shape functions
    // 0 : hypothenuse edge, 1 : vertical edge, 2 : horizontal edge
    std::vector<element_type> u_vec( Xh->nLocalDof() );

    std::string shape_name = "shape_functions";
    export_ptrtype exporter_shape( export_type::New( this->vm(),
                                   ( boost::format( "%1%-%2%-%3%" )
                                     % this->about().appName()
                                     % mesh_path.stem()
                                     % shape_name ).str() ) );

    exporter_shape->step( 0 )->setMesh( mesh );

    for ( size_type i = 0; i < Xh->nLocalDof(); ++i )
    {
        // U_ref corresponds to shape function (on reference element)
        U_ref.zero();
        U_ref( i ) = 1;

        u_vec[i] = U_ref;

        std::ostringstream ostr;
        ostr <<  mesh_path.stem() << "-" << Xh->basis()->familyName() << "-" << i;
        exporter_shape->step( 0 )->add( ostr.str(), U_ref );
    }

    exporter_shape->save();

    auto F = M_backend->newVector( Xh );

    //// *********************** Check  alpha_i(N_j) evaluations  on reference element (with idv keyword) ********////
    int check_size = Xh->nLocalDof()*Xh->nLocalDof();
    std::vector<double> checkidv( check_size );
    std::vector<double> checkform1( check_size );
    std::vector<std::string> edges = boost::assign::list_of( "hypo" )( "vert" )( "hor" );

    std::vector<double> checkStokesidv( 2*Xh->nLocalDof() );
    std::vector<double> checkStokesform1( 2*Xh->nLocalDof() );

    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        int edgeid = 0;
        BOOST_FOREACH( std::string edge, edges )
        {
            auto int_u_t = integrate( markedfaces( oneelement_mesh, edge ), trans( T() )*idv( u_vec[i] ) ).evaluate()( 0,0 );

            if ( edgeid == i )
                BOOST_CHECK_CLOSE( int_u_t, 1, 1e-13 );

            else
                BOOST_CHECK_SMALL( int_u_t, 1e-13 );

            checkidv[3*i+edgeid] = int_u_t;

            form1( _test=Xh, _vector=F, _init=true ) = integrate( markedfaces( oneelement_mesh, edge ), trans( T() )*id( V_ref ) );
            auto form_v_t = inner_product( u_vec[i], *F );

            if ( edgeid == i )
                BOOST_CHECK_CLOSE( form_v_t, 1, 1e-13 );
            else
                BOOST_CHECK_SMALL( form_v_t, 1e-13 );

            checkform1[3*i+edgeid] = form_v_t;

            ++edgeid;
        }

        // curl2D(u1,u2) = d(u2)/dx1 - d(u1)/dx2
        // 2D case : only curlx is initialized to curl2D(u1,u2) (curl is a vector with only one component initialized)
        auto int_curlxv = integrate( elements( oneelement_mesh ), curlxv( u_vec[i] ) ).evaluate()( 0,0 );
        auto int_vt = integrate( boundaryfaces( oneelement_mesh ), trans( T() )*idv( u_vec[i] ) ).evaluate()( 0,0 );

        BOOST_CHECK_CLOSE( int_vt, 1, 1e-13 );
        BOOST_CHECK_CLOSE( int_curlxv, int_vt, 1e-13 );
        checkStokesidv[i] = int_curlxv;
        checkStokesidv[i + Xh->nLocalDof()] = int_vt;

        // curl2D(u1,u2) = d(u2)/dx1 - d(u1)/dx2
        // 2D case : only curlx is initialized to curl2D(u1,u2) (curl is a vector with only one component initialized)
        form1( _test=Xh, _vector=F, _init=true ) = integrate( elements( oneelement_mesh ), curlx( V_ref ) );
        auto form_curlxv = inner_product( u_vec[i], *F );
        form1( _test=Xh, _vector=F, _init=true ) = integrate( boundaryfaces( oneelement_mesh ),trans( T() )*id( V_ref ) );
        auto form_vt = inner_product( u_vec[i], *F );

        BOOST_CHECK_CLOSE( form_vt, 1, 1e-13 );
        BOOST_CHECK_CLOSE( form_curlxv, form_vt, 1e-13 );
        checkStokesform1[i] = form_curlxv;
        checkStokesform1[i + Xh->nLocalDof()] = form_vt;
    }

    BOOST_TEST_MESSAGE( " ********** Values of alpha_i (N_j ) = delta_{i,j}  ********** \n"
                        << "\n"
                        << " ********** Using idv keyword ********************* "
                        << "\n" );

    for ( int i = 0; i < 3; ++i )
    {
        int edgeid = 0;
        BOOST_FOREACH( std::string edge, edges )
        {
            BOOST_TEST_MESSAGE( " *** dof N_"<< i << " (associated with " << edge << " edge) *** \n"
                                << "alpha_"<< edgeid << "(N_"<<i<<") = " << checkidv[3*i+edgeid] << "\n" );
            ++edgeid;
        }
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    BOOST_TEST_MESSAGE( " ********** Values of alpha_i (N_j ) = delta_{i,j}  ********** \n"
                        << "\n"
                        << " ********** Using form1 keyword ********************* "
                        << "\n" );

    for ( int i = 0; i < 3; ++i )
    {
        int edgeid = 0;
        BOOST_FOREACH( std::string edge, edges )
        {
            BOOST_TEST_MESSAGE( " *** dof N_"<< i << " (associated with " << edge << " edge) *** \n"
                                << "alpha_"<< edgeid << "(N_"<<i<<") = " << checkform1[3*i+edgeid] << "\n" );
            ++edgeid;
        }
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    BOOST_TEST_MESSAGE( " ********** Stokes theorem  ********** \n"
                        << "\n"
                        << " ********** Using idv keyword ********************* "
                        << "\n" );
    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        BOOST_TEST_MESSAGE("int(Omega) curl(N_" << i << ") = " << checkStokesidv[i] );
        BOOST_TEST_MESSAGE("int(boundary) N_" << i << ".t = " << checkStokesidv[i+Xh->nLocalDof()] );
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    BOOST_TEST_MESSAGE( " ********** Stokes theorem  ********** \n"
                        << "\n"
                        << " ********** Using form1 keyword ********************* "
                        << "\n" );
    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        BOOST_TEST_MESSAGE("int(Omega) curl(N_" << i << ") = " << checkStokesform1[i] );
        BOOST_TEST_MESSAGE("int(boundary) N_" << i << ".t = " << checkStokesform1[i+Xh->nLocalDof()] );
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    //// ************************************************************************************ ////
}

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( HCURL )

BOOST_AUTO_TEST_CASE( test_hcurl_N0_shape_func )
{
    TestHCurlOneElt t;
    std::vector<std::string> mygeoms = vsoption(_name="meshes");//.template as< std::vector<std::string> >();
    for(std::string geo : mygeoms)
        {
            std::cout << "*** shape functions on " << geo << " *** \n";
            t.shape_functions( geo );
            std::cout << "*** projections on " << geo << " *** \n";
            t.testProjector( geo );
        }
}

BOOST_AUTO_TEST_SUITE_END()
#else

int
main( int argc, char* argv[] )
{
    Feel::Environment env( argc,argv,
                           makeAbout(), makeOptions() );

    Feel::TestHCurlOneElt app_hcurl;

    std::vector<std::string> mygeoms = option(_name="meshes").template as< std::vector<std::string> >();
    for(std::string geo : mygeoms)
        {
            app_hcurl.shape_functions(geo);
            app_hcurl.testProjector(geo);

        }
}

#endif
