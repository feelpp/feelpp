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
   \file test_hdiv3D.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2014-01-29
 */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE H_curl3D approximation
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>
#include <feel/feel.hpp>
#include <feel/feelvf/print.hpp>
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
    po::options_description testHcurl3DOptions( "test h_div 3D options" );
    testHcurl3DOptions.add_options()
        ( "meshes", po::value< std::vector<std::string> >(), "vector containing mesh names" )
        ;
    return testHcurl3DOptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_hcurl3D_oneelt" ,
                     "test_hcurl3D_oneelt" ,
                     "0.1",
                     "Test for h_curl space",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

using namespace Feel;

class TestHCurl3DOneElt
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
    typedef Simplex<3,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    typedef Mesh< Simplex<1,1,3> > submesh1d_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef boost::shared_ptr<submesh1d_type> submesh1d_ptrtype;

    //! the basis type of our approximation space
    typedef bases<Nedelec<0,NedelecKind::NED1> > basis_type;
    typedef bases<Lagrange<1,Vectorial> > lagrange_basis_v_type; //P1 vectorial space
    typedef bases<Lagrange<1,Scalar> > lagrange_basis_s_type; //P1 scalar space
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
    TestHCurl3DOneElt()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        std::cout << "[TestHCurl3DOneElt]\n";

        this->changeRepository( boost::format( "%1%" ) % this->about().appName());
    }

    /**
     * run the application
     */
    void shape_functions( std::string one_element_mesh );
    void testProjector(std::string one_element_mesh );

private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! exporter factory
    export_ptrtype exporter;

}; //TestHDiv

void
TestHCurl3DOneElt::testProjector(std::string one_element_mesh )
{
    std::string mesh_name = one_element_mesh + ".msh";
    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type,
                                  _filename=mesh_name );

    space_ptrtype Nh = space_type::New( mesh );
    lagrange_space_v_ptrtype Yh_v = lagrange_space_v_type::New( mesh ); //lagrange vectorial space
    lagrange_space_s_ptrtype Yh_s = lagrange_space_s_type::New( mesh ); //lagrange scalar space

    auto E = (Py()+Pz())*unitX() + (Px()+Pz())*unitY() + (Px()+Py())*unitZ(); //E : (y+z; x+z; x+y)
    auto f = cst(0.)*unitX() + cst(0.)*unitY() + cst(0.)*unitZ(); //curld[3d](E) = f

    // L2 projection (Lagrange)
    auto l2_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=L2 ); //l2 vectorial proj
    auto E_pL2_lag = l2_lagV->project( _expr= E );
    auto error_pL2_lag = l2_lagV->project( _expr=curlv(E_pL2_lag) - f );

    // H1 projection (Lagrange)
    auto h1_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=H1 ); //h1 vectorial proj
    auto E_pH1_lag = h1_lagV->project( _expr= E, _grad_expr=mat<3,3>(cst(0.),cst(1.),cst(1.),
                                                                     cst(1.),cst(0.),cst(1.),
                                                                     cst(1.),cst(1.),cst(0.) ) );
    auto error_pH1_lag = l2_lagV->project( _expr=curlv(E_pH1_lag) - f );

    // HCURL projection (Lagrange)
    auto hcurl_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=HCURL );
    auto E_pHCURL_lag = hcurl_lagV->project( _expr= E /*, _curl_expr=f*/ );
    auto error_pHCURL_lag = l2_lagV->project( _expr=curlv(E_pHCURL_lag) - f );

    // L2 projection (Nedelec)
    auto l2_ned = opProjection( _domainSpace=Nh, _imageSpace=Nh, _type=L2 );
    auto E_pL2_ned = l2_ned->project( _expr= E );
    auto error_pL2_ned = l2_lagV->project( _expr=curlv(E_pL2_lag) - f );

    BOOST_TEST_MESSAGE("L2 projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagV->energy( error_pL2_lag, error_pL2_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagV->energy( error_pL2_lag, error_pL2_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("H1 projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagV->energy( error_pH1_lag, error_pH1_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagV->energy( error_pH1_lag, error_pH1_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("HCURL projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagV->energy( error_pHCURL_lag, error_pHCURL_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagV->energy( error_pHCURL_lag, error_pHCURL_lag ) ), 1e-13 );

    BOOST_TEST_MESSAGE("L2 projection [NED]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagV->energy( error_pL2_ned, error_pL2_ned ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagV->energy( error_pL2_ned, error_pL2_ned ) ), 1e-13 );
#if 0
    BOOST_TEST_MESSAGE("H1 projection [NED]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagV->energy( error_pH1_ned, error_pH1_ned ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagV->energy( error_pH1_ned, error_pH1_ned ) ), 1e-13 );
    BOOST_TEST_MESSAGE("HCURL projection [NED]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagV->energy( error_pHCURL_ned, error_pHCURL_ned ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagV->energy( error_pHCURL_ned, error_pHCURL_ned ) ), 1e-13 );
#endif

    std::string proj_name = "projection";
    export_ptrtype exporter_proj( export_type::New( this->vm(),
                                  ( boost::format( "%1%-%2%-%3%" )
                                    % this->about().appName()
                                    % ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str()
                                    % proj_name ).str() ) );

    exporter_proj->step( 0 )->setMesh( mesh );
    exporter_proj->step( 0 )->add( "proj_L2_E[Lagrange]", E_pL2_lag );
    exporter_proj->step( 0 )->add( "proj_H1_E[Lagrange]", E_pH1_lag );
    exporter_proj->step( 0 )->add( "proj_HDiv_E[Lagrange]", E_pHCURL_lag );
    exporter_proj->step( 0 )->add( "proj_L2_E[NED]", E_pL2_ned );
    exporter_proj->save();
}

void
TestHCurl3DOneElt::shape_functions( std::string one_element_mesh )
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
                                     % mesh_path.stem().string()
                                     % shape_name ).str() ) );

    exporter_shape->step( 0 )->setMesh( mesh );

    for ( size_type i = 0; i < Xh->nLocalDof(); ++i )
    {
        // U_ref corresponds to shape function (on reference element)
        U_ref.zero();
        U_ref( i ) = 1;

        u_vec[i] = U_ref;

        std::ostringstream ostr;
        ostr << mesh_path.stem().string() << "-" << Xh->basis()->familyName() << "-" << i;
        exporter_shape->step( 0 )->add( ostr.str(), U_ref );
    }

    exporter_shape->save();

    auto F = M_backend->newVector( Xh );

    //// *********************** Check  alpha_i(N_j) evaluations  on reference element (with idv keyword) ********////
    int check_size = Xh->nLocalDof()*Xh->nLocalDof();
    std::vector<double> checkidv( check_size );
    std::vector<double> checkform1( check_size );
    std::vector<std::string> edges = {"zAxis","yAxis","yzAxis","xyAxis","xzAxis","xAxis"};

    submesh1d_ptrtype edgeMesh( new submesh1d_type );
    edgeMesh = createSubmesh(oneelement_mesh, boundaryedges(oneelement_mesh) ); //submesh of edges

    double div1sqrt2 = 1/sqrt(2.);

    // Tangents (normalized) on ref element
    auto t0 = vec(cst(0.),cst(0.),cst(-1.));
    auto t1 = vec(cst(0.),cst(1.),cst(0.));
    auto t2 = vec(cst(0.),cst(-div1sqrt2),cst(div1sqrt2));
    auto t3 = vec(cst(div1sqrt2),cst(-div1sqrt2),cst(0.));
    auto t4 = vec(cst(div1sqrt2),cst(0.),cst(-div1sqrt2));
    auto t5 = vec(cst(1.),cst(0.),cst(0.));
    std::vector< decltype(t0) > tangentRef{t0,t1,t2,t3,t4,t5};

    // Jacobian of geometrical transforms
    std::string jac;
    std::string jac_name = "jac_"+ mesh_path.stem().string();
    if(mesh_path.stem().string() == "one-elt-ref-3d" || mesh_path.stem().string() == "one-elt-real-homo-3d" )
        jac = "{1,0,0,0,1,0,0,0,1}:x:y:z";
    else if(mesh_path.stem().string() == "one-elt-real-rotx" )
        jac = "{1,0,0,0,0,-1,0,1,0}:x:y:z";
    else if(mesh_path.stem().string() == "one-elt-real-roty" )
        jac = "{0,0,1,0,1,0,-1,0,0}:x:y:z";
    else if(mesh_path.stem().string() == "one-elt-real-rotz" )
        jac = "{0,-1,0,1,0,0,0,0,1}:x:y:z";

    std::vector<double> checkStokesidv( 2*Xh->nLocalDof() );
    std::vector<double> checkStokesform1( 2*Xh->nLocalDof() );

    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        int edgeid = 0;
        BOOST_FOREACH( std::string edge, edges )
        {
            CHECK( oneelement_mesh->hasMarkers({edge}) );
            auto int_u_t = integrate( markedelements(edgeMesh, edge),
                                      trans(expr<3,3>(jac,jac_name)*tangentRef[edgeid])*idv( u_vec[i] ) ).evaluate()(0,0);

            if ( edgeid == i )
                BOOST_CHECK_CLOSE( int_u_t, 1, 1e-13 );
            else
                BOOST_CHECK_SMALL( int_u_t, 1e-13 );

            checkidv[Xh->nLocalDof()*i+edgeid] = int_u_t;

            form1( _test=Xh, _vector=F, _init=true ) = integrate( markedelements(edgeMesh, edge),
                                                                  trans(expr<3,3>(jac,jac_name)*tangentRef[edgeid])*id( V_ref ) );
            auto form_v_t = inner_product( u_vec[i], *F );

            if ( edgeid == i )
                BOOST_CHECK_CLOSE( form_v_t, 1, 1e-13 );
            else
                BOOST_CHECK_SMALL( form_v_t, 1e-13 );

            checkform1[Xh->nLocalDof()*i+edgeid] = form_v_t;

            ++edgeid;
        }

        //Stokes Theorem for curl :
        auto int_curlv = integrate( elements( oneelement_mesh ), curlv( u_vec[i] ) ).evaluate()(0,0);
        int_curlv += integrate( elements( oneelement_mesh ), curlv( u_vec[i] ) ).evaluate()(1,0);
        int_curlv += integrate( elements( oneelement_mesh ), curlv( u_vec[i] ) ).evaluate()(2,0);

        form1( _test=Xh, _vector=F, _init=true ) = integrate( elements( oneelement_mesh ), trans(one())*curl( V_ref ) );
        auto form_curlv = inner_product( u_vec[i], *F );

        auto int_vn = integrate(boundaryfaces(oneelement_mesh), -cross( idv(u_vec[i]),N() )).evaluate()(0,0);
        int_vn += integrate(boundaryfaces(oneelement_mesh), -cross( idv(u_vec[i]),N() )).evaluate()(1,0);
        int_vn += integrate(boundaryfaces(oneelement_mesh), -cross( idv(u_vec[i]),N() )).evaluate()(2,0);

        form1( _test=Xh, _vector=F, _init=true ) = integrate( boundaryfaces( oneelement_mesh ),
                                                              -trans(one())*cross( idv(u_vec[i]),N() ) );
        auto form_vn = inner_product( u_vec[i], *F );

        BOOST_CHECK_SMALL( int_curlv - int_vn, 1e-13 );
        checkStokesidv[i] = int_curlv;
        checkStokesidv[i + Xh->nLocalDof()] = int_vn;

        BOOST_CHECK_SMALL( form_curlv - form_vn, 1e-13 );
        checkStokesform1[i] = form_curlv;
        checkStokesform1[i + Xh->nLocalDof()] = form_vn;
    }

    BOOST_TEST_MESSAGE( " ********** Values of alpha_i (N_j ) = delta_{i,j}  ********** \n"
                        << "\n"
                        << " ********** Using idv keyword ********************* "
                        << "\n" );

    for ( int i = 0; i < Xh->nLocalDof(); ++i )
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

    for ( int i = 0; i < Xh->nLocalDof(); ++i )
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
        BOOST_TEST_MESSAGE("int(boundary) N_" << i << "^n = " << checkStokesidv[i+Xh->nLocalDof()] );
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    BOOST_TEST_MESSAGE( " ********** Stokes theorem  ********** \n"
                        << "\n"
                        << " ********** Using form1 keyword ********************* "
                        << "\n" );
    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        BOOST_TEST_MESSAGE("int(Omega) curl(N_" << i << ") = " << checkStokesform1[i] );
        BOOST_TEST_MESSAGE("int(boundary) N_" << i << "^n = " << checkStokesform1[i+Xh->nLocalDof()] );
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    //// ************************************************************************************ ////
}

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( HCURL3D_oneelt )

BOOST_AUTO_TEST_CASE( test_hcurl3D_N0 )
{
    TestHCurl3DOneElt t;
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
    Environment env( argc,argv,
                           makeAbout(), makeOptions() );
    TestHCurl3DOneElt app_hcurl;
    std::vector<std::string> mygeoms = option(_name="meshes").template as< std::vector<std::string> >();
    for(std::string geo : mygeoms)
        {
            app_hcurl.shape_functions(geo);
            //app_hcurl.testProjector(geo);

        }
}

#endif
