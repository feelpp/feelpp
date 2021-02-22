/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE H_curl approximation
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/pchv.hpp>

#include <feel/feelfilters/loadmesh.hpp>

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
    po::options_description testhcurloptions( "test h_curl options" );
    testhcurloptions.add_options()
        //( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        //( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
        //( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
        //( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
        ( "tol", po::value<double>()->default_value( 2.5e-1 ), "tolerance for the solution" )
        ( "u_exact_2d", po::value<std::string>()->default_value("{1-y^2,1-x^2}:x:y"), "exact solution" )
        ( "f_exact_2d", po::value<std::string>()->default_value("{3-y^2,3-x^2}:x:y"), "exact rhs" )
        ( "u_exact_3d", po::value<std::string>()->default_value("{(1-y^2)*(1-z^2),(1-x^2)*(1-z^2),(1-x^2)*(1-y^2)}:x:y:z"), "exact solution" )
        ( "f_exact_3d", po::value<std::string>()->default_value("{y^2*z^2-3*y^2-3*z^2+5,x^2*z^2-3*x^2-3*z^2+5,x^2*y^2-3*x^2-3*y^2+5}:x:y:z"), "exact rhs" )
    ;
    return testhcurloptions.add( Feel::feel_options() );
}
#if 0
inline
AboutData
makeAbout()
{
    AboutData about( "test_hcurl" ,
                     "test_hcurl" ,
                     "0.1",
                     "Test for h_curl space (Dim=2 Order=1)",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "cecile.daversin@lncmi.cnrs.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}
#endif
using namespace Feel;

template<int Dim>
class TestHCurl
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
    typedef typename std::shared_ptr<backend_type> backend_ptrtype ;
    //! sparse matrix type associated with backend
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    //! sparse matrix type associated with backend (shared_ptr<> type)
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    //! vector type associated with backend
    typedef typename backend_type::vector_type vector_type;
    //! vector type associated with backend (shared_ptr<> type)
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

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
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef std::shared_ptr<lagrange_space_s_type> lagrange_space_s_ptrtype;
    typedef std::shared_ptr<lagrange_space_v_type> lagrange_space_v_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef std::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    TestHCurl()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        meshSize( doption("gmsh.hsize") ),
        M_penaldir(1e5)
    {
        Feel::cout << "[TestHCurl]\n";
        Feel::cout << "hsize = " << meshSize << std::endl;

        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % doption("gmsh.hsize") );

        mesh = loadMesh(_mesh = new mesh_type);
        M_Xh = space_type::New(mesh);

        M_a = M_backend->newMatrix(_test=M_Xh,_trial=M_Xh);
        M_f = M_backend->newVector(_test=M_Xh);
    }

    /**
     * run the application
     */
    //void twoElementsMesh();
    //void eightElementsMesh();
    void testProjector();
    void exampleProblem1();

private:
    mesh_ptrtype mesh;
    space_ptrtype M_Xh;

    //! linear algebra backend
    backend_ptrtype M_backend;
    vector_ptrtype M_f;
    sparse_matrix_ptrtype M_a;

    //! mesh characteristic size
    double meshSize;

    double M_penaldir;

}; //TestHCurl

// Resolve problem curl(curl(u)) + u = f with cross_prod(u,n) = 0 on boundary
template<int Dim>
void
TestHCurl<Dim>::exampleProblem1()
{
    // M_Xh : space build with Nedelec elements
    auto Vh = Pchv<1>( mesh );
    auto u = M_Xh->element();
    auto phi = M_Xh->element();

    //auto u_exact = expr<2,1>( "{1-y*y,1-x*x}:x:y" );
    //auto v = Vh->element();

    auto a = form2( _test=M_Xh, _trial=M_Xh, _matrix=M_a );
    auto l = form1( _test=M_Xh, _vector=M_f );
    //variationnal formulation : curl(curl(u)) + u = f

    std::string u_opt = "u_exact_"+std::to_string(Dim)+"d";
    std::string f_opt = "f_exact_"+std::to_string(Dim)+"d";
    auto u_exact = expr<Dim,1>( soption(u_opt) );
    auto f = expr<Dim,1>( soption(f_opt) );

    l = integrate( _range=elements(mesh), _expr=inner(f,id(phi)) );
    //Dirichlet bc on weak form
    l += integrate( _range=boundaryfaces(mesh),
                    _expr=- inner(curl(phi),cross(u_exact,N()))
                    + M_penaldir*trans( cross(u_exact,N()) )*cross(id(phi),N())/hFace() );

    a = integrate(_range=elements(mesh),
                  _expr=inner(curlt(u),curl(phi)) + inner(idt(u),id(phi)) );
    //a += on( _range=boundaryfaces( mesh ),_element=u, _rhs=l, _expr=cst(0.) );

    //Dirichlet bc on weak form
    a += integrate(_range=boundaryfaces(mesh),
                   _expr=-inner(curlt(u),cross(id(phi),N()))
                   - inner(curl(phi),cross(idt(u),N()))
                   + M_penaldir*trans( cross(idt(u),N()) )*cross(id(phi),N())/hFace() );

    a.solve( _solution=u, _rhs=l, _rebuild=true);

    // L2 projection of solution u
    auto u_L2proj = M_Xh->element();
    auto phi_L2proj = M_Xh->element();

    auto a_L2proj = form2( _test=M_Xh, _trial=M_Xh );
    a_L2proj = integrate( _range=elements(mesh), _expr = trans(idt(u_L2proj))*id(phi_L2proj) );
    auto f_L2proj = form1( _test=M_Xh );
    f_L2proj = integrate( _range=elements(mesh), _expr = trans(idv(u))*id(phi_L2proj) );
    a_L2proj.solve( _solution=u_L2proj, _rhs=f_L2proj, _rebuild=true);

    double errorU = normL2(_range=elements(mesh), _expr= idv(u)-u_exact);
    auto errorProj = normL2(_range=elements(mesh), _expr=idv(u)-idv(u_L2proj) );

    BOOST_CHECK_SMALL( errorProj, 1e-13 );
    BOOST_CHECK_SMALL( errorU, doption("tol") );

#if 0
    auto l2proj = opProjection( _domainSpace=M_Xh, _imageSpace=M_Xh, _type=L2 );
    auto error_hcurl = l2proj->project( _expr=trans( u_exact - idv(u) ) );
#endif

    std::string pro1_name = "problem1";
    std::string exporterName = ( boost::format( "%1%-%2%-%3%" )
                                 % this->about().appName()
                                 % ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str()
                                 % pro1_name ).str();
    auto exporter_pro1 = exporter(_mesh=mesh,_name=exporterName );
    exporter_pro1->addRegions();
    exporter_pro1->add( "u", u );
    exporter_pro1->add( "proj_L2_u", u_L2proj );
    exporter_pro1->add( "u_exact", u_exact );
    exporter_pro1->save();

}

#if 0
void
TestHCurl::testProjector()
{
    mesh_ptrtype mesh = loadMesh(_mesh = new mesh_type);
    auto Nh = Ned1h<0>( mesh );
    lagrange_space_v_ptrtype Yh_v = lagrange_space_v_type::New( mesh ); //lagrange vectorial space
    lagrange_space_s_ptrtype Yh_s = lagrange_space_s_type::New( mesh ); //lagrange scalar space

    auto E = Py()*unitX() + Px()*unitY();
    auto f = cst(0.); //curld[2d](E) = f

    // curl2D(u1,u2) = d(u2)/dx1 - d(u1)/dx2
    // 2D case : only curlx is initialized to curl2D(u1,u2) (curl is a vector with only one component initialized)

    // L2 projection (Lagrange)
    auto l2_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=L2 ); //l2 vectorial proj
    auto l2_lagS = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=L2 ); //l2 scalar proj
    auto E_pL2_lag = l2_lagV->project( _expr= trans(E) );
    auto error_pL2_lag = l2_lagS->project( _expr=curlxv(E_pL2_lag) - f );

    // L2 projection (Nedelec)
    auto l2_ned = opProjection( _domainSpace=Nh, _imageSpace=Nh, _type=L2 );
    auto E_pL2_ned = l2_ned->project( _expr= trans(E) );
    auto error_pL2_ned = l2_lagS->project( _expr=curlxv(E_pL2_lag) - f );

    // H1 projection (Lagrange)
    auto h1_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=H1 ); //h1 vectorial proj
    auto h1_lagS = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=H1 ); //h1 scalar proj
    auto E_pH1_lag = h1_lagV->project( _expr= trans(E), _grad_expr=mat<2,2>(cst(0.),cst(1.),cst(1.),cst(0.)) );
    auto error_pH1_lag = l2_lagS->project( _expr=curlxv(E_pH1_lag) - f );

    // H1 projection (Nedelec)
    auto h1_ned = opProjection( _domainSpace=Nh, _imageSpace=Nh, _type=H1 ); //h1 vectorial proj
    auto E_pH1_ned = h1_ned->project( _expr= trans(E), _grad_expr=mat<2,2>(cst(0.),cst(1.),cst(1.),cst(0.)) );
    auto error_pH1_ned = l2_lagS->project( _expr=curlxv(E_pH1_ned) - f );

    // HCURL projection (Lagrange)
    auto hcurl_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=HCURL );
    auto hcurl_lagS = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=HCURL );
    auto E_pHCURL_lag = hcurl_lagV->project( _expr= trans(E), _div_expr=cst(0.) );
    auto error_pHCURL_lag = l2_lagS->project( _expr=curlxv(E_pHCURL_lag) - f );

    // HCURL projection (Nedelec)
    auto hcurl = opProjection( _domainSpace=Nh, _imageSpace=Nh, _type=HCURL ); //hdiv proj (RT elts)
    auto E_pHCURL_ned = hcurl->project( _expr= trans(E), _div_expr=cst(0.) );
    auto error_pHCURL_ned = l2_lagS->project( _expr=curlxv(E_pHCURL_ned) - f );

    BOOST_TEST_MESSAGE("L2 projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pL2_lag, error_pL2_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pL2_lag, error_pL2_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("L2 projection [NED]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pL2_ned, error_pL2_ned ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pL2_ned, error_pL2_ned ) ), 1e-13 );
    BOOST_TEST_MESSAGE("H1 projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pH1_lag, error_pH1_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pH1_lag, error_pH1_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("H1 projection [NED]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pH1_ned, error_pH1_ned ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pH1_ned, error_pH1_ned ) ), 1e-13 );
    BOOST_TEST_MESSAGE("HDIV projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pHCURL_lag, error_pHCURL_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pHCURL_lag, error_pHCURL_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("HDIV projection [NED]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pHCURL_ned, error_pHCURL_ned ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pHCURL_ned, error_pHCURL_ned ) ), 1e-13 );

    std::string proj_name = "projection";
    std::string exporterName = ( boost::format( "%1%-%2%-%3%" )
                                 % this->about().appName()
                                 % ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str()
                                 % proj_name ).str();
    auto exporter_proj = exporter(_mesh=mesh,_name=exporterName );
    exporter_proj->step( 0 )->add( "proj_L2_E[Lagrange]", E_pL2_lag );
    exporter_proj->step( 0 )->add( "proj_L2_E[NED]", E_pL2_ned );
    exporter_proj->step( 0 )->add( "proj_H1_E[Lagrange]", E_pH1_lag );
    exporter_proj->step( 0 )->add( "proj_H1_E[NED]", E_pH1_ned );
    exporter_proj->step( 0 )->add( "proj_HDiv_E[Lagrange]", E_pHCURL_lag );
    exporter_proj->step( 0 )->add( "proj_HDiv_E[NED]", E_pHCURL_ned );
    exporter_proj->save();
}
#endif
}


FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault("test_hcurl"), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( space )

#if 0
BOOST_AUTO_TEST_CASE( test_hcurl_projection )
{
    BOOST_TEST_MESSAGE( "test_hcurl_N0 : projection on one real element" );
    Feel::TestHCurl t;
    t.testProjector();
    BOOST_TEST_MESSAGE( "test_hcurl_N0 : projection on one real element done" );
}
#endif

BOOST_AUTO_TEST_CASE( test_hcurl_example_2d )
{
    BOOST_TEST_MESSAGE( "test_hcurl on example 1" );
    Feel::TestHCurl<2> t;
    t.exampleProblem1();
    BOOST_TEST_MESSAGE( "test_hcurl_N0 on example 1 done" );
}

BOOST_AUTO_TEST_CASE( test_hcurl_example_3d )
{
    BOOST_TEST_MESSAGE( "test_hcurl on example 1" );
    Feel::TestHCurl<3> t;
    t.exampleProblem1();
    BOOST_TEST_MESSAGE( "test_hcurl_N0 on example 1 done" );
}

BOOST_AUTO_TEST_SUITE_END()
