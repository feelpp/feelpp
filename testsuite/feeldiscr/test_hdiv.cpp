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
   \file test_hdiv.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2014-01-29
 */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE H_div approximation
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <feel/feel.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

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
    po::options_description testhdivoptions( "test h_div options" );
    testhdivoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
    ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
    ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
    ;
    return testhdivoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_hdiv" ,
                     "test_hdiv" ,
                     "0.1",
                     "Test for h_div space",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "cecile.daversin@lncmi.cnrs.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

using namespace Feel;

class TestHDiv
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
    typedef bases<RaviartThomas<0> > basis_type;
    typedef bases<Lagrange<1,Vectorial> > lagrange_basis_v_type; //P1 vectorial space
    typedef bases<Lagrange<1,Scalar> > lagrange_basis_s_type; //P1 scalar space
    typedef bases< RaviartThomas<0>, Lagrange<1,Scalar> > prod_basis_type; //For Darcy : (u,p) (\in H_div x L2)
    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef FunctionSpace<mesh_type, lagrange_basis_s_type> lagrange_space_s_type;
    typedef FunctionSpace<mesh_type, lagrange_basis_v_type> lagrange_space_v_type;
    typedef FunctionSpace<mesh_type, prod_basis_type> prod_space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<lagrange_space_s_type> lagrange_space_s_ptrtype;
    typedef boost::shared_ptr<lagrange_space_v_type> lagrange_space_v_ptrtype;
    typedef boost::shared_ptr<prod_space_type> prod_space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;
    typedef typename prod_space_type::element_type prod_element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    TestHDiv()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        std::cout << "[TestHDiv]\n";

        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % this->vm()["hsize"].as<double>()
                              );
    }

    /**
     * run the application
     */
    inline double hSize(){return meshSize;}
    //void shape_functions( gmsh_ptrtype ( *one_element_mesh )( double ));
    void testProjector();
    void exampleProblem1();

private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! exporter factory
    export_ptrtype exporter;

}; //TestHDiv

void
TestHDiv::exampleProblem1()
{
    mesh_ptrtype mesh = loadMesh(_mesh = new mesh_type);

    //auto K = ones<2,2>(); // Hydraulic conductivity tensor
    auto K = mat<2,2>(cst(2.),cst(1.),cst(1.),cst(2.));
    //auto Lambda = ones<2,2>(); // Hydraulic resistivity tensor
    auto Lambda = (1.0/3.0)*mat<2,2>(cst(2.),cst(-1.),cst(-1.),cst(2.)); // Lambda = inv(K)
    auto f = Px()+Py(); // int_omega f = 0
    auto epsilon = 1e-7;

    // ****** Primal formulation - with Lagrange ******
    lagrange_space_s_ptrtype Xh = lagrange_space_s_type::New( mesh );
    lagrange_space_v_ptrtype Xhvec = lagrange_space_v_type::New( mesh );
    auto p_l = Xh->element( "p" );
    auto u_l = Xhvec->element( "u" );
    auto q_l = Xh->element( "q" );

    auto F_l = M_backend->newVector( Xh );
    auto darcyL_rhs = form1( _test=Xh, _vector=F_l );
    // fq
    darcyL_rhs += integrate( _range=elements(mesh), _expr=f*id(q_l) );
    F_l->close();

    auto M_l = M_backend->newMatrix( Xh, Xh );
    auto darcyL = form2( _test=Xh, _trial=Xh, _matrix=M_l);
    // K \grap p \grad q
    darcyL += integrate( _range=elements(mesh), _expr=gradt(p_l)*K*trans(grad(q_l)) );
    M_l->close();

    darcyL += on( boundaryfaces(mesh), _element=p_l, _rhs=darcyL_rhs, _expr=cst(0.) );

    // Solve problem (p)
    M_backend->solve( _matrix=M_l, _solution=p_l, _rhs=F_l );
    // Deduce u :
    u_l = vf::project( Xhvec, elements(mesh), -K*trans(gradv(p_l)) );

    std::cout << "[Darcy] Lagrange solve done" << std::endl;

    // ****** Dual-mixed solving - with Raviart Thomas ******
    prod_space_ptrtype Yh = prod_space_type::New( mesh );
    auto U_rt = Yh->element( "(u,p)" ); //trial
    auto V_rt = Yh->element( "(v,q)" ); //test

    auto u_rt = U_rt.element<0>( "u" ); //velocity field
    auto v_rt = V_rt.element<0>( "v" ); // potential field
    auto p_rt = U_rt.element<1>( "p" );
    auto q_rt = V_rt.element<1>( "q" );

    auto F_rt = M_backend->newVector( Yh );
    auto darcyRT_rhs = form1( _test=Yh, _vector=F_rt );
    // fq
    darcyRT_rhs += integrate( _range=elements(mesh), _expr=f*id(q_rt) );

    auto M_rt = M_backend->newMatrix( Yh, Yh );
    auto darcyRT = form2( _test=Yh, _trial=Yh, _matrix=M_rt);
    // Lambda u v
    darcyRT += integrate( _range=elements(mesh), _expr = -trans(Lambda*idt(u_rt))*id(v_rt) );
    // p div(v)
    darcyRT += integrate( _range=elements(mesh), _expr = idt(p_rt)*div(v_rt) );
    // div(u) q
    darcyRT += integrate( _range=elements(mesh), _expr = divt(u_rt)*id(q_rt) );

    // Solve problem
    backend(_rebuild=true)->solve( _matrix=M_rt, _solution=U_rt, _rhs=F_rt );

    std::cout << "[Darcy] RT solve done" << std::endl;

    // ****** Compute error ******
    auto l2err_u = normL2( _range=elements(mesh), _expr=idv(u_l) - idv(u_rt) );
    auto l2err_p = normL2( _range=elements(mesh), _expr=idv(p_l) - idv(p_rt) );

    if( Environment::isMasterRank() )
        {
            std::cout << "||u(primal) - u(dual-mixed)|| = " << l2err_u << std::endl;
            std::cout << "||p(primal) - p(dual-mixed)|| = " << l2err_p << std::endl;
        }

    BOOST_CHECK_SMALL(l2err_u, 1.0 );
    BOOST_CHECK_SMALL(l2err_p, 1.0 );

    // ****** Export results ******
    export_ptrtype exporter_pro1( export_type::New( this->vm(),
                                  ( boost::format( "%1%-%2%-%3%" )
                                    % this->about().appName()
                                    % ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str()
                                    % "darcy" ).str() ) );

    exporter_pro1->step( 0 )->setMesh( mesh );
    exporter_pro1->step( 0 )->add( "velocity_L", u_l );
    exporter_pro1->step( 0 )->add( "potential_L", p_l );
    exporter_pro1->step( 0 )->add( "velocity_RT", u_rt );
    exporter_pro1->step( 0 )->add( "potential_RT", p_rt );
    exporter_pro1->save();
}

void
TestHDiv::testProjector()
{
    mesh_ptrtype mesh = loadMesh(_mesh = new mesh_type);

    auto RTh = Dh<0>( mesh );
    lagrange_space_v_ptrtype Yh_v = lagrange_space_v_type::New( mesh ); //lagrange vectorial space
    lagrange_space_s_ptrtype Yh_s = lagrange_space_s_type::New( mesh ); //lagrange scalar space

    auto E = Py()*unitX() + Px()*unitY();
    auto f = cst(0.); //div(E) = f

    // L2 projection (Lagrange)
    auto l2_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=L2 ); //l2 vectorial proj
    auto l2_lagS = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=L2 ); //l2 scalar proj
    auto E_pL2_lag = l2_lagV->project( _expr= E );
    auto error_pL2_lag = l2_lagS->project( _expr=divv(E_pL2_lag) - f );

    // H1 projection (Lagrange)
    auto h1_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=H1 ); //h1 vectorial proj
    auto h1_lagS = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=H1 ); //h1 scalar proj
    auto E_pH1_lag = h1_lagV->project( _expr= E, _grad_expr=mat<2,2>(cst(0.),cst(1.),cst(1.),cst(0.)) );
    auto error_pH1_lag = l2_lagS->project( _expr=divv(E_pH1_lag) - f );

    // HDIV projection (Lagrange)
    auto hdiv_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=HDIV );
    auto E_pHDIV_lag = hdiv_lagV->project( _expr= E, _div_expr=cst(0.) );
    auto error_pHDIV_lag = l2_lagS->project( _expr=divv(E_pHDIV_lag) - f );

    // L2 projection (RT)
    auto l2_rt = opProjection( _domainSpace=RTh, _imageSpace=RTh, _type=L2 );
    auto E_pL2_rt = l2_rt->project( _expr= E );
    auto error_pL2_rt = l2_lagS->project( _expr=divv(E_pL2_lag) - f );

    BOOST_TEST_MESSAGE("L2 projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pL2_lag, error_pL2_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pL2_lag, error_pL2_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("H1 projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pH1_lag, error_pH1_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pH1_lag, error_pH1_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("HDIV projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pHDIV_lag, error_pHDIV_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pHDIV_lag, error_pHDIV_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("L2 projection [RT]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pL2_rt, error_pL2_rt ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pL2_rt, error_pL2_rt ) ), 1e-13 );

    std::string proj_name = "projection";
    export_ptrtype exporter_proj( export_type::New( this->vm(),
                                  ( boost::format( "%1%-%2%-%3%" )
                                    % this->about().appName()
                                    % ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str()
                                    % proj_name ).str() ) );

    exporter_proj->step( 0 )->setMesh( mesh );
    exporter_proj->step( 0 )->add( "proj_L2_E[Lagrange]", E_pL2_lag );
    exporter_proj->step( 0 )->add( "proj_H1_E[Lagrange]", E_pH1_lag );
    exporter_proj->step( 0 )->add( "proj_HDiv_E[Lagrange]", E_pHDIV_lag );
    exporter_proj->step( 0 )->add( "proj_L2_E[RT]", E_pL2_rt );
    exporter_proj->save();
}

}
#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( HDIV )

BOOST_AUTO_TEST_CASE( test_hdiv_projection_ref )
{
    BOOST_TEST_MESSAGE( "*** projection on square ***" );
    Feel::TestHDiv t;
    Feel::Environment::changeRepository( boost::format( "/%1%/test_projection/" )
                                         % Feel::Environment::about().appName() );
    t.testProjector();
}

BOOST_AUTO_TEST_CASE( test_hdiv_example_1 )
{
    BOOST_TEST_MESSAGE( "*** resolution of Darcy problem ***" );
    Feel::TestHDiv t;

    Feel::Environment::changeRepository( boost::format( "/%1%/test_Darcy/h_%2%/" )
                                         % Feel::Environment::about().appName()
                                         % t.hSize() );

    t.exampleProblem1();
}

BOOST_AUTO_TEST_SUITE_END()
#else

int
main( int argc, char* argv[] )
{
    Feel::Environment env( argc,argv,
                           makeAbout(), makeOptions() );
    Feel::TestHDiv app_hdiv;
    app_hdiv.testProjector();
    app_hdiv.exampleProblem1();
}

#endif
