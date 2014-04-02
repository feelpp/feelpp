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
#define BOOST_TEST_MODULE H_div approximation
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>

#include <feel/feel.hpp>
#include <feel/feelpoly/raviartthomas.hpp>

namespace Feel
{
/// Geometry for one-element meshes
std::string
oneelement_geometry_ref()
{
    std::string name = "one-elt-ref";

    if(!fs::exists( name+".msh" ))
        {
            std::ofstream costr(name+".msh");
            costr << "$MeshFormat\n"
                  << "2.2 0 8\n"
                  << "$EndMeshFormat\n"
                  << "$PhysicalNames\n"
                  << "5\n"
                  << "2 1 \"xyFace\"\n"
                  << "2 2 \"yzFace\"\n"
                  << "2 3 \"xzFace\"\n"
                  << "2 4 \"xyzFace\"\n"
                  << "3 5 \"volume\"\n"
                  << "$EndPhysicalNames\n"
                  << "$Nodes\n"
                  << "4\n"
                  << "1 -1 -1 -1\n"
                  << "2 1 -1 -1\n"
                  << "3 -1 1 -1\n"
                  << "4 -1 -1 1\n"
                  << "$EndNodes\n"
                  << "$Elements\n"
                  << "5\n"
                  << "1 2 2 1 11 1 2 3\n"
                  << "2 2 2 2 12 1 4 3\n"
                  << "3 2 2 3 13 1 2 4\n"
                  << "4 2 2 4 14 2 3 4\n"
                  << "5 4 2 5 20 3 4 1 2\n"
                  << "$EndElements\n";
            costr.close();
        }

        return name;
}

// homothetic transformation of reference element (center 0, rate 2)
std::string
oneelement_geometry_real_1()
{
    std::string name = "one-elt-real-homo";

    if(!fs::exists( name+".msh" ))
        {
            std::ofstream costr(name+".msh");
            costr << "$MeshFormat\n"
                  << "2.2 0 8\n"
                  << "$EndMeshFormat\n"
                  << "$PhysicalNames\n"
                  << "5\n"
                  << "2 1 \"xyFace\"\n"
                  << "2 2 \"yzFace\"\n"
                  << "2 3 \"xzFace\"\n"
                  << "2 4 \"xyzFace\"\n"
                  << "3 5 \"volume\"\n"
                  << "$EndPhysicalNames\n"
                  << "$Nodes\n"
                  << "4\n"
                  << "1 -2 -2 -2\n"
                  << "2 2 -2 -2\n"
                  << "3 -2 2 -2\n"
                  << "4 -2 -2 2\n"
                  << "$EndNodes\n"
                  << "$Elements\n"
                  << "5\n"
                  << "1 2 2 1 11 1 2 3\n"
                  << "2 2 2 2 12 1 4 3\n"
                  << "3 2 2 3 13 1 2 4\n"
                  << "4 2 2 4 14 2 3 4\n"
                  << "5 4 2 5 20 3 4 1 2\n"
                  << "$EndElements\n";
            costr.close();
        }

    return name;
}

// Rotation of angle (pi/2) around x axis
std::string
oneelement_geometry_real_2()
{
    std::string name = "one-elt-real-rotx";
    std::ofstream costr(name+".msh");

    costr << "$MeshFormat\n"
          << "2.2 0 8\n"
          << "$EndMeshFormat\n"
          << "$PhysicalNames\n"
          << "5\n"
          << "2 1 \"xyFace\"\n"
          << "2 2 \"yzFace\"\n"
          << "2 3 \"xzFace\"\n"
          << "2 4 \"xyzFace\"\n"
          << "3 5 \"volume\"\n"
          << "$EndPhysicalNames\n"
          << "$Nodes\n"
          << "4\n"
          << "1 -1 1 -1\n"
          << "2 1 1 -1\n"
          << "3 -1 1 1\n"
          << "4 -1 -1 -1\n"
          << "$EndNodes\n"
          << "$Elements\n"
          << "5\n"
          << "1 2 2 1 11 1 2 3\n"
          << "2 2 2 2 12 1 4 3\n"
          << "3 2 2 3 13 1 2 4\n"
          << "4 2 2 4 14 2 3 4\n"
          << "5 4 2 5 20 3 4 1 2\n"
          << "$EndElements\n";

    costr.close();
    return name;
}

// Rotation of angle (pi/2) around y axis
std::string
oneelement_geometry_real_3()
{
    std::string name = "one-elt-real-roty";

    if(!fs::exists( name+".msh" ))
        {
            std::ofstream costr(name+".msh");
            costr << "$MeshFormat\n"
                  << "2.2 0 8\n"
                  << "$EndMeshFormat\n"
                  << "$PhysicalNames\n"
                  << "5\n"
                  << "2 1 \"xyFace\"\n"
                  << "2 2 \"yzFace\"\n"
                  << "2 3 \"xzFace\"\n"
                  << "2 4 \"xyzFace\"\n"
                  << "3 5 \"volume\"\n"
                  << "$EndPhysicalNames\n"
                  << "$Nodes\n"
                  << "4\n"
                  << "1 -1 -1 1\n"
                  << "2 -1 -1 -1\n"
                  << "3 -1 1 1\n"
                  << "4 1 -1 -1\n"
                  << "$EndNodes\n"
                  << "$Elements\n"
                  << "5\n"
                  << "1 2 2 1 11 1 2 3\n"
                  << "2 2 2 2 12 1 4 3\n"
                  << "3 2 2 3 13 1 2 4\n"
                  << "4 2 2 4 14 2 3 4\n"
                  << "5 4 2 5 20 3 4 1 2\n"
                  << "$EndElements\n";
            costr.close();
        }

    return name;
}

// Rotation of angle (pi/2) around z axis
std::string
oneelement_geometry_real_4()
{
    std::string name = "one-elt-real-rotz";

    if(!fs::exists( name+".msh" ))
        {
            std::ofstream costr(name+".msh");
            costr << "$MeshFormat\n"
                  << "2.2 0 8\n"
                  << "$EndMeshFormat\n"
                  << "$PhysicalNames\n"
                  << "5\n"
                  << "2 1 \"xyFace\"\n"
                  << "2 2 \"yzFace\"\n"
                  << "2 3 \"xzFace\"\n"
                  << "2 4 \"xyzFace\"\n"
                  << "3 5 \"volume\"\n"
                  << "$EndPhysicalNames\n"
                  << "$Nodes\n"
                  << "4\n"
                  << "1 1 -1 -1\n"
                  << "2 1 1 -1\n"
                  << "3 -1 -1 -1\n"
                  << "4 1 -1 1\n"
                  << "$EndNodes\n"
                  << "$Elements\n"
                  << "5\n"
                  << "1 2 2 1 11 1 2 3\n"
                  << "2 2 2 2 12 1 4 3\n"
                  << "3 2 2 3 13 1 2 4\n"
                  << "4 2 2 4 14 2 3 4\n"
                  << "5 4 2 5 20 3 4 1 2\n"
                  << "$EndElements\n";
            costr.close();
        }

    return name;
}

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
    po::options_description testhdiv3Doptions( "test hdiv3D options" );
    testhdiv3Doptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
    ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
    ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
    ;
    return testhdiv3Doptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_hdiv3D" ,
                     "test_hdiv3D" ,
                     "0.1",
                     "Test for h_div space",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "cecile.daversin@lncmi.cnrs.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

using namespace Feel;

class TestHDiv3D
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
    TestHDiv3D()
        :
        super(),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        std::cout << "[TestHDiv3D]\n";

        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % this->vm()["hsize"].as<double>()
                              );
    }

    /**
     * run the application
     */
    inline double hSize(){return meshSize;}
    void shape_functions( std::string ( *one_element_mesh )() );

    void testProjector(std::string ( *one_element_mesh_desc_fun )() );
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
TestHDiv3D::exampleProblem1()
{
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str() ,
                                                      _shape="hypercube",
                                                      _usenames=true,
                                                      _dim=3,
                                                      _h=meshSize,
                                                      _xmin=-1,_xmax=1,
                                                      _ymin=-1,_ymax=1,
                                                      _zmin=-1,_zmax=1) );

    //auto K = ones<2,2>(); // Hydraulic conductivity tensor
    auto K = mat<3,3>(cst(2.),cst(1.),cst(1.),
                      cst(1.),cst(2.),cst(1.),
                      cst(1.),cst(1.),cst(2.));
    //auto Lambda = ones<2,2>(); // Hydraulic resistivity tensor
    auto Lambda = (1.0/4.0)*mat<3,3>(cst(3.),cst(-1.),cst(-1.),
                                     cst(-1.),cst(3.),cst(-1.),
                                     cst(-1.),cst(-1.),cst(3.)); // Lambda = inv(K)
    auto f = Px()+Py()+Pz(); // int_omega f = 0
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

    std::cout << "||u(primal) - u(dual-mixed)|| = " << l2err_u << std::endl;
    std::cout << "||p(primal) - p(dual-mixed)|| = " << l2err_p << std::endl;

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
TestHDiv3D::testProjector(std::string ( *one_element_mesh_desc_fun )())
{
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str() ,
                                                      _shape="hypercube",
                                                      _usenames=true,
                                                      _dim=3,
                                                      _h=meshSize,
                                                      _xmin=-1,_xmax=1,
                                                      _ymin=-1,_ymax=1,
                                                      _zmin=-1,_zmax=1 ) );

    // Only one element in the mesh - TEMPORARLY
    // mesh_ptrtype mesh = loadGMSHMesh( _mesh=new mesh_type,
    //                                   _filename=one_element_mesh_desc_fun() );

    auto RTh = Dh<0>( mesh );
    lagrange_space_v_ptrtype Yh_v = lagrange_space_v_type::New( mesh ); //lagrange vectorial space
    lagrange_space_s_ptrtype Yh_s = lagrange_space_s_type::New( mesh ); //lagrange scalar space

    auto E = Px()*unitX() + Py()*unitY() + Pz()*unitZ(); //(x,y,z)
    auto f = cst(3.); //div(E) = f

    // L2 projection (Lagrange)
    auto l2_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=L2 ); //l2 vectorial proj
    auto l2_lagS = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=L2 ); //l2 scalar proj
    auto E_pL2_lag = l2_lagV->project( _expr= trans(E) );
    auto error_pL2_lag = l2_lagS->project( _expr=divv(E_pL2_lag) - f );

    // L2 projection (RT)
    auto l2_rt = opProjection( _domainSpace=RTh, _imageSpace=RTh, _type=L2 );
    auto E_pL2_rt = l2_rt->project( _expr= trans(E) );
    auto error_pL2_rt = l2_lagS->project( _expr=divv(E_pL2_lag) - f );

    // H1 projection (Lagrange)
    auto h1_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=H1 ); //h1 vectorial proj
    auto h1_lagS = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=H1 ); //h1 scalar proj
    auto E_pH1_lag = h1_lagV->project( _expr= trans(E), _grad_expr=eye<3>() );
    auto error_pH1_lag = l2_lagS->project( _expr=divv(E_pH1_lag) - f );

    // H1 projection (RT)
    auto h1_rt = opProjection( _domainSpace=RTh, _imageSpace=RTh, _type=H1 ); //h1 vectorial proj
    auto E_pH1_rt = h1_rt->project( _expr= trans(E), _grad_expr=eye<3>() );
    auto error_pH1_rt = l2_lagS->project( _expr=divv(E_pH1_rt) - f );

    // HDIV projection (Lagrange)
    auto hdiv_lagV = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=HDIV );
    auto hdiv_lagS = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=HDIV );
    auto E_pHDIV_lag = hdiv_lagV->project( _expr= trans(E), _div_expr=cst(0.) );
    auto error_pHDIV_lag = l2_lagS->project( _expr=divv(E_pHDIV_lag) - f );

    // HDIV projection (RT)
    auto hdiv = opProjection( _domainSpace=RTh, _imageSpace=RTh, _type=HDIV ); //hdiv proj (RT elts)
    auto E_pHDIV_rt = hdiv->project( _expr= trans(E), _div_expr=cst(0.) );
    auto error_pHDIV_rt = l2_lagS->project( _expr=divv(E_pHDIV_rt) - f );

    BOOST_TEST_MESSAGE("L2 projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pL2_lag, error_pL2_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pL2_lag, error_pL2_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("L2 projection [RT]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pL2_rt, error_pL2_rt ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pL2_rt, error_pL2_rt ) ), 1e-13 );
    BOOST_TEST_MESSAGE("H1 projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pH1_lag, error_pH1_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pH1_lag, error_pH1_lag ) ), 1e-13 );
    BOOST_TEST_MESSAGE("H1 projection [RT]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pH1_rt, error_pH1_rt ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pH1_rt, error_pH1_rt ) ), 1e-13 );
    BOOST_TEST_MESSAGE("HDIV projection [Lagrange]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pHDIV_lag, error_pHDIV_lag ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pHDIV_lag, error_pHDIV_rt ) ), 1e-13 );
    BOOST_TEST_MESSAGE("HDIV projection [RT]: error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_lagS->energy( error_pHDIV_rt, error_pHDIV_rt ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_lagS->energy( error_pHDIV_rt, error_pHDIV_rt ) ), 1e-13 );

    std::string proj_name = "projection";
    export_ptrtype exporter_proj( export_type::New( this->vm(),
                                  ( boost::format( "%1%-%2%-%3%" )
                                    % this->about().appName()
                                    % ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str()
                                    % proj_name ).str() ) );

    exporter_proj->step( 0 )->setMesh( mesh );
    exporter_proj->step( 0 )->add( "proj_L2_E[Lagrange]", E_pL2_lag );
    exporter_proj->step( 0 )->add( "proj_L2_E[RT]", E_pL2_rt );
    exporter_proj->step( 0 )->add( "proj_H1_E[Lagrange]", E_pH1_lag );
    exporter_proj->step( 0 )->add( "proj_H1_E[RT]", E_pH1_rt );
    exporter_proj->step( 0 )->add( "proj_HDiv_E[Lagrange]", E_pHDIV_lag );
    exporter_proj->step( 0 )->add( "proj_HDiv_E[RT]", E_pHDIV_rt );
    exporter_proj->save();
}


void
//TestHDiv3D::shape_functions( gmsh_ptrtype ( *one_element_mesh_desc_fun )(double) )
TestHDiv3D::shape_functions( std::string ( *one_element_mesh_desc_fun )() )
{

    auto mesh_name = one_element_mesh_desc_fun()+".msh"; //create the mesh and load it
    mesh_ptrtype oneelement_mesh = loadGMSHMesh( _mesh=new mesh_type,
                                                 _filename=mesh_name );

    auto refine_level = std::floor(1 - math::log( meshSize )); //Deduce refine level from meshSize (option)
    mesh_ptrtype mesh = loadGMSHMesh( _mesh=new mesh_type,
                                      _filename=mesh_name ,
                                      _refine=( int )refine_level );


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
                                     % one_element_mesh_desc_fun()//->prefix()
                                     % shape_name ).str() ) );

    exporter_shape->step( 0 )->setMesh( mesh );

    for ( size_type i = 0; i < Xh->nLocalDof(); ++i )
    {
        U_ref.zero();
        U_ref( i ) = 1;

        // U_vec[i] = shape function associated with dof i
        u_vec[i] = U_ref;

        std::ostringstream ostr;
        ostr <<  one_element_mesh_desc_fun() << "-" << Xh->basis()->familyName() << "-" << i;
        exporter_shape->step( 0 )->add( ostr.str(), U_ref );
    }

    exporter_shape->save();

    auto F = M_backend->newVector( Xh );

    std::cout << "Xh->nLocalDof() = " << Xh->nLocalDof() << std::endl;
    std::cout << "Xh->dof()... = " << Xh->dof()->nLocalDofWithGhost() << std::endl;

    int check_size = Xh->nLocalDof()*Xh->nLocalDof();
    std::vector<double> checkidv( check_size );
    std::vector<double> checkform1( check_size );
    // Faces have to be in the correct order
    std::vector<std::string> faces = boost::assign::list_of( "xzFace" )( "xyFace" )( "xyzFace" )( "yzFace" );

    std::vector<double> checkStokesidv( 2*Xh->nLocalDof() );
    std::vector<double> checkStokesform1( 2*Xh->nLocalDof() );

    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        // ****** Check shape function property : alpha_i(N_j) = \delta_{i,j} \forall i,j ***** //
        // alpha_i(N_j) = \int_face(i) u.n
        // Piola transformation : u -> (J/detJ)u_ref
        int faceid = 0;
        BOOST_FOREACH( std::string face, faces )
        {
            auto int_u_n = integrate( markedfaces( oneelement_mesh, face ), trans(N())*idv( u_vec[i]) ).evaluate()( 0,0 );

            if ( faceid == i )
                BOOST_CHECK_CLOSE( int_u_n, 1, 1e-13 );
            else
                BOOST_CHECK_SMALL( int_u_n, 1e-13 );
            checkidv[3*i+faceid] = int_u_n;

            form1( _test=Xh, _vector=F, _init=true ) = integrate( markedfaces( oneelement_mesh, face), trans( N() )*id( V_ref ) );
            auto form_v_n = inner_product( u_vec[i], *F );

            if ( faceid == i )
                BOOST_CHECK_CLOSE( form_v_n, 1, 1e-13 );
            else
                BOOST_CHECK_SMALL( form_v_n, 1e-13 );
            checkform1[3*i+faceid] = form_v_n;

            ++faceid;
        }

        // ****** Check Stokes theorem : \int_\Omega div(x)w = \int_\delta\Omega (x.n)w \forall w ***** //
        // Test : \int_\Omega div(N_i) = \int_\delta\Omega (N_i.n) = 1 [RT : \int_\delta\Omega (N_i.n) = \sum \face (N_i.n) = 1]
        // Piola transformation : (ref elt -> real elt) : u -> ( J/detJ )u [div u -> (1/detJ)trace(J grad u)]

        auto int_divx = integrate( elements( oneelement_mesh ), divv(u_vec[i])  ).evaluate()( 0,0 );
        auto int_xn = integrate(boundaryfaces(oneelement_mesh), trans(N())*idv(u_vec[i]) ).evaluate()( 0,0 );

        BOOST_CHECK_CLOSE( int_divx, 1, 1e-13 );
        BOOST_CHECK_CLOSE( int_divx, int_xn, 1e-13 );
        checkStokesidv[i] = int_divx;
        checkStokesidv[i + Xh->nLocalDof()] = int_xn;

        form1( _test=Xh, _vector=F, _init=true ) = integrate( elements( oneelement_mesh ), div(V_ref) );
        auto form_divx = inner_product( u_vec[i], *F );
        form1( _test=Xh, _vector=F, _init=true ) = integrate( boundaryfaces( oneelement_mesh ), trans(N())*id(V_ref) );
        auto form_xn = inner_product( u_vec[i], *F );

        BOOST_CHECK_CLOSE( form_divx, 1, 1e-13 );
        BOOST_CHECK_CLOSE( form_divx, form_xn, 1e-13 );
        checkStokesform1[i] = form_divx;
        checkStokesform1[i + Xh->nLocalDof()] = form_xn;

        // **********************************************************************************************************//
    }

    BOOST_TEST_MESSAGE( " ********** Values of alpha_i(N_j ) = delta_{i,j}  ********** \n"
                        << "\n"
                        << " ********** Using idv keyword ********************* "
                        << "\n" );

    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        int faceid = 0;
        BOOST_FOREACH( std::string face, faces )
        {
            BOOST_TEST_MESSAGE( " *** dof N_"<< i << " (associated with " << faces[i] << " face) *** \n"
                                << "alpha_"<< faceid << "(N_"<<i<<") = " << checkidv[3*i+faceid] );
            ++faceid;
        }
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    BOOST_TEST_MESSAGE( " ********** Values of alpha_i (N_j ) = delta_{i,j}  ********** \n"
                        << "\n"
                        << " ********** Using form1 keyword ********************* "
                        << "\n" );

    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        int faceid = 0;
        BOOST_FOREACH( std::string face, faces )
        {
            BOOST_TEST_MESSAGE( " *** dof N_"<< i << " (associated with " << faces[i] << " edge) *** \n"
                                << "alpha_"<< faceid << "(N_"<<i<<") = " << checkform1[3*i+faceid] );
            ++faceid;
        }
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    BOOST_TEST_MESSAGE( " ********** Stokes theorem  ********** \n"
                        << "\n"
                        << " ********** Using idv keyword ********************* "
                        << "\n" );
    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        BOOST_TEST_MESSAGE("int(Omega) div(N_" << i << ") = " << checkStokesidv[i] );
        BOOST_TEST_MESSAGE("int(boundary) N_" << i << ".n = " << checkStokesidv[i+Xh->nLocalDof()] );
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    BOOST_TEST_MESSAGE( " ********** Stokes theorem  ********** \n"
                        << "\n"
                        << " ********** Using form1 keyword ********************* "
                        << "\n" );
    for ( int i = 0; i < Xh->nLocalDof(); ++i )
    {
        BOOST_TEST_MESSAGE("int(Omega) div(N_" << i << ") = " << checkStokesform1[i] );
        BOOST_TEST_MESSAGE("int(boundary) N_" << i << ".n = " << checkStokesform1[i+Xh->nLocalDof()] );
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }

    //// ************************************************************************************ ////
}

}
#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( HDIV3D )

BOOST_AUTO_TEST_CASE( test_hdiv3D_N0_ref )
{
    BOOST_TEST_MESSAGE( "*** shape functions on reference element (1 elt) ***" );
    Feel::TestHDiv3D t;
    t.shape_functions( &Feel::oneelement_geometry_ref );
}
BOOST_AUTO_TEST_CASE( test_hdiv_N0_real1 )
{
    BOOST_TEST_MESSAGE( "*** shape functions on real element - homothetic transfo (1 elt) ***" );
    Feel::TestHDiv3D t;
    t.shape_functions( &Feel::oneelement_geometry_real_1 );
}
BOOST_AUTO_TEST_CASE( test_hdiv_N0_real2 )
{
    BOOST_TEST_MESSAGE( "*** shape functions on real element - rotation pi/2 - x axis (1 elt) ***" );
    Feel::TestHDiv3D t;
    t.shape_functions( &Feel::oneelement_geometry_real_2 );
}
BOOST_AUTO_TEST_CASE( test_hdiv_N0_real3 )
{
    BOOST_TEST_MESSAGE( "*** shape functions on real element - rotation pi/2 - y axis (1 elt) ***" );
    Feel::TestHDiv3D t;
    t.shape_functions( &Feel::oneelement_geometry_real_3 );
}
BOOST_AUTO_TEST_CASE( test_hdiv_N0_real4 )
{
    BOOST_TEST_MESSAGE( "*** shape functions on real element - rotation pi/2 - z axis (1 elt) ***" );
    Feel::TestHDiv3D t;
    t.shape_functions( &Feel::oneelement_geometry_real_4 );
}

BOOST_AUTO_TEST_CASE( test_hdiv_projection_ref )
{
    BOOST_TEST_MESSAGE( "*** projection on cube ***" );
    Feel::TestHDiv3D t;
    Feel::Environment::changeRepository( boost::format( "/%1%/test_projection/" )
                                         % Feel::Environment::about().appName() );
    t.testProjector(&Feel::oneelement_geometry_ref);
}

BOOST_AUTO_TEST_CASE( test_hdiv_example_1 )
{
    BOOST_TEST_MESSAGE( "*** resolution of Darcy problem ***" );
    Feel::TestHDiv3D t;

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
    Feel::TestHDiv3D app_hdiv;
    app_hdiv.shape_functions();
    //app_hdiv.exampleProblem1();
}

#endif
