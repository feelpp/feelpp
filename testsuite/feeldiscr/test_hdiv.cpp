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
/// Geometry for one-element meshes
gmsh_ptrtype
oneelement_geometry_ref( double h = 2 )
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n" //
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          <<"h=" << h << ";\n"
          <<"Point(1) = {-1,-1,0,h};\n"
          <<"Point(2) = {1,-1,0,h};\n"
          <<"Point(3) = {-1,1,0,h};\n"
          <<"Line(1) = {1,2};\n"
          <<"Line(2) = {2,3};\n"
          <<"Line(3) = {3,1};\n";

    if ( std::abs( h - 2 ) < 1e-10 )
        costr <<"Transfinite Line{1} = 1;\n"
              <<"Transfinite Line{2} = 1;\n"
              <<"Transfinite Line{3} = 1;\n";

    costr <<"Line Loop(4) = {3,1,2};\n"
          <<"Plane Surface(5) = {4};\n"
          <<"Physical Line(\"hor\") = {1};\n"
          <<"Physical Line(\"hypo\") = {2};\n"
          <<"Physical Line(\"vert\") = {3};\n"
          <<"Physical Surface(9) = {5};\n";

    std::ostringstream nameStr;

    if ( std::abs( h - 2 ) < 1e-10 )
        nameStr << "one-elt-ref";
    else
        nameStr << "one-elt-mesh-ref";

    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

// homothetic transformation of reference element (center 0, rate 2)
gmsh_ptrtype
oneelement_geometry_real_1( double h = 2 )
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n"
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          <<"h=" << h << ";\n"
          <<"Point(1) = {-2,-2,0,h};\n"
          <<"Point(2) = {2,-2,0,h};\n"
          <<"Point(3) = {-2,2,0,h};\n"
          <<"Line(1) = {1,2};\n"
          <<"Line(2) = {2,3};\n"
          <<"Line(3) = {3,1};\n";

    if ( std::abs( h - 2 ) < 1e-10 )
        costr <<"Transfinite Line{1} = 1;\n"
              <<"Transfinite Line{2} = 1;\n"
              <<"Transfinite Line{3} = 1;\n";

    costr <<"Line Loop(4) = {3,1,2};\n"
          <<"Plane Surface(5) = {4};\n"
          <<"Physical Line(\"hor\") = {1};\n"
          <<"Physical Line(\"hypo\") = {2};\n"
          <<"Physical Line(\"vert\") = {3};\n"
          <<"Physical Surface(9) = {5};\n";

    std::ostringstream nameStr;

    if ( std::abs( h - 2 ) < 1e-10 )
        nameStr << "one-elt-real-homo";
    else
        nameStr << "one-elt-mesh-homo";

    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

// Rotation of angle (pi/2)
gmsh_ptrtype
oneelement_geometry_real_2( double h = 2 )
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n"
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          <<"h=" << h << ";\n"
          <<"Point(1) = {1,-1,0,h};\n"
          <<"Point(2) = {1,1,0,h};\n"
          <<"Point(3) = {-1,-1,0,h};\n"
          <<"Line(1) = {1,2};\n"
          <<"Line(2) = {2,3};\n"
          <<"Line(3) = {3,1};\n";

    if ( std::abs( h - 2 ) < 1e-10 )
        costr <<"Transfinite Line{1} = 1;\n"
              <<"Transfinite Line{2} = 1;\n"
              <<"Transfinite Line{3} = 1;\n";

    costr <<"Line Loop(4) = {3,1,2};\n"
          <<"Plane Surface(5) = {4};\n"
          <<"Physical Line(\"vert\") = {3};\n"
          <<"Physical Line(\"hypo\") = {2};\n"
          <<"Physical Line(\"hor\") = {1};\n"
          <<"Physical Surface(9) = {5};\n";

    std::ostringstream nameStr;

    if ( std::abs( h - 2 ) < 1e-10 )
        nameStr << "one-elt-real-rot";
    else
        nameStr << "one-elt-mesh-rot";

    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

// Rotation of angle (-pi/2)
gmsh_ptrtype
oneelement_geometry_real_3( double h = 2 )
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n"
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          <<"h=" << h << ";\n"
          <<"Point(1) = {-1,1,0,h};\n"
          <<"Point(2) = {-1,-1,0,h};\n"
          <<"Point(3) = {1,1,0,h};\n"
          <<"Line(1) = {1,2};\n"
          <<"Line(2) = {2,3};\n"
          <<"Line(3) = {3,1};\n";

    if ( std::abs( h - 2 ) < 1e-10 )
        costr <<"Transfinite Line{1} = 1;\n"
              <<"Transfinite Line{2} = 1;\n"
              <<"Transfinite Line{3} = 1;\n";

    costr <<"Line Loop(4) = {3,1,2};\n"
          <<"Plane Surface(5) = {4};\n"
          <<"Physical Line(\"vert\") = {3};\n"
          <<"Physical Line(\"hypo\") = {2};\n"
          <<"Physical Line(\"hor\") = {1};\n"
          <<"Physical Surface(9) = {5};\n";

    std::ostringstream nameStr;

    if ( std::abs( h - 2 ) < 1e-10 )
        nameStr << "one-elt-real-rot2";
    else
        nameStr << "one-elt-mesh-rot2";

    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

gmsh_ptrtype
twoelement_geometry_( double h = 2 )
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n"
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          << "    h=" << h << ";\n"
          << "Point(1) = {0, 0, 0, h};\n"
          << "Point(2) = {1, 0, 0, h};\n"
          << "Point(3) = {1, 1, 0, h};\n"
          << "Point(4) = {2, 0.5, 0, 1.0};\n"
          << "Line(1) = {1, 2};\n"
          << "Line(2) = {2, 3};\n"
          << "Line(3) = {3, 1};\n"
          << "Line(4) = {2, 4};\n"
          << "Line(5) = {4, 3};\n"
          << "Line Loop(6) = {3, 1, 2};\n"
          << "Line Loop(8) = {5, -2, 4};\n";

    if ( std::abs( h - 2 ) < 1e-10 )
        costr << "Transfinite Line(1) = 1;\n"
              << "Transfinite Line(2) = 1;\n"
              << "Transfinite Line(3) = 1;\n"
              << "Transfinite Line(4) = 1;\n"
              << "Transfinite Line(5) = 1;\n";

    costr << "Plane Surface(11) = {6};\n"
          << "Plane Surface(12) = {8};\n";
    std::ostringstream nameStr;

    if ( std::abs( h - 2 ) < 1e-10 )
        nameStr << "two-elt-mesh";
    else
        nameStr << "two-elt-mesh-fine";

    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
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
        M_backend( backend_type::build( this->vm() ) ),
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
    void shape_functions( gmsh_ptrtype ( *one_element_mesh )( double ));
    void testProjector(gmsh_ptrtype ( *one_element_mesh_desc_fun )( double ));
    void exampleProblem1();

private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! exporter factory
    export_ptrtype exporter;

}; //TestHDiv

// Resolve problem curl(curl(u)) + u = f with cross_prod(u,n) = 0 on boundary

void
TestHDiv::exampleProblem1()
{
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str() ,
                                                _shape="hypercube",
                                                _usenames=true,
                                                _dim=2,
                                                _h=meshSize,
                                                _xmin=-1,_xmax=1,
                                                _ymin=-1,_ymax=1 ) );

    auto K = eye<2>(); // Hydraulic conductivity tensor
    auto Lambda = eye<2>(); // Hydraulic resistivity tensor
    auto f = Px()+Py(); // int_omega f = 0

    // ****** Primal formulation - with Lagrange ******
    lagrange_space_s_ptrtype Xh = lagrange_space_s_type::New( mesh );
    lagrange_space_v_ptrtype Xhvec = lagrange_space_v_type::New( mesh );
    auto p_l = Xh->element( "p" );
    auto u_l = Xhvec->element( "u" );
    auto q_l = Xh->element( "q" );

    auto F_l = M_backend->newVector( Xh );
    auto darcyL_rhs = form1( _test=Xh, _vector=F_l );
    // fq
    darcyL_rhs = integrate( _range=elements(mesh), _expr=f*id(q_l) );

    auto M_l = M_backend->newMatrix( Xh, Xh );
    auto darcyL = form2( _test=Xh, _trial=Xh, _matrix=M_l);
    darcyL = integrate( _range=elements(mesh), _expr=grad(q_l)*K*trans(gradt(p_l)) );

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

    auto epsilon = 1e-7;

    auto F_rt = M_backend->newVector( Yh );
    auto darcyRT_rhs = form1( _test=Yh, _vector=F_rt );
    // fq
    darcyRT_rhs += integrate( _range=elements(mesh), _expr=f*id(q_rt) );

    auto M_rt = M_backend->newMatrix( Yh, Yh );
    auto darcyRT = form2( _test=Yh, _trial=Yh, _matrix=M_rt);
    // Lambda u v
    darcyRT = integrate( _range=elements(mesh), _expr = -trans(idt(u_rt))*Lambda*id(v_rt) );
    // p div(v)
    darcyRT += integrate( _range=elements(mesh), _expr = idt(p_rt)*div(v_rt) );
    // div(u) q
    darcyRT += integrate( _range=elements(mesh), _expr = divt(u_rt)*id(q_rt) );

    //espilon pq
    darcyRT += integrate( _range=elements(mesh), _expr=epsilon*idt(p_rt)*id(q_rt) );

    //boundary condition
    //darcyRT += on( _range=elements(mesh), _element=u_rt, _rhs=darcyRT_rhs,_expr=cst(0.) );

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
TestHDiv::testProjector(gmsh_ptrtype ( *one_element_mesh_desc_fun )( double ))
{
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str() ,
                                                _shape="hypercube",
                                                _usenames=true,
                                                _dim=2,
                                                _h=meshSize,
                                                _xmin=-1,_xmax=1,
                                                _ymin=-1,_ymax=1 ) );

    // Only one element in the mesh - TEMPORARLY
    // mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
    //                                     _desc = one_element_mesh_desc_fun( 2 ) );

    space_ptrtype Xh = space_type::New( mesh );
    lagrange_space_v_ptrtype Yh_v = lagrange_space_v_type::New( mesh ); //lagrange vectorial space
    lagrange_space_s_ptrtype Yh_s = lagrange_space_s_type::New( mesh ); //lagrange scalar space

    auto E = Py()*unitX() + Px()*unitY();
    auto f = cst(0.); //div(E) = f

    // L2 projection
    auto l2_v = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=L2 ); //l2 vectorial proj
    auto l2_s = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=L2 ); //l2 scalar proj
    auto E_l2 = l2_v->project( _expr= trans(E) );
    auto f_l2 = l2_s->project( _expr= f );
    auto error_l2 = l2_s->project( _expr=divv(E_l2) - idv(f_l2) );

    // H1 projection
    auto h1_v = opProjection( _domainSpace=Yh_v, _imageSpace=Yh_v, _type=H1 ); //h1 vectorial proj
    auto h1_s = opProjection( _domainSpace=Yh_s, _imageSpace=Yh_s, _type=H1 ); //h1 scalar proj
    auto E_h1 = h1_v->project( _expr= trans(E), _grad_expr=mat<2,2>(cst(0.),cst(1.),cst(1.),cst(0.)) );
    auto f_h1 = h1_s->project( _expr= f );
    auto error_h1 = h1_s->project( _expr=divv(E_h1) - idv(f_h1) );

    auto hdiv = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=HDIV ); //hdiv proj (RT elts)
    auto E_hdiv = hdiv->project( _expr= trans(E), _div_expr=cst(0.) );
    // Piola transfo (ref elt -> real elt ) : div( (1/detJ)*J*u ) = (1/detJ)*div(Ju) = (1/detJ)*trace(J*grad(u))
    auto l2norm_div = normL2( _range=elements(mesh), _expr=(1/detJ())*trace(J()*gradv(E_hdiv)) );

    BOOST_TEST_MESSAGE("L2 projection : error[div(E)-f]");
    std::cout << "error L2: " << math::sqrt( l2_s->energy( error_l2, error_l2 ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( l2_s->energy( error_l2, error_l2 ) ), 1e-13 );
    BOOST_TEST_MESSAGE("H1 projection : error[div(E)-f]");
    std::cout << "error H1: " << math::sqrt( h1_s->energy( error_h1, error_h1 ) ) << "\n";
    BOOST_CHECK_SMALL( math::sqrt( h1_s->energy( error_h1, error_h1 ) ), 1e-13 );
    BOOST_TEST_MESSAGE("HDIV projection : error[div(E)-f]");
    std::cout << "error div(f): " << l2norm_div << "\n";
    BOOST_CHECK_SMALL( l2norm_div, 1e-13 );

    std::string proj_name = "projection";
    export_ptrtype exporter_proj( export_type::New( this->vm(),
                                  ( boost::format( "%1%-%2%-%3%" )
                                    % this->about().appName()
                                    % ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % 2 % 1 ).str()
                                    % proj_name ).str() ) );

    exporter_proj->step( 0 )->setMesh( mesh );
    exporter_proj->step( 0 )->add( "proj_L2_E", E_l2 );
    exporter_proj->step( 0 )->add( "proj_H1_E", E_h1 );
    exporter_proj->step( 0 )->add( "proj_HDiv_E", E_hdiv );
    exporter_proj->save();
}

void
TestHDiv::shape_functions( gmsh_ptrtype ( *one_element_mesh_desc_fun )( double ) )
{
    //    using namespace Feel::vf;

    mesh_ptrtype oneelement_mesh = createGMSHMesh( _mesh=new mesh_type,
                                   _desc = one_element_mesh_desc_fun( 2 ) );

    // then a fine mesh which we use to export the basis function to  visualize them
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=one_element_mesh_desc_fun( meshSize ) );

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
                                     % one_element_mesh_desc_fun( 2 )->prefix()
                                     % shape_name ).str() ) );

    exporter_shape->step( 0 )->setMesh( mesh );

    for ( size_type i = 0; i < Xh->nLocalDof(); ++i )
    {
        U_ref.zero();
        U_ref( i ) = 1;

        // U_vec[i] = shape function associated with dof i
        u_vec[i] = U_ref;

        std::ostringstream ostr;
        ostr <<  one_element_mesh_desc_fun( 2 )->prefix()<< "-" << Xh->basis()->familyName() << "-" << i;
        exporter_shape->step( 0 )->add( ostr.str(), U_ref );
    }

    exporter_shape->save();

    auto F = M_backend->newVector( Xh );

    int check_size = Xh->nLocalDof()*Xh->nLocalDof();
    std::vector<double> checkidv( check_size );
    std::vector<double> checkform1( check_size );
    std::vector<std::string> faces = boost::assign::list_of( "hypo" )( "vert" )( "hor" );

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
            auto int_u_n = integrate( markedfaces( oneelement_mesh, face ), trans(N())*(1/detJ())*J()*idv( u_vec[i] ) ).evaluate()( 0,0 );
            if ( faceid == i )
                BOOST_CHECK_CLOSE( int_u_n, 1, 1e-13 );
            else
                BOOST_CHECK_SMALL( int_u_n, 1e-13 );
            checkidv[3*i+faceid] = int_u_n;

            form1( _test=Xh, _vector=F, _init=true ) = integrate( markedfaces( oneelement_mesh, face), trans( N() )*(1/detJ())*J()*id( V_ref ) );
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

        //BOOST_TEST_MESSAGE( "*** Stokes theorem - Integration ***");
        auto int_divx = integrate( elements( oneelement_mesh ), (1/detJ())*trace(J()*gradv(u_vec[i]))  ).evaluate()( 0,0 );
        auto int_xn = integrate(boundaryfaces(oneelement_mesh), trans(N())*(1/detJ())*J()*idv(u_vec[i])).evaluate()( 0,0 );

        BOOST_CHECK_CLOSE( int_divx, 1, 1e-13 );
        BOOST_CHECK_CLOSE( int_divx, int_xn, 1e-13 );
        checkStokesidv[i] = int_divx;
        checkStokesidv[i + Xh->nLocalDof()] = int_xn;

        //BOOST_TEST_MESSAGE( "*** Stokes theorem - Linear form ***");
        form1( _test=Xh, _vector=F, _init=true ) = integrate( elements( oneelement_mesh ), (1/detJ())*trace(J()*grad(V_ref)) );
        auto form_divx = inner_product( u_vec[i], *F );
        form1( _test=Xh, _vector=F, _init=true ) = integrate( boundaryfaces( oneelement_mesh ), trans(N())*(1/detJ())*J()*id(V_ref) );
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
            BOOST_TEST_MESSAGE( " *** dof N_"<< i << " (associated with " << face << " face) *** \n"
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
            BOOST_TEST_MESSAGE( " *** dof N_"<< i << " (associated with " << face << " edge) *** \n"
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

BOOST_AUTO_TEST_SUITE( HDIV )
BOOST_AUTO_TEST_CASE( test_hdiv_N0_ref )
{
    BOOST_TEST_MESSAGE( "*** shape functions on reference element (1 elt) ***" );
    Feel::TestHDiv t;
    t.shape_functions( &Feel::oneelement_geometry_ref );
}
BOOST_AUTO_TEST_CASE( test_hdiv_N0_real1 )
{
    BOOST_TEST_MESSAGE( "*** shape functions on real element - homothetic transfo (1 elt) ***" );
    Feel::TestHDiv t;
    t.shape_functions( &Feel::oneelement_geometry_real_1 );
}
BOOST_AUTO_TEST_CASE( test_hdiv_N0_real2 )
{
    BOOST_TEST_MESSAGE( "*** shape functions on real element - rotation pi/2 (1 elt) ***" );
    Feel::TestHDiv t;
    t.shape_functions( &Feel::oneelement_geometry_real_2 );
}

BOOST_AUTO_TEST_CASE( test_hdiv_N0_real3 )
{
    BOOST_TEST_MESSAGE( "*** shape functions on real element - rotation -pi/2 (1 elt) ***" );
    Feel::TestHDiv t;
    t.shape_functions( &Feel::oneelement_geometry_real_3 );
}

BOOST_AUTO_TEST_CASE( test_hdiv_projection_ref )
{
    BOOST_TEST_MESSAGE( "*** projection on reference element ***" );
    Feel::TestHDiv t;
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/test_projection/%2%/" )
                                         % Feel::Environment::about().appName()
                                         % "ref" );
    t.testProjector(&Feel::oneelement_geometry_ref);
}
BOOST_AUTO_TEST_CASE( test_hdiv_projection_real1 )
{
    BOOST_TEST_MESSAGE( "*** projection on real element - homothetic transfo***" );
    Feel::TestHDiv t;
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/test_projection/%2%/" )
                                         % Feel::Environment::about().appName()
                                         % "real_1" );
    t.testProjector(&Feel::oneelement_geometry_real_1);
}
BOOST_AUTO_TEST_CASE( test_hdiv_projection_real2 )
{
    BOOST_TEST_MESSAGE( "*** projection on real element - rotation pi/2 ***" );
    Feel::TestHDiv t;
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/test_projection/%2%/" )
                                         % Feel::Environment::about().appName()
                                         % "real_2" );
    t.testProjector(&Feel::oneelement_geometry_real_2);
}
BOOST_AUTO_TEST_CASE( test_hdiv_projection_real3 )
{
    BOOST_TEST_MESSAGE( "*** projection on real element - rotation -pi/2 ***" );
    Feel::TestHDiv t;
    Feel::Environment::changeRepository( boost::format( "/testsuite/feeldiscr/%1%/test_projection/%2%/" )
                                         % Feel::Environment::about().appName()
                                         % "real_2" );
    t.testProjector(&Feel::oneelement_geometry_real_3);
}

BOOST_AUTO_TEST_CASE( test_hdiv_example_1 )
{
    BOOST_TEST_MESSAGE( "*** resolution of Darcy problem ***" );
    Feel::TestHDiv t;
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
    //app_hdiv.shape_functions();
    //app_hdiv.exampleProblem1();
}

#endif
