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
    po::options_description testhcurloptions( "test h_curl_assembly options" );
    testhcurloptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
    ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
    ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
    ;
    return testhcurloptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_hcurl_assembly" ,
                     "test_hcurl_assembly" ,
                     "0.1",
                     "Test for h_curl space (Dim=2 Order=1)",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "cecile.daversin@lncmi.cnrs.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

using namespace Feel;

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
    TestHCurl()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        meshSize( doption("gmsh.hsize") ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        std::cout << "[TestHCurl]\n";
        std::cout << "hsize = " << meshSize << std::endl;

        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % doption("gmsh.hsize") );
    }

    /**
     * run the application
     */
    void twoElementsMesh();
    void eightElementsMesh();

private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! exporter factory
    export_ptrtype exporter;

}; //TestHCurl

void
TestHCurl::twoElementsMesh()
{
    auto geo_name = "two_elements_mesh.geo"; //create the mesh and load it
    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type,
                                  _filename=geo_name);

    auto Xh = Ned1h<0>( mesh );
    auto u = Xh->element();
    auto phi = Xh->element();
    for( auto const& dof : Xh->dof()->localDof() )
        {
            LOG(INFO) << "test local dof element " << dof.first.elementId() << " id:" << dof.first.localDof()
                      << " global dof : " << dof.second.index() << " pts: " << Xh->dof()->dofPoint( dof.second.index() ).get<0>() << std::endl;
        }

    // assembly curl(curl(u)) + u
    auto a = form2( _test=Xh, _trial=Xh);
    a = integrate(elements( mesh ), curlxt(u)*curlx(phi) + trans(idt(u))*id(phi) );
    a.matrix().printMatlab( "mass_assembly.m" );

    //Check each component of the mass matrix
    BOOST_CHECK_CLOSE( a.matrix()(0,0), 4./3., 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(1,0), 0.5, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(2,0), 0.5, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(3,0), -0.5, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(4,0), -0.5, 1e-10 );

    BOOST_CHECK_CLOSE( a.matrix()(0,1), 0.5, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(1,1), 5./6., 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(2,1), 1./3., 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(3,1), 0, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(4,1), 0, 1e-10 );

    BOOST_CHECK_CLOSE( a.matrix()(0,2), 0.5, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(1,2), 1./3., 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(2,2), 5./6., 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(3,2), 0, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(4,2), 0, 1e-10 );

    BOOST_CHECK_CLOSE( a.matrix()(0,3), -0.5, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(1,3), 0, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(2,3), 0, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(3,3), 5./6., 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(4,3), 1./3., 1e-10 );

    BOOST_CHECK_CLOSE( a.matrix()(0,4), -0.5, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(1,4), 0, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(2,4), 0, 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(3,4), 1./3., 1e-10 );
    BOOST_CHECK_CLOSE( a.matrix()(4,4), 5./6., 1e-10 );

    // Find u = l2_proj( (1,1) ) on Hcurl space
    // \int u_proj*phi = \int u_ex*phi \forall phi \in \Hcurl
    auto u_cst = Xh->element();
    phi.zero();
    auto u_cst_exact = vec( cst(1.), cst(1.) );
    auto a_cst = form2( _test=Xh, _trial=Xh );
    a_cst = integrate( _range=elements(mesh), _expr = trans(idt(u_cst))*id(phi) );
    auto f_cst = form1( _test=Xh );
    f_cst = integrate( _range=elements(mesh), _expr = trans(u_cst_exact)*id(phi) );
    a_cst.solve( _solution=u_cst, _rhs=f_cst );

    BOOST_CHECK_SMALL( u_cst(0), 1e-10 );
    BOOST_CHECK_CLOSE( u_cst(1), -2, 1e-10 );
    BOOST_CHECK_CLOSE( u_cst(2), 2, 1e-10 );
    BOOST_CHECK_CLOSE( u_cst(3), -2, 1e-10 );
    BOOST_CHECK_CLOSE( u_cst(4), 2, 1e-10 );

    // Find u = l2_proj( (-y,x) ) on Hcurl space
    // \int u_proj*phi = \int u_ex*phi \forall phi \in \Hcurl
    auto u_yx = Xh->element();
    phi.zero();
    auto u_yx_exact = vec( -Py(), Px() );
    auto a_yx = form2( _test=Xh, _trial=Xh );
    a_yx = integrate( _range=elements(mesh), _expr = trans(idt(u_yx))*id(phi) );
    auto f_yx = form1( _test=Xh );
    f_yx = integrate( _range=elements(mesh), _expr = trans(u_yx_exact)*id(phi) );
    a_yx.solve( _solution=u_yx, _rhs=f_yx );

    BOOST_CHECK_SMALL( u_yx(0), 1e-10 );
    BOOST_CHECK_CLOSE( u_yx(1), 2, 1e-10 );
    BOOST_CHECK_CLOSE( u_yx(2), 2, 1e-10 );
    BOOST_CHECK_CLOSE( u_yx(3), 2, 1e-10 );
    BOOST_CHECK_CLOSE( u_yx(4), 2, 1e-10 );

    // Find u = l2_proj( (1-y^2,1-x^2) ) on Hcurl space
    // \int u_proj*phi = \int u_ex*phi \forall phi \in \Hcurl
    // Coarse mesh (2 elements, 5 dofs) : interpolant is not exact (u_exact is 2nd degree)
    // => Projection of u_exact is the constant vector (2/3, 2/3)
    auto u_pb1=Xh->element();
    phi.zero();
    auto u_pb1_exact = vec( 1-Py()*Py(), 1-Px()*Px() );
    auto a_pb1 = form2( _test=Xh, _trial=Xh );
    a_pb1 = integrate( _range=elements(mesh), _expr = trans(idt(u_pb1))*id(phi) );
    auto f_pb1 = form1( _test=Xh );
    f_pb1 = integrate( _range=elements(mesh), _expr = trans(u_pb1_exact)*id(phi) );
    a_pb1.solve( _solution=u_pb1, _rhs=f_pb1 );

    BOOST_CHECK_SMALL( f_pb1.vector()(0), 1e-10 );
    BOOST_CHECK_CLOSE( f_pb1.vector()(1), -2./3., 1e-10 );
    BOOST_CHECK_CLOSE( f_pb1.vector()(2), 2./3., 1e-10 );
    BOOST_CHECK_CLOSE( f_pb1.vector()(3), -2./3., 1e-10 );
    BOOST_CHECK_CLOSE( f_pb1.vector()(4), 2./3., 1e-10 );

    BOOST_CHECK_SMALL( u_pb1(0), 1e-10 );
    BOOST_CHECK_CLOSE( u_pb1(1), -4./3., 1e-10 );
    BOOST_CHECK_CLOSE( u_pb1(2), 4./3., 1e-10 );
    BOOST_CHECK_CLOSE( u_pb1(3), -4./3., 1e-10 );
    BOOST_CHECK_CLOSE( u_pb1(4), 4./3., 1e-10 );

    // Finer mesh, interpolant is not exact but has to be more precise
    mesh_ptrtype refined_mesh = loadMesh(_mesh = new mesh_type);
    auto Xhf = Ned1h<0>( refined_mesh );
    auto u_pb1f = Xhf->element();
    auto phi_pb1f = Xhf->element();

    auto a_pb1f = form2( _test=Xhf, _trial=Xhf );
    a_pb1f = integrate( _range=elements(refined_mesh), _expr = trans(idt(u_pb1f))*id(phi_pb1f) );
    auto f_pb1f = form1( _test=Xhf );
    f_pb1f = integrate( _range=elements(refined_mesh), _expr = trans(u_pb1_exact)*id(phi_pb1f) );
    a_pb1f.solve( _solution=u_pb1f, _rhs=f_pb1f, _rebuild=true);

    // Check exportation of these Hcurl elements
    export_ptrtype exporterProj( export_type::New( this->vm(),std::string("test_L2proj_2elts")) );
    exporterProj->step( 0 )->setMesh( refined_mesh );
    exporterProj->step( 0 )->add( "vec_cst_proj_2elts", u_cst );
    exporterProj->step( 0 )->add( "vec_yx_proj_2elts", u_yx );
    exporterProj->step( 0 )->add( "vec_pb1_coarse_proj_2elts", u_pb1 );
    exporterProj->step( 0 )->add( "vec_pb1_fine_proj_2elts", u_pb1f );
    exporterProj->save();
}

void
TestHCurl::eightElementsMesh()
{
    auto msh_name = "eight_elements_mesh.msh"; //create the mesh and load it
    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type,
                                  _filename=msh_name);

    auto Xh = Ned1h<0>( mesh );
    auto u = Xh->element();
    auto phi = Xh->element();
    for( auto const& dof : Xh->dof()->localDof() )
        {
            LOG(INFO) << "test local dof element " << dof.first.elementId() << " id:" << dof.first.localDof()
                      << " global dof : " << dof.second.index() << " pts: " << Xh->dof()->dofPoint( dof.second.index() ).get<0>() << std::endl;
        }

    auto u_cst = Xh->element();
    phi.zero();
    auto u_cst_exact = vec( cst(1.), cst(1.) );
    auto a_cst = form2( _test=Xh, _trial=Xh );
    a_cst = integrate( _range=elements(mesh), _expr = trans(idt(u_cst))*id(phi) );
    auto f_cst = form1( _test=Xh );
    f_cst = integrate( _range=elements(mesh), _expr = trans(u_cst_exact)*id(phi) );
    a_cst.solve( _solution=u_cst, _rhs=f_cst, _rebuild=true );

    //Check matrix assembly (non zero components)
    std::vector<int> diag_48 = {0,1,2,3,5,6,7,9,12,13,14,15};
    std::vector<int> diag_24 = {4,8,10,11};
    BOOST_FOREACH( int comp, diag_48 )
        {
            BOOST_TEST_MESSAGE( "a(N" << comp << ",N" << comp << ") = " << a_cst.matrix()(comp,comp) << "\n");
            BOOST_CHECK_CLOSE( a_cst.matrix()(comp,comp), 1./3., 1e-10 );
        }
    BOOST_FOREACH( int comp, diag_24 )
        {
            BOOST_TEST_MESSAGE( "a(N" << comp << ",N" << comp << ") = " << a_cst.matrix()(comp,comp) << "\n");
            BOOST_CHECK_CLOSE( a_cst.matrix()(comp,comp), 2./3., 1e-10 );
        }

    typedef std::pair<int,int> map_value_type;
    std::map<int,int> other_comp_inf1 ={{1,2},{5,4},{8,4},{8,7},{11,10},{15,14}};
    std::map<int,int> other_comp_inf2 ={{11,12},{10,13}};
    BOOST_FOREACH( map_value_type comp, other_comp_inf1)
        {
            BOOST_TEST_MESSAGE( "a(N" << comp.first << ",N" << comp.second << ") = " << a_cst.matrix()(comp.first,comp.second) << "\n");
            BOOST_CHECK_CLOSE( a_cst.matrix()(comp.first,comp.second), -1./6., 1e-10 );
            BOOST_CHECK_CLOSE( a_cst.matrix()(comp.second,comp.first), -1./6., 1e-10 );
        }
    BOOST_FOREACH( map_value_type comp, other_comp_inf2)
        {
            BOOST_TEST_MESSAGE( "a(N" << comp.first << ",N" << comp.second << ") = " << a_cst.matrix()(comp.first,comp.second) << "\n");
            BOOST_CHECK_CLOSE( a_cst.matrix()(comp.first,comp.second), 1./6., 1e-10 );
            BOOST_CHECK_CLOSE( a_cst.matrix()(comp.second,comp.first), 1./6., 1e-10 );
        }

    auto u_yx = Xh->element();
    phi.zero();
    auto u_yx_exact = vec( -Py(), Px() );
    auto a_yx = form2( _test=Xh, _trial=Xh );
    a_yx = integrate( _range=elements(mesh), _expr = trans(idt(u_yx))*id(phi) );
    auto f_yx = form1( _test=Xh );
    f_yx = integrate( _range=elements(mesh), _expr = trans(u_yx_exact)*id(phi) );
    a_yx.solve( _solution=u_yx, _rhs=f_yx, _rebuild=true);

    auto u_pb1=Xh->element();
    phi.zero();
    auto u_pb1_exact = vec( 1-Py()*Py(), 1-Px()*Px() );
    auto a_pb1 = form2( _test=Xh, _trial=Xh );
    a_pb1 = integrate( _range=elements(mesh), _expr = trans(idt(u_pb1))*id(phi) );
    auto f_pb1 = form1( _test=Xh );
    f_pb1 = integrate( _range=elements(mesh), _expr = trans(u_pb1_exact)*id(phi) );
    a_pb1.solve( _solution=u_pb1, _rhs=f_pb1, _rebuild=true);

    mesh_ptrtype refined_mesh = loadMesh(_mesh = new mesh_type);
    auto Xhf = Ned1h<0>( refined_mesh );
    auto u_pb1f = Xhf->element();
    auto phi_pb1f = Xhf->element();

    auto a_pb1f = form2( _test=Xhf, _trial=Xhf );
    a_pb1f = integrate( _range=elements(refined_mesh), _expr = trans(idt(u_pb1f))*id(phi_pb1f) );
    auto f_pb1f = form1( _test=Xhf );
    f_pb1f = integrate( _range=elements(refined_mesh), _expr = trans(u_pb1_exact)*id(phi_pb1f) );
    a_pb1f.solve( _solution=u_pb1f, _rhs=f_pb1f, _rebuild=true);

    // Check exportation of these Hcurl elements
    export_ptrtype exporterProj( export_type::New( this->vm(),std::string("test_L2proj_8elts")) );
    exporterProj->step( 0 )->setMesh( refined_mesh );
    exporterProj->step( 0 )->add( "vec_cst_proj_8elts", u_cst );
    exporterProj->step( 0 )->add( "vec_yx_proj_8elts", u_yx );
    exporterProj->step( 0 )->add( "vec_pb1_coarse_proj_8elts", u_pb1 );
    exporterProj->step( 0 )->add( "vec_pb1_fine_proj_8elts", u_pb1f );
    exporterProj->save();
}
}


#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( space )

BOOST_AUTO_TEST_CASE( test_hcurl_two_elements )
{
    BOOST_TEST_MESSAGE( "test_hcurl_N0  : assembly on two elements mesh" );
    Feel::TestHCurl t;
    t.twoElementsMesh();
    BOOST_TEST_MESSAGE( "test_hcurl_N0  : assembly on two elements mesh done" );
}

BOOST_AUTO_TEST_CASE( test_hcurl_eight_elements )
{
    BOOST_TEST_MESSAGE( "test_hcurl_N0  : assembly on eight elements mesh" );
    Feel::TestHCurl t;
    t.eightElementsMesh();
    BOOST_TEST_MESSAGE( "test_hcurl_N0  : assembly on eight elements mesh done" );
}

BOOST_AUTO_TEST_SUITE_END()
#else

int
main( int argc, char* argv[] )
{
    Feel::Environment env( argc,argv,
                           makeAbout(), makeOptions() );

    Feel::TestHCurl app_hcurl;

    app_hcurl.twoElementsMesh();
    app_hcurl.eightElementsMesh();
}

#endif
