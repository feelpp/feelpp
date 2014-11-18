/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
Date: 2007-12-19

Copyright (C) 2007-2012 Université Joseph Fourier (Grenoble I)

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
  \file test_interpolation.cpp
  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  \date 2007-12-19
  */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE interpolation testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <testsuite/testsuite.hpp>
#include <feel/feel.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelvf/print.hpp>

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
    po::options_description testHcurlInterpolationOptions( "test h_div options" );
    testHcurlInterpolationOptions.add_options()
        ( "meshes-3d", po::value< std::vector<std::string> >(), "vector containing mesh names" )
        ;
    return testHcurlInterpolationOptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_interpolation_Nedelec3d" ,
                     "test_interpolation_Nedelec3d" ,
                     "0.1",
                     "Test for interpolation with h_curl space",
                     AboutData::License_GPL,
                     "Copyright (c) 2014 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

    //using namespace Feel;

class TestInterpolationHCurl3D
    :
public Application
{
    typedef Application super;

public :
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
    typedef Mesh< Simplex<2,1,3> > submesh2d_type;
    typedef Mesh< Simplex<1,1,3> > submesh1d_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef boost::shared_ptr<submesh2d_type> submesh2d_ptrtype;
    typedef boost::shared_ptr<submesh1d_type> submesh1d_ptrtype;

    //! the basis type of our approximation space
    typedef bases< Nedelec<0,NedelecKind::NED1> > basis_type;
    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef FunctionSpace<submesh2d_type, Nedelec<0,NedelecKind::NED1> > subspace2d_type;
    typedef FunctionSpace<submesh1d_type, Lagrange<1,Scalar> > subspace1d_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef boost::shared_ptr<subspace2d_type> subspace2d_ptrtype;
    typedef boost::shared_ptr<subspace1d_type> subspace1d_ptrtype;
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    typedef Exporter<submesh2d_type> export2d_type;
    typedef Exporter<submesh1d_type> export1d_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef boost::shared_ptr<export2d_type> export2d_ptrtype;
    typedef boost::shared_ptr<export1d_type> export1d_ptrtype;

    /**
     * Constructor
     */
    TestInterpolationHCurl3D()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        this->changeRepository( boost::format( "%1%/" )
                                % this->about().appName()
                              );
    }

    void testInterpolation( std::string one_element_mesh );

private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! exporter factory
    export_ptrtype exporter;

};

void
TestInterpolationHCurl3D::testInterpolation( std::string one_element_mesh )
{
    //auto myexpr = unitX() + unitY() + unitZ() ; //(1,1,1)
    auto myexpr = vec( cst(1.), cst(1.), cst(1.));

    // one element mesh
    auto mesh_name = one_element_mesh + ".msh"; //create the mesh and load it
    fs::path mesh_path( mesh_name );

    mesh_ptrtype oneelement_mesh = loadMesh( _mesh=new mesh_type,
                                             _filename=mesh_name);

    // refined mesh (export)
    auto refine_level = std::floor(1 - math::log( 0.1 )); //Deduce refine level from meshSize (option)
    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type,
                                  _filename=mesh_name,
                                  _refine=( int )refine_level);

    space_ptrtype Xh = space_type::New( oneelement_mesh );

    std::vector<std::string> faces = {"yzFace","xyzFace","xyFace"};
    std::vector<std::string> edges = {"zAxis","yAxis","yzAxis","xyAxis","xzAxis","xAxis"};

    element_type U_h_int = Xh->element();
    element_type U_h_on = Xh->element();
    element_type U_h_on_boundary = Xh->element();

    submesh1d_ptrtype edgeMesh( new submesh1d_type );
    edgeMesh = createSubmesh(oneelement_mesh, boundaryedges(oneelement_mesh) ); //submesh of edges

    // Tangents on ref element
    auto t0 = vec(cst(0.),cst(0.),cst(-2.));
    auto t1 = vec(cst(0.),cst(2.),cst(0.));
    auto t2 = vec(cst(0.),cst(-2.),cst(2.));
    auto t3 = vec(cst(2.),cst(-2.),cst(0.));
    auto t4 = vec(cst(2.),cst(0.),cst(-2.));
    auto t5 = vec(cst(2.),cst(0.),cst(0.));

    // Jacobian of geometrical transforms
    std::string jac;
    if(mesh_path.stem().string() == "one-elt-ref-3d" || mesh_path.stem().string() == "one-elt-real-homo-3d" )
        jac = "{1,0,0,0,1,0,0,0,1}:x:y:z";
    else if(mesh_path.stem().string() == "one-elt-real-rotx" )
        jac = "{1,0,0,0,0,-1,0,1,0}:x:y:z";
    else if(mesh_path.stem().string() == "one-elt-real-roty" )
        jac = "{0,0,1,0,1,0,-1,0,0}:x:y:z";
    else if(mesh_path.stem().string() == "one-elt-real-rotz" )
        jac = "{0,-1,0,1,0,0,0,0,1}:x:y:z";

    U_h_int(0) = integrate( markedelements(edgeMesh, edges[0]), trans(expr<3,3>(jac)*t0)*myexpr ).evaluate()(0,0);
    U_h_int(1) = integrate( markedelements(edgeMesh, edges[1]), trans(expr<3,3>(jac)*t1)*myexpr ).evaluate()(0,0);
    U_h_int(2) = integrate( markedelements(edgeMesh, edges[2]), trans(expr<3,3>(jac)*t2)*myexpr ).evaluate()(0,0);
    U_h_int(3) = integrate( markedelements(edgeMesh, edges[3]), trans(expr<3,3>(jac)*t3)*myexpr ).evaluate()(0,0);
    U_h_int(4) = integrate( markedelements(edgeMesh, edges[4]), trans(expr<3,3>(jac)*t4)*myexpr ).evaluate()(0,0);
    U_h_int(5) = integrate( markedelements(edgeMesh, edges[5]), trans(expr<3,3>(jac)*t5)*myexpr ).evaluate()(0,0);

    for(int i=0; i<edges.size(); i++)
        {
            double edgeLength = integrate( markedelements(edgeMesh, edges[i]), cst(1.) ).evaluate()(0,0);
            U_h_int(i) /= edgeLength;
        }

#if 0 //Doesn't work for now
    for(int i=0; i<Xh->nLocalDof(); i++)
        {
            CHECK( edgeMesh->hasMarkers( {edges[i]} ) );
            U_h_int(i) = integrate( markedelements(edgeMesh, edges[i]), trans( print(T(),"T=") )*myexpr ).evaluate()(0,0);
            std::cout << "U_h_int(" << i << ")= " << U_h_int(i) << std::endl;
        }
#endif

    // nedelec interpolant using on keyword
    // interpolate on element
    U_h_on.zero();
    U_h_on.on(_range=elements(oneelement_mesh), _expr=myexpr);
    U_h_on_boundary.on(_range=boundaryfaces(oneelement_mesh), _expr=myexpr);

    export_ptrtype exporter_proj( export_type::New( this->vm(),
                                  ( boost::format( "%1%" ) % this->about().appName() ).str() ) );

    exporter_proj->step( 0 )->setMesh( mesh );
    exporter_proj->step( 0 )->add( "U_interpolation_handly_"+mesh_path.stem().string(), U_h_int );
    exporter_proj->step( 0 )->add( "U_interpolation_on_"+mesh_path.stem().string(), U_h_on );
    exporter_proj->save();

    // print coefficient only for reference element
    U_h_int.printMatlab( "U_h_int_" + mesh_path.stem().string() + ".m" );
    U_h_on.printMatlab( "U_h_on_" + mesh_path.stem().string() + ".m" );
    U_h_on_boundary.printMatlab( "U_h_on_boundary_" + mesh_path.stem().string() + ".m" );

    //L2 norm of error
    auto error = vf::project(_space=Xh, _range=elements(oneelement_mesh), _expr=idv(U_h_int) - idv(U_h_on) );
    double L2error = error.l2Norm();
    std::cout << "L2 error (elements) = " << L2error << std::endl;

    auto error_boundary = vf::project(_space=Xh, _range=elements(oneelement_mesh), _expr=idv(U_h_int) - idv(U_h_on_boundary) );
    double L2error_boundary = error_boundary.l2Norm();
    std::cout << "L2 error (boundary) = " << L2error_boundary << std::endl;
    BOOST_CHECK_SMALL( L2error_boundary - L2error, 1e-13 );
}
}// namespace Feel

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( HCURL_INTERPOLANT )

BOOST_AUTO_TEST_CASE( test_hcurl_interpolant_1 )
{
    using namespace Feel;

    TestInterpolationHCurl3D t3;
    std::vector<std::string> mygeoms3d = vsoption("meshes-3d"); //option(_name="meshes-3d").template as< std::vector<std::string> >();
    for(std::string geo3d : mygeoms3d)
        {
            std::cout << "*** interpolant on " << geo3d << " *** \n";
            t3.testInterpolation(geo3d);
        }
}
BOOST_AUTO_TEST_SUITE_END()
#else

int
main( int argc, char* argv[] )
{
    Feel::Environment env( argc,argv,
                           makeAbout(), makeOptions() );
    Feel::TestInterpolationHCurl3D app_hcurl;
    std::vector<std::string> mygeoms = option(_name="meshes-3d").template as< std::vector<std::string> >();
    for(std::string geo : mygeoms)
        {
            app_hcurl.testInterpolant(geo);
        }
}

#endif
