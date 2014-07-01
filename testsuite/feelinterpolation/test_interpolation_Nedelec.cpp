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
        ( "meshes", po::value< std::vector<std::string> >(), "vector containing mesh names" )
        ;
    return testHcurlInterpolationOptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_interpolation_hcurl" ,
                     "test_interpolation_hcurl" ,
                     "0.1",
                     "Test for interpolation with h_curl space",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

    //using namespace Feel;

template<int Dim>
class TestInterpolationHCurl
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
    typedef Simplex<Dim,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    typedef bases< Nedelec<0,NedelecKind::NED1> > basis_type;
    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    TestInterpolationHCurl()
        :
        super(),
        M_backend( backend_type::build( this->vm() ) ),
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

template<int Dim>
void
TestInterpolationHCurl<Dim>::testInterpolation( std::string one_element_mesh )
{
    // expr to interpolate
    auto myexpr = unitX(); //(1,0)

    // one element mesh
    auto mesh_name = one_element_mesh + ".msh"; //create the mesh and load it
    mesh_ptrtype oneelement_mesh = loadMesh( _mesh=new mesh_type,
                                             _filename=mesh_name);

    // refined mesh (export)
    auto refine_level = std::floor(1 - math::log( 0.1 )); //Deduce refine level from meshSize (option)
    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type,
                                  _filename=mesh_name,
                                  _refine=( int )refine_level);

    space_ptrtype Xh = space_type::New( oneelement_mesh );
    std::cout << "nb dof = " << Xh->nDof() << std::endl;
    std::vector<std::string> edges{ "hypo","vert","hor"}; //list of faces

    element_type U_h_int = Xh->element();
    element_type U_h_on = Xh->element();

    // handly computed interpolant coeff (in hdiv basis)
    for ( int i = 0; i < Xh->nLocalDof(); ++i )
        U_h_int(i) = integrate( markedfaces( oneelement_mesh, edges[i] ), trans(print(T()))*myexpr ).evaluate()(0,0);

    // raviart-thomas interpolant using on
    U_h_on.zero();
    U_h_on.on(_range=elements(oneelement_mesh), _expr=myexpr);

    export_ptrtype exporter_proj( export_type::New( this->vm(),
                                  ( boost::format( "%1%" ) % this->about().appName() ).str() ) );

    exporter_proj->step( 0 )->setMesh( mesh );
    exporter_proj->step( 0 )->add( "U_interpolation_handly"+one_element_mesh, U_h_int );
    exporter_proj->step( 0 )->add( "U_interpolation_on"+one_element_mesh, U_h_on );
    exporter_proj->save();

    std::cout << "one_elt_mesh = " << one_element_mesh << std::endl;
    if(one_element_mesh == "$datadir/gmsh/one-elt-meshes/one-elt-ref")
        {
            U_h_int.printMatlab( "U_h_int.m" );
            U_h_on.printMatlab( "U_h_on.m" );
        }

    //L2 norm of error
    auto error = vf::project(_space=Xh, _range=elements(oneelement_mesh), _expr=idv(U_h_int) - idv(U_h_on) );
    double L2error = error.l2Norm();
    std::cout << "L2 error  = " << L2error << std::endl;
}
}// namespace Feel

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( HCURL_INTERPOLANT )

BOOST_AUTO_TEST_CASE( test_hcurl_interpolant_1 )
{
    using namespace Feel;
    TestInterpolationHCurl<2> t;
    std::vector<std::string> mygeoms = option(_name="meshes").template as< std::vector<std::string> >();
    for(std::string geo : mygeoms)
        {
            BOOST_TEST_MESSAGE( "*** interpolant on " << geo << " ***" );
            t.testInterpolation(geo);
        }
}
BOOST_AUTO_TEST_SUITE_END()
#else

int
main( int argc, char* argv[] )
{
    Feel::Environment env( argc,argv,
                           makeAbout(), makeOptions() );
    Feel::TestInterpolationHCurl app_hcurl;
    std::vector<std::string> mygeoms = option(_name="meshes").template as< std::vector<std::string> >();
    for(std::string geo : mygeoms)
        {
            app_hcurl.testInterpolant(geo);
        }
}

#endif
