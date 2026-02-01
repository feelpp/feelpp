/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 21 août 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#define BOOST_TEST_MODULE test_interpolation_grad
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/testsuite.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/dh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

namespace bdata = boost::unit_test::data;

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Options" );
    options.add_options()
        ( "grad",po::value<std::vector<std::string>>()->multitoken(),
         "list of functions to test the Grad interpolation operator (default : {\"1\"}" )
        ( "curl",po::value<std::vector<std::string>>()->multitoken(),
          "list of functions to test the Curl interpolation operator (default : {\"{1,1,1}\"}" )
        ( "div",po::value<std::vector<std::string>>()->multitoken(),
          "list of functions to test the Div interpolation operator (default : {\"{1,1,1}\"}" )

        ;
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_interpolation_grad" ,
                     "test_interpolation_grad" ,
                     "0.2",
                     "test_interpolation_grad",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );
    
    return about;
}

template<int Dim>
struct TestFixture
{
    using mesh_type = Mesh<Simplex<Dim>>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
  
    //! Hcurl space
    using curl_space_type = Ned1h_type<mesh_type,0>;
    using curl_space_ptrtype = Ned1h_ptrtype<mesh_type,0>;

    //! DT space
    using rt_space_type = Dh_type<mesh_type,0>;
    using rt_space_ptrtype = Dh_ptrtype<mesh_type,0>;

    //! Pch space
    using lag_space_type = Pch_type<mesh_type,1>;
    using lag_space_ptrtype = Pch_ptrtype<mesh_type,1>;

    //! Pch 0 space
    using lag_0_space_type = Pdh_type<mesh_type,0>;
    using lag_0_space_ptrtype = Pdh_ptrtype<mesh_type,0>;

    //! Projection operators
    using i_type = I_t<lag_space_type, lag_space_type>;
    using grad_type = Grad_t<lag_space_type, curl_space_type>;
    using curl_type = Curl_t<curl_space_type, rt_space_type>;
    using div_type = Div_t<rt_space_type, lag_0_space_type>;

    /// Mesh
    mesh_ptrtype mesh;
  
    /// Spaces
    lag_space_ptrtype Xh;
    curl_space_ptrtype Gh;
    rt_space_ptrtype Ch;
    lag_0_space_ptrtype P0h;

    /// Projections - using shared_ptr to avoid default construction issues
    std::shared_ptr<grad_type> Igrad;
    std::shared_ptr<curl_type> Icurl;
    std::shared_ptr<div_type> Idiv;
    
    TestFixture()
    {
        mesh = loadMesh( _mesh=new mesh_type );
        Xh = Pch<1>(mesh);
        Gh = Ned1h<0>(mesh);
        Ch = Dh<0>(mesh);
        P0h = Pdh<0>(mesh);
        Igrad = std::make_shared<grad_type>( Grad( _domainSpace = Xh, _imageSpace=Gh ) );
        Icurl = std::make_shared<curl_type>( Curl( _domainSpace = Gh, _imageSpace=Ch ) );
        Idiv = std::make_shared<div_type>( Div( _domainSpace = Ch, _imageSpace=P0h ) );
    }
};

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

// Global fixture shared across all test cases to avoid re-creating mesh
struct GlobalFixture 
{
    GlobalFixture() 
    {
        BOOST_TEST_MESSAGE("Setting up global test fixture");
        fixture = std::make_shared<TestFixture<3>>();
    }
    
    ~GlobalFixture()
    {
        BOOST_TEST_MESSAGE("Tearing down global test fixture");
        // Explicitly reset the fixture BEFORE Feel++ tears down MPI
        // to avoid MPI operations after MPI_Finalize
        fixture.reset();
    }
    
    static std::shared_ptr<TestFixture<3>> fixture;
};

std::shared_ptr<TestFixture<3>> GlobalFixture::fixture;

BOOST_GLOBAL_FIXTURE( GlobalFixture );

BOOST_AUTO_TEST_SUITE( test_interpolation_grad )

// Test data for grad operator
auto grad_test_data = bdata::make( 
    Environment::vm().count("grad") ? vsoption(_name="grad") : std::vector<std::string>{"1"} 
);

// Test data for curl operator  
auto curl_test_data = bdata::make(
    Environment::vm().count("curl") ? vsoption(_name="curl") : std::vector<std::string>{"{1,1,1}:x:y:z"}
);

// Test data for div operator
auto div_test_data = bdata::make(
    Environment::vm().count("div") ? vsoption(_name="div") : std::vector<std::string>{"{1,1,1}:x:y:z"}
);

BOOST_DATA_TEST_CASE( test_grad_operator, grad_test_data, test_function )
{
    auto& fixture = *GlobalFixture::fixture;
    BOOST_TEST_MESSAGE( "Testing grad operator with function: " << test_function );
    
    auto u = fixture.Xh->element( expr(test_function), "u_grad", test_function );
    auto w = (*fixture.Igrad)(u);
    
    auto v = fixture.Gh->element( trans(grad<3>(expr(test_function))), "grad_u", 
                         str(grad<3>(expr(test_function))) );
    
    auto const errL2 = normL2( _range=elements(fixture.mesh), _expr=idv(w)-idv(v) );
    
    BOOST_TEST_MESSAGE( "errL2( pi_h grad(u) = grad(pi_h(u)) ): " << errL2 );
    BOOST_CHECK_SMALL( errL2, 1e-12 );
}

BOOST_DATA_TEST_CASE( test_curl_operator, curl_test_data, test_function )
{
    auto& fixture = *GlobalFixture::fixture;
    BOOST_TEST_MESSAGE( "Testing curl operator with function: " << test_function );
    
    auto u = fixture.Gh->element( expr<3,1>(test_function), "u_curl", test_function );
    auto w = (*fixture.Icurl)(u);
    
    auto v = fixture.Ch->element( curl(expr<3,1>(test_function)), "curl_u",
                         str(curl(expr<3,1>(test_function))) );
    
    auto const errL2 = normL2( _range=elements(fixture.mesh), _expr=idv(w)-idv(v) );
    
    BOOST_TEST_MESSAGE( "errL2( pi_h curl(u) = curl(pi_h(u)) ): " << errL2 );
    BOOST_CHECK_SMALL( errL2, 1e-12 );
}

BOOST_DATA_TEST_CASE( test_div_operator, div_test_data, test_function )
{
    auto& fixture = *GlobalFixture::fixture;
    BOOST_TEST_MESSAGE( "Testing div operator with function: " << test_function );
    
    auto u = fixture.Ch->element( expr<3,1>(test_function), "u_div", test_function );
    auto w = (*fixture.Idiv)(u);
    
    auto v = fixture.P0h->element( div(expr<3,1>(test_function)), "div_u",
                          str(div(expr<3,1>(test_function))) );
    
    auto const errL2 = normL2( _range=elements(fixture.mesh), _expr=idv(w)-idv(v) );
    
    BOOST_TEST_MESSAGE( "errL2( pi_h div(u) = div(pi_h(u)) ): " << errL2 );
    BOOST_CHECK_SMALL( errL2, 1e-12 );
}
BOOST_AUTO_TEST_SUITE_END()


