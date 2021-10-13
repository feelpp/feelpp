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

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Options" );
    options.add_options()
        ( "grad",po::value<std::vector<std::string>>()->default_value( {"1"} ),
         "list of functions to test the Grad interpolation operator (default : {\"1\"}" )
        ( "curl",po::value<std::vector<std::string>>()->default_value( {"{1,1,1}:x:y:z"} ),
          "list of functions to test the Curl interpolation operator (default : {\"1,1,1\"}" )
        ( "div",po::value<std::vector<std::string>>()->default_value( {"{1,1,1}:x:y:z"} ),
          "list of functions to test the Div interpolation operator (default : {\"1,1,1\"}" )

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
class Test:
    public Simget
{
    typedef Mesh<Simplex<Dim>> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
  
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
    

    //! Pchv space
    using lag_v_space_type = Pchv_type<mesh_type,1>;
    using lag_v_space_ptrtype = Pchv_ptrtype<mesh_type,1>;
  
    //! Projection 
    //Id 
    typedef I_t<lag_space_type, lag_space_type> i_type;
    typedef I_ptr_t<lag_space_type, lag_space_type> i_ptrtype;
    //Grad 
    typedef Grad_t<lag_space_type, curl_space_type> grad_type;
    typedef Grad_ptr_t<lag_space_type, curl_space_type> grad_ptrtype;
    //Curl
    typedef Curl_t<curl_space_type, rt_space_type> curl_type;
    typedef Curl_ptr_t<curl_space_type, rt_space_type> curl_ptrtype;
    //Div 
    typedef Div_t<rt_space_type, lag_0_space_type> div_type;
    typedef Div_ptr_t<rt_space_type, lag_0_space_type> div_ptrtype;
private:
  
  /// Mesh
  mesh_ptrtype mesh;
  
  /// Spaces
  lag_space_ptrtype Xh;
  curl_space_ptrtype Gh;
  rt_space_ptrtype Ch;
  lag_0_space_ptrtype P0h;

  /// Projections
  i_type Ih;
  grad_type Igrad;
  curl_type Icurl;
  div_type Idiv;

public :
    
    void run()
        {
            mesh = loadMesh( _mesh=new mesh_type );
            Xh = Pch<1>(mesh);
            Gh = Ned1h<0>(mesh);
            Ch = Dh<0>(mesh);
            P0h = Pdh<0>(mesh);
            Ih = I( _domainSpace = Xh, _imageSpace=Xh );
            Igrad = Grad( _domainSpace = Xh, _imageSpace=Gh );
            Icurl = Curl( _domainSpace = Gh, _imageSpace=Ch );
            Idiv = Div( _domainSpace = Ch, _imageSpace=P0h );
            auto e = exporter(_mesh=mesh);
            
            int i = 0;
            for( auto f : vsoption( _name="grad" ) )
            {
                
                std::string n { str(format("u%1%")%i) };
                auto u = Xh->element( expr(f), n, f );

                std::string Ign { str(format("Igrad_u%1%")%i) };
                e->add(n,u);
                auto w = Igrad(u);
                
                if ( Environment::isSequential() )
                {
                    Igrad.matPtr()->printMatlab("Igrad.m");
                    w.printMatlab("w.m");
                }
                e->add(Ign,w);
                
                std::string gn = str(boost::format("grad_u%1%")%i);
                auto v = Gh->element(trans(grad<Dim>(expr(f))), gn, str(grad<Dim>(expr(f))));
                auto errL2 = normL2( _range=elements(mesh), _expr=idv(w)-idv(v) );
                BOOST_TEST_MESSAGE( "errL2( pi_h grad(u)=grad(pi_h(u)):" << errL2 );
                e->add(gn,v);
                ++i;
            }
            for( auto f : vsoption( _name="curl" ) )
            {
                
                std::string n { str(format("u%1%")%i) };
                auto u = Gh->element( expr<Dim,1>(f), n, f );

                std::string Idn { str(format("Icurl_u%1%")%i) };
                e->add(n,u);
                auto w = Icurl(u);
                
                if ( Environment::isSequential() )
                {
                    Icurl.matPtr()->printMatlab("Icurl.m");
                    w.printMatlab("wcurl.m");
                    u.printMatlab("ucurl.m");
                }
                e->add(Idn,w);
                
                std::string dn = str(boost::format("curl_u%1%")%i);
                auto v = Ch->element(curl(expr<Dim,1>(f)), dn, str(curl(expr<Dim,1>(f))));
                auto errL2 = normL2( _range=elements(mesh), _expr=idv(w)-idv(v) );
                BOOST_TEST_MESSAGE( "errL2( pi_h curl(u)=curl(pi_h(u)):" << errL2 );
                //BOOST_CHECK_SMALL( errL2, 1e-12 );
                e->add(dn,v);
                ++i;
            }
            for( auto f : vsoption( _name="div" ) )
            {
                
                std::string n { str(format("u%1%")%i) };
                auto u = Ch->element( expr<Dim,1>(f), n, f );

                std::string Idn { str(format("Idiv_u%1%")%i) };
                e->add(n,u);
                auto w = Idiv(u);
                
                if ( Environment::isSequential() )
                {
                    Idiv.matPtr()->printMatlab("Idiv.m");
                    w.printMatlab("wdiv.m");
                    u.printMatlab("udiv.m");
                }
                e->add(Idn,w);
                
                std::string dn = str(boost::format("div_u%1%")%i);
                auto v = P0h->element(div(expr<Dim,1>(f)), dn, str(div(expr<Dim,1>(f))));
                auto errL2 = normL2( _range=elements(mesh), _expr=idv(w)-idv(v) );
                BOOST_TEST_MESSAGE( "errL2( pi_h div(u)=div(pi_h(u)):" << errL2 );
                //BOOST_CHECK_SMALL( errL2, 1e-12 );
                e->add(dn,v);
                ++i;
            }
            e->save();
        }
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( test_interpolation_grad )


BOOST_AUTO_TEST_CASE( test )
{
    Test<3> test;
    test.run();
}


BOOST_AUTO_TEST_SUITE_END()


