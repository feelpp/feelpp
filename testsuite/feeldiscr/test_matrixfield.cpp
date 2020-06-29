/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 18 Apr 2015

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

#define BOOST_TEST_MODULE test_matrixfield
#include <feel/feelcore/testsuite.hpp>


#include <feel/options.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/redux.hpp>
#include <feel/feelvf/eig.hpp>
#include <feel/feelfilters/geotool.hpp>


#include <iostream>

using namespace Feel;

namespace test_matrixfield
{

typedef Application Application_type;
typedef std::shared_ptr<Application_type> Application_ptrtype;

/*_________________________________________________*
 * Options
 *_________________________________________________*/

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_matrixfield options" );
    desc_options.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ;
    return desc_options.add( Feel::feel_options() );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "Test_Matrixfield" ,
                     "Test_Matrixfield" ,
                     "0.1",
                     "test matrix fields",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}


}

FEELPP_ENVIRONMENT_WITH_OPTIONS( test_matrixfield::makeAbout(), test_matrixfield::makeOptions() );

BOOST_AUTO_TEST_SUITE( interp_matrixfield )

typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( interp_matrixfield, T, dim_types )
{
    typedef Mesh<Simplex<T::value,1,T::value> > mesh_type;
    typedef std::shared_ptr<  mesh_type > mesh_ptrtype;

    double meshSize = doption("gmsh.hsize");

    auto create_mesh = [&](){
                           if constexpr ( T::value == 2 )
                           {
                               GeoTool::Node x1( 0,0 );
                               GeoTool::Node x2( 4,1 );
                               GeoTool::Rectangle R( meshSize,"OMEGA",x1,x2 );
                               return R;
                           }
                           else
                           {
                               GeoTool::Node x1( 0,0,0 );
                               GeoTool::Node x2( 1,2,1);
                               GeoTool::Cube R( meshSize,"OMEGA",x1,x2 );
                               return R;
                           }
                       };
    auto R = create_mesh();
    auto mesh = R.createMesh(_mesh=new mesh_type,_name= "domain" );
    auto Mh = Pchm<2>( mesh) ;
    auto Xh = Pch<2>( mesh) ;
    auto v = Mh->element();
    auto u = Xh->element(), w = Xh->element();
    //

    auto check_interp = [&]( auto const& e ){
                            v.on(_range=elements(mesh), _expr=e ) ;
                            auto i0 = integrate( _range=elements(mesh), _expr=idv(v) ).evaluate();
                            auto i0_res = integrate( _range=elements(mesh), _expr=e ).evaluate();
                            BOOST_TEST_MESSAGE( "i0 = " << i0  );
                            BOOST_TEST_MESSAGE( "i0_res = " << i0_res );
                            BOOST_CHECK_SMALL( (i0-i0_res).norm(), 1e-13 );
                        };
    if constexpr( T::value == 2 )
    {
        check_interp( mat<2,2>( cst(1.), cst(1.), cst(1.), cst(1.)  ) );
        check_interp( mat<2,2>( cst(1.), cst(2.), cst(3.), cst(4.)  ) );
    }
    if constexpr( T::value == 3 )
    {
        check_interp( ones<3,3>() );
        check_interp( P()*trans(P()) );
    }


    auto check_sum_trace = [&]( auto const& e, auto const& expected_res ){
                               v.on(_range=elements(mesh), _expr=e ) ;
                               u.on(_range=elements(mesh), _expr=sum(eig(idv(v))) ) ;
                               w.on(_range=elements(mesh), _expr=trace(idv(v))  );
                               double l2error  = normL2(_range=elements(mesh), _expr=idv(u)-idv(w));
                               BOOST_TEST_MESSAGE( "l2error = " << l2error  );
                               BOOST_CHECK_SMALL( l2error, 1e-13 );
                               u.on(_range=elements(mesh), _expr=sum(trans(eig(idv(v)))) ) ;
                               double l2error_trans  = normL2(_range=elements(mesh), _expr=idv(u)-idv(w));
                               BOOST_TEST_MESSAGE( "l2error trans = " << l2error_trans  );
                               BOOST_CHECK_SMALL( l2error_trans, 1e-13 );
                               double int_u = integrate( _range=elements( mesh ), _expr= idv(u)  ).evaluate()( 0, 0 );
                               double int_exp = integrate( _range=elements( mesh ), _expr= expected_res  ).evaluate()( 0, 0 );
                               BOOST_CHECK_CLOSE( int_u, int_exp, 1e-13 );

                           };
    if constexpr ( T::value == 2 )
    {
        check_sum_trace( mat<2,2>( cst(1.), cst(0.), cst(0.), cst(4.) ), cst(5.) );
        check_sum_trace( mat<2,2>( Px(), cst(0.), cst(0.), Py()  ), sum(P()) );
    }
    else
    {
        check_sum_trace( eye<3,3>(), cst(3.) );
        check_sum_trace( mat<3,3>( Px(), cst(0.), cst(2.),
                                   cst(0.),  Py(), cst(0.),
                                   cst(2.),  cst(0.), Pz()
                                   ), sum(P()) );
    }   
}
BOOST_AUTO_TEST_SUITE_END()

