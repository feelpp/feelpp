/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Goncalo Pena <gpena@mat.uc.pt>
       Date: 2011-01-16

  Copyright (C) 2011 University of Coimbra

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
/**
   \file test_elementminmax.cpp
   \author Goncalo Pena <gpena@mat.uc.pt>
   \date 2011-01-16
 */
#define USE_BOOST_TEST 1
// Boost.Test

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE elementminmax testsuite
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN


#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/minmax.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( elementminmaxsuite )

template<typename mesh_type>
void check( std::shared_ptr<mesh_type> const& mesh, std::map<std::string,double>&& minmax )
{
    typedef FunctionSpace<mesh_type,bases<Lagrange<1, Scalar> > > space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    space_ptrtype Xh = space_type::New( mesh );
    auto u = Xh->element();
    u = vf::project( _space=Xh, _range=elements( mesh ), _expr=vf::sqrt( Px()*Px() + Py()*Py() + Pz()*Pz() ) );

    double maximum = u.max();
    double minimum = u.min();
    auto [min_r,arg_min_r] = minelt(_range=elements(mesh), _element=u);
    auto [max_r,arg_max_r] = maxelt(_range=elements(mesh), _element=u);
    auto minmax_r = minmaxelt(_range=elements(mesh), _element=u);
    BOOST_CHECK_CLOSE( maximum, minmax["max"], 2e-1 );
    BOOST_CHECK_CLOSE( max_r, minmax["max"], 2e-1 );
    BOOST_CHECK_CLOSE( std::get<0>(minmax_r[1]), minmax["max"], 2e-1 );
    if ( std::abs(minmax["min"]) < 1e-10 )
    {
        BOOST_CHECK_SMALL( minimum, 2e-1 );
        BOOST_CHECK_SMALL( min_r, 2e-1 );
        BOOST_CHECK_SMALL( std::get<0>(minmax_r[0]), 2e-1 );
    }
    else
    {
        BOOST_CHECK_CLOSE( minimum, minmax["min"], 2e-1 );
        BOOST_CHECK_CLOSE( min_r, minmax["min"], 2e-1 );
        BOOST_CHECK_CLOSE( std::get<0>( minmax_r[0] ), minmax["min"], 2e-1 );
    }
    auto check_with_marker = [&]( std::string m )
    {
        if ( mesh->hasMarker( m ) && mesh->markerDim( m ) == mesh_type::nRealDim )
        {
            auto [min_r, arg_min_r] = minelt( _range = markedelements( mesh, m ), _element = u );
            auto minmax_r = minmaxelt( _range = markedelements( mesh, m ), _element = u );
            if ( std::abs( minmax[fmt::format("min{}",m)] ) < 1e-10 )
            {
                BOOST_CHECK_SMALL( min_r, 2e-1  );
                BOOST_CHECK_SMALL( std::get<0>(minmax_r[0]), 2e-1 );
            }
            else
            {
                BOOST_CHECK_CLOSE( min_r, minmax[fmt::format("min{}",m)], 2e-1 );
                BOOST_CHECK_CLOSE( std::get<0>(minmax_r[0]), minmax[fmt::format("min{}",m)], 2e-1 );
            }
            auto [max_r, arg_max_r] = maxelt( _range = markedelements( mesh,m ), _element = u );
            BOOST_CHECK_CLOSE( max_r, minmax[fmt::format( "max{}", m )] , 2e-1 );
            BOOST_CHECK_CLOSE( std::get<0>(minmax_r[1]), minmax[fmt::format( "max{}", m )] , 2e-1 );
        }
        // facet tests
        if ( mesh->hasMarker( m ) && mesh->markerDim( m ) == mesh_type::nRealDim-1 )
        {
            BOOST_TEST_MESSAGE( "Checking facet minmax for marker " << m );
            auto minmax_r = minmaxelt( _range = markedfaces( mesh, m ), _element = u );
            if ( std::abs( minmax[fmt::format("min{}",m)] ) < 1e-10 )
            {
                //BOOST_CHECK_SMALL( min_r, 2e-1  );
                BOOST_CHECK_SMALL( std::get<0>(minmax_r[0]), 2e-1 );
            }
            else
            {
                //BOOST_CHECK_CLOSE( min_r, minmax[fmt::format("min{}",m)], 2e-1 );
                BOOST_CHECK_CLOSE( std::get<0>(minmax_r[0]), minmax[fmt::format("min{}",m)], 2e-1 );
            }
            //auto [max_r, arg_max_r] = maxelt( _range = markedelements( mesh,m ), _element = u );
            //BOOST_CHECK_CLOSE( max_r, minmax[fmt::format( "max{}", m )] , 2e-1 );
            BOOST_CHECK_CLOSE( std::get<0>(minmax_r[1]), minmax[fmt::format( "max{}", m )] , 2e-1 );
        }
    };
    BOOST_TEST_MESSAGE( "check with marker Inner" );
    check_with_marker( "Inner" );   // check the inner part of the mesh
    BOOST_TEST_MESSAGE( "check with marker Outer" );
    check_with_marker( "Outter" );  // check the outter marker
    BOOST_TEST_MESSAGE( "check with marker InnerBoundary" );
    check_with_marker( "InnerBoundary" );  // check the outter marker
    BOOST_TEST_MESSAGE( "check with marker OutterBoundary" );
    check_with_marker( "OutterBoundary" );  // check the outter marker

}

typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( elementminmax1, T, dim_types )
{
    BOOST_TEST_MESSAGE( "check function min and max for dim = " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;


    // The mesh created is either the [0,1] interval,
    // the triangle obtained by halving [0,1]^2,
    // or the tetrahedron obtained in a similar way from the cube [0,1]^3
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=fmt::format( "elementminmax1-{}", T::value ),
                                                      _shape="simplex",
                                                      _h=doption("gmsh.hsize") ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );

    check<mesh_type>( mesh, {{"min",0.},{"max",1.}} );
}
using dim_t2 = boost::mpl::list<boost::mpl::int_<2>, boost::mpl::int_<3>>;
BOOST_AUTO_TEST_CASE_TEMPLATE( elementminmax2, T, dim_t2 )
{
    BOOST_TEST_MESSAGE( "check function min and max for dim = " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;


    // The mesh created is either the [0,1] interval,
    // the triangle obtained by halving [0,1]^2,
    // or the tetrahedron obtained in a similar way from the cube [0,1]^3
    mesh_ptrtype mesh = loadMesh(_mesh=new mesh_type,
                                 _filename=fmt::format("test_elementminmax{}d.geo", T::value) );

    double max_v = ( T::value == 3 )?0.25:std::sqrt( 2 * 0.5 * 0.5);
    double maxInner_v = 0.25;
    double minInnerBoundary = maxInner_v;
    double maxInnerBoundary = maxInner_v;
    double minOutterBoundary = std::sqrt(0.5*0.5);
    double maxOutterBoundary = max_v;
    check<mesh_type>( mesh, { { "min", 0. }, { "max", max_v }, { "minInner", 0. }, { "maxInner", maxInner_v }, { "minOutter", maxInner_v }, { "maxOutter", max_v }, 
                              {"minInnerBoundary",minInnerBoundary}, {"maxInnerBoundary",maxInnerBoundary},
                              {"minOutterBoundary",minOutterBoundary}, {"maxOutterBoundary",maxOutterBoundary} } );
}
BOOST_AUTO_TEST_CASE_TEMPLATE( elementminmax3, T, dim_types )
{
    BOOST_TEST_MESSAGE( "check broken function min and max for dim = " << T::value << "\n" );
    typedef Mesh<Simplex<T::value,1> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;


    // The mesh created is either the [0,1] interval,
    // the triangle obtained by halving [0,1]^2,
    // or the tetrahedron obtained in a similar way from the cube [0,1]^3
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "elementminmax3-%1%" )  % T::value ).str() ,
                                                _shape="simplex",
                                                _dim=T::value,
                                                _h=0.2 ),
                                        _update=MESH_CHECK|MESH_UPDATE_EDGES|MESH_UPDATE_FACES );


    typedef FunctionSpace<mesh_type,bases<Lagrange<4, Scalar> > > space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    space_ptrtype Xh = space_type::New( mesh );
    auto u = Xh->element();
    u = vf::project( _space=Xh, _range=elements( mesh ), _expr=vf::sqrt( Px()*Px() + Py()*Py() + Pz()*Pz() ) );

    typedef FunctionSpace<mesh_type,bases<Lagrange<0, Scalar, Discontinuous> > > p0_space_type;
    typedef std::shared_ptr<p0_space_type> p0_space_ptrtype;
    p0_space_ptrtype P0h = p0_space_type::New( mesh );

    auto maxPerElem = P0h->element();
    auto minPerElem = P0h->element();
    maxPerElem = u.max( P0h );
    minPerElem = u.min( P0h );

    double maximum = maxPerElem.max();
    double minimum = minPerElem.min();

    const double eps = 1000*Feel::type_traits<double>::epsilon();

#if defined(USE_BOOST_TEST)
    BOOST_CHECK_CLOSE( maximum, 1.0, 2e-1 );
    BOOST_CHECK_SMALL( minimum, 2e-1 );
#endif

    std::cout << "\n";
}


BOOST_AUTO_TEST_SUITE_END()


