/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Jan 2016

 Copyright (C) 2016 Feel++ Consortium

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
#define BOOST_TEST_MODULE test_forms
#include <testsuite/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feeldiscr/pdh.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test Forms  Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_forms" ,
                     "test_forms" ,
                     "0.2",
                     "nD(n=2,3) test forms",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}



FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( forms_suite )

//using dim_t = boost::mpl::list<boost::mpl::int_<2>, boost::mpl::int_<3> >;
using dim_t = boost::mpl::list<boost::mpl::int_<2>>;


BOOST_AUTO_TEST_CASE_TEMPLATE( test_form2_faces, T, dim_t )
{
    BOOST_MESSAGE( "test_form2_faces starts for dim=" << T::value);
    auto meshnd = unitHypercube<T::value>();
    auto mesh = createSubmesh( meshnd, faces(meshnd), EXTRACTION_KEEP_MESH_RELATION, 0 );
    size_type nFaceInParallelMeshnd = nelements(faces(meshnd),true) - nelements(interprocessfaces(meshnd),true)/2;
    BOOST_CHECK_EQUAL( nelements(elements(mesh),true), nFaceInParallelMeshnd  );
    auto Vh=Pdhv<1>(meshnd);
    auto u = Vh->element();
    auto Wh=Pdh<1>(meshnd);
    auto p = Wh->element();
    p.on(_range=elements(meshnd),_expr=cst(1.));
    auto Mh=Pdh<1>(mesh);
    auto l=Mh->element();
    l.on(_range=elements(meshnd),_expr=cst(1.));

    auto a = form2(_test=Mh, _trial=Wh );
    a = integrate( _range=internalfaces(mesh), _expr=id(l)*(leftfacet(idt(p))+rightfacet(idt(p))));
    a.close();

    if ( Environment::isMasterRank() )
        std::cout << "a(1,1) = " << a.energy(l,p) << std::endl;

    BOOST_MESSAGE( "test_form2_faces ends for dim=" << T::value);
}


BOOST_AUTO_TEST_SUITE_END()
