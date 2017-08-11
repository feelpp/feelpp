//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 03 Feb 2017
//! @copyright 2017 Feel++ Consortium
//!
//!
#define BOOST_TEST_MODULE test_space_cst
#include <testsuite/testsuite.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

/** use Feel namespace */
using namespace Feel;

inline
po::options_description makeOptions()
{
    po::options_description options( "Test space with constant functions Options" );
    return options;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_space_cst" ,
                     "test_space_cst" ,
                     "0.2",
                     "test spaces of constant functions",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2017 Feel++ Consortium" );

    about.addAuthor( "C Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

 FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
 BOOST_AUTO_TEST_SUITE( productspace_suite )


 BOOST_AUTO_TEST_CASE( test1 )
 {
     using namespace Feel;
     using Feel::cout;

     auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
     //std::cout << "nelts mesh: " << nelements(elements(mesh),true) << std::endl;
     //std::cout << "nelts mesh/proc: " << nelements(elements(mesh)) << std::endl;
     auto Ch=Pch<0>(mesh);
     auto Xh=Pch<1>(mesh);
     auto u = Xh->element(cst(1.),"u");

     auto face_mesh = createSubmesh( mesh, boundaryfaces(mesh), EXTRACTION_KEEP_MESH_RELATION, 0 );
     auto Cbh= Pch<0>(face_mesh);
     //std::cout << "nelts faces_mesh: " << nelements(elements(face_mesh),true) << std::endl;
     //std::cout << "nelts faces_mesh/proc: " << nelements(elements(face_mesh)) << std::endl;
     auto ub=Cbh->element(cst(1.),"ub");

     auto l = form1(_test=Cbh);
     l = integrate(_range=elements(face_mesh),_expr=id(ub));
     auto lb = l(ub);
     BOOST_CHECK_CLOSE( lb, 4, 1e-12 );
     BOOST_TEST_MESSAGE( "l(ub)=" << lb  );
     //cout << "l(u)=" << l(ub) << std::endl;

     auto a = form2(_trial=Cbh, _test=Cbh);
     a = integrate(_range=elements(face_mesh),_expr=id(ub)*idt(ub));
     auto ea = a( ub, ub );
     BOOST_CHECK_CLOSE( ea, 4, 1e-12 );
     BOOST_TEST_MESSAGE( "a(u,u)=" << ea  );

     auto b = form2(_trial=Xh, _test=Cbh);
     b = integrate(_range=elements(face_mesh),_expr=idt(u)*id(ub));
     auto eb = b( ub, u );
     BOOST_TEST_MESSAGE( "b(ub,u)=" << eb  );
     BOOST_CHECK_CLOSE( eb, 4, 1e-12 );


}

BOOST_AUTO_TEST_SUITE_END()
