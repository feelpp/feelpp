/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 27 sept. 2015

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

#define BOOST_TEST_MODULE test_removecomments
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/ptreetools.hpp>



FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( removecomments )

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace Feel;
    std::string s = "int main() /* toto */ { // a comment \n int i = 0; /* tutu */ }\n";
    auto res = removeComments(s);
    BOOST_TEST_MESSAGE( "remove comment of \n" << s << "\n gives\n" << res );
    BOOST_CHECK_EQUAL( res, "int main()  { \n int i = 0;  }\n" );
}

BOOST_AUTO_TEST_SUITE_END()
