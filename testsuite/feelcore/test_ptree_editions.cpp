/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): JB Wahl <wahl.jb@gmail.com>
 Date: 27 sept. 2015

 Copyright (C) 2018 Feel++ Consortium

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

#define BOOST_TEST_MODULE test_ptree_editions
#include <feel/feelcore/testsuite.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault("test_ptree_editions"),
                                 ptree_options() )
BOOST_AUTO_TEST_SUITE( ptree_editions )

BOOST_AUTO_TEST_CASE( test_0 )
{
    pt::ptree p;

    std::string filename = Environment::expand("$cfgdir/test_ptree_editions.json");

    auto json_str_wo_comments = removeComments(readFromFile(filename));
    std::istringstream istr( json_str_wo_comments );
    pt::read_json(istr, p);

    std::string my_value = p.get<std::string>( "Field1.subfield.myvariable" );
    BOOST_CHECK_EQUAL( my_value, "value1" );

    editPtreeFromOptions( p );
    my_value = p.get<std::string>( "Field1.subfield.myvariable" );
    BOOST_CHECK_EQUAL( my_value, "value2" );
}


BOOST_AUTO_TEST_SUITE_END()
