/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Guillaume Dolle <dolle.guillaume@gmail.com>
       Date: 2018-08-21

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
#define BOOST_TEST_MODULE mongo test
#include <testsuite.hpp>
#include <feel/feelcore/mongocxx.hpp>

BOOST_AUTO_TEST_CASE( test_mongo )
{
    using namespace Feel;

#if defined(FEELPP_HAS_MONGOCXX)
    MongoCxx::instance(); // Instance created.
    MongoCxx::instance(); // No instance created.

    // Note: mongo gives an error, if two different instance are created.
#endif
}


