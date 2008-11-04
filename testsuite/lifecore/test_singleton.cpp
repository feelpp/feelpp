/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <prudhomm@mit.edu>
       Date: 2004-09-10

  Copyright (C) 2004 EPFL

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
   \file test_singleton.cpp
   \author Christophe Prud'homme <prudhomm@mit.edu>
   \date 2004-09-10
 */
#define BOOST_TEST_MODULE singleton test
#include <boost/test/unit_test.hpp>

#include <iostream>

#include <life/lifecore/singleton.hpp>




class A
{
public:
    A()
        :
        _M_value( -1 )
        {}
    double value() const { return _M_value; }

private:
    double _M_value;
};

typedef Life::Singleton<A> aSingleton;

BOOST_AUTO_TEST_CASE( singleton_test )
{

    BOOST_CHECK_EQUAL( aSingleton::instance().value(), -1 );

}
