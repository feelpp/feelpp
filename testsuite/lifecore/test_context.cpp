/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-04-07

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file test_context.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-04-07
 */
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
#define BOOST_TEST_MODULE context test
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <life/lifecore/life.hpp>
#include <life/lifecore/context.hpp>
#include <life/lifealg/enums.hpp>



BOOST_AUTO_TEST_CASE( test_context )
{
    using namespace Life;
    Context ctx( ON_ELIMINATION | ON_ELIMINATION_KEEP_DIAGONAL );
    BOOST_CHECK_EQUAL( ctx.test( ON_PENALISATION ), false );
    BOOST_CHECK_EQUAL( ctx.test( ON_ELIMINATION ), true );
    BOOST_CHECK_EQUAL( ctx.test( ON_ELIMINATION_KEEP_DIAGONAL ), true );
    BOOST_CHECK_EQUAL( ctx.test( ON_ELIMINATION_SYMMETRIC ), false );


    ctx = ON_ELIMINATION | ON_ELIMINATION_SYMMETRIC;
    BOOST_CHECK_EQUAL( ctx.test( ON_PENALISATION ), false );
    BOOST_CHECK_EQUAL( ctx.test( ON_ELIMINATION ), true );
    BOOST_CHECK_EQUAL( ctx.test( ON_ELIMINATION_KEEP_DIAGONAL ), false );
    BOOST_CHECK_EQUAL( ctx.test( ON_ELIMINATION_SYMMETRIC ), true );

#if 0
    detail::Flags<on_context_type> oct;

    //LIFE_DECLARE_FLAGS( onct, on_context_type );


    onct |= (ON_ELIMINATION|ON_ELIMINATION_SYMMETRIC);
    BOOST_CHECK_EQUAL( onct.testFlag(ON_ELIMINATION), true );
    BOOST_CHECK_EQUAL( onct.testFlag(ON_ELIMINATION_KEEP_DIAGONAL), false );
    BOOST_CHECK_EQUAL( onct.testFlag(ON_ELIMINATION_SYMMETRIC), true );
#endif
}

