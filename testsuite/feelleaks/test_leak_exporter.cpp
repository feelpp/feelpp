/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-07-30

  Copyright (C) 2013 Université de Strasbourg

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
   \file test_leak_exporter.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-07-30
 */

#include <boost/timer.hpp>

#define BOOST_TEST_MODULE evaluate testsuite
#include <testsuite/testsuite.hpp>

#include <boost/mpl/list.hpp>

#if defined(FEELPP_HAS_GPERFTOOLS)
#include <gperftools/heap-checker.h>
#endif /* FEELPP_HAS_GPERFTOOLS */

#include <feel/feel.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( exportersuite )
//typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
typedef boost::mpl::list<boost::mpl::int_<2> > dim_types;
//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3>,boost::mpl::int_<1> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( leak1, T, dim_types )
{
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker checker("checker for exporter leaks");
#endif /* FEELPP_HAS_GPERFTOOLS */

    {
        auto mesh = unitHypercube<T::value>();
        auto Xh = Pch<2>( mesh );

        auto ex = exporter( _mesh=mesh, _name="leak1" );
        for( int i = 0; i < 10; ++i )
        {
            auto v = project( _space=Xh, _range=elements(mesh),
                              _expr=sin(pi*Px()/i));
            ex->step(i)->add( "v", v );
            ex->save();
        }
        Environment::clearSomeMemory();
    }
#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(checker.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */

}

BOOST_AUTO_TEST_SUITE_END()
