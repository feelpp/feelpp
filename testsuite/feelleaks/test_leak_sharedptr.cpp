/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-11-23

  Copyright (C) 2006-2009 Université Joseph Fourier (Grenoble I)

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file laplacian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-11-23
 */
#include <feel/feel.hpp>
#include <boost/smart_ptr/make_shared.hpp>

#if defined(FEELPP_HAS_GPERFTOOLS)
#include <gperftools/heap-checker.h>
#endif /* FEELPP_HAS_GPERFTOOLS */

class A
{
public:
    ~A()
        {
            std::cerr << "Destructor\n";
        }
};
//boost::shared_ptr<A> f() { return new A; }
boost::shared_ptr<A> g() { return boost::shared_ptr<A>(new A); }

int main(int argc, char**argv )
{
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check0("checker 0");
#endif /* FEELPP_HAS_GPERFTOOLS */
    {
    }
#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(check0.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */


#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check3("checker 3");
#endif /* FEELPP_HAS_GPERFTOOLS */
#if 1
    {
    boost::shared_ptr<A> m( new A );
    CHECK( m.use_count() == 1 ) << "Invalid shared_ptr, count: " << m.use_count();
    std::cout << "before\n";
    m.reset();
    std::cout << "after\n";
    CHECK( m.use_count() == 0 ) << "Invalid shared_ptr, count: " << m.use_count();
    }
    {
    boost::shared_ptr<A> m( g() );
    CHECK( m.use_count() == 1 ) << "Invalid shared_ptr, count: " << m.use_count();
    m.reset();
    CHECK( m.use_count() == 0 ) << "Invalid shared_ptr, count: " << m.use_count();
    }
    {
    boost::shared_ptr<A> m( g() );
    CHECK( m.use_count() == 1 ) << "Invalid shared_ptr, count: " << m.use_count();
    m.reset();
    CHECK( m.use_count() == 0 ) << "Invalid shared_ptr, count: " << m.use_count();
    }
#endif
#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(check3.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */
}
