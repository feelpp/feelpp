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

int main(int argc, char**argv )
{
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check("checker -1");
#endif /* FEELPP_HAS_GPERFTOOLS */
    {
    }
#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(check.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */

    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="test_leak_mesh",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    boost::shared_ptr<Mesh<Simplex<2>>> aMesh;
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check0("checker 0");
#endif /* FEELPP_HAS_GPERFTOOLS */
    {
    }
#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(check0.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */


#if 0
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check1("checker 1");
#endif /* FEELPP_HAS_GPERFTOOLS */
    {

        for(int i = 0; i < 2; ++i )
        {
            auto mesh = loadMesh(_mesh=boost::make_shared<>(new Mesh<Simplex<2>>));
            mesh->clear();
        }
    }
#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(check1.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */


#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check2("checker 2");
#endif /* FEELPP_HAS_GPERFTOOLS */
    {
        aMesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
        aMesh->clear();
        aMesh.reset();
    }
#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(check2.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */
#endif
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check3("checker 3");
#endif /* FEELPP_HAS_GPERFTOOLS */
    {
        //boost::shared_ptr<Mesh<Simplex<2>>> m( new Mesh<Simplex<2>> );
        aMesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
        CHECK( aMesh.use_count() == 1 ) << "Invalid shared_ptr, count: " << aMesh.use_count();
        aMesh.reset();
        CHECK( aMesh.use_count() == 0 ) << "Invalid shared_ptr, count: " << aMesh.use_count();
    }
#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(check3.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */
    CHECK( aMesh.use_count() == 0 ) << "Invalid shared_ptr, count: " << aMesh.use_count();
}
