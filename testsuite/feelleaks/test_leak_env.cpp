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

#if defined(FEELPP_HAS_GPERFTOOLS)
#include <gperftools/heap-checker.h>
#endif /* FEELPP_HAS_GPERFTOOLS */

int main(int argc, char**argv )
{


#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker checkere("checker");
#endif /* FEELPP_HAS_GPERFTOOLS */
    {
        using namespace Feel;
        Environment env( _argc=argc, _argv=argv,
                         _desc=feel_options(),
                         _about=about(_name="test_leak_mesh",
                                      _author="Feel++ Consortium",
                                      _email="feelpp-devel@feelpp.org"));
    }
#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(checkere.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */

}
