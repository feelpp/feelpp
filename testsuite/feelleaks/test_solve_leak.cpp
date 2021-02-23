/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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

#define BOOST_TEST_MODULE leak solve

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>

#if defined(FEELPP_HAS_GPERFTOOLS)
#include <gperftools/heap-checker.h>
#endif /* FEELPP_HAS_GPERFTOOLS */

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault("test_solve_leak"), Feel::backend_options("toto") )

BOOST_AUTO_TEST_SUITE( solve_leak )

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace Feel;

#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker checkere("checker");
#endif /* FEELPP_HAS_GPERFTOOLS */

    std::shared_ptr < MatrixSparse<double> > ptr_m;
    {
        for(int i = 0; i < 2; ++i )
        {
            auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
            auto Vh = Pch<1>( mesh );
            auto u = Vh->element();
            auto v = Vh->element();

            auto g = expr( soption(_name="functions.g") );
            auto syms = symbols<2>();
            auto laplacian_g = laplacian( g  );

            auto l = form1( _test=Vh );
            l = integrate(_range=elements(mesh),
                          _expr=-laplacian_g*id(v));

            auto a = form2( _trial=Vh, _test=Vh);
            a = integrate(_range=elements(mesh),
                          _expr=gradt(u)*trans(grad(v)) );
            a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr= g );
            a.solve(_rhs=l,_solution=u,_rebuild=true);

            LOG(INFO) << "pointing on matrix from bilinear form = "<<a.matrixPtr().use_count()<<std::endl;

            LOG(INFO) << " 1- L2 error norm : " << normL2( _range=elements(mesh), _expr=idv(u)-g );
            backend(_name="toto",_rebuild=true)->solve(_matrix=a.matrixPtr(),_rhs=l.vectorPtr(),_solution=u);

            LOG(INFO) << "pointing on matrix after using backend->solve() = "<<a.matrixPtr().use_count()<<std::endl;

            LOG(INFO) << " 2- L2 error norm : " << normL2( _range=elements(mesh), _expr=idv(u)-g );

            ptr_m = a.matrixPtr();
            LOG(INFO) << "pointing on ptr_m = "<<ptr_m.use_count()<<std::endl;

        }
    }

    LOG(INFO) << "pointing on ptr_m before clearSomeMemory = "<<ptr_m.use_count()<<std::endl;
    Environment::clearSomeMemory();
    LOG(INFO) << "pointing on ptr_m after clearSomeMemory = "<<ptr_m.use_count()<<std::endl;

#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(checkere.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */
}

BOOST_AUTO_TEST_SUITE_END()

