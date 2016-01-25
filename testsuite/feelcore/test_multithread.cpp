/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Alexandre Ancel <alexandre.ancel@cemosis.fr>
       Date: 2015-04-14

  Copyright (C) 2014 Cemosis, Université de Strasbourg

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
/**
   \file test_multithread.cpp
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2015-04-14
 */

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/ginac.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#if defined(FEELPP_HAS_OPENMP)
#include <omp.h>
#endif

#if defined(FEELPP_HAS_HARTS)
#include "RunTimeSystem/Model/RunTimeSysEnv.h"
#include "RunTimeSystem/DataMng/DataHandler.h"
#include "RunTimeSystem/TaskMng/TaskMng.h"
#include "RunTimeSystem/TaskMng/AsynchTask.h"
#include "RunTimeSystem/TaskMng/StdScheduler.h"
#include "RunTimeSystem/TaskMng/StdDriver.h"
#include "RunTimeSystem/TaskMng/TBBScheduler.h"
#include "RunTimeSystem/TaskMng/TBBDriver.h"
#include "RunTimeSystem/TaskMng/PTHScheduler.h"
#include "RunTimeSystem/TaskMng/PTHDriver.h"
#include "RunTimeSystem/DataMng/DataArgs.h"
#include "Utils/PerfTools/PerfCounterMng.h"
#endif

using namespace Feel;


typedef Mesh<Simplex<3>> MeshType;

double getElapsedTime(struct timespec start, struct timespec end)
{
    double ts, te;
    ts = (double)(start.tv_sec) + (double)(start.tv_nsec) / (double)(1e9); 
    te = (double)(end.tv_sec) + (double)(end.tv_nsec) / (double)(1e9); 

    return (te - ts);
}

#if defined(FEELPP_HAS_OPENMP)
template<typename MeshTypePtr>
void omp_integrate(MeshTypePtr mesh, int nCores)
{
    double seqRes, seqTime;
    struct timespec ts1, ts2;
    double * coarseParRes, * coarseParTime;
    double fineParRes, fineParTime;

    /* store option values */
    bool parCPUEnable = boption( _name="parallel.cpu.enable" );
    std::string parCPUImpl = soption( _name="parallel.cpu.impl" );
    int parCPURestrict = ioption( _name="parallel.cpu.restrict" );
    
    /* Set the number of threads to use */
    omp_set_num_threads(nCores);

    coarseParRes = new double[nCores];
    coarseParTime = new double[nCores];

    /* sequential with a constant */
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);

    auto intfCst = integrate( _range = elements( mesh ), _expr = cst(1.0) ).evaluate(false);
    seqRes = intfCst(0 , 0);

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
    seqTime = getElapsedTime(ts1, ts2);

    /* parallel with constant */
#pragma omp parallel
    {
        struct timespec tts1, tts2;
        clock_gettime(CLOCK_MONOTONIC_RAW, &tts1);

        auto intf = integrate( _range = elements( mesh ), _expr = cst(1.0) ).evaluate(false);
        coarseParRes[omp_get_thread_num()] = intf(0, 0);

        clock_gettime(CLOCK_MONOTONIC_RAW, &tts2);
        coarseParTime[omp_get_thread_num()] = getElapsedTime(tts1, tts2);
    }

    for(int i = 0; i < nCores; i++)
    { CHECK( coarseParRes[i] == seqRes ) << "Test failed (" << __FUNCTION__ << ") in thread " << i << "(" << coarseParRes[i] << " != " << seqRes << ")" << std::endl; }

    /* paralel implementation in Feel++ with a constant */
    Environment::setOptionValue("parallel.cpu.enable", bool(true));
    Environment::setOptionValue("parallel.cpu.impl", std::string("openmp"));
    Environment::setOptionValue("parallel.cpu.restrict", int(nCores));

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);

    auto intfCstPar = integrate( _range = elements( mesh ), _expr = cst(1.0) ).evaluate(false);
    fineParRes = intfCstPar(0 , 0);

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
    fineParTime = getElapsedTime(ts1, ts2);

    Environment::setOptionValue("parallel.cpu.enable", parCPUEnable);
    Environment::setOptionValue("parallel.cpu.impl", parCPUImpl);
    Environment::setOptionValue("parallel.cpu.restrict", parCPURestrict);

    CHECK( fineParRes  == seqRes ) << "Test failed (" << __FUNCTION__ << ") (" << fineParRes << " != " << seqRes << ")" << std::endl;

    std::cout << "Test OK: seqTime=" << seqTime << ", coarseParTime=";
    if(nCores > 0)
    {
        std::cout << coarseParTime[0];
        for(int i = 1; i < nCores; i++)
        { std::cout << ", " << coarseParTime[i]; }
    }
    std::cout << "), fineParTime=" << fineParTime << std::endl;

    /* sequential with an expression */
    // our function to integrate
    auto g = expr( "sin(pi*x)*cos(pi*y):x:y" );

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);

    auto intfExpr = integrate( _range = elements( mesh ), _expr = g ).evaluate(false);
    seqRes = intfExpr(0 , 0);

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
    seqTime = getElapsedTime(ts1, ts2);

#pragma omp parallel
    {
        struct timespec tts1, tts2;
        clock_gettime(CLOCK_MONOTONIC_RAW, &tts1);

        // our function to integrate
        auto g = expr( "sin(pi*x)*cos(pi*y):x:y" );
        auto intf = integrate( _range = elements( mesh ), _expr = g ).evaluate(false);
        coarseParRes[omp_get_thread_num()] = intf(0, 0);

        clock_gettime(CLOCK_MONOTONIC_RAW, &tts2);
        coarseParTime[omp_get_thread_num()] = getElapsedTime(tts1, tts2);
    }

    for(int i = 0; i < nCores; i++)
    { CHECK( coarseParRes[i] == seqRes ) << "Test failed (" << __FUNCTION__ << ") in thread " << i << "(" << coarseParRes[i] << " != " << seqRes << ")\n"; }

    std::cout << "Test OK: seqTime=" << seqTime << ", parTime=(";
    if(nCores > 0)
    {
        std::cout << coarseParTime[0];
        for(int i = 1; i < nCores; i++)
        { std::cout << ", " << coarseParTime[i]; }
    }
    std::cout << ")" << std::endl;

    delete[] coarseParRes;
    coarseParRes = nullptr;
    delete[] coarseParTime;
    coarseParTime = nullptr;
}
#endif

#if defined(FEELPP_HAS_HARTS)
template<typename MeshTypePtr>
void harts_integrate(MeshTypePtr mesh, int nCores)
{
    typedef typename Feel::detail::quadptlocrangetype< decltype(elements(mesh)) >::type range_type;
    typedef decltype(cst(1.0)) expr_type;
    typedef _Q< ExpressionOrder<range_type,expr_type>::value > quad_type;
    typedef _Q< ExpressionOrder<range_type,expr_type>::value_1 > quad1_type;
    typedef Integrator<range_type, quad_type, expr_type, quad1_type> int_type;
}
#endif

int main(int argc, char ** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="test_multithread",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    // create the mesh (specify the dimension of geometric entity)
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<3>> );

    int nCores = 4;//omp_get_num_threads();
    std::cout << "Using " << nCores << " threads" << std::endl;

#ifdef FEELPP_HAS_OPENMP
    omp_integrate(mesh, nCores);
#else
    std::cout << "OpenMP is not activated. Skipping associated tests." << std::endl;
#endif 

#if 0
#ifdef FEELPP_HAS_HARTS
    harts_integrate(mesh, nCores);
#else
    std::cout << "HARTS is not activated. Skipping associated tests." << std::endl;
#endif 
#endif

    return 0;
}
