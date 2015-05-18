/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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

#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE multithread
#include <testsuite/testsuite.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

#if 0
template<value_type>
class Task
{
    public:
        void computeCPU(DataArgsType& args)
        {

        }

        value_type result() const
        {
            return M_ret;
        }
    private:
        value_type M_ret;
}
#endif

void reduceAndCompare( double * array, int n, double refval )
{
    double total = 0.0;

    for( int i = 0; i < n; i++ )
    {
        total += array[i];
    }

    std::cout << total << " " << refval << std::endl;
    BOOST_CHECK_CLOSE_FRACTION(total, refval, 0.00000001);
}

#if 1
BOOST_AUTO_TEST_SUITE(HARTS)

BOOST_AUTO_TEST_CASE(integrate_coarse)
{
    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()

#endif

#if defined(FEELPP_HAS_OPENMP)
BOOST_AUTO_TEST_SUITE(OpenMP)

BOOST_AUTO_TEST_CASE(integrate_coarse)
{
    using namespace Feel;

    typedef Mesh<Simplex<3>> MeshType;
    // create the mesh (specify the dimension of geometric entity)
    auto mesh = loadMesh( _mesh=new MeshType );

    typedef typename Feel::detail::quadptlocrangetype< decltype(elements(mesh)) >::type range_type;
    typedef decltype(cst(1.0)) expr_type;
    typedef _Q< ExpressionOrder<range_type,expr_type>::value > quad_type;
    typedef _Q< ExpressionOrder<range_type,expr_type>::value_1 > quad1_type;
    typedef Integrator<range_type, quad_type, expr_type, quad1_type> int_type;

    //int nbcores = 4;
    //omp_set_num_threads(nbcores);
    int nbcores = omp_get_num_threads();
    std::cout << "Using " << nbcores << " threads" << std::endl;

    //std::cout << "Using " << nbcores << " threads" << std::endl;
    //std::cout << boption( _name="parallel.cpu.enable" ) << " " << soption( _name="parallel.cpu.impl") << " " << ioption( _name="parallel.cpu.restrict") << std::endl;

    struct timespec ts1;
    struct timespec ts2;
    double t1, t2;

    /* sequential with a constant */
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);

    auto intf1 = integrate( _range = elements( mesh ),
            _expr = cst(1.0) ).evaluate(false);
    std::cout << "sequential res= " << intf1( 0, 0 ) << std::cout;

#if 1
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
    t1 = (double)(ts1.tv_sec) + (double)(ts1.tv_nsec) / (1000000000.0); 
    t2 = (double)(ts2.tv_sec) + (double)(ts2.tv_nsec) / (1000000000.0); 
    std::cout << Environment::worldComm().globalRank() << " sequential (MONOTONIC_RAW)=" <<  (t2 - t1) << std::endl;

    /* parallel with constant */
    double * res = new double[nbcores];

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);

#pragma omp parallel
    {
        //std::cout << omp_get_thread_num() << std::endl;
        auto intf = integrate( _range = elements( mesh ),
                                     _expr = cst(1.0) ).evaluate(false);
        res[omp_get_thread_num()] = intf(0, 0);
        std::cout << omp_get_thread_num() << " res: " << intf( 0, 0) << std::endl;
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
    t1 = (double)(ts1.tv_sec) + (double)(ts1.tv_nsec) / (1000000000.0); 
    t2 = (double)(ts2.tv_sec) + (double)(ts2.tv_nsec) / (1000000000.0); 
    std::cout << Environment::worldComm().globalRank() << " parallel max (MONOTONIC_RAW)=" <<  (t2 - t1) << std::endl;

    //reduceAndCompare( res, nbcores, nbcores * 1.0 );
#endif
#if 0

    /* sequential with an expression */
    // our function to integrate
    auto g = expr( "sin(pi*x)*cos(pi*y):x:y" );

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);
    auto intf2 = integrate( _range = elements( mesh ),
                                 _expr = g ).evaluate(false);
    std::cout << "sequential res= " << intf2( 0, 0 ) << std::cout;

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
    t1 = (double)(ts1.tv_sec) + (double)(ts1.tv_nsec) / (1000000000.0); 
    t2 = (double)(ts2.tv_sec) + (double)(ts2.tv_nsec) / (1000000000.0); 
    std::cout << Environment::worldComm().globalRank() << " sequential (MONOTONIC_RAW)=" <<  (t2 - t1) << std::endl;

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);

#pragma omp parallel
    {
        // our function to integrate
        auto g = expr( "sin(pi*x)*cos(pi*y):x:y" );
        //std::cout << omp_get_thread_num() << std::endl;
        auto intf = integrate( _range = elements( mesh ),
                                     _expr = g ).evaluate(false);
        res[omp_get_thread_num()] = intf(0, 0);
        std::cout << omp_get_thread_num() << " res: " << intf(0, 0) << std::endl;
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &ts2);
    t1 = (double)(ts1.tv_sec) + (double)(ts1.tv_nsec) / (1000000000.0); 
    t2 = (double)(ts2.tv_sec) + (double)(ts2.tv_nsec) / (1000000000.0); 
    std::cout << Environment::worldComm().globalRank() << " elapsed (MONOTONIC_RAW)=" <<  (t2 - t1) << std::endl;
#endif

    delete[] res;
    res = nullptr;
}

BOOST_AUTO_TEST_CASE(integrate_fine)
{
    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
#endif
