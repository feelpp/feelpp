/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-20

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2010-2014 Feel++ Consortium

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
   \file integrator.hpp
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2014-07-17
 */

#ifndef FEELPP_DETAIL_INTEGRATOR_HPP
#define FEELPP_DETAIL_INTEGRATOR_HPP 1

#include <feel/feelmesh/elements.hpp>
#include <feel/feelvf/integrator.hpp>

#if defined(FEELPP_HAS_HARTS)
#include "RunTimeSystem/Model/RunTimeSysEnv.h"
#include "RunTimeSystem/DataMng/DataHandler.h"
#include "RunTimeSystem/TaskMng/TaskMng.h"
#include "RunTimeSystem/TaskMng/AsynchTask.h"
#include "RunTimeSystem/TaskMng/StdScheduler.h"
#include "RunTimeSystem/TaskMng/StdDriver.h"
#include "RunTimeSystem/TaskMng/TBBScheduler.h"
#include "RunTimeSystem/TaskMng/TBBDriver.h"
#include "RunTimeSystem/TaskMng/PTHDriver.h"

#include "RunTimeSystem/DataMng/DataArgs.h"

#include "Utils/PerfTools/PerfCounterMng.h"

#include <hwloc.h>
#endif //defined(FEELPP_HAS_HARTS)

namespace Feel
{
namespace vf
{

namespace integrator
{

namespace parallel
{
    template<typename ExprType, typename IMType, typename ElementIterator, typename Eval>
    class HartsContextEvaluate
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        //
        // some typedefs
        //

        typedef Eval eval;
        typedef typename eval::the_element_type EltType;
        typedef ElementIterator element_iterator;

        typedef ExprType expression_type;
        typedef typename eval::gm_type gm_type;
        typedef boost::shared_ptr<gm_type> gm_ptrtype;
        typedef typename eval::gmc_type gmc_type;
        typedef typename eval::gmpc_type gmpc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        typedef boost::shared_ptr<gmpc_type> gmpc_ptrtype;
        typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
        typedef typename expression_type::template tensor<map_gmc_type> eval_expr_type;
        typedef typename eval_expr_type::shape shape;
        typedef IMType im_type;
        typedef typename eval::matrix_type value_type;

#if defined(FEELPP_HAS_HARTS)
        typedef RunTimeSystem::DataHandler                               DataHandlerType ;
        typedef RunTimeSystem::DataArgs<DataHandlerType>                 DataArgsType ;
#endif

        HartsContextEvaluate( int threadId, ExprType const& _expr,
                         IMType const& _im,
                         EltType const& _elt )
            :
            M_threadId( threadId ),
            M_gm( new gm_type( *_elt.gm() ) ),
            M_geopc( new gmpc_type( M_gm, _im.points() ) ),
            M_c( new gmc_type( M_gm, _elt, M_geopc ) ),
            M_expr( _expr, map_gmc_type( fusion::make_pair<vf::detail::gmc<0> >( M_c ) ) ),
            M_im( _im ),
            M_ret( eval::matrix_type::Zero() ),
            M_cpuTime( 0.0 )
        { }

        HartsContextEvaluate( HartsContextEvaluate const& c )
            :
            M_gm( new gm_type( *c.M_gm ) ),
            //M_geopc( new gmpc_type( M_gm, c.M_im.points() ) ),
            M_geopc( c.M_geopc ),
            M_c( new gmc_type( M_gm, c.M_c->element(), M_geopc ) ),
            M_expr( c.M_expr ),
            M_im( c.M_im ),
            M_ret( c.M_ret ),
            M_cpuTime( 0.0 )
        { }

#if defined(FEELPP_HAS_HARTS)
        void computeCPU(DataArgsType& args)
        {
            char * a;
            int cid;
            hwloc_cpuset_t set;
            std::ostringstream oss;

            /* get a cpuset object */
            set = hwloc_bitmap_alloc();

            /* Get the cpu thread affinity info of the current process/thread */
            hwloc_get_cpubind(Environment::getHwlocTopology(), set, 0);
            hwloc_bitmap_asprintf(&a, set);
            oss << a;
            free(a); 
            
            cid = hwloc_bitmap_first(set);
            oss << "(";
            while(cid != -1)
            {
                oss << cid << " ";
                cid = hwloc_bitmap_next(set, cid);
            }
            oss << ")|";
            std::cout << Environment::worldComm().rank() << "|" << M_threadId << " " << oss.str() << std::endl;

            /* Get the latest core location of the current process/thread */
            hwloc_get_last_cpu_location(Environment::getHwlocTopology(), set, 0);
            hwloc_bitmap_asprintf(&a, set);
            oss << a;
            free(a);

            cid = hwloc_bitmap_first(set);
            oss << "(";
            while(cid != -1)
            {
                oss << cid << " ";
                cid = hwloc_bitmap_next(set, cid);
            }
            oss << ");";
            std::cout << Environment::worldComm().rank() << "|" << M_threadId << " " << oss.str() << std::endl;

            /* free memory */
            hwloc_bitmap_free(set);

            perf_mng.init("cpu") ;
            perf_mng.start("cpu") ;

            // DEFINE the range to be iterated on
            std::vector<std::pair<element_iterator, element_iterator> > * r =
                args.get("r")->get<std::vector<std::pair<element_iterator, element_iterator> > >();

            /*
            perf_mng.init("section1") ;
            perf_mng.init("section2") ;
            perf_mng.init("loop") ;
            */

            for (int i = 0; i < r->size(); i++)
            {
                //std::cout << Environment::worldComm().rank() <<  " nbItems: " << r->size() << " nbElts " << std::distance(r->at(i), r->at(i+1)) << std::endl;
                for ( auto _elt = r->at(i).first; _elt != r->at(i).second; ++_elt )
                {
                    //perf_mng.init("section1");
                    //perf_mng.start("section1");
                    M_c->update( *_elt );
                    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( M_c ) );
                    //perf_mng.stop("section1");

                    //perf_mng.init("section2");
                    //perf_mng.start("section2");
                    M_expr.update( mapgmc );
                    M_im.update( *M_c );
                    //perf_mng.stop("section2");

                    //perf_mng.init("loop");
                    //perf_mng.start("loop");
                    for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    {
                        for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        {
                            M_ret( c1,c2 ) += M_im( M_expr, c1, c2 );
                        }
                    }
                    //perf_mng.stop("loop");
                }
            }

            perf_mng.stop("cpu") ;
            M_cpuTime = perf_mng.getValueInSeconds("cpu");
        }
#endif

        void computeCPUOMP(std::vector<std::pair<element_iterator, element_iterator> > * veit)
        {
#if defined(FEELPP_HAS_HARTS)
            perf_mng.init("cpu") ;
            perf_mng.start("cpu") ;
#endif

            for (int i = 0; i < veit->size(); i+=2)
            {
                /*
                std::cout << Environment::worldComm().rank() <<  " nbItems: " << veit->size() 
                          << " nbElts " << std::distance(veit->at(i), veit->at(i+1))
                          << " 1st id " << veit->at(i)->id() << std::endl;
                */

                //std::cout << Environment::worldComm().rank() << "|" << theadId << " fid=" veit.at(i).first.id() << std::endl;
                for ( auto _elt = veit->at(i).first; _elt != veit->at(i).second; ++_elt )
                {
                    M_c->update( *_elt );
                    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( M_c ) );

                    M_expr.update( mapgmc );
                    M_im.update( *M_c );

                    for ( uint16_type c1 = 0; c1 < eval::shape::M; ++c1 )
                    {
                        for ( uint16_type c2 = 0; c2 < eval::shape::N; ++c2 )
                        {
                            M_ret( c1,c2 ) += M_im( M_expr, c1, c2 );
                        }
                    }
                }
            }

#if defined(FEELPP_HAS_HARTS)
            perf_mng.stop("cpu") ;
            M_cpuTime = perf_mng.getValueInSeconds("cpu");
#endif
        }

        value_type result() const
        {
            return M_ret;
        }

        double elapsed()
        {
            return M_cpuTime;
        }

#if defined(FEELPP_HAS_HARTS)
        void printPerfInfo()
        {
            perf_mng.printInfo();
        }
#endif

        int M_threadId;

        gm_ptrtype M_gm;
        gmpc_ptrtype M_geopc;
        gmc_ptrtype M_c;
        eval_expr_type M_expr;
        im_type M_im;
        value_type M_ret;

#if defined(FEELPP_HAS_HARTS)
        RunTimeSystem::PerfCounterMng<std::string> perf_mng ;
#endif

        double M_cpuTime;
    };
} //parallel
    
} // integrator
} // vf
} // Feel


#endif //FEELPP_DETAIL_INTEGRATOR_HPP
