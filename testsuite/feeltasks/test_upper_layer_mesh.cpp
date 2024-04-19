

#define BOOST_TEST_MODULE test_specx

#include <feel/feelmesh/ranges.hpp>

#include <feel/feelcore/enumerate.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <fmt/chrono.h>

/*
#include <specx/Data/SpDataAccessMode.hpp>
#include <specx/Legacy/SpRuntime.hpp>
#include <specx/Task/SpPriority.hpp>
#include <specx/Utils/SpArrayView.hpp>
*/


#include <specx/Data/SpDataAccessMode.hpp>
#include <specx/Legacy/SpRuntime.hpp>
#include <specx/Task/SpPriority.hpp>
#include <specx/Task/SpProbability.hpp>
#include <specx/Utils/SpArrayView.hpp>
#include <specx/Utils/SpTimer.hpp>
#include <specx/Utils/small_vector.hpp>
#include <specx/Utils/SpBufferDataView.hpp>
#include <specx/Utils/SpBufferDataView.hpp>
#include <specx/Utils/SpHeapBuffer.hpp>
#include <specx/Utils/SpUtils.hpp>
#include <specx/Utils/SpConsumerThread.hpp>
#include <specx/Legacy/SpRuntime.hpp>

#include "Taskflow_HPC.hpp"

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( testspecx_suite )


BOOST_AUTO_TEST_CASE( test_upper_layer_mesh_2 )
{
    using namespace Feel;
    using mesh_t = Mesh<Simplex<3, 1, 3>>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    auto mesh = loadMesh( _mesh = new mesh_t );
    auto therange = elements( mesh );
    using range_t = decay_type<decltype(therange)>;
    double area = 1;

    for(int NumThreads : {1,2,4,8,16})
    {
        LEM::Task TsK(NumThreads,3);
        TsK.setSave(false);
        TsK.setInfo(false);
        TsK.setFileName("./Test_Mesh2");
        TsK.getInformation();
        
        double area = 0;
        auto ranges = parTaskflow_HPCionRange( therange, NumThreads );
        tic();

        Eigen::VectorXd area_vec( NumThreads );
        area_vec.setZero();
        for ( auto [i,r] : enumerate(ranges) )
        {
            int threadid = i;
            auto const& ra = r;
            
            TsK.add(_parameters=Frontend::parameters(area_vec(threadid )),
                _task= [ra,NumThreads,threadid]( double& a )
                {
                    for( element_t<range_t> const& e : ra )
                    {
                        a += e.measure();
                    }
                }
            );
        }
        
        TsK.run();
        BOOST_MESSAGE( fmt::format( "[waitall] areas: {}", area_vec ) );
        area = area_vec.sum();
        double t = toc(fmt::format( "measure-area-{}", NumThreads ) );
        TsK.close();
        TsK.debriefingTasks();
        
        BOOST_MESSAGE( fmt::format( "time for {} threads: {}, area: {}", NumThreads, t, area ) );
        BOOST_CHECK_CLOSE( area, 1, 1e-10 );
    }
}




BOOST_AUTO_TEST_CASE( test_upper_layer_integrate_2 )
{
    using namespace Feel;
    using mesh_t = Mesh<Simplex<3, 1, 3>>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    auto mesh = loadMesh( _mesh = new mesh_t );
    auto therange = elements( mesh );
    using range_t = decay_type<decltype( therange )>;
    double area = 1;

    for ( int NumThreads : { 1, 2, 4, 8, 16 } )
    {
        LEM::Task TsK(NumThreads,3);
        TsK.setSave(false);
        TsK.setInfo(false);
        TsK.setFileName("./Test_Integrate2");
        TsK.getInformation();

        double area = 0;
        auto ranges = parTaskflow_HPCionRange( therange, NumThreads );
        tic();

        Eigen::VectorXd area_vec( NumThreads );
        area_vec.setZero();
        for ( auto [i, r] : enumerate( ranges ) )
        {
            int threadid = i;
            auto const& ra = r;
            
            TsK.add(_parameters=Frontend::parameters(area_vec(threadid )),
                _task=[ra, NumThreads, threadid]( double& a )
                {
                    a = integrate( _range=ra, _expr=expr("1") ).evaluate()(0,0);
                }
            );
        }
        TsK.run();
        BOOST_MESSAGE( fmt::format( "[waitall] areas: {}", area_vec ) );
        area = area_vec.sum();
        double t = toc( fmt::format( "int-area-{}", NumThreads ) );
        TsK.close();
        TsK.debriefingTasks();
        BOOST_MESSAGE( fmt::format( "time for {} threads: {}, area: {}", NumThreads, t, area ) );
        BOOST_CHECK_CLOSE( area, 1, 1e-10 );
    }
}




BOOST_AUTO_TEST_SUITE_END()
