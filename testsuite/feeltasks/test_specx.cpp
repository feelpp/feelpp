#define BOOST_TEST_MODULE test_specx
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

#include <specx/Data/SpDataAccessMode.hpp>
#include <specx/Legacy/SpRuntime.hpp>
#include <specx/Task/SpPriority.hpp>
#include <specx/Utils/SpArrayView.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( testspecx_suite )

BOOST_AUTO_TEST_CASE( test_specx_1 )
{
    SpRuntime<> runtime( 2 );
    runtime.task([]()
                    {
                        std::cout << "Hello World!\n";
                    }).setTaskName("hello");
    runtime.task([]()
                    {
                        std::cout << "Howdi!\n";
                    }).setTaskName("howdi");
    runtime.waitAllTasks();
    runtime.stopAllThreads();
}
BOOST_AUTO_TEST_SUITE_END()
