#define BOOST_TEST_MODULE test_specx
#include <feel/feelcore/testsuite.hpp>

#include <specx/Data/SpDataAccessMode.hpp>
#include <specx/Legacy/SpRuntime.hpp>
#include <specx/Task/SpPriority.hpp>
#include <specx/Utils/SpArrayView.hpp>


#include <specx/Utils/SpTimer.hpp>
#include <specx/Utils/small_vector.hpp>
#include <specx/Utils/SpBufferDataView.hpp>
#include <specx/Utils/SpBufferDataView.hpp>
#include <specx/Utils/SpHeapBuffer.hpp>

#include <specx/Legacy/SpRuntime.hpp>



FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( testspecx_suite )

BOOST_AUTO_TEST_CASE( test_specx_1 )
{
    SpRuntime runtime( 2 );
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


BOOST_AUTO_TEST_CASE( test_specx_2 )
{
    const int NumThreads = SpUtils::DefaultNumThreads();
    SpRuntime runtime(NumThreads);
    auto start_time= std::chrono::steady_clock::now();
    const int initVal = 100;
    int writeVal = 0;
    // Create a task with lambda function
    runtime.task(SpRead(initVal), SpWrite(writeVal),
        [](const int& initValParam, int& writeValParam){
        writeValParam += initValParam;
    });

    // Create a task with lambda function (that returns a bool)
    auto returnValue = runtime.task(SpRead(initVal), SpWrite(writeVal),
        [](const int& initValParam, int& writeValParam) -> bool {
        writeValParam += initValParam;
        return true;
    });

    // Wait completion of a single task
    returnValue.wait();
    // Get the value of the task
    const bool res = returnValue.getValue();
    //We are waiting for all tasks to complete.
    runtime.waitAllTasks();
    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    std::cout<<"[INFO SPECX] : CTRL VALUES initval="<<initVal<<" writeVal="<<writeVal<<" res="<<res<<"\n";
    
    //We calculate the time frame with the “SPECX” clock in ms
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	std::cout<<"[INFO SPECX] : Execution Time in ms since start :"<<run_time.count()<<"\n";
}



BOOST_AUTO_TEST_CASE( test_specx_3 )
{
    const int NumThreads = 6;
    SpRuntime runtime(NumThreads);
    auto start_time= std::chrono::steady_clock::now();
    SpHeapBuffer<small_vector<int>> heapBuffer;
    int valueN=0; 
    for(int idx = 0 ; idx < 6 ; ++idx)
    {
        
        auto vectorBuffer = heapBuffer.getNewBuffer();
        runtime.task( SpWrite(vectorBuffer.getDataDep() ) ,
            [&](SpDataBuffer<small_vector<int>> ) mutable 
            {
                valueN=idx;
                usleep(1000); //<== simply to see the task boxes better on the graph.
            }
            ).setTaskName("Write Vector Buffer");

        for(int idxSub = 0 ; idxSub < 2 ; ++idxSub)
        { 
            runtime.task( SpRead( vectorBuffer.getDataDep() ),
            [=] (const SpDataBuffer<small_vector<int>>) 
            {
                 usleep(2000);  //<== simply to see the task boxes better on the graph.
            }
            ).setTaskName("Read Vector Buffer");
        }
        
    }
    //We are waiting for all tasks to complete.
    runtime.waitAllTasks(); 
    //We stop the completion of all tasks.
    runtime.stopAllThreads(); 

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	std::cout<<"[INFO SPECX] : Execution Time in ms since start :"<<run_time.count()<<"\n";

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_3.dot");
    // We generate an Svg trace of the execution in ms
    runtime.generateTrace("test_specx_3.svg");
}


BOOST_AUTO_TEST_CASE( test_specx_4 )
{
    const int NumThreads = SpUtils::DefaultNumThreads();
    SpRuntime runtime(NumThreads);
    auto start_time= std::chrono::steady_clock::now();
    const int baseTime = 10000;
    for(int idx = 0 ; idx < 5 ; ++idx)
    {
        SpBufferDataView<small_vector<int>> vectorBuffer;
        runtime.task(SpWrite(vectorBuffer.getDataDep()),         
            []([[maybe_unused]] SpDataBuffer<small_vector<int>> vector)
            {
                usleep(1000); //<== simply to see the task boxes better on the graph.            
            }
        );

        usleep(baseTime);

        for(int idxSub = 0 ; idxSub < 3 ; ++idxSub)
        {
            runtime.task(SpRead(vectorBuffer.getDataDep()),       
                []([[maybe_unused]] const SpDataBuffer<small_vector<int>> vector)
                {
                    usleep(2000); //<== simply to see the task boxes better on the graph.
                }
            );
        }
    }

    //We are waiting for all tasks to complete
    runtime.waitAllTasks();
    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	std::cout<<"[INFO SPECX] : Execution Time in ms since start :"<<run_time.count()<<"\n";

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_4.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_4.dot");

}

BOOST_AUTO_TEST_CASE( test_specx_5 )
{

    SpRuntime<SpSpeculativeModel::SP_MODEL_2> runtime;
    auto start_time= std::chrono::steady_clock::now();

    runtime.setSpeculationTest([](const int,const SpProbability&) -> bool{
        return true; // <== Always speculate
    });

    int val = 0;
    std::promise<int> promise1;

    runtime.task(SpRead(val), [&promise1](const int&){
        promise1.get_future().get();
    }).setTaskName("First-task");
        
    //standard task
    for(int idx = 0; idx < 2; idx++) {
        runtime.task(SpWrite(val), [](int&)  
        {
        }).setTaskName("Certain task -- " + std::to_string(idx));
    }
        
    const int nbUncertainTasks = 3;

    for(int idx = 0 ; idx < nbUncertainTasks ; ++idx)
    {
        runtime.task(SpPotentialWrite(val), [](int&) -> bool 
        {
            return true;
        }).setTaskName("Uncertain task -- " + std::to_string(idx));
     }

    for(int idx = 2; idx < 4; idx++) 
    {
        runtime.task(SpWrite(val), [](int&)  
        {
            //...
        }).setTaskName("Certain task -- " + std::to_string(idx));
    }
        
    runtime.task(SpWrite(val), []([[maybe_unused]] int& valParam)
    {
        //...
    }).setTaskName("Last-task");

        
    promise1.set_value(0);

    runtime.waitAllTasks();
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	std::cout<<"[INFO SPECX] : Execution Time in ms since start :"<<run_time.count()<<"\n";

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_4.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_4.dot");

       
}


BOOST_AUTO_TEST_SUITE_END()



