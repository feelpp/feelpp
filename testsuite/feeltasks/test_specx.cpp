#define BOOST_TEST_MODULE test_specx
#include <feel/feelcore/testsuite.hpp>


//#include <feel/feelmesh/ranges.hpp>

///nvme0/lemoinep/feelpp/feelpp/feel/feelcore/range.hpp

//#include <feel/feelcore/enumerate.hpp>
//#include <feel/feelcore/environment.hpp>
//#include <feel/feeldiscr/mesh.hpp>
//#include <feel/feelfilters/loadmesh.hpp>



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


//#include "UTester.hpp"
//#include "utestUtils.hpp"

#define DEVICE_FUNC __device__ __host__
//#define USE_GPU_HIP



#include "Taskflow_HPC.hpp"

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( testspecx_suite )

BOOST_AUTO_TEST_CASE( test_specx_1 )
{
    SpRuntime runtime( 2 );
    runtime.task([]()
                    {
                       BOOST_MESSAGE("Hello World!\n");
                    }).setTaskName("hello");
    runtime.task([]()
                    {
                       BOOST_MESSAGE("Howdi!\n");
                    }).setTaskName("howdi");
    runtime.waitAllTasks();
    runtime.stopAllThreads();
}


BOOST_AUTO_TEST_CASE( test_specx_2 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T2>\n");
    //const int NumThreads = SpUtils::DefaultNumThreads();
    const int NumThreads = 6;
    SpRuntime runtime(NumThreads);
    auto start_time= std::chrono::steady_clock::now();
    const int initVal = 100;
    int writeVal = 0;
    // Create a task with lambda function
    runtime.task(SpRead(initVal), SpWrite(writeVal),
        [](const int& initValParam, int& writeValParam){
        writeValParam += initValParam;
        usleep(10000); //<== simply to see the task boxes better on the graph.
    });

    // Create a task with lambda function (that returns a bool)
    auto returnValue = runtime.task(SpRead(initVal), SpWrite(writeVal),
        [](const int& initValParam, int& writeValParam) -> bool {
        writeValParam += initValParam;
        usleep(20000); //<== simply to see the task boxes better on the graph.
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

    BOOST_MESSAGE("[INFO SPECX] : CTRL VALUES initval="<<initVal<<" writeVal="<<writeVal<<" res="<<res<<"\n");
    
    //We calculate the time frame with the “SPECX” clock in ms
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T2> in ms since start :"<<run_time.count()<<"\n");

        // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_2.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_2.svg");   
}



BOOST_AUTO_TEST_CASE( test_specx_3 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T3>\n");
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
    BOOST_MESSAGE("[INFO SPECX] : Execution Time <T3> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_3.dot");
    // We generate an Svg trace of the execution in ms
    runtime.generateTrace("test_specx_3.svg");
}


BOOST_AUTO_TEST_CASE( test_specx_4 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T4>\n");
    //const int NumThreads = SpUtils::DefaultNumThreads();
    const int NumThreads = 6;
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
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T4> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_4.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_4.svg");

}

BOOST_AUTO_TEST_CASE( test_specx_5 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T5>\n");
    // We instantiate a runtime object and we specify that the
    // runtime should use speculation model 2.
    SpRuntime<SpSpeculativeModel::SP_MODEL_2> runtime;
    auto start_time= std::chrono::steady_clock::now();

    // Next we set a predicate that will be called by the runtime each
    // time a speculative task becomes ready to run. It is used to
    // decide if the speculative task should be allowed to run.
    runtime.setSpeculationTest([](const int,const SpProbability&) -> bool{
        return true; // <== Always speculate
    });

    int val = 0;
    std::promise<int> promise1;

    runtime.task(SpRead(val), [&promise1](const int&) {
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
            usleep(1000);  //<== simply to see the task boxes better on the graph.
            return true;
        }).setTaskName("Uncertain task -- " + std::to_string(idx));
     }

    for(int idx = 2; idx < 4; idx++) 
    {
        runtime.task(SpWrite(val), [](int&)  
        {
            usleep(2000);  //<== simply to see the task boxes better on the graph.
            //...
        }).setTaskName("Certain task -- " + std::to_string(idx));
    }
        
    runtime.task(SpWrite(val), []([[maybe_unused]] int& valParam)
    {
        //...
    }).setTaskName("Last-task");

        
    promise1.set_value(0);

    //We are waiting for all tasks to complete
    runtime.waitAllTasks();
    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T5> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_5.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_5.svg");
 
}


BOOST_AUTO_TEST_CASE( test_specx_6 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T6>\n");
    //const int numThreads = SpUtils::DefaultNumThreads();
    int numThreads = std::min(10,SpUtils::DefaultNumThreads());
    SpRuntime runtime(numThreads);
    auto start_time= std::chrono::steady_clock::now();
    const int nbTasksToSubmit = runtime.getNbThreads();
    static const int nbLoops = 10;
    int initVal = 1;

    small_vector<SpAbstractTaskWithReturn<double>::SpTaskViewer> elapsed;
    elapsed.reserve(nbTasksToSubmit*nbLoops);


    for(int idxLoop = 0 ; idxLoop < nbLoops ; ++idxLoop){
            for(int idx = 0 ; idx < nbTasksToSubmit ; ++idx){
                SpTimer timerTask;
                elapsed.emplace_back( 
                    runtime.task(SpRead(initVal),
                        [timerTask](const int&) mutable -> double {
                        timerTask.stop();
                        usleep(1000);  //<== simply to see the task boxes better on the graph.
                        return timerTask.getElapsed();
                    }
                ) 
                );
            }
            runtime.waitAllTasks();
        }

    double averageToExecute = 0;
    for(const SpAbstractTaskWithReturn<double>::SpTaskViewer& viewer : elapsed){
            averageToExecute += viewer.getValue()/double(nbTasksToSubmit*nbLoops);
    }

    BOOST_MESSAGE("[INFO SPECX] : Average time <T6> for a task to be executed without pressure = " << averageToExecute << "s\n");

    //We are waiting for all tasks to complete
    runtime.waitAllTasks();
    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T6> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_6.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_6.svg");
}



BOOST_AUTO_TEST_CASE( test_specx_7 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T7>\n");
    [[maybe_unused]] const size_t seedSpeculationSuccess = 42;
    [[maybe_unused]] const size_t seedSpeculationFailure = 0;
    const size_t seed = seedSpeculationSuccess;
    // We instantiate a runtime object and we specify that the
    // runtime should use speculation model 2.
    SpRuntime<SpSpeculativeModel::SP_MODEL_2> runtime;
    auto start_time= std::chrono::steady_clock::now();
    
    // Next we set a predicate that will be called by the runtime each
    // time a speculative task becomes ready to run. It is used to
    // decide if the speculative task should be allowed to run.
    runtime.setSpeculationTest(
    []([[maybe_unused]] const int nbReadyTasks,
    [[maybe_unused]] const SpProbability& meanProbability) -> bool {
        usleep(1000);  //<== simply to see the task boxes better on the graph.
        return true; // Here we always return true, this basically means
                     // that we always allow speculative tasks to run
                     // regardless of runtime conditions.
    });

    int a = 41, b = 0, c = 0;
    auto task1 = runtime.task(SpRead(a), [](const int& inA) -> int {
        usleep(2000);  //<== simply to see the task boxes better on the graph.
        return inA + 1;
    });
    task1.setTaskName("First-task");
    b = task1.getValue();
    
    // Next we create a potential task, i.e. a task which might write to some data.
    // In this case the task may write to "a" with a probability of 0.5.
    // Subsequent tasks will be allowed to speculate over this task.
    // The task returns a boolean to inform the runtime of whether or 
    // not it has written to its maybe-write data dependency a.
    std::mt19937_64 mtEngine(seed);
    std::uniform_real_distribution<double> dis01(0,1);
    
    runtime.task(SpPriority(0), SpProbability(0.5), SpRead(b),
    SpPotentialWrite(a),
    [dis01, mtEngine] (const int &inB, int &inA) mutable -> bool {
        double val = dis01(mtEngine);
        
        if(inB == 42  && val < 0.5) {
            inA = 43;
            return true;
        }
        usleep(3000);  //<== simply to see the task boxes better on the graph.
        return false;
        
    }).setTaskName("Second-task");
    
    // We create a final normal task that reads from a and writes to c.
    // The task reads from a so there should be a strict write -> read
    // dependency between the second and the final task but since the
    // second task may not always write to a, the runtime will try to
    // execute a speculative version of the final task in parallel
    // with the second task in case the second task doesn't write to a.
    runtime.task(SpRead(a), SpWrite(c), [] (const int &inA, int &inC) {
        if(inA == 41) {
            inC = 1;
        } else {
            inC = 2;
        }
        usleep(4000);  //<== simply to see the task boxes better on the graph.
    }).setTaskName("Final-task");

    // We wait for all tasks to finish
    runtime.waitAllTasks();
    
    // We make all runtime threads exit
    runtime.stopAllThreads();
    
    assert((a == 41 || a == 43) && b == 42 && (c == 1 || c == 2)
            && "Try again!");
    
    //We are waiting for all tasks to complete
    runtime.waitAllTasks();
    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T7> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_7.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_7.svg");
}


BOOST_AUTO_TEST_CASE( test_specx_8 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T8>\n");

    /*
    //Small Thread Race
    std::array<unsigned int,3> SleepTimes{0, 500,1000};
    int const NumThreads = 10;
    SpTimer T0;

    auto start_time= std::chrono::steady_clock::now();

    for(auto SleepTime : SleepTimes){   
        SpRuntime runtime(NumThreads);      
        runtime.setSpeculationTest(
            [](const int, const SpProbability&) -> bool
            {
                return true;
            });

        const int arraySize = 6;
        int val[arraySize] = {0};

        UTestRaceChecker counterAccess;
        std::string ChInfo,ChInfoPlus;
        ChInfo = "Runtime Array Process";

        runtime.task(SpReadArray(val,SpArrayView(arraySize)), 
        [](SpArrayAccessor<const int>&){}).setTaskName("Small Race");
 
            for(int idx = 0 ; idx < arraySize ; ++idx){
                ChInfoPlus = ChInfo+std::to_string(idx);

                runtime.task(SpWrite(val[idx]),
                                      SpReadArray(val,SpArrayView(arraySize).removeItem(idx)),
                                      [SleepTime,idx,&counterAccess]
                                      (int& valParam, const SpArrayAccessor<const int>& valArray) -> bool {
                    {
                        counterAccess.lock();
                        counterAccess.addWrite(&valParam);
                        for(int idxTest = 0 ; idxTest < valArray.getSize() ; ++idxTest){
                            counterAccess.addRead(&valArray.getAt(idxTest));
                        }
                        counterAccess.unlock();
                    }

                    std::string ChNumidx=std::to_string(idx);
                    if (1==0) {
                        std::cout<<"   [INFO SPECX] : Values of TimeSleep="<<SleepTime<<" Idx="<<idx<<" valParam1="<<valParam;
                    }
                    if (1==1) {
                         if(idx == 1){
                            valParam += 2;
                        }
                        if(idx == 3){
                            valParam += 1;
                        }
                        if(idx == 5){
                            valParam += 10;
                        }
                    }
                    if (1==0) {
                        std::cout<<" valParam2="<<valParam<<"\n";
                    }

                    usleep(SleepTime);

                    {
                        counterAccess.lock();
                        counterAccess.releaseWrite(&valParam);
                        for(int idxTest = 0 ; idxTest < valArray.getSize() ; ++idxTest){
                            counterAccess.releaseRead(&valArray.getAt(idxTest));
                        }
                        counterAccess.unlock();
                    }

                    return (idx == 3 || idx == 5);
                }).setTaskName(ChInfoPlus); //END runtime;
            }
        //We are waiting for all tasks to complete
        runtime.waitAllTasks();
        //We stop the completion of all tasks.
        runtime.stopAllThreads();

        // We generate the task graph corresponding to the execution 
        runtime.generateDot("test_specx_8_SleepTime"+std::to_string(SleepTime)+".dot", true);
    
        // We generate an Svg trace of the execution
        runtime.generateTrace("test_specx_8_SleepTime"+std::to_string(SleepTime)+".svg");   
    }
    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);

	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T8> in ms since start :"<<run_time.count()<<"\n"); 

    */
}   



BOOST_AUTO_TEST_CASE( test_specx_9 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T9>\n");
    SpRuntime runtime(2);
    auto start_time= std::chrono::steady_clock::now();

        //std::cout << "Task Block1\n";   
        {
            std::promise<int> promise1;
            std::promise<int> promise2;

            const int initVal = 1;

            runtime.task(SpRead(initVal),
                         [&](const int& /*initValParam*/){
                promise1.set_value(0);
                promise2.get_future().wait();
            });
            
            promise1.get_future().wait();
            
            runtime.task(SpRead(initVal),
                         [&](const int& /*initValParam*/){
                promise2.set_value(0);
            });        
            
            runtime.waitAllTasks();
        }

        //std::cout << "Task Block2\n";   
        {
            std::promise<int> promise1;
            std::promise<int> promise2;

            const int initVal = 1;
            int writeVal = 0;

            runtime.task(SpRead(initVal), SpWrite(writeVal),
                         [&](const int& /*initValParam*/, int& /*writeValParam*/){
                promise1.set_value(0);
                promise2.get_future().wait();
            });
            
            promise1.get_future().wait();
            
            runtime.task(SpRead(initVal),
                         [&](const int& /*initValParam*/){
                promise2.set_value(0);
            });        
            
            runtime.waitAllTasks();
        }

        //std::cout << "Task Block3\n";       
        {
            std::promise<int> promise1;
            std::promise<int> promise2;

            const int initVal[10] = {0};

            runtime.task(SpReadArray(initVal, SpArrayView(10)),
                         [&](const SpArrayAccessor<const int>& /*initValParam*/){
                promise1.set_value(0);
                promise2.get_future().wait();
            });
            
            promise1.get_future().wait();
            
            runtime.task(SpReadArray(initVal, SpArrayView(10)),
                         [&](const SpArrayAccessor<const int>& /*initValParam*/){
                promise2.set_value(0);
            });
            
            runtime.waitAllTasks();
        }

        //std::cout << "Task Block4\n";   
        {
            std::promise<int> promise1;
            std::promise<int> promise2;

            const int initVal[10] = {0};

            runtime.task(SpReadArray(initVal, SpArrayView(10).removeItems(5,9)),
                         [&](const SpArrayAccessor<const int>& /*initValParam*/){
                promise1.set_value(0);
                promise2.get_future().wait();
            });
            
            promise1.get_future().wait();
            
            runtime.task(SpRead(initVal[0]),
                         [&](const int& /*initValParam*/){
                promise2.set_value(0);
            });        
            
            runtime.waitAllTasks();
        }

        //std::cout << "Task Block5\n";   
        {
            std::promise<int> promise1;
            std::promise<int> promise2;

            int initVal[10] = {0};

            runtime.task(SpReadArray(initVal, SpArrayView(10).removeItems(1)),
                         [&](const SpArrayAccessor<const int>& /*initValParam*/){
                promise1.set_value(0);
                promise2.get_future().wait();
            });
            
            promise1.get_future().wait();
            
            runtime.task(SpWrite(initVal[1]),
                         [&](int& /*initValParam*/){
                promise2.set_value(0);
            });        
            
            runtime.waitAllTasks();
        }

        //std::cout << "Task Block6\n";   
        {
            std::promise<int> promise1;
            std::promise<int> promise2;

            int initVal[10] = {0};

            runtime.task(SpReadArray(initVal, SpArrayView(10).removeItems(0)),
                         [&](const SpArrayAccessor<const int>& /*initValParam*/){
                promise1.set_value(0);
                promise2.get_future().wait();
            });

            promise1.get_future().wait();

            runtime.task(SpWrite(initVal[0]),
                         [&](int& /*initValParam*/){
                promise2.set_value(0);
            });

            //We are waiting for all tasks to complete
            runtime.waitAllTasks();
        }

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T9> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_9.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_9.svg");   
}



BOOST_AUTO_TEST_CASE( test_specx_10 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T10>\n");
    auto start_time= std::chrono::steady_clock::now();
    int  nbObjects = 12;
    

    int numThreads = std::min(5,std::min(nbObjects,SpUtils::DefaultNumThreads()));
    SpRuntime runtime(numThreads);
    
    int nbTasksToSubmit = runtime.getNbThreads();
    int nbLoops=std::max(1,nbObjects/nbTasksToSubmit);
    int nbCoor=round((float(nbObjects)/float(nbTasksToSubmit)-float(nbObjects/nbTasksToSubmit))*float(nbTasksToSubmit));
    int nbidx=nbTasksToSubmit;
    small_vector<SpAbstractTaskWithReturn<double>::SpTaskViewer> elapsed;
    elapsed.reserve(nbTasksToSubmit*nbLoops);
    
    bool qWaitTask=true;
    int  initVal = 1;
    int  it=0;
    for(int idxLoop = 0 ; idxLoop < nbLoops ; ++idxLoop)
    {
        if (idxLoop==nbLoops-1) 
        { 
            nbidx=nbidx+nbCoor;
        }

        for(int idx = 0 ; idx < nbidx ; ++idx){
            SpTimer timerTask;
            elapsed.emplace_back( 
                runtime.task(SpRead(initVal),
                    [timerTask](const int&) mutable -> double {
                        timerTask.stop();
                        usleep(10000); //<== simply to see the task boxes better on the graph.
                        return timerTask.getElapsed();
                    }
                )   
            );
            it++;
        }
        if (qWaitTask) { runtime.waitAllTasks(); }
    }
    //std::cout<<"CTRL Nb="<<it<<"\n";

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T10> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_10.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_10.svg");   
}


BOOST_AUTO_TEST_CASE( test_specx_11 )
{
    //PARALLEL_WRITE WITH ATOMIC
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T11>\n");
    auto start_time= std::chrono::steady_clock::now();
    int NumThreads = std::min(6,SpUtils::DefaultNumThreads());
    SpRuntime runtime(NumThreads);
    std::atomic<int> initVal(0);
    int nbTh=4;
    int val0=5;
    for(int idxTh = 0 ; idxTh < nbTh; ++idxTh){
        runtime.task(SpParallelWrite(initVal),
            [&](std::atomic<int>& initValParam){
                initValParam += 1;
                //val0++;
                //std::cout<<idxTh<<" "<<initValParam<<"...."<<val0<<"\n";
                while(initValParam != nbTh){
                    usleep(1000); //<== simply to see the task boxes better on the graph.
                } 
               usleep(2000);  //<== simply to see the task boxes better on the graph.
            }).setTaskName("OpPW("+std::to_string(idxTh)+")");

        runtime.task(SpParallelWrite(val0),
            [&](int& valParam0)
            {
                valParam0++;
                usleep(valParam0*1000); //<== simply to see the task boxes better on the graph.
                //...
                //std::cout<<"- val0="<<valParam0<<"\n";
            }).setTaskName("OpPW("+std::to_string(idxTh)+")");
    }

    //We are waiting for all tasks to complete
    runtime.waitAllTasks();

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T11> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_11.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_11.svg");   
}


BOOST_AUTO_TEST_CASE( test_specx_12 )
{
    //PARALLEL_WRITE WITH PROMISE
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T12>\n");
    auto start_time= std::chrono::steady_clock::now();
    int NumThreads = std::min(6,SpUtils::DefaultNumThreads());
    SpRuntime runtime(NumThreads);
    int nbTh=3;
    int val0=30;
    std::promise<long int> promises[5];
    int dumbVal = 0;

    for(int idxTh = 0 ; idxTh < nbTh ; ++idxTh){
        runtime.task(SpParallelWrite(dumbVal),
            [&,idxTh](int&){
            promises[idxTh].set_value(idxTh);
            const long int res = promises[(idxTh+1)%nbTh].get_future().get();
        }).setTaskName("OpPW("+std::to_string(idxTh)+")");

        runtime.task(SpRead(val0),
            [&](const int& valParam0){
                usleep(40000); //<== simply to see the task boxes better on the graph.
                //...
                //std::cout<<"- val0="<<valParam0<<"\n";
        }).setTaskName("OpR("+std::to_string(idxTh)+")");
    }
  
    //We are waiting for all tasks to complete
    runtime.waitAllTasks();

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T12> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_12.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_12.svg");   

    
}

BOOST_AUTO_TEST_CASE( test_specx_13 )
{
    //MERGE DEMO
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T12>\n");
    auto start_time= std::chrono::steady_clock::now();
    
    int NumThreads = std::min(6,SpUtils::DefaultNumThreads());
    SpRuntime runtime(NumThreads);

    runtime.setSpeculationTest([](const int, const SpProbability&) -> bool {
        return true;
    });
            
    int a=0, b=0, c=0, d=0;

    std::promise<bool> promise1;

    runtime.task(SpWrite(a), [&promise1](int& param_a){
        param_a = 1;
        usleep(100000);
        promise1.get_future().get();
    });

    runtime.task(SpRead(a), SpPotentialWrite(b), 
        []([[maybe_unused]] const int& a_param, [[maybe_unused]] int&) -> bool {
            usleep(200000);
            return false;
    });

    runtime.task(SpRead(b), SpPotentialWrite(c), 
        [](const int& param_b, [[maybe_unused]] int&) -> bool {
        bool res = false;
        if(param_b != 0) {
            res = true;
        }
        usleep(300000);
        return res;
    });
            
    runtime.task(SpRead(c), SpWrite(d), [](const int& param_c, int&param_d){
        if (param_c == 0) { param_d = 1; }
        usleep(400000);
    });
    
    promise1.set_value(true);
    
    //We are waiting for all tasks to complete
    runtime.waitAllTasks();

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T13> in ms since start :"<<run_time.count()<<"\n");

    // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_13.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_13.svg");   
}


BOOST_AUTO_TEST_CASE( test_specx_14 )
{
    
    //Calculation of PI value estimate
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T14>\n");
    auto start_time= std::chrono::steady_clock::now();

    int nbThreads = std::min(6,SpUtils::DefaultNumThreads());
    SpRuntime runtime(nbThreads);

    long int nbN=1000000;
    int sizeBlock=nbN/nbThreads;
    int diffBlock=nbN-sizeBlock*nbThreads;
    double h=1.0/double(nbN);
    double integralValue=0.0;
    std::vector<double> valuesVec(nbThreads,0.0);
 
    for(int k1 = 0 ; k1 < nbThreads ; ++k1){
        int vkBegin=k1*sizeBlock;
        int vkEnd=(k1+1)*sizeBlock;
        if ((k1==nbThreads-1) && (k1>0) && (diffBlock>0)) { vkEnd=vkBegin+diffBlock; }
        int threadid = k1;
        runtime.task(
            SpWrite(valuesVec.at(threadid)),
                [h,vkBegin,vkEnd](double& s) -> bool {
                    double sum=0.0; double x;
                    for(int j=vkBegin;j<vkEnd;j++)
                    {
                        x=h*double(j);
                        sum+=4.0/(1.0+x*x);
                    }
                    s=sum;
                    return true;
                }
            ).setTaskName("Op("+std::to_string(k1)+")");
    }

    //We are waiting for all tasks to complete
    runtime.waitAllTasks();

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_14.dot", true);
    
    //We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_14.svg");   

    //Sum of vector elements
    integralValue=h*std::reduce(valuesVec.begin(),valuesVec.end());
    double DeltaError=std::abs(M_PI-integralValue);   
    //std::cout<<"PI Value= "<<integralValue<<"\n";
    BOOST_MESSAGE("[INFO SPECX] : PI Value= "<<integralValue<<"Err="<<DeltaError<<"\n");

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T14> in ms since start :"<<run_time.count()<<"\n");
}

BOOST_AUTO_TEST_CASE( test_specx_15 )
{
    //calculation of PI value estimate with range
    
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T15>\n");
    auto start_time= std::chrono::steady_clock::now();

    int nbThreads = std::min(6,SpUtils::DefaultNumThreads());
    SpRuntime runtime(nbThreads);
    
    long int nbN=1000000;
    int sizeBlock=nbN/nbThreads;
    int diffBlock=nbN-sizeBlock*nbThreads;
    double h=1.0/double(nbN);
    double integralValue=0.0;
    
    //We build a list of numbers. This can be a list of floats.
    std::vector<int> v(boost::counting_iterator<int>(0), boost::counting_iterator<int>(nbN));
    //auto ranges = parTaskflow_HPCionRange(v,nbThreads);

    //We construct the ranges by cutting the list of numbers
    int lenClusters = v.size()/nbThreads;
    int size = (v.size() - 1) / lenClusters + 1;
    std::vector<int> ranges[size];
    for (int k = 0; k < size; ++k)
    {
        auto start_itr = std::next(v.cbegin(), k*lenClusters);
        auto end_itr   = std::next(v.cbegin(), k*lenClusters + lenClusters);
        ranges[k].resize(lenClusters);
        if (k*lenClusters + lenClusters > v.size())
        {
            end_itr = v.cend(); ranges[k].resize(v.size() - k*lenClusters);
        }
        std::copy(start_itr, end_itr,ranges[k].begin());
    }

    //We calculate the integral
    std::vector<double> valuesVec(size,0.0);
    for(int k1 = 0 ; k1 < size ; ++k1){
        int threadid = k1;
        auto const& ra = ranges[k1];
        runtime.task(
            SpWrite(valuesVec.at(threadid)),
                [h,ra](double& s) -> bool {
                    double sum=0.0; double x;
                    for (int i=0;i<ra.size();i++) 
                    { 
                        x=h*double(ra.at(i));
                        sum+=4.0/(1.0+x*x);
                    }
                    s=sum;
                    return true;
                }
            ).setTaskName("Op("+std::to_string(k1)+")");
    }

    //We are waiting for all tasks to complete
    runtime.waitAllTasks();

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_15.dot", true);
    
    //We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_15.svg");   

    //Sum of vector elements
    integralValue=h*std::reduce(valuesVec.begin(),valuesVec.end());
    double DeltaError=std::abs(M_PI-integralValue);   
    //std::cout<<"PI Value= "<<integralValue<<"\n";
    BOOST_MESSAGE("[INFO SPECX] : PI Value= "<<integralValue<<"Err="<<DeltaError<<"\n");

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T15> in ms since start :"<<run_time.count()<<"\n");

}


BOOST_AUTO_TEST_CASE( test_specx_16 )
{
    //An example to see the difference, whether qWaitTask is activated or not.
    //look at the graph "test_specx_16.svg" to understand the difference...
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T16>\n");
    auto start_time= std::chrono::steady_clock::now();

    int nbThreads = std::min(6,SpUtils::DefaultNumThreads());
    SpRuntime runtime(nbThreads);
    int  nbObjects = 12;
    
    int nbTasksToSubmit = runtime.getNbThreads();
    int nbLoops=std::max(1,nbObjects/nbTasksToSubmit);
    int nbCoor=round((float(nbObjects)/float(nbTasksToSubmit)-float(nbObjects/nbTasksToSubmit))*float(nbTasksToSubmit));
    int nbidx=nbTasksToSubmit;
    

    bool qWaitTask=true;

    for(int qi = 0 ; qi < 2 ; ++qi)
    {
        qWaitTask=qi;
        int  initVal = 1;
        int  it=0;
        for(int idxLoop = 0 ; idxLoop < nbLoops ; ++idxLoop)
        {
            if (idxLoop==nbLoops-1) 
            { 
                nbidx=nbidx+nbCoor;
            }

            for(int idx = 0 ; idx < nbidx ; ++idx){
                    auto returnValue=runtime.task( SpRead(initVal),
                        [](const int& value)-> double {
                            usleep(10000); //<== simply to see the task boxes better on the graph.
                            return true;
                        }
                    );
                it++;
                if (qWaitTask) { returnValue.wait(); }
            }
        }

        usleep(100000);
    }

    //We are waiting for all tasks to complete
    runtime.waitAllTasks();

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_16.dot", true);
    
    //We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_16.svg");   

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T16> in ms since start :"<<run_time.count()<<"\n");

}


BOOST_AUTO_TEST_CASE( test_specx_17 )
{
    //Discovering another function.
    //An example using the consum function
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T17>\n");
    auto start_time= std::chrono::steady_clock::now();
    SpConsumerThread cs;
    int Sum=0;
    cs.submitJob([&](){
        Sum++; 
        BOOST_MESSAGE(" 1 Sum="<<Sum<<"\n");
    });

    cs.submitJob([&](){ 
        Sum++; 
        BOOST_MESSAGE(" 2 Sum2="<<Sum<<"\n");
    });

    cs.submitJobAndWait([&](){
        Sum++; 
        BOOST_MESSAGE(" 3 Sum2="<<Sum<<"\n");
    });

    cs.submitJob([&](){
        Sum++; 
        BOOST_MESSAGE(" 4 Su2="<<Sum<<"\n");
    });

    cs.submitJobAndWait([&](){
        Sum++;  
        BOOST_MESSAGE(" 5 Sum="<<Sum<<"\n");
    });
    cs.stop();
    
    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T16> in ms since start :"<<run_time.count()<<"\n");

}


BOOST_AUTO_TEST_CASE( test_specx_18 )
{
    //Speculative write
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T18>\n");
    auto start_time= std::chrono::steady_clock::now();

    int nbThreads = std::min(6,SpUtils::DefaultNumThreads());
    SpRuntime runtime(nbThreads);


    runtime.setSpeculationTest([](const int /*inNbReadyTasks*/, const SpProbability& /*inProbability*/) -> bool {
        return true;
    });

    int val = 0;
    std::promise<int> promise1;

    runtime.task(SpRead(val), [](const int& /*valParam*/){ });
    // val is 0

    runtime.task(SpRead(val), [&promise1](const int& /*valParam*/){
        promise1.get_future().get();
    });
    // val is 0

    runtime.task(SpPotentialWrite(val), [](int& /*valParam*/) -> bool {
        BOOST_MESSAGE("Maybe task will return false\n");
        //std::cout.flush();
        return false;
    });
    // val is 0

    std::atomic_int counterFirstSpec(0);

    runtime.task(SpWrite(val), [&val,&counterFirstSpec](int& valParam){
        BOOST_MESSAGE("Speculative task, valParam is "<< valParam<<" at "<< &valParam <<"\n");
        BOOST_MESSAGE("Speculative task, val is "<<val<<" at "<<&val<<"\n");
        valParam += 1;
        BOOST_MESSAGE("Speculative task, valParam is "<< valParam<<" at "<<&valParam <<"\n");
        BOOST_MESSAGE("Speculative task, val is "<<val<<" at "<<&val<<"\n");
        //std::cout.flush();
        counterFirstSpec += 1;
    });
    // val is 1


    runtime.task(SpPotentialWrite(val), [](int& valParam) -> bool {
        // valParam should be 1
        BOOST_MESSAGE("Maybe task 2, valParam is "<< valParam<<" at "<< &valParam<<"\n");
        //std::cout.flush();
        valParam += 2;
        BOOST_MESSAGE("Maybe task 2, return true with valParam is "<<valParam<<" at "<<&valParam <<"\n");
        return true;
    });
    // val is 3

    std::atomic_int counterSecondSpec(0);

    runtime.task(SpWrite(val), [&val,&counterSecondSpec](int& valParam){
        BOOST_MESSAGE("Speculative last write, valParam is "<< valParam<< " at "<<&valParam<<"\n");
        BOOST_MESSAGE("Speculative last write, val is "<<val<< " at "<<&val<<"\n");
        valParam *= 2;
        BOOST_MESSAGE("Speculative last write, valParam is "<< valParam << " at "<< &valParam<<"\n");
        BOOST_MESSAGE("Speculative last write, val is "<< val<<" at "<< &val<<"\n");
        //std::cout.flush();
        counterSecondSpec += 1;
    });
    // val is 6

    promise1.set_value(0);


    //We are waiting for all tasks to complete
    runtime.waitAllTasks();

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_18.dot", true);
    
    //We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_18.svg");   

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T18> in ms since start :"<<run_time.count()<<"\n");
}


BOOST_AUTO_TEST_CASE( test_specx_19 )
{
    //Demo to see communications between tasks and see the dependencies
    
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T19>\n");
    auto start_time= std::chrono::steady_clock::now();

    int nbThreads = std::min(6,SpUtils::DefaultNumThreads());
    SpRuntime runtime(nbThreads);

    double valueA=1.0;
    double valueS=0;
    double valueS1=0;
    double valueS2=0;
    int size=6;
    std::vector<double> valuesVec(size,0.0);
    int deltaTimeSleepBetweenTasks=200000; //break time between tasks in order to better see the graphic representation.

    runtime.task(
            SpWrite(valueS),
                [valueA](double& s) -> bool {
                    s=valueA;
                    usleep(100000);
                    return true;
                }
            ).setTaskName("Op(1)");

    usleep(deltaTimeSleepBetweenTasks); //break time 
 

    for(int k1 = 0 ; k1 < size ; ++k1){
        int threadid = k1;
        runtime.task(SpRead(valueS),SpWrite(valuesVec.at(threadid)), 
            [size](const double& a_param,double& writeValParam) -> bool {
                writeValParam = (1.0/size)*a_param;
                usleep(writeValParam*100000);
                return true;
        }).setTaskName("OpExp("+std::to_string(k1)+")");
    }

    runtime.waitAllTasks();
    usleep(deltaTimeSleepBetweenTasks); //break time 


    runtime.task(SpRead(valuesVec),SpWrite(valueS1), 
        [](const std::vector<double> v,double& writeValParam) -> bool {
                writeValParam = std::reduce(v.begin(),v.begin()+v.size()/2);
                usleep(writeValParam*100000);
                return true;
    }).setTaskName("Op(1)");

     runtime.task(SpRead(valuesVec),SpWrite(valueS2), 
        [](const std::vector<double> v,double& writeValParam) -> bool {
                writeValParam = std::reduce(v.begin()+v.size()/2,v.end());
                usleep(writeValParam*100000);
                return true;
    }).setTaskName("Op(2)");

    runtime.waitAllTasks();
    usleep(deltaTimeSleepBetweenTasks); //break time 

    runtime.task(SpRead(valueS1),SpRead(valueS2),SpWrite(valueS), 
            [](const double& a_param1,const double& a_param2,double& writeValParam) -> bool {
                writeValParam = a_param1+a_param2;
                usleep(writeValParam*100000);
                return true;
    }).setTaskName("OpSum(1)");


    //We are waiting for all tasks to complete
    runtime.waitAllTasks();

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_19.dot", true);
    
    //We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_19.svg");   
    
    BOOST_MESSAGE("[INFO SPECX] : Value S= "<<valueS<<"\n"); 

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T19> in ms since start :"<<run_time.count()<<"\n");

}




BOOST_AUTO_TEST_CASE( test_specx_20 )
{
    //Demo to see the average time to insert a task with pressure 
    
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T20>\n");
    auto start_time= std::chrono::steady_clock::now();

    //int nbThreads = std::min(96,SpUtils::DefaultNumThreads());
    //int nbThreads = 6;
    int nbThreads = SpUtils::DefaultNumThreads();
    SpRuntime runtime(nbThreads);
    std::vector<int> valuesVec(nbThreads,0);
    
    int initVal = 1;
    int nVal=0;
    int nbLoops=1;
    int deltaTimeSleepBetweenTasks=10000; //break time between tasks in order to better see the graphic representation.
    bool QViewInfo=false;
    
    small_vector<SpAbstractTaskWithReturn<double>::SpTaskViewer> elapsed;
    elapsed.reserve(nbThreads*nbLoops);
    for(int idxLoop = 0 ; idxLoop < nbLoops ; ++idxLoop){
        for(int idx = 0 ; idx < nbThreads; ++idx){
            SpTimer timerTask;
            elapsed.emplace_back(     
                runtime.task(SpRead(idx),SpWrite(valuesVec.at(idx)),
                    [timerTask](const int& id,int& writeValParam) mutable -> double {
                        timerTask.stop();
                        usleep(1000);
                        writeValParam=sched_getcpu(); //returns the processor number on which the calling process is currently running.
                        return timerTask.getElapsed();
                    }
                )
            );
        }
        runtime.waitAllTasks();
        usleep(deltaTimeSleepBetweenTasks); 
    }

    double MeanTime = 0.0;
    int k=0;
    for(const SpAbstractTaskWithReturn<double>::SpTaskViewer& dt : elapsed) {
        if (QViewInfo) { std::cout<<"[INFO SPECX] : ["<<k<<"] dtime = "<<dt.getValue()<<" Num CPU="<<valuesVec[k]<<"\n"; }
        k++;
        MeanTime += dt.getValue();
    }
    MeanTime=MeanTime/double(nbThreads*nbLoops);
    
    double Variance=0.0;  
    for(const SpAbstractTaskWithReturn<double>::SpTaskViewer& dt : elapsed) {
        Variance+=(dt.getValue()-MeanTime)*(dt.getValue()-MeanTime);    
    }
    Variance=Variance/double(nbThreads*nbLoops);

    //We are waiting for all tasks to complete
    runtime.waitAllTasks();

    //We stop the completion of all tasks.
    runtime.stopAllThreads();

    //We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_20.dot", true);
    
    //We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_20.svg");   
    
    BOOST_MESSAGE("[INFO SPECX] : Value= "<<MeanTime<<"\n");    
    BOOST_MESSAGE("[INFO SPECX] : Variance= "<<Variance<<"\n");   

    if (QViewInfo) {
        std::cout<<"[INFO SPECX] : nbThreads= "<<nbThreads<<"\n";
        std::cout<<"[INFO SPECX] : Mean= "<<MeanTime<<"\n";
        std::cout<<"[INFO SPECX] : Variance= "<<Variance<<"\n";
        std::cout<<"[INFO SPECX] : Ecart-Type= "<<std::sqrt(Variance)<<"\n";
    }

    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T20> in ms since start :"<<run_time.count()<<"\n");

}



BOOST_AUTO_TEST_SUITE_END()



