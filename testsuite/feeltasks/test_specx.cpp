#define BOOST_TEST_MODULE test_specx
#include <feel/feelcore/testsuite.hpp>

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

#include <specx/Legacy/SpRuntime.hpp>


#include "UTester.hpp"
#include "utestUtils.hpp"




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

    BOOST_MESSAGE("[INFO SPECX] : CTRL VALUES initval="<<initVal<<" writeVal="<<writeVal<<" res="<<res<<"\n");
    
    //We calculate the time frame with the “SPECX” clock in ms
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T2> in ms since start :"<<run_time.count()<<"\n");

        // We generate the task graph corresponding to the execution 
    runtime.generateDot("test_specx_2.dot", true);
    
    // We generate an Svg trace of the execution
    runtime.generateTrace("test_specx_2.dot");   
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
    const int NumThreads = SpUtils::DefaultNumThreads();
    SpRuntime runtime(NumThreads);
    auto start_time= std::chrono::steady_clock::now();
    const int NbTasksToSubmit = runtime.getNbThreads();
    static const int NbLoops = 1;
    int initVal = 1;

    small_vector<SpAbstractTaskWithReturn<double>::SpTaskViewer> elapsed;
    elapsed.reserve(NbTasksToSubmit*NbLoops);


    for(int idxLoop = 0 ; idxLoop < NbLoops ; ++idxLoop){
            for(int idx = 0 ; idx < NbTasksToSubmit ; ++idx){
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
            averageToExecute += viewer.getValue()/double(NbTasksToSubmit*NbLoops);
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
        return true; // Here we always return true, this basically means
                     // that we always allow speculative tasks to run
                     // regardless of runtime conditions.
    });

    int a = 41, b = 0, c = 0;
    auto task1 = runtime.task(SpRead(a), [](const int& inA) -> int {
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
    [dis01, mtEngine] (const int &inB, int &inA) mutable -> bool{
        double val = dis01(mtEngine);
        
        if(inB == 42  && val < 0.5) {
            inA = 43;
            return true;
        }
        
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
        runtime.generateDot("test_specx_8.dot", true);
    
        // We generate an Svg trace of the execution
        runtime.generateTrace("test_specx_8.svg");   
    }
    //We calculate the time frame with the “SPECX” clock
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T8> in ms since start :"<<run_time.count()<<"\n"); 
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
    runtime.generateTrace("test_specx_10.dot");   
}


BOOST_AUTO_TEST_CASE( test_specx_11 )
{
    //PARALLEL_WRITE WITH ATOMIC
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T11>\n");
    auto start_time= std::chrono::steady_clock::now();
    int NumThreads = std::min(6,SpUtils::DefaultNumThreads());
    SpRuntime runtime(NumThreads);
    std::atomic<int> initVal(0);
    int nbTh=3;
    int val0=20;
    for(int idxTh = 0 ; idxTh < nbTh; ++idxTh){
        runtime.task(SpParallelWrite(initVal),
            [&](std::atomic<int>& initValParam){
                initValParam += 1;
                val0++;
                //std::cout<<idxTh<<" "<<initValParam<<"...."<<val0<<"\n";
                while(initValParam != nbTh){
                    usleep(10000); //<== simply to see the task boxes better on the graph.
                } 
                //std::cout<<"<OK>\n";
            }).setTaskName("OpPW("+std::to_string(idxTh)+")");

        runtime.task(SpRead(val0),
            [&](const int& valParam0)
            {
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
    //std::cout<<"Start 13\n";
    //MERGE
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T12>\n");
    auto start_time= std::chrono::steady_clock::now();
    
    int NumThreads = std::min(25,SpUtils::DefaultNumThreads());
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


BOOST_AUTO_TEST_SUITE_END()



