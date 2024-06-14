
 
//=======================================================================================================================
// Meta function tools allowing you to process an expression defined in Task
//=======================================================================================================================


constexpr auto& _parameters = NA::identifier<struct parameters_tag>;
constexpr auto& _tasks = NA::identifier<struct task_tag>;



namespace SUBTASK{

    template<typename ...T, size_t... I>
    auto extractParametersAsTuple( std::tuple<T...> && t, std::index_sequence<I...>)
    {
        return std::forward_as_tuple( std::get<I>(t).getValue()...);
    }

    struct Runtime{
        template <typename ... Ts>
        void task(Ts && ... ts ) {
            auto t = std::make_tuple( std::forward<Ts>(ts)... );
            auto callback = std::get<sizeof...(Ts) - 1>(t);
            auto parameters = extractParametersAsTuple( std::move(t), std::make_index_sequence<sizeof...(Ts)-1>{} );
            std::apply( callback, std::move(parameters) );
        }
    };

    template <typename T,bool b>
    class SpData
    {
        static_assert(std::is_reference<T>::value,
                    "The given type must be a reference");
    public:
        using value_type = T;
        static constexpr bool isWrite = b;

        template <typename U, typename = std::enable_if_t<std::is_convertible_v<U,T>> >
        constexpr explicit SpData( U && u ) : M_val( std::forward<U>(u) ) {}

        constexpr value_type getValue() { return M_val; }
    private:
        value_type M_val;
    };

    template <typename T>
    auto spRead( T && t )
    {
        return SpData<T,false>{ std::forward<T>( t ) };
    }
    template <typename T>
    auto spWrite( T && t )
    {
        return SpData<T,true>{ std::forward<T>( t ) };
    }

    template<typename T>
    auto toSpData( T && t )
    {
        if constexpr ( std::is_const_v<std::remove_reference_t<T>> )
            return spRead( std::forward<T>( t ) );
        else
            return spWrite( std::forward<T>( t ) );
    }

    template<typename ...T, size_t... I>
    auto makeSpDataHelper( std::tuple<T...>& t, std::index_sequence<I...>)
    {
        return std::make_tuple( toSpData(std::get<I>(t))...);
    }
    template<typename ...T>
    auto makeSpData( std::tuple<T...>& t ){
        return makeSpDataHelper<T...>(t, std::make_index_sequence<sizeof...(T)>{});
    }

    template<typename T>
    auto toSpDataSpecx( T && t )
    {
        if constexpr ( std::is_const_v<std::remove_reference_t<T>> )
            return SpRead(std::forward<T>( t ));
        else
            return SpWrite(std::forward<T>( t ));
    }

    template<typename ...T, size_t... I>
    auto makeSpDataHelperSpecx( std::tuple<T...>& t, std::index_sequence<I...>)
    {
        return std::make_tuple( toSpDataSpecx(std::get<I>(t))...);
    }
    template<typename ...T>
    auto makeSpDataSpecx( std::tuple<T...>& t ){
        return makeSpDataHelperSpecx<T...>(t, std::make_index_sequence<sizeof...(T)>{});
    }



    template<typename T>
    auto toSpDataSpecxGPU( T && t )
    {
        if constexpr ( std::is_const_v<std::remove_reference_t<T>> )
            return SpRead(std::forward<T>( t ));
        else
            //return SpWrite(std::forward<T>( t ));
            return SpCommutativeWrite(std::forward<T>( t ));
    }

    template<typename ...T, size_t... I>
    auto makeSpDataHelperSpecxGPU( std::tuple<T...>& t, std::index_sequence<I...>)
    {
        return std::make_tuple( toSpDataSpecxGPU(std::get<I>(t))...);
    }
    template<typename ...T>
    auto makeSpDataSpecxGPU( std::tuple<T...>& t ){
        return makeSpDataHelperSpecxGPU<T...>(t, std::make_index_sequence<sizeof...(T)>{});
    }


    template <typename ... Ts>
    auto parameters(Ts && ... ts)
    {
        return std::forward_as_tuple( std::forward<Ts>(ts)... );
    }
}


/*
namespace Frontend
{
    template <typename ... Ts>
    auto parameters(Ts && ... ts)
    {
        return std::forward_as_tuple( std::forward<Ts>(ts)... );
    }
}
*/


//================================================================================================================================
// CLASS Task: Provide a family of multithreaded functions...
//================================================================================================================================

// Nota: The objective is to provide a range of tools in the case of using a single variable in multithreading.
// In the case of work with several variables use the class TasksDispatchComplex.

//#define USE_jthread



namespace TASK {


void *workerNumCPU(void *arg) {
    // Function used to run a task on a given CPU number.
    std::function<void()> *func = (std::function<void()>*)arg;
    (*func)();
    pthread_exit(NULL);
}


class Task
{
    private:
        int M_nbThTotal;                                       // variable indicator of the total number of threads available
        std::string M_FileName;                                // name of the recording file that will be used for the debriefing
        template <typename ... Ts>
        auto parameters(Ts && ... ts);                         // meta function to use multiple variables
        bool M_qEmptyTask;                                     // variable indicator if there is no task
        bool M_qFlagDetachAlert;                               // flag internal variable indicator for activation of the detach module
        
        //SpTaskGraph<SpSpeculativeModel::SP_NO_SPEC> M_mytg;    // Specx TaskGraph function 
        SpTaskGraph<SpSpeculativeModel::SP_MODEL_1> M_mytg;    // Specx TaskGraph function 
        SpComputeEngine M_myce;                                // Specx engine function 

        std::vector<int>               M_idType;                // This vector saves the thread type according to the task id number
        std::vector<int>               M_numTaskStatus;         // This vector saves the thread type status according to the task id number
        std::vector<std::future<bool>> M_myfutures;             // Vector table of std::future
        std::vector<std::future<bool>> M_myfuturesdetach;       // Vector table of std::future detach
        std::vector<std::thread>       M_mythreads;             // Vector table of std::thread
        std::vector<pthread_t>         M_mypthread_t;           // Vector table of pthread_t

        #ifdef COMPILE_WITH_CXX_20
            std::vector<std::jthread>      myjthreads;              // Vector table of <std::jthread
        #endif

        //pthread_t M_mypthread_t[100]; 
        pthread_attr_t                 M_mypthread_attr_t;      // An attribute set to supply to pthread_create()
        std::vector<int>               M_mypthread_cpu;         // Vector table of thread with
        std::mutex                     M_mymtx;                 // Mutex
        //std::promise<int> promise0;
        std::chrono::steady_clock::time_point M_t_begin,M_t_end;

        //BEGIN::GPU part

        //END::GPU part


        // variables indicatrices utilisées par les fonctions set et get. Voir les descriptions ci-dessous.
        long int M_t_laps;
        bool M_qFirstTask;
        int  M_idk;
        int  M_numLevelAction;
        int  M_nbThreadDetach;
        bool M_qReady;

        bool M_qCUDA;
        bool M_qHIP;

        int  M_nbTh;
        int  M_numTypeTh;

        bool M_qDetach;
        bool M_qYield;
        bool M_qDeferred;
        bool M_qUseIndex;

        bool M_qViewChrono;
        bool M_qInfo;
        bool M_qSave;

        //BEGIN::GPU part
        int  M_numBlocksGPU;
        int  M_nThPerBckGPU;
        int  M_numOpGPU;
        //END::GPU part

        template <typename ... Ts> 
            auto common(Ts && ... ts);

        void init(); // Function that initializes all variables we need
        

    public:

        //BEGIN::No copy and no move
        Task(const Task&) = delete;
        Task(Task&&) = delete;
        Task& operator=(const Task&) = delete;
        Task& operator=(Task&&) = delete;
        //END::No copy and no move

        //BEGIN::Small functions and variables to manage initialization parameters
        
        void setDetach           (bool b)  {  M_qDetach     = b; } // This function allows you to indicate to the next task whether it should be detached.
        void setYield            (bool b)  {  M_qYield      = b; } // This function allows you to indicate to the next task whether it should be yield.
        void setDeferred         (bool b)  {  M_qDeferred   = b; } // This function allows you to indicate to the next task whether it should be deferred. Used only for std::async part
        void setUseIndex         (bool b)  {  M_qUseIndex   = b; } // 
        void setSave             (bool b)  {  M_qSave       = b; } // Indicates whether to save debriefing data
        void setInfo             (bool b)  {  M_qInfo       = b; } // Allows Indicates whether we must display progress information during the execution of the class.
        void setViewChrono       (bool b)  {  M_qViewChrono = b; } // Allows Indicates whether to display time measurement information from the stopwatch.
        


        bool isDetach            () const  {  return(M_qDetach);     } // Indicates the state of the boolean variable isDetach
        bool isYield             () const  {  return(M_qYield);      } // Indicates the state of the boolean variable isYield 
        bool isDeferred          () const  {  return(M_qDeferred);   } // Indicates the state of the boolean variable isDeferred
        bool isSave              () const  {  return(M_qSave);       } // Indicates the state of the boolean variable isSave
        bool isInfo              () const  {  return(M_qInfo);       } // Indicates the state of the boolean variable isInfo 
        bool isViewChrono        () const  {  return(M_qViewChrono); } // Indicates the state of the boolean variable isViewChrono


        void setNbThread         (int v)   { M_nbTh=std::min(v,M_nbThTotal); } // Fix the desired number of threads
        int  getNbMaxThread      ()        { M_nbThTotal=std::thread::hardware_concurrency(); return(M_nbThTotal); } // Gives the maximum number of threads
        int  getNbThreads        () const  { int val=M_nbTh; if (M_numTypeTh==3) { val=static_cast<int>(M_myce.getCurrentNbOfWorkers()); } return val; } // Gives the number of threads used
        int  getNbCpuWorkers     () const  { int val=M_nbTh; if (M_numTypeTh==3) { val=static_cast<int>(M_myce.getNbCpuWorkers()); } return val; } //Gives the number of CpuWorkers used

        auto getIdThread         (int i);  // Gives the memory address of the thread used
        int  getNbThreadPerBlock (int i);  // Gives the number of threads used per GPU block

        void setNumOpGPU         (int v)  {  M_numOpGPU = v; } 


        long int  getTimeLaps ()        { return M_t_laps; } // Gives the time laps simulation


        #ifdef COMPILE_WITH_CUDA
            int getNbCudaWorkers() const {
                return static_cast<int>(M_myce.getNbCudaWorkers());
                } // Returns the total number of Cuda cards
        #endif
        #ifdef COMPILE_WITH_HIP
            int getNbHipWorkers() const {
                return static_cast<int>(M_myce.getNbHipWorkers());
            } // Returns the total number of Hip cards
        #endif
            
        void setFileName(std::string s) { M_FileName=s; } 
        //END::Small functions and variables to manage initialization parameters


        void subTask(const int nbThread,const int nbBlocks,int M_numTypeThread);
           

        Task(void);  // Constructor, initializes the default variables that are defined in the init() function.
        ~Task(void); // Destructor is invoked automatically whenever an object is going to be destroyed. Closes the elements properly


        #if defined(COMPILE_WITH_HIP) || defined(COMPILE_WITH_CUDA)
            #ifdef COMPILE_WITH_HIP
                explicit Task(const int nbThread,const int nbBlocks,int M_numTypeThread):M_mytg(),M_myce(SpWorkerTeamBuilder::TeamOfCpuHipWorkers()) // Class constructor in classic mode
                {
                    std::cout<<"[INFO]: WELCOME TO GPU:HIP"<<"\n";
                    subTask(nbThread,nbBlocks,M_numTypeThread);
                }


                explicit Task(const int nbThread,int M_numTypeThread):M_mytg(),M_myce(SpWorkerTeamBuilder::TeamOfCpuHipWorkers()) // Class constructor in classic mode
                {
                    std::cout<<"[INFO]: WELCOME TO GPU:HIP"<<"\n";
                    subTask(nbThread,1,M_numTypeThread);
                }


            #endif

            #ifdef COMPILE_WITH_CUDA
                explicit Task(const int nbThread,const int nbBlocks,int M_numTypeThread):M_mytg(),M_myce(SpWorkerTeamBuilder::TeamOfCpuCudaWorkers()) // Class constructor in classic mode
                {
                    std::cout<<"[INFO]: WELCOME TO GPU:CUDA"<<"\n";
                    subTask(nbThread,nbBlocks,M_numTypeThread);
                }
            #endif
        #else
            explicit Task(const int nbThread,int M_numTypeThread):M_mytg(),M_myce(SpWorkerTeamBuilder::TeamOfCpuWorkers(nbThread)) // Class constructor in classic mode
            {
                std::cout<<"[INFO]: WELCOME TO CPU"<<"\n";
                subTask(nbThread,1,M_numTypeThread);
            }
        #endif

        template <class ClassFunc> void execOnWorkers(ClassFunc&& func) { M_myce.execOnWorkers(std::forward<ClassFunc>(func)); } //Execute a ClassFunc on workers
  
        void setSpeculationTest(std::function<bool(int,const SpProbability&)> inFormula){
            if (M_numTypeTh==3) { M_mytg.setSpeculationTest(std::move(inFormula)); }
        }

        template <typename ... Ts>
        void add( Ts && ... ts ); // This main function allows you to add a task

            template <typename ... Ts>
                void addTaskSimple( Ts && ... ts ); // This subfunction allows you to add a simple task

            template <typename ... Ts>
                void addTaskSpecx( Ts && ... ts ); // This subfunction allows you to add a specx task

            template <typename ... Ts>
                void addTaskAsync( Ts && ... ts ); // This subfunction allows you to add a std::async task

            template <typename ... Ts>
                void addTaskMultithread( Ts && ... ts ); // This subfunction allows you to add a multithread task

            
            #ifdef COMPILE_WITH_CXX_20
            template <typename ... Ts>
                void addTaskjthread( Ts && ... ts ); // This subfunction allows you to add a jthread task. Only works under C++20
            #endif

        #ifdef COMPILE_WITH_HIP
            template <typename ... Ts>
                void add_GPU( Ts && ... ts); // This main function allows you to add a task O:CPU Normal  1:SpHip  2:SpCuda  3:GPU to CPU
        #endif

        template <class InputIterator,typename ... Ts>
            void for_each(InputIterator first, InputIterator last,Ts && ... ts); // This function allows you to apply the same treatment to a set of elements of a task.

        
        template <typename ... Ts>
            void add(int numCPU,Ts && ... ts);  // Add a task on specific CPU number

        template<typename FctDetach>
            auto add_detach(FctDetach&& func) -> std::future<decltype(func())>; // Add a detach thread task

        template <typename ... Ts>
            void runInCPUs(const std::vector<int> & numCPU,Ts && ... ts); // Execute all tasks on specific CPU number
        
        void run();                     // Execution of all added tasks. 
        void close();                   // Memory cleanup of all variables used before closing the class.
        void debriefingTasks();         // Wrote a report on execution times and generated .dot .svg files regarding Specx and .csv to save the times.
        void getInformation();          // Provides all information regarding graphics cards (CUDA and HIP)
        void getGPUInformationError();  // Provides error types from the GPU (CUDA and HIP)

        
        //GPU-AMD-CUDA

        #ifdef USE_GPU_HIP
        // set of functions allowing you to use eigen under Hip gpu.
        template<typename Kernel, typename Input, typename Output>
            void run_gpu_1D(const Kernel& kernel_function,dim3 blocks,int n,const Input& in,Output& out);
        template<typename Kernel, typename Input, typename Output>
            void run_gpu_2D(const Kernel& kernel_function,dim3 blocks,int n,const Input& in,Output& out);
        template<typename Kernel, typename Input, typename Output>
            void run_gpu_3D(const Kernel& kernel_function,dim3 blocks,int n,const Input& in,Output& out);
        // set of functions allowing you to use eigen under Hip cpu. Juste to control the results
        template<typename Kernel, typename Input, typename Output>
            void run_cpu_1D(const Kernel& kernel_function, int n, const Input& in, Output& out);
        #endif

        #ifdef UseHIP

        #endif

        #if defined(COMPILE_WITH_HIP) || defined(COMPILE_WITH_CUDA)

        template <class... ParamsTy>
            void addTaskSpecxPureGPU(ParamsTy&&...params); 

        template <typename ... Ts>
            void addTaskSpecxGPU( Ts && ... ts); // This subfunction allows you to add a specx task

        #endif

};



Task::Task()
{
    init();
}

Task::~Task()
{
    //Add somes   
    if ((M_numTypeTh== 3) && (M_numLevelAction==3)) {  M_myce.stopIfNotAlreadyStopped(); } 
    if ((M_numTypeTh==33) && (M_numLevelAction==3)) {  M_myce.stopIfNotAlreadyStopped(); } // <== [ ] see if we really need it
    //Specx
}


void Task::init()
{
    M_nbThTotal=std::thread::hardware_concurrency();
    M_nbTh             = M_nbThTotal;
    M_qInfo            = true;
    M_qSave            = false;
    M_qDeferred        = false;
    M_numTypeTh        = 0;
    M_qUseIndex        = false;
    M_FileName         ="NoName";
    M_qFirstTask       = true;
    M_idk              = 0;
    M_numLevelAction   = 0;
    M_qCUDA            = false;
    M_qHIP             = false;
    M_qEmptyTask       = true;
    M_qFlagDetachAlert = false;
    M_nbThreadDetach   = 0;
    M_qReady           = false;
    M_qYield           = false;
    M_numBlocksGPU     = 128; //<- see after

    M_numOpGPU         = 0;
    
    M_mythreads.clear();
    M_myfutures.clear();
    
    #ifdef COMPILE_WITH_CXX_20
    myjthreads.clear();
    #endif
}


void Task::subTask(const int nbThread,const int nbBlocks,int M_numTypeThread)
{
    M_numBlocksGPU=nbBlocks;
    M_nbTh=nbThread;
    M_numTypeTh = M_numTypeThread;
    M_idk=0; M_numTaskStatus.clear(); M_qDetach=0; M_nbThreadDetach=0;
    M_idType.clear(); 
    M_t_begin   = std::chrono::steady_clock::now();
    M_qDeferred = false;
    M_qReady    = false;
                
    if (M_numTypeTh==0) { } // No Thread
    if (M_numTypeTh==1) { } // multithread
    if (M_numTypeTh==2) { } // std::async
    if (M_numTypeTh==3) { M_mytg.computeOn(M_myce); } // Specx

    if (M_numTypeTh==10) { pthread_attr_init(&M_mypthread_attr_t); } //pthread

    if (M_numTypeTh==33) { 
        #ifdef COMPILE_WITH_HIP
            //static_assert(SpDeviceDataView<std::vector<int>>::MoveType == SpDeviceDataUtils::DeviceMovableType::STDVEC,"should be stdvec"); //will see after
            M_mytg.computeOn(M_myce);
        #endif
    } // Specx GPU
}

int Task::getNbThreadPerBlock(int i)
{
    int numDevices=0;
    #ifdef COMPILE_WITH_HIP
        hipGetDeviceCount(&numDevices);
        hipDeviceProp_t devProp;
        if ((i>=0) && (i<numDevices))
        {
            HIP_CHECK(hipGetDeviceProperties(&devProp,i));
            return(devProp.maxThreadsPerBlock);
        }
    #endif

    #ifdef COMPILE_WITH_CUDA
        cudaGetDeviceCount(&numDevices);
        cudaDeviceProp_t devProp;
        if ((i>=0) && (i<numDevices))
        {
            cudaGetDeviceProperties(&devProp,i);
            return(devProp.maxThreadsPerBlock);
        }
    #endif
}



void Task::getInformation()
{
    // Provides all information regarding graphics cards CUDA and AMD
    if (M_qInfo)
    {
        if (M_numTypeTh== 0) { std::cout<<"[INFO]: Mode No Thread\n"; }
        if (M_numTypeTh== 1) { std::cout<<"[INFO]: Mode Multithread\n"; }
        if (M_numTypeTh== 2) { std::cout<<"[INFO]: Mode Std::async\n"; }
        if (M_numTypeTh== 3) { std::cout<<"[INFO]: Mode Specx\n"; }

        if (M_numTypeTh==10) { std::cout<<"[INFO]: Mode Thread in CPU\n"; }
        M_nbThTotal=getNbMaxThread(); 
        std::cout<<"[INFO]: Nb max Thread="<<M_nbThTotal<<"\n";
        HARD::getInformationSystem();
    }
}


void Task::getGPUInformationError()
{
    #ifdef COMPILE_WITH_CUDA
        cudaError_t num_err = cudaGetLastError();
        if (num_err != cudaSuccess) {
            std::cout<<"[INFO]: hip Error Name ="<<cudaGetErrorName(num_err)<<" "<<cudaGetErrorString(num_err)<<"\n";
        }
        num_err = cudaDeviceSynchronize();
        if (num_err != cudaSuccess) {
            std::cout<<"[INFO]: hip Error Name="<<cudaGetErrorName(err)<<" "<<cudaGetErrorString(num_err)<<"\n";
        }
    #endif

    #ifdef COMPILE_WITH_HIP
        hipError_t num_err = hipGetLastError();
        if (num_err != hipSuccess) {
            std::cout<<"[INFO]: hip Error Name ="<<hipGetErrorName(num_err)<<" "<<hipGetErrorString(num_err)<<"\n";
        }
        num_err = hipDeviceSynchronize();
        if (num_err != hipSuccess) {
            std::cout<<"[INFO]: hip Error Name="<<hipGetErrorName(num_err)<<" "<<hipGetErrorString(num_err)<<"\n";
        }
    #endif
}


#ifdef USE_GPU_HIP
template<typename Kernel, typename Input, typename Output>
void Task::run_gpu_1D(const Kernel& kernel_function,dim3 blocks,int n, const Input& in, Output& out)
{
    bool qInfo_GPU_process=true;
    if (M_qFirstTask) { M_t_begin = std::chrono::steady_clock::now(); M_qFirstTask=false;}
    std::chrono::steady_clock::time_point M_t_begin_all_GPU_process,M_t_end_all_GPU_process;
    std::chrono::steady_clock::time_point M_t_begin_inside,M_t_end_inside;

    M_t_begin_all_GPU_process = std::chrono::steady_clock::now();
        typename Input::Scalar*  d_in;
        typename Output::Scalar* d_out;
        std::ptrdiff_t in_bytes  = in.size()  * sizeof(typename Input::Scalar);
        std::ptrdiff_t out_bytes = out.size() * sizeof(typename Output::Scalar);
        
        HIP_ASSERT(hipMalloc((void**)(&d_in),  in_bytes));
        HIP_ASSERT(hipMalloc((void**)(&d_out), out_bytes));
        
        HIP_ASSERT(hipMemcpy(d_in,  in.data(),  in_bytes,  hipMemcpyHostToDevice));
        HIP_ASSERT(hipMemcpy(d_out, out.data(), out_bytes, hipMemcpyHostToDevice));
        
        //dim3 blocks(128);
        //dim3 grids( (n+int(blocks.x)-1)/int(blocks.x) );
        int num_blocks = (n + blocks.x - 1) / blocks.x;
	    dim3 thread_block(blocks.x, 1, 1);
	    dim3 grid(num_blocks,1,1);

        hipDeviceSynchronize();
            M_t_begin_inside = std::chrono::steady_clock::now();
            hipLaunchKernelGGL(HIP_KERNEL_NAME(
                OP_IN_KERNEL_GPU_1D<Kernel,typename std::decay<decltype(*d_in)>::type,typename std::decay<decltype(*d_out)>::type>), 
                        grid, thread_block, 0, 0, kernel_function, n, d_in, d_out);

            M_t_end_inside = std::chrono::steady_clock::now();
        getGPUInformationError();

        hipMemcpy(const_cast<typename Input::Scalar*>(in.data()),  d_in,  in_bytes,  hipMemcpyDeviceToHost);
        hipMemcpy(out.data(), d_out, out_bytes, hipMemcpyDeviceToHost);
        HIP_ASSERT(hipFree(d_in));
        HIP_ASSERT(hipFree(d_out));    
    M_t_end_all_GPU_process = std::chrono::steady_clock::now();

    long int M_t_laps3= std::chrono::duration_cast<std::chrono::microseconds>(M_t_end_inside - M_t_begin_inside).count();
    long int M_t_laps2= std::chrono::duration_cast<std::chrono::microseconds>(M_t_end_all_GPU_process - M_t_begin_all_GPU_process).count();
    if (qInfo_GPU_process) {
        std::cout << "[INFO]: Elapsed microseconds inside: "<<M_t_laps3<< " us\n";    
        std::cout << "[INFO]: Elapsed microseconds inside + memory copy: "<<M_t_laps2<< " us\n";
        std::cout << "[INFO]: nb grids: "<<grid.x<<" "<<grid.y<<" "<<grid.z<< "\n";
        std::cout << "[INFO]: nb block: "<<thread_block.x<<" "<<thread_block.y<<" "<<thread_block.z<< "\n";
    }
}

template<typename Kernel, typename Input, typename Output>
void Task::run_gpu_2D(const Kernel& kernel_function,dim3 blocks,int n, const Input& in, Output& out)
{
    if (M_qFirstTask) { M_t_begin = std::chrono::steady_clock::now(); M_qFirstTask=false;}
}

template<typename Kernel, typename Input, typename Output>
void Task::run_gpu_3D(const Kernel& kernel_function,dim3 blocks,int n, const Input& in, Output& out)
{
    if (M_qFirstTask) { M_t_begin = std::chrono::steady_clock::now(); M_qFirstTask=false;}
}


template<typename Kernel, typename Input, typename Output>
void Task::run_cpu_1D(const Kernel& kernel_function, int n, const Input& in, Output& out)
{
  if (M_qFirstTask) { M_t_begin = std::chrono::steady_clock::now(); M_qFirstTask=false;}
  for(int i=0; i<n; i++)
    kernel_function(i, in.data(), out.data());
}
#endif







#if defined(COMPILE_WITH_HIP) || defined(COMPILE_WITH_CUDA)

template <class... ParamsTy>
void Task::addTaskSpecxPureGPU(ParamsTy&&...params)
{
      M_mytg.task(std::forward<ParamsTy>(params)...);
}


template <typename ... Ts>
void Task::addTaskSpecxGPU( Ts && ... ts)
{
    //=33
    auto args = NA::make_arguments( std::forward<Ts>(ts)... );
    auto && task = args.get(_tasks);
    auto && parameters = args.get_else(_parameters,std::make_tuple());
    
    //std::cout<<"GET Task="<<std::get<0>(args)<<"\n";
    //std::cout<<"GET Parameters="<<std::get<0>(parameters)<<"\n";



    auto tpSUBTASK=SUBTASK::makeSpDataSpecxGPU( parameters );
    int NbtpSUBTASK=std::tuple_size<decltype(tpSUBTASK)>::value;
    std::cout <<"Size Parameters="<<NbtpSUBTASK<< std::endl;

    auto LambdaExpression=std::make_tuple(task);
    int NbLambdaExpression=std::tuple_size<decltype(LambdaExpression)>::value;
    std::cout <<"Size Tuple Task="<<NbLambdaExpression<< std::endl;

    //auto tpHip=std::tuple_cat(tpSUBTASK,LambdaExpression);  

    //M_mytg.task(std::get<0>(tpSUBTASK),std::get<1>(tpSUBTASK),SpHip(task));


    //std::apply([&](auto &&... args) { M_mytg.task(args...).setTaskName("Op("+std::to_string(M_idk)+")"); },tpHip);    
    //addTaskSpecxPureGPU(tpHip);   



    if (!M_qDetach)
    {
        if (M_numOpGPU==0) { 
            
            auto tp=std::tuple_cat( SUBTASK::makeSpDataSpecxGPU( parameters ), std::make_tuple( task ));  
            std::apply([&](auto &&... args) { M_mytg.task(args...).setTaskName("Op("+std::to_string(M_idk)+")"); },tp);       
        } //CPU Normal


        if (M_numOpGPU==1)
        { 
            //auto tp=std::tuple_cat( SUBTASK::makeSpDataSpecxGPU( parameters ), SpHip(std::make_tuple(task))); 
            //std::apply([&](auto &&... args) { M_mytg.task(args...).setTaskName("Op("+std::to_string(M_idk)+")"); },tp);
        } //op in Hip


/*

        if (M_numOpGPU==2) { 
            auto tp=std::tuple_cat( SUBTASK::makeSpDataSpecxGPU( parameters ), SpCuda(std::make_tuple( task )));   
            std::apply([&](auto &&... args) { M_mytg.task(args...).setTaskName("Op("+std::to_string(M_idk)+")"); },tp);  
        } //op in Cuda

        if (M_numOpGPU==3) { 
            auto tp=std::tuple_cat( SUBTASK::makeSpDataSpecx( parameters ), SpCPU(std::make_tuple( task )));
            std::apply([&](auto &&... args) { M_mytg.task(args...).setTaskName("Op("+std::to_string(M_idk)+")"); },tp);
        } //op in GPU to CPU
*/
        usleep(0); std::atomic_int counter(0);
    }
    else
    {
        addTaskAsync(std::forward<Ts>(ts)...);
    }
}
#endif





template <typename ... Ts>
auto Task::parameters(Ts && ... ts)
{
    return std::forward_as_tuple( std::forward<Ts>(ts)... );
}

template <typename ... Ts>
    auto Task::common(Ts && ... ts)
{
    auto args = NA::make_arguments( std::forward<Ts>(ts)... );
    auto && task = args.get(_tasks);
    auto && parameters = args.get_else(_parameters,std::make_tuple());
    auto tp=std::tuple_cat(  SUBTASK::makeSpData( parameters ), std::make_tuple( task ) );
    return(tp);
}


template <typename ... Ts>
void Task::addTaskSimple( Ts && ... ts )
{
    auto tp=common(std::forward<Ts>(ts)...);
    SUBTASK::Runtime runtime;
    std::apply( [&runtime](auto... args){ runtime.task(args...); }, tp );
}


template <typename ... Ts>
void Task::addTaskMultithread( Ts && ... ts )
{
    auto tp=common(std::forward<Ts>(ts)...);
    SUBTASK::Runtime runtime;
	auto LamdaTransfert = [&]() {
			std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
            return true; 
	};

    if (!M_qDetach)
    {
        std::thread th(LamdaTransfert);
        M_mythreads.push_back(std::move(th));
    }
    else
    {
        if (M_qInfo) { std::cout<<"[INFO]: detach in process...\n"; }
        M_myfuturesdetach.emplace_back(add_detach(LamdaTransfert)); 
    }
    usleep(1);
}


template <typename ... Ts>
void Task::addTaskAsync( Ts && ... ts )
{
    auto tp=common(std::forward<Ts>(ts)...);
    SUBTASK::Runtime runtime;
    auto LamdaTransfert = [&]() {
			std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
            return true; 
		};

    if (!M_qDetach)
    {
        if (M_qDeferred) { M_myfutures.emplace_back(std::async(std::launch::deferred,LamdaTransfert));}
        else             { M_myfutures.emplace_back(std::async(std::launch::async,LamdaTransfert)); }
    }
    else
    {
        if (M_qInfo) { std::cout<<"[INFO]: detach in process...\n"; }
        M_myfuturesdetach.emplace_back(add_detach(LamdaTransfert)); 
    }

    usleep(1);
}



template <typename ... Ts>
void Task::addTaskSpecx( Ts && ... ts )
{
    auto args = NA::make_arguments( std::forward<Ts>(ts)... );
    auto && task = args.get(_tasks);
    auto && parameters = args.get_else(_parameters,std::make_tuple());
    auto tp=std::tuple_cat( 
					SUBTASK::makeSpDataSpecx( parameters ), 
					std::make_tuple( task ) 
				);
    if (!M_qDetach)
    {
        std::apply([&](auto &&... args) { M_mytg.task(args...).setTaskName("Op("+std::to_string(M_idk)+")"); },tp);
        usleep(0); std::atomic_int counter(0);
    }
    else
    {
        addTaskAsync(std::forward<Ts>(ts)...);
    }
}

#ifdef COMPILE_WITH_CXX_20
template <typename ... Ts>
void Task::addTaskjthread( Ts && ... ts )
{
    auto tp=common(std::forward<Ts>(ts)...);
    SUBTASK::Runtime runtime;
	auto LamdaTransfert = [&]() {
			std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
            return true; 
	};

    if (!M_qDetach)
    {
        std::jthread th(LamdaTransfert);
        myjthreads.push_back(std::move(th));
    }
    else
    {
        if (M_qInfo) { std::cout<<"[INFO]: detach in process...\n"; }
        M_myfuturesdetach.emplace_back(add_detach(LamdaTransfert)); 
    }
    usleep(1);
}
#endif


template <typename ... Ts>
void Task::add( Ts && ... ts )
{
    M_numLevelAction=1;
    M_qEmptyTask=false;
    M_idk++; M_idType.push_back(M_numTypeTh); M_numTaskStatus.push_back(M_qDetach);
    if (M_qDetach) { M_qFlagDetachAlert=true; M_nbThreadDetach++; }
    if (M_qFirstTask) { M_t_begin = std::chrono::steady_clock::now(); M_qFirstTask=false;}
    //std::cout<<"M_numTypeTh="<<M_numTypeTh<<"\n";
    switch(M_numTypeTh) {
        case 1: addTaskMultithread(std::forward<Ts>(ts)...); //multithread
        break;
        case 2: addTaskAsync(std::forward<Ts>(ts)...); //std::async
        break;
        case 3: addTaskSpecx(std::forward<Ts>(ts)...); //Specx
        break;
        #ifdef COMPILE_WITH_CXX_20
            case 4: addTaskjthread(std::forward<Ts>(ts)...); //std::jthread
            break;
        #endif
        default: addTaskSimple(std::forward<Ts>(ts)...); //No Thread
    }
    M_qDetach=false;
    //std::cout<<"[INFO] :ADD OK"<<std::endl;
}

#ifdef COMPILE_WITH_HIP
template <typename ... Ts>
void Task::add_GPU( Ts && ... ts)
{
    M_numLevelAction=1;
    M_qEmptyTask=false;
    M_idk++; M_idType.push_back(M_numTypeTh); M_numTaskStatus.push_back(M_qDetach);
    if (M_qDetach) { M_qFlagDetachAlert=true; M_nbThreadDetach++; }
    if (M_qFirstTask) { M_t_begin = std::chrono::steady_clock::now(); M_qFirstTask=false;}
    switch(M_numTypeTh) {
        case 33: addTaskSpecxGPU(std::forward<Ts>(ts)...); //Specx GPU
        break;
    }
    M_qDetach=false;
}
#endif


auto Task::getIdThread(int i)
{
    if ((i>0) && (i<M_idType.size()))
    {
        int nb=M_idType.size();
        int numM_idType=M_idType[i];
        int k=0; int ki=-1; bool qOn=true;
        while (qOn)
        {
            if (M_idType[k]==numM_idType) { ki++; }
            k++; if (k>i) { qOn=false; }
        }

        switch(numM_idType) {
            case 1: return(M_mythreads[ki].get_id());//multithread
            break;
            case 2:  //std::async
            break;
            case 3: //Specx
            break;
            #ifdef COMPILE_WITH_CXX_20
            case 4: return(myjthreads[ki].get_id()); //std::jtread
            break;
            #endif

        }
    }
    return(M_mythreads[0].get_id());
}


template <class InputIterator,typename ... Ts>
    void Task::for_each(InputIterator first, InputIterator last,Ts && ... ts)
{
    M_qUseIndex=true; //Iterator used
    M_numLevelAction=1;
    M_qEmptyTask=false;
    if (M_qFirstTask) { M_t_begin = std::chrono::steady_clock::now(); M_qFirstTask=false;}
    for ( ; first!=last; ++first )
    {
        M_idk++; M_idType.push_back(M_numTypeTh);
        auto const& ivdk = *first;
        //std::cout <<ivdk;

        auto args = NA::make_arguments( std::forward<Ts>(ts)... );
        auto && task = args.get(_tasks);
        auto && parameters = args.get_else(_parameters,std::make_tuple());
        if (M_qUseIndex) { std::get<0>(parameters)=std::cref(ivdk); } 
        auto tp=std::tuple_cat( 
					SUBTASK::makeSpData( parameters ), 
					std::make_tuple( task ) 
				);

        SUBTASK::Runtime runtime;

        if (M_numTypeTh==0) {
            std::apply( [&runtime](auto... args){ runtime.task(args...); }, tp );
        }

        if (M_numTypeTh==1) {
            auto LamdaTransfert = [&]() {
			    std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
            return true; 
            };
            std::thread th(LamdaTransfert);
            M_mythreads.push_back(std::move(th));
            usleep(1);
        }

        if (M_numTypeTh==2) {
            auto LamdaTransfert = [&]() {
			    std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
            return true; 
		    };

            if (M_qDeferred) { M_myfutures.emplace_back(std::async(std::launch::deferred,LamdaTransfert));}
            else             { M_myfutures.emplace_back(std::async(std::launch::async,LamdaTransfert)); }
            usleep(1);
        }

        if (M_numTypeTh==3) {
            auto tp=std::tuple_cat( 
					SUBTASK::makeSpDataSpecx( parameters ), 
					std::make_tuple( task ) 
			);
            std::apply([&](auto &&... args) { M_mytg.task(args...).setTaskName("Op("+std::to_string(M_idk)+")"); },tp);
            usleep(0); std::atomic_int counter(0);
        }

        #ifdef COMPILE_WITH_CXX_20
        if (M_numTypeTh==4) {
            auto LamdaTransfert = [&]() {
			    std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
            return true; 
            };
            std::jthread th(LamdaTransfert);
            myjthreads.push_back(std::move(th));
            usleep(1);
        }
        #endif

    }
}


void Task::run()
{
    if (M_qEmptyTask) { std::cout<<"[INFO]: Run failed empty task\n"; exit(0); }
    M_numLevelAction=2;
    if (M_qInfo) { std::cout<<"[INFO]: Run\n"; }
    if (M_numTypeTh==0) { } //No Thread
    if (M_numTypeTh==1) { 
        for (std::thread &t : M_mythreads) { t.join();} 
    } //multithread

    if (M_numTypeTh==2) { 
        for( auto& r : M_myfutures){ auto a =  r.get(); }; 
    } //std::async

    if (M_numTypeTh==3) { 
        /*promise0.set_value(0);*/ 
        M_mytg.waitAllTasks();
    } //Specx

    #ifdef COMPILE_WITH_CXX_20
    if (M_numTypeTh==4) { 
        for (std::jthread &t : myjthreads) { t.join();} 
    } //std::jthread
    #endif

    if (M_numTypeTh==10) { 
        for (int i = 0; i < M_mypthread_cpu.size(); i++) { 
            std::cout<<"[INFO]: Joint "<<i<<"\n";
            pthread_join(M_mypthread_t[i], NULL); 
        }
        //pthread_attr_destroy(&M_mypthread_attr_t);
    } //In CPU

     if (M_numTypeTh==33) { 
        #ifdef COMPILE_WITH_HIP
            M_mytg.waitAllTasks();
        #endif
    } //Specx GPU


    M_mythreads.clear();
    M_myfutures.clear();
    #ifdef COMPILE_WITH_CXX_20
    myjthreads.clear();
    #endif

    M_numTaskStatus.clear();
    M_idType.clear();

    if (M_qInfo) { std::cout<<"[INFO]: All Tasks Accomplished\n"; }
    M_qEmptyTask=true;
}

void Task::close()
{
    if (M_qFlagDetachAlert) {
        if (M_qInfo) { std::cout<<"[INFO]: Detach processes are still running...\n"; }
        if (M_qInfo) { std::cout<<"[INFO]: Please wait before closing.\n"; }

        if ((M_numTypeTh>0) && (M_numTypeTh<10))
        {
            if (M_myfuturesdetach.size()>0)
            {
                for( auto& r : M_myfuturesdetach) {
                    //std::cout <<"Detach Status=" <<r.valid() << '\n';
                    r.wait(); 
                }
            }
        }      
    }

    M_numLevelAction=3;
    if (M_numTypeTh== 0) {  } //No Thread
    if (M_numTypeTh== 1) {  } //multithread
    if (M_numTypeTh== 2) {  } //std::async
    if (M_numTypeTh== 3) { M_myce.stopIfNotAlreadyStopped(); } //Specx
    if (M_numTypeTh==33) { M_myce.stopIfNotAlreadyStopped(); } //Specx GPU <== see if we really need it

    if (!M_qFirstTask) { M_t_end = std::chrono::steady_clock::now(); M_qFirstTask=true;}

    if (M_myfuturesdetach.size()>0) { M_myfuturesdetach.clear(); }

    if (M_qInfo) { std::cout<<"[INFO]: Close All Tasks and process\n"; }
    M_idk=0;
}

void Task::debriefingTasks()
{
    M_t_laps=std::chrono::duration_cast<std::chrono::microseconds>(M_t_end - M_t_begin).count();

    if (M_qViewChrono) {  
        if (M_qInfo) { std::cout << "[INFO]: Elapsed microseconds: "<<M_t_laps<< " us\n"; }
    }

    if (M_qSave)
    {
        if ((M_numTypeTh==3) || (M_numTypeTh==33)) { 
            std::cout << "[INFO]: Save "<<M_FileName<< "\n";
            M_mytg.generateDot(M_FileName+".dot",true); 
            M_mytg.generateTrace(M_FileName+".svg",true); 
        }

        std::ofstream myfile;
        myfile.open (M_FileName+".csv");
        myfile << "Elapsed microseconds,"<< M_t_laps<<"\n";
        myfile << "Nb max Thread,"<< M_nbThTotal<<"\n";
        myfile << "Nb Thread used,"<< M_nbTh<<"\n";
        myfile << "Nb Thread Detach used,"<<M_nbThreadDetach<<"\n";
        myfile << "Mode,"<< M_numTypeTh<<"\n";
        myfile.close();
    }
} 


template<typename FctDetach>
auto Task::add_detach(FctDetach&& func) -> std::future<decltype(func())>
{
    auto task   = std::packaged_task<decltype(func())()>(std::forward<FctDetach>(func));
    auto future = task.get_future();
    std::thread(std::move(task)).detach();
    return std::move(future);
}

template <typename ... Ts>
void Task::add(int numCPU,Ts && ... ts)
{
    //Run in num CPU
    M_mypthread_cpu.push_back(numCPU);
    pthread_t new_thread;
    

    auto args = NA::make_arguments( std::forward<Ts>(ts)... );
    auto && task = args.get(_tasks);
    auto && parameters = args.get_else(_parameters,std::make_tuple());
    SUBTASK::Runtime runtime;
    auto tp=std::tuple_cat( 
            SUBTASK::makeSpData( parameters ), 
                      std::make_tuple( task ) 
    );
    auto LamdaTransfert = [&]() {
                std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
                return true; 
    };


    std::function<void()> func =LamdaTransfert;
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(numCPU, &cpuset);
    pthread_attr_setaffinity_np(&M_mypthread_attr_t, sizeof(cpuset), &cpuset);
    int ret=pthread_create(&new_thread,&M_mypthread_attr_t,workerNumCPU,&func);
    if (ret) { std::cerr << "Error in creating thread" << std::endl; }

    M_mypthread_t.push_back(new_thread);

    //pthread_attr_destroy(&M_mypthread_attr_t);

    //vectorOfThreads.resize(NUM_THREADS);
}




template <typename ... Ts>
void Task::runInCPUs(const std::vector<int> & numCPU,Ts && ... ts)
{
    int M_nbTh=numCPU.size();
    pthread_t thread_array[M_nbTh];
    pthread_attr_t pta_array[M_nbTh];

    auto args = NA::make_arguments( std::forward<Ts>(ts)... );
    auto && task = args.get(_tasks);
    auto && parameters = args.get_else(_parameters,std::make_tuple());
    SUBTASK::Runtime runtime;
    M_qUseIndex=true;
    
    for (int i = 0; i < M_nbTh; i++) {
        int const& M_idk = i;
        if (M_qUseIndex) { std::get<0>(parameters)=M_idk; }
        auto tp=std::tuple_cat( 
                SUBTASK::makeSpData( parameters ), 
                        std::make_tuple( task ) 
        );

        auto LamdaTransfert = [&]() {
                    std::apply([&runtime](auto... args){ runtime.task(args...); }, tp);
                    return true; 
        };
        std::function<void()> func =LamdaTransfert;
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(numCPU[i], &cpuset);
        std::cout<<"Num CPU="<< numCPU[i] <<" activated"<<std::endl;
        pthread_attr_init(&pta_array[i]);
        pthread_attr_setaffinity_np(&pta_array[i], sizeof(cpuset), &cpuset);
        if (pthread_create(&thread_array[i],&pta_array[i],workerNumCPU,&func)) { std::cerr << "Error in creating thread" << std::endl; }
    }

    for (int i = 0; i < M_nbTh; i++) {
            pthread_join(thread_array[i], NULL);
    }

    for (int i = 0; i < M_nbTh; i++) {
            pthread_attr_destroy(&pta_array[i]);
    }
}

} //End namespace TASK
//--------------------------------------------------------------------------------------------------------------------------------


//================================================================================================================================
// THE END.
//================================================================================================================================
