namespace Feel
{



//================================================================================================================================
// Tools to manage memories between CPU and GPU
//================================================================================================================================

#ifdef UseHIP
struct VectorBuffer {
      unsigned int dimension; 
      unsigned int dimensionSizeof; 
 	  unsigned int pitch; 
 	  double* data;

      /////////////////////////////////////////////////////////////

      class DataDescr {
        std::size_t size;
        public:
            explicit DataDescr(const std::size_t inSize = 0) : size(inSize){}

            auto getSize() const{
                return size;
            }
        };

        using DataDescriptor = DataDescr;

        std::size_t memmovNeededSize() const{
            return dimensionSizeof;
        }

        template <class DeviceMemmov>
            auto memmovHostToDevice(DeviceMemmov& mover, void* devicePtr,[[maybe_unused]] std::size_t size){
                assert(size == dimensionSizeof);
                double* doubleDevicePtr = reinterpret_cast<double*>(devicePtr);
                mover.copyHostToDevice(doubleDevicePtr, data, dimensionSizeof);
                return DataDescr(dimension);
            }

        template <class DeviceMemmov>
            void memmovDeviceToHost(DeviceMemmov& mover, void* devicePtr,[[maybe_unused]] std::size_t size, const DataDescr& /*inDataDescr*/){
                assert(size == dimensionSizeof);
                double* doubleDevicePtr = reinterpret_cast<double*>(devicePtr);
                mover.copyDeviceToHost(data, doubleDevicePtr, dimensionSizeof);
            }
};

 #endif


//================================================================================================================================
// 
//================================================================================================================================


#ifdef UseHIP
template<typename ptrtype>
    struct bufferGraphHIP {
        unsigned int size; 
        ptrtype* data;
        ptrtype* deviceBuffer;

    void memoryInit(int dim)
    {
        size=dim; data=(ptrtype *)malloc(sizeof(ptrtype) * size);
    }

    void memmovHostToDevice()
    {
        hipMalloc((void **) &deviceBuffer, sizeof(ptrtype) * size);
        hipMemcpy(deviceBuffer,data,sizeof(ptrtype) * size, hipMemcpyHostToDevice);
    };

    void memmovDeviceToHost()
    {
        hipMemcpy(data,deviceBuffer,sizeof(ptrtype) * size, hipMemcpyDeviceToHost);
        hipFree(deviceBuffer);
    }

};
 #endif

 #ifdef UseCUDA
template<typename ptrtype>
    struct bufferGraphCUDA {
        unsigned int size; 
        ptrtype* data;
        ptrtype* deviceBuffer;

    void memoryInit(int dim)
    {
        size=dim; data=(ptrtype *)malloc(sizeof(ptrtype) * size);
    }

    void memmovHostToDevice()
    {
        cudaMalloc((void **) &deviceBuffer, sizeof(ptrtype) * size);
        cudaMemcpy(deviceBuffer,data,sizeof(ptrtype) * size, cudaMemcpyHostToDevice);
    };

    void memmovDeviceToHost()
    {
        cudaMemcpy(data,deviceBuffer,sizeof(ptrtype) * size, cudaMemcpyDeviceToHost);
        cudaFree(deviceBuffer);
    }
};
 #endif


#ifdef UseHIP
template<typename ptrtype>
    struct bufferGraphUnifiedHIP {
        unsigned int size; 
        ptrtype* data;

    void memoryInit(int dim) {  size=dim; hipMallocManaged(&data,sizeof(ptrtype) * size); }
    void memoryInit(int dim,ptrtype v) { 
      size=dim; hipMallocManaged(&data,sizeof(ptrtype) * size); 
      for (long int i = 0; i < dim; i++) { data[i]=v; }
    }
};
#endif

#ifdef UseCUDA
template<typename ptrtype>
    struct bufferGraphUnifiedCUDA {
        unsigned int size; 
        ptrtype* data;

    void memoryInit(int dim) { size=dim; cudaMallocManaged(&data,sizeof(ptrtype) * size); }
    void memoryInit(int dim,ptrtype v) { 
      size=dim; cudaMallocManaged(&data,sizeof(ptrtype) * size); 
      for (long int i = 0; i < dim; i++) { data[i]=v; }
    }
};
#endif

 
//================================================================================================================================
// GPU KERNEL FUNCTIONS
//================================================================================================================================


 #if defined(UseHIP) || defined(UseCUDA)
    template<typename Kernel, typename Input, typename Output>
    __global__ void OP_IN_KERNEL_GPU_1D(const Kernel kernel_function, int n,Input* in, Output* out)
    {
        int i = threadIdx.x + blockIdx.x*blockDim.x;
        if(i<n) {
            kernel_function(i, in, out);
        }
    }

    template<typename Kernel,typename Input>
    __global__ void OP_IN_KERNEL_LAMBDA_GPU_1D(Kernel op,Input *A, int nb) 
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x; 
        if (idx < nb)
            op(idx,A);
        //__syncthreads();
    }

  
     template<typename Kernel,typename Input>
    __global__ void OP_IN_KERNEL_LAMBDA_GPU_1D_3I(Kernel op,Input *R,Input *A,Input *B, int nb) 
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x; 
        if (idx < nb)
            op(idx,R,A,B);
        //__syncthreads();
    }
 
    template<typename Kernel,typename Input>
    __global__ void OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D(Kernel op,Input *A,int iBegin,int iEnd) 
    {
        int idx = blockIdx.x * blockDim.x + threadIdx.x; 
        if ((idx>=iBegin) && (idx<iEnd))
            op(idx,A);
        //__syncthreads();
    }


    template<typename Kernel,typename Input>
    __global__ void OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D(Kernel op,Input *A,int offset) 
    {
        int idx =  offset + blockIdx.x * blockDim.x + threadIdx.x; 
        op(idx,A);
        //__syncthreads();
    }

    template<typename Kernel,typename Input>
    __global__ void OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I(Kernel op,Input *R,Input *A,Input *B, int offset) 
    {
        int idx =  offset + blockIdx.x * blockDim.x + threadIdx.x; 
        op(idx,R,A,B);
    }
#endif




//================================================================================================================================
// CLASS Task: Provide to manipulate graph Hip/Cuda hybrid system Explicit and Implicit
//================================================================================================================================

// GRAPH EXPLICIT

namespace Taskgpu {

class SingleTask
{
    private:
        std::string M_FileName;
        bool        M_qSave;            
        int         M_numType;  // 1:HIP 2:CUDA
        int         M_nbStreams; 
        bool        M_qViewInfo; 
        int         M_block_size; 
        int         M_numModel; // 1:Serial 2:Stream     
        bool        M_qUsedUnifiedMemory;  
        int         M_solveLevel; //will be used in the future
        long int    M_time_laps;        
       
        std::chrono::steady_clock::time_point M_t_begin,M_t_end;


        template<typename Kernel, typename Input>
            void serial_hip(const Kernel& kernel_function,
                                        int numElems,Input* buffer);

        template<typename Kernel, typename Input>
            void serial_cuda(const Kernel& kernel_function,
                                        int numElems,Input* buffer);

        template<typename Kernel, typename Input>
            void stream_hip(const Kernel& kernel_function,
                                        int numElems,Input* buffer);

        template<typename Kernel, typename Input>
            void stream_cuda(const Kernel& kernel_function,
                                        int numElems,Input* buffer);


        template<typename Kernel, typename Input>
            void serial_hip_3I(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB);

        template<typename Kernel, typename Input>
            void serial_cuda_3I(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB);

        template<typename Kernel, typename Input>
            void stream_hip_3I(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB);

        template<typename Kernel, typename Input>
            void stream_cuda_3I(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB);


    public:
        void setNumType       (int v)         {  M_numType      = v; }
        void setNumModel      (int v)         {  M_numModel     = v; }
        void setNbStreams     (int v)         {  M_nbStreams    = v; }
        void setViewInfo      (bool b)        {  M_qViewInfo    = b; }
        void setFileName      (std::string s) {  M_FileName=s;       }
        void setSave          (bool b)        {  M_qSave        = b; }
        void setNbBlock       (int v)         {  M_block_size   = v; }
        void setDevice        (int v);
        void setUnifiedMemory (bool b)        {  M_qUsedUnifiedMemory = b; }
        void setSolveLevel    (int v)         {  M_solveLevel   = v; }
    
        int  getNumType       () const        {  return (M_numType);   }
        int  getNumModel      () const        {  return (M_numModel);  }
        int  getNbStreams     () const        {  return (M_nbStreams); }
        int  getSolveLevel    () const        {  return (M_solveLevel); }
        bool isSave           () const        {  return (M_qSave);     }
        bool isUnifiedMemory  () const        {  return (M_qUsedUnifiedMemory); }

        long int  getTimeLaps     () const        {  return (M_time_laps);   }

        SingleTask();
        ~SingleTask();

        template<typename Kernel, typename Input>
            void run(const Kernel& kernel_function,
                                int numElems,
                                Input* buffer);   

        template<typename Kernel, typename Input>
            void run(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB);               
        void debriefing();
        void close();
};


SingleTask::SingleTask()
{
    M_nbStreams          = 3;
    M_numType            = 1; // 1:HIP 2:CUDA
    M_qViewInfo          = true;
    M_time_laps          = 0;
    M_FileName           = "NoName";
    M_qSave              = false;
    M_block_size         = 512;
    M_numModel           = 1; // 1:Serial 2:Stream
    M_qUsedUnifiedMemory = false;
    M_solveLevel         = 1; 
}

SingleTask::~SingleTask()
{
}

void SingleTask::debriefing()
{
    if (M_qViewInfo) {
        std::cout<<"<=====================================================================>"<<"\n";
        std::cout<<"[INFO]: Debriefing"<<"\n";
        std::cout<<"[INFO]: Elapsed microseconds : "<<M_time_laps<< " us\n";
    }
    if (M_qSave)
    {
        if (M_qViewInfo) { std::cout<<"[INFO]: Save Informations"<<"\n"; }
        std::ofstream myfile;
        myfile.open (M_FileName+".csv");
        myfile << "Elapsed microseconds,"<<M_time_laps<<"\n";
        myfile <<"\n";
        myfile.close();
    }
}

void SingleTask::close()
{

}

void SingleTask::setDevice(int v)
{
    if (M_numType==1) { 
        #ifdef UseHIP  
            int numDevices=0; hipGetDeviceCount(&numDevices);
            if ((v>=0) && (v<=numDevices)) { hipSetDevice(v); }
        #endif
    }

    if (M_numType==2) { 
        #ifdef UseCUDA  
            int numDevices=0; cudaGetDeviceCount(&numDevices);
        if ((v>=0) && (v<=numDevices)) { cudaSetDevice(v); }
        #endif
    }
}


template<typename Kernel, typename Input>
    void SingleTask::run(const Kernel& kernel_function,
                                int numElems,
                                Input* buffer)
{   
    M_t_begin = std::chrono::steady_clock::now();
    if (M_numModel==1) { //SERIAL MODEL
        if (M_numType==1) { 
            #ifdef UseHIP  
                serial_hip(kernel_function,numElems,buffer);
            #endif
        }

        if (M_numType==2) { 
            #ifdef UseCUDA  
                serial_cuda(kernel_function,numElems,buffer);
            #endif
        }
    }

    if (M_numModel==2) {  //STREAM MODEL
        if (M_numType==1) { 
            #ifdef UseHIP  
                stream_hip(kernel_function,numElems,buffer);
            #endif
        }

        if (M_numType==2) { 
            #ifdef UseCUDA  
                stream_cuda(kernel_function,numElems,buffer);
            #endif
        }
    }
    M_t_end = std::chrono::steady_clock::now();
    M_time_laps= std::chrono::duration_cast<std::chrono::microseconds>(M_t_end - M_t_begin).count();
}


template<typename Kernel, typename Input>
    void SingleTask::serial_hip(const Kernel& kernel_function,
                                int numElems,
                                Input* buffer)
{   
    #ifdef UseHIP  
        int sz=sizeof(decltype(buffer)) * numElems;
        int iEnd=numElems; 
        int iBegin=0; 
        int num_blocks = (numElems + M_block_size - 1) / M_block_size;
        dim3 thread_block(M_block_size, 1, 1);
        dim3 grid(num_blocks, 1);
        if (M_qUsedUnifiedMemory)
        {
            hipLaunchKernelGGL(OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D,grid,thread_block,0,0,kernel_function,buffer,iBegin,iEnd);
            hipDeviceSynchronize();
        }
        else 
        {
            Input *deviceBuffer;
            hipMalloc((void **) &deviceBuffer,sz);
            hipMemcpy(deviceBuffer,buffer,sz, hipMemcpyHostToDevice);
            hipLaunchKernelGGL(OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D,grid,thread_block,0,0,kernel_function,deviceBuffer,iBegin,iEnd);
            hipMemcpy(buffer,deviceBuffer,sz, hipMemcpyDeviceToHost);
            hipFree(deviceBuffer);
        }
    #endif
}

template<typename Kernel, typename Input>
    void SingleTask::serial_cuda(const Kernel& kernel_function,
                                int numElems,
                                Input* buffer)
{   
    #ifdef UseCUDA 
        int sz=sizeof(decltype(buffer)) * numElems;
        int iEnd=numElems; 
        int iBegin=0; 
        int num_blocks = (numElems + M_block_size - 1) / M_block_size;
        dim3 thread_block(M_block_size, 1, 1);
        dim3 grid(num_blocks, 1);   
        if (M_qUsedUnifiedMemory)
        {
            OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<<<gridthread_block>>>(kernel_function,buffer,iBegin,iEnd);
            cudaDeviceSynchronize();
        }     
        else
        {
            Input *deviceBuffer;
            cudaMalloc((void **) &deviceBuffer,sz);
            cudaMemcpy(deviceBuffer,buffer,sz, hipMemcpyHostToDevice);
            OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<<<gridthread_block>>>(kernel_function,deviceBuffer,iBegin,iEnd);
            cudaMemcpy(buffer,deviceBuffer,sz, hipMemcpyDeviceToHost);
            cudaFree(deviceBuffer);
        }
    #endif
}




template<typename Kernel, typename Input>
    void SingleTask::run(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB)
{   
    M_t_begin = std::chrono::steady_clock::now();
    if (M_numModel==1) { //SERIAL MODEL
        if (M_numType==1) { 
            #ifdef UseHIP  
                serial_hip_3I(kernel_function,numElems,bufferR,bufferA,bufferB);
            #endif
        }

        if (M_numType==2) { 
            #ifdef UseCUDA  
                serial_cuda_3I(kernel_function,numElems,bufferR,bufferA,bufferB);
            #endif
        }
    }

    if (M_numModel==2) {  //STREAM MODEL
        if (M_numType==1) { 
            #ifdef UseHIP  
                stream_hip_3I(kernel_function,numElems,bufferR,bufferA,bufferB);
            #endif
        }

        if (M_numType==2) { 
            #ifdef UseCUDA  
                stream_cuda_3I(kernel_function,numElems,bufferR,bufferA,bufferB);
            #endif
        }
    }
    M_t_end = std::chrono::steady_clock::now();
    M_time_laps= std::chrono::duration_cast<std::chrono::microseconds>(M_t_end - M_t_begin).count();
}


template<typename Kernel, typename Input>
    void SingleTask::serial_hip_3I(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB)
{   
    #ifdef UseHIP  
        
        int sz=sizeof(decltype(bufferR)) * numElems;
        int iEnd=numElems; 
        int iBegin=0; 
        int num_blocks = (numElems + M_block_size - 1) / M_block_size;
        dim3 thread_block(M_block_size, 1, 1);
        dim3 grid(num_blocks, 1);
        if (M_qUsedUnifiedMemory)
        {
            hipLaunchKernelGGL(OP_IN_KERNEL_LAMBDA_GPU_1D_3I,grid,thread_block,0,0,kernel_function,
                bufferR,bufferA,bufferB,iBegin,iEnd);
            hipDeviceSynchronize();
        }
        else 
        {
            Input *deviceBufferR,*deviceBufferA,*deviceBufferB;
            hipMalloc((void **) &deviceBufferR,sz);
            hipMalloc((void **) &deviceBufferA,sz);
            hipMalloc((void **) &deviceBufferB,sz);
            hipMemcpy(deviceBufferR,bufferR,sz, hipMemcpyHostToDevice);
            hipMemcpy(deviceBufferA,bufferA,sz, hipMemcpyHostToDevice);
            hipMemcpy(deviceBufferB,bufferB,sz, hipMemcpyHostToDevice);
            hipLaunchKernelGGL(OP_IN_KERNEL_LAMBDA_GPU_1D_3I,grid,thread_block,0,0,kernel_function,
                deviceBufferR,deviceBufferA,deviceBufferB,iBegin,iEnd);
            hipMemcpy(bufferR,deviceBufferR,sz, hipMemcpyDeviceToHost);
            hipFree(deviceBufferR);
            hipFree(deviceBufferA);
            hipFree(deviceBufferB);
        }
    #endif
}


template<typename Kernel, typename Input>
    void SingleTask::serial_cuda_3I(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB)
{   
    #ifdef UseCUDA  
        
        int sz=sizeof(decltype(bufferR)) * numElems;
        int iEnd=numElems; 
        int iBegin=0; 
        int num_blocks = (numElems + M_block_size - 1) / M_block_size;
        dim3 thread_block(M_block_size, 1, 1);
        dim3 grid(num_blocks, 1);
        if (M_qUsedUnifiedMemory)
        {
            OP_IN_KERNEL_LAMBDA_GPU_1D_3I<<<gridthread_block>>>(kernel_function,
                bufferR,bufferA,bufferB,iBegin,iEnd);
            cudaDeviceSynchronize();
        }
        else 
        {
            Input *deviceBufferR,*deviceBufferA,*deviceBufferB;
            cudaMalloc((void **) &deviceBufferR,sz);
            cudaMalloc((void **) &deviceBufferA,sz);
            cudaMalloc((void **) &deviceBufferB,sz);
            cudaMemcpy(deviceBufferR,bufferR,sz, cudaMemcpyHostToDevice);
            cudaMemcpy(deviceBufferA,bufferA,sz, cudaMemcpyHostToDevice);
            cudaMemcpy(deviceBufferB,bufferB,sz, cudaMemcpyHostToDevice);
            OP_IN_KERNEL_LAMBDA_GPU_1D_3I<<<gridthread_block>>>(kernel_function,
                deviceBufferR,deviceBufferA,deviceBufferB,iBegin,iEnd);
            cudaMemcpy(bufferR,deviceBufferR,sz, cudaMemcpyDeviceToHost);
            cudaFree(deviceBufferR);
            cudaFree(deviceBufferA);
            cudaFree(deviceBufferB);
        }
    #endif
}


template<typename Kernel, typename Input>
    void SingleTask::stream_hip(const Kernel& kernel_function,
                                int numElems,
                                Input* buffer)
{   
    #ifdef UseHIP  
        int   numVersion=1;
        float ms;
        int   i;

        //const int blockSize = 20;  //const int blockSize M_block_size;
        const int blockSize = M_block_size;
        
        const int streamSize = numElems / M_nbStreams;
        const int streamBytes = streamSize * sizeof(decltype(buffer));
        const int sz = numElems * sizeof(decltype(buffer));
        int grid=streamSize/blockSize;
        if (M_qViewInfo) {  
            std::cout<<"[INFO]: Block size :"<<M_block_size<<"\n"; 
            std::cout<<"[INFO]: nb Streams :"<<M_nbStreams<<"\n"; 
        }
        if ((grid==0) && (M_qViewInfo)) {  std::cout<<"[INFO]: Grid size error : ==> auto correction"<<"\n"; }
        grid=max(streamSize/blockSize,1);
        Input *d_buffer; hipMalloc((void **) &d_buffer,sz);
        hipEvent_t startEvent, stopEvent, dummyEvent;
        hipStream_t stream[M_nbStreams];
        hipEventCreate(&startEvent) ;
        hipEventCreate(&stopEvent) ;
        hipEventCreate(&dummyEvent) ;
        for (i = 0; i < M_nbStreams; ++i) { hipStreamCreate(&stream[i]) ; }
        hipEventRecord(startEvent,0) ;

        if (numVersion==1) {
            for (i = 0; i < M_nbStreams; ++i) {
                int offset = i * streamSize;
                hipMemcpyAsync(&d_buffer[offset], &buffer[offset],streamBytes, hipMemcpyHostToDevice,stream[i]) ;
                hipLaunchKernelGGL(OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D,grid,blockSize,0,stream[i],kernel_function,d_buffer,offset);
                hipMemcpyAsync(&buffer[offset], &d_buffer[offset],streamBytes, hipMemcpyDeviceToHost,stream[i]);
            }
            hipEventRecord(stopEvent, 0);
            hipEventSynchronize(stopEvent);
            //hipEventElapsedTime(&ms, startEvent, stopEvent);
            //printf("Time (ms): %f\n", ms);
        }

        if (numVersion==2) {    
            hipEventRecord(startEvent,0);
            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                hipMemcpyAsync(&d_buffer[offset], &buffer[offset], streamBytes, hipMemcpyHostToDevice,stream[i]) ;
            }

            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                hipLaunchKernelGGL(OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D,grid,blockSize,0,stream[i],kernel_function,d_buffer,offset);
            }

            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                hipMemcpyAsync(&buffer[offset], &d_buffer[offset],streamBytes, hipMemcpyDeviceToHost,stream[i]) ;
            }
            hipEventRecord(stopEvent,0);
            hipEventSynchronize(stopEvent);
            //hipEventElapsedTime(&ms, startEvent, stopEvent);
            //printf("Time (ms): %f\n", ms);
        }

       // cleanup
        hipEventDestroy(startEvent) ;
        hipEventDestroy(stopEvent) ;
        hipEventDestroy(dummyEvent);
        for (int i = 0; i < M_nbStreams; ++i) { hipStreamDestroy(stream[i]); }
        hipFree(d_buffer);
    #endif
}


template<typename Kernel, typename Input>
    void SingleTask::stream_cuda(const Kernel& kernel_function,
                                int numElems,
                                Input* buffer)
{   
    #ifdef UseCUDA  
        int   numVersion=1;
        float ms;
        int   i;

        //const int blockSize = 20;  //const int blockSize M_block_size;
        const int blockSize = M_block_size;
        
        const int streamSize = numElems / M_nbStreams;
        const int streamBytes = streamSize * sizeof(decltype(buffer));
        const int sz = numElems * sizeof(decltype(buffer));
        int grid=streamSize/blockSize;
        if (M_qViewInfo) {  
            std::cout<<"[INFO]: Block size :"<<M_block_size<<"\n"; 
            std::cout<<"[INFO]: nb Streams :"<<M_nbStreams<<"\n"; 
        }
        if ((grid==0) && (M_qViewInfo)) {  std::cout<<"[INFO]: Grid size error : ==> auto correction"<<"\n"; }
        grid=max(streamSize/blockSize,1);
        Input *d_buffer; cudaMalloc((void **) &d_buffer,sz);
        cudaEvent_t startEvent, stopEvent, dummyEvent;
        cudaStream_t stream[M_nbStreams];
        cudaEventCreate(&startEvent) ;
        cudaEventCreate(&stopEvent) ;
        cudaEventCreate(&dummyEvent) ;
        for (i = 0; i < M_nbStreams; ++i) { cudaStreamCreate(&stream[i]) ; }
        cudaEventRecord(startEvent,0) ;

        if (numVersion==1) {
            for (i = 0; i < M_nbStreams; ++i) {
                int offset = i * streamSize;
                cudaMemcpyAsync(&d_buffer[offset], &buffer[offset],streamBytes, cudaMemcpyHostToDevice,stream[i]) ;
                OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D<<<grid,blockSize,0,stream[i]>>>(kernel_function,d_buffer,offset);
                cudaMemcpyAsync(&buffer[offset], &d_buffer[offset],streamBytes, cudaMemcpyDeviceToHost,stream[i]);
            }
            cudaEventRecord(stopEvent, 0);
            cudaEventSynchronize(stopEvent);
            //cudaEventElapsedTime(&ms, startEvent, stopEvent);
            //printf("Time (ms): %f\n", ms);
        }

        if (numVersion==2) {    
            cudaEventRecord(startEvent,0);
            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                cudaMemcpyAsync(&d_buffer[offset], &buffer[offset], streamBytes, cudaMemcpyHostToDevice,stream[i]) ;
            }

            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D<<<grid,blockSize,0,stream[i]>>>(kernel_function,d_buffer,offset);
            }

            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                cudaMemcpyAsync(&buffer[offset], &d_buffer[offset],streamBytes, cudaMemcpyDeviceToHost,stream[i]) ;
            }
            cudaEventRecord(stopEvent,0);
            cudaEventSynchronize(stopEvent);
            //cudaEventElapsedTime(&ms, startEvent, stopEvent);
            //printf("Time (ms): %f\n", ms);
        }

       // cleanup
        cudaEventDestroy(startEvent) ;
        cudaEventDestroy(stopEvent) ;
        cudaEventDestroy(dummyEvent);
        for (int i = 0; i < M_nbStreams; ++i) { cudaStreamDestroy(stream[i]); }
        cudaFree(d_buffer);
       
    #endif
}



template<typename Kernel, typename Input>
    void SingleTask::stream_hip_3I(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB)
{   
    #ifdef UseHIP  
        int   numVersion=1;
        float ms;
        int   i;

        //const int blockSize = 20;  //const int blockSize M_block_size;
        const int blockSize = M_block_size;
        
        const int streamSize = numElems / M_nbStreams;
        const int streamBytes = streamSize * sizeof(decltype(bufferR));
        const int sz = numElems * sizeof(decltype(bufferR));
        int grid=streamSize/blockSize;
        if (M_qViewInfo) {  
            std::cout<<"[INFO]: Block size :"<<M_block_size<<"\n"; 
            std::cout<<"[INFO]: nb Streams :"<<M_nbStreams<<"\n"; 
        }
        if ((grid==0) && (M_qViewInfo)) {  std::cout<<"[INFO]: Grid size error : ==> auto correction"<<"\n"; }
        grid=max(streamSize/blockSize,1);
        Input *d_bufferR; hipMalloc((void **) &d_bufferR,sz);
        Input *d_bufferA; hipMalloc((void **) &d_bufferA,sz);
        Input *d_bufferB; hipMalloc((void **) &d_bufferB,sz);
        hipEvent_t startEvent, stopEvent, dummyEvent;
        hipStream_t stream[M_nbStreams];
        hipEventCreate(&startEvent) ;
        hipEventCreate(&stopEvent) ;
        hipEventCreate(&dummyEvent) ;
        for (i = 0; i < M_nbStreams; ++i) { hipStreamCreate(&stream[i]) ; }
        hipEventRecord(startEvent,0) ;

        if (numVersion==1) {
            for (i = 0; i < M_nbStreams; ++i) {
                int offset = i * streamSize;
                hipMemcpyAsync(&d_bufferR[offset], &bufferR[offset],streamBytes, hipMemcpyHostToDevice,stream[i]) ;
                hipMemcpyAsync(&d_bufferA[offset], &bufferA[offset],streamBytes, hipMemcpyHostToDevice,stream[i]) ;
                hipMemcpyAsync(&d_bufferB[offset], &bufferB[offset],streamBytes, hipMemcpyHostToDevice,stream[i]) ;
                hipLaunchKernelGGL(OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I,grid,blockSize,0,stream[i],kernel_function,
                    d_bufferR,d_bufferA,d_bufferB,offset);
                hipMemcpyAsync(&bufferR[offset], &d_bufferR[offset],streamBytes, hipMemcpyDeviceToHost,stream[i]);
            }
            hipEventRecord(stopEvent, 0);
            hipEventSynchronize(stopEvent);
        }

         if (numVersion==2) {    
            hipEventRecord(startEvent,0);
            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                hipMemcpyAsync(&d_bufferR[offset], &bufferR[offset], streamBytes, hipMemcpyHostToDevice,stream[i]) ;
                hipMemcpyAsync(&d_bufferA[offset], &bufferA[offset], streamBytes, hipMemcpyHostToDevice,stream[i]) ;
                hipMemcpyAsync(&d_bufferB[offset], &bufferB[offset], streamBytes, hipMemcpyHostToDevice,stream[i]) ;
            }

            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                hipLaunchKernelGGL(OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I,grid,blockSize,0,stream[i],kernel_function,
                    d_bufferR,d_bufferA,d_bufferB,offset);
            }

            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                hipMemcpyAsync(&bufferR[offset], &d_bufferR[offset],streamBytes, hipMemcpyDeviceToHost,stream[i]) ;
            }
            hipEventRecord(stopEvent,0);
            hipEventSynchronize(stopEvent);
            //hipEventElapsedTime(&ms, startEvent, stopEvent);
            //printf("Time (ms): %f\n", ms);
        }
        // cleanup
        hipEventDestroy(startEvent) ;
        hipEventDestroy(stopEvent) ;
        hipEventDestroy(dummyEvent);
        for (int i = 0; i < M_nbStreams; ++i) { hipStreamDestroy(stream[i]); }
        hipFree(d_bufferR); hipFree(d_bufferA); hipFree(d_bufferB);
    #endif
}


template<typename Kernel, typename Input>
    void SingleTask::stream_cuda_3I(const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR,Input* bufferA,Input* bufferB)
{   
    #ifdef UseCUDA  
        int   numVersion=1;
        float ms;
        int   i;

        //const int blockSize = 20;  //const int blockSize M_block_size;
        const int blockSize = M_block_size;
        
        const int streamSize = numElems / M_nbStreams;
        const int streamBytes = streamSize * sizeof(decltype(bufferR));
        const int sz = numElems * sizeof(decltype(bufferR));
        int grid=streamSize/blockSize;
        if (M_qViewInfo) {  
            std::cout<<"[INFO]: Block size :"<<M_block_size<<"\n"; 
            std::cout<<"[INFO]: nb Streams :"<<M_nbStreams<<"\n"; 
        }
        if ((grid==0) && (M_qViewInfo)) {  std::cout<<"[INFO]: Grid size error : ==> auto correction"<<"\n"; }
        grid=max(streamSize/blockSize,1);
        Input *d_bufferR; cudaMalloc((void **) &d_bufferR,sz);
        Input *d_bufferA; cudaMalloc((void **) &d_bufferA,sz);
        Input *d_bufferB; cudaMalloc((void **) &d_bufferB,sz);
        cudaEvent_t startEvent, stopEvent, dummyEvent;
        cudaStream_t stream[M_nbStreams];
        cudaEventCreate(&startEvent) ;
        cudaEventCreate(&stopEvent) ;
        cudaEventCreate(&dummyEvent) ;
        for (i = 0; i < M_nbStreams; ++i) { cudaStreamCreate(&stream[i]) ; }
        cudaEventRecord(startEvent,0) ;

        if (numVersion==1) {
            for (i = 0; i < M_nbStreams; ++i) {
                int offset = i * streamSize;
                cudaMemcpyAsync(&d_bufferR[offset], &bufferR[offset],streamBytes, cudaMemcpyHostToDevice,stream[i]) ;
                cudaMemcpyAsync(&d_bufferA[offset], &bufferA[offset],streamBytes, cudaMemcpyHostToDevice,stream[i]) ;
                cudaMemcpyAsync(&d_bufferB[offset], &bufferB[offset],streamBytes, cudaMemcpyHostToDevice,stream[i]) ;
                OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I<<<grid,blockSize,0,stream[i]>>>(kernel_function,
                    d_bufferR,d_bufferA,d_bufferB,offset);
                cudaMemcpyAsync(&bufferR[offset], &d_bufferR[offset],streamBytes, cudaMemcpyDeviceToHost,stream[i]);
            }
            cudaEventRecord(stopEvent, 0);
            cudaEventSynchronize(stopEvent);
        }

         if (numVersion==2) {    
            cudaEventRecord(startEvent,0);
            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                cudaMemcpyAsync(&d_bufferR[offset], &bufferR[offset], streamBytes, cudaMemcpyHostToDevice,stream[i]) ;
                cudaMemcpyAsync(&d_bufferA[offset], &bufferA[offset], streamBytes, cudaMemcpyHostToDevice,stream[i]) ;
                cudaMemcpyAsync(&d_bufferB[offset], &bufferB[offset], streamBytes, cudaMemcpyHostToDevice,stream[i]) ;
            }

            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I<<<grid,blockSize,0,stream[i]>>>(kernel_function,
                    d_bufferR,d_bufferA,d_bufferB,offset);
            }

            for (i = 0; i < M_nbStreams; ++i)
            {
                int offset = i * streamSize;
                cudaMemcpyAsync(&bufferR[offset], &d_bufferR[offset],streamBytes, cudaMemcpyDeviceToHost,stream[i]) ;
            }
            cudaEventRecord(stopEvent,0);
            cudaEventSynchronize(stopEvent);
            //cudaEventElapsedTime(&ms, startEvent, stopEvent);
            //printf("Time (ms): %f\n", ms);
        }
        // cleanup
        cudaEventDestroy(startEvent) ;
        cudaEventDestroy(stopEvent) ;
        cudaEventDestroy(dummyEvent);
        for (int i = 0; i < M_nbStreams; ++i) { cudaStreamDestroy(stream[i]); }
        cudaFree(d_bufferR); cudaFree(d_bufferA); cudaFree(d_bufferB);
    #endif
}



//--------------------------------------------------------------------------------------------------------------------------------

class Task
{
    private:
        std::string M_FileName;
        int         M_nbTh;
        int         M_numBlocksGPU;
        int         M_nThPerBckGPU;
        bool        M_q_graph;
        bool        M_qViewInfo;
        bool        M_qSave;
        long int    M_time_laps;
        bool        M_qDeviceReset;
        int         M_numType;

        std::vector<int>  M_ListGraphDependencies;

        
        std::chrono::steady_clock::time_point M_t_begin,M_t_end;

        //#ifdef COMPILE_WITH_HIP && UseHIP
        //#if defined(COMPILE_WITH_HIP) && defined(UseHIP)
        #ifdef UseHIP
            hipGraph_t                          hip_graph;
            hipGraphExec_t                      hip_graphExec;
            hipStream_t                         hip_graphStream;
            hipKernelNodeParams                 hip_nodeParams;
            std::vector<hipGraphNode_t>         M_hipGraphNode_t;              
        #endif

        #ifdef UseCUDA
            cudaGraph_t                         cuda_graph;
            cudaGraphExec_t                     cuda_graphExec;
            cudaStream_t                        cuda_graphStream;
            cudaKernelNodeParams                cuda_nodeParams;
            std::vector<cudaGraphNode_t>        M_cudaGraphNode_t;   
        #endif


    public:
        
        Task();
        ~Task();

        void setSave         (bool b)        {  M_qSave        = b; }
        void setViewInfo     (bool b)        {  M_qViewInfo    = b; }
        void setDeviceReset  (bool b)        {  M_qDeviceReset = b; }
        void setNumType      (int v)         {  M_numType      = v; }
    
        void setDeviceHIP    (int v);
        void setDeviceCUDA   (int v);
        void setFileName     (std::string s) {  M_FileName=s;           }
        int  getNumType      () const        {  return(M_numType);      }
        bool isSave          () const        {  return(M_qSave);        }
        bool isDeviceReset   () const        {  return(M_qDeviceReset); }

        void open(int nbBlock,int NbTh);

        template<typename Kernel, typename Input, typename Output>
            void add_hip(const Kernel& kernel_function,
                     int numElems,
                     int iBegin,int iEnd,
                     Input* buffer,
                     Output* hostbuffer,
                     std::vector<int> links);

        template<typename Kernel, typename Input, typename Output>
            void add_cuda(const Kernel& kernel_function,
                     int numElems,
                     int iBegin,int iEnd,
                     Input* buffer,
                     Output* hostbuffer,
                     std::vector<int> links);
        void run();
        void close();
        void debriefing(); 
};

Task::Task()
{
    M_numType         = 1;  // 1: HIP  2: CUDA
    M_nbTh            = 1;
    M_numBlocksGPU    = 1;
    M_q_graph         = false;
    M_qViewInfo       = true;
    M_time_laps       = 0;
    M_FileName        = "NoName";
    M_qSave           = true;
    M_qDeviceReset    = false;
    M_ListGraphDependencies.clear();
}

Task::~Task()
{
    M_ListGraphDependencies.clear();
}

void Task::debriefing()
{
    if (M_qViewInfo) { 
        std::cout<<"<=====================================================================>"<<"\n"; 
        std::cout<<"[INFO]: Debriefing"<<"\n";  
        std::cout<<"[INFO]: Elapsed microseconds : "<<M_time_laps<< " us\n";
        std::cout<<"[INFO]: List Graph Dependencie  >>> [";
        for (int i = 0; i < M_ListGraphDependencies.size(); i++) { std::cout<<M_ListGraphDependencies[i];}
        std::cout<<"] <<<\n";
    }

    if (M_qSave)
    {
        if (M_qViewInfo) { std::cout<<"[INFO]: Save Informations"<<"\n"; }   
        std::ofstream myfile;
        myfile.open (M_FileName+".csv");
        myfile << "Elapsed microseconds,"<<M_time_laps<<"\n";
        myfile << "Nb Thread,"<< M_nbTh<<"\n";
        myfile << "Nb Block used,"<<M_numBlocksGPU<<"\n";
        myfile << "Nb Th/Block,"<<M_nThPerBckGPU<<"\n";
        myfile << "List Graph Dependencie,";
        for (int i = 0; i < M_ListGraphDependencies.size(); i++) { myfile <<M_ListGraphDependencies[i];}
        myfile <<"\n";
        myfile.close();
   }
   if (M_qViewInfo) { std::cout<<"<=====================================================================>"<<"\n"; }
}

void Task::setDeviceHIP(int v) 
{
    #ifdef UseHIP
    int numDevices=0; hipGetDeviceCount(&numDevices);
    if ((v>=0) && (v<=numDevices)) { hipSetDevice(v); }
    #endif
}

void Task::setDeviceCUDA(int v) 
{
    #ifdef UseCUDA
    int numDevices=0; cudaGetDeviceCount(&numDevices);
    if ((v>=0) && (v<=numDevices)) { cudaSetDevice(v); }
    #endif
}

void Task::open(int nbBlock,int NbTh)
{
    M_q_graph=false;
    M_nbTh= NbTh;
    M_numBlocksGPU=nbBlock;
    M_nThPerBckGPU=M_nbTh/M_numBlocksGPU;
    if (M_qViewInfo)
    {
        std::cout<<"<=====================================================================>"<<"\n";
        std::cout<<"[INFO]: Open Graph"<<"\n";
        std::cout<<"[INFO]: nb Thread        : "<<M_nbTh<<"\n";
        std::cout<<"[INFO]: nb Block         : "<<M_numBlocksGPU<<"\n";
        std::cout<<"[INFO]: Thread Per Block : "<<M_nThPerBckGPU<<"\n";
        std::cout<<"<=====================================================================>"<<"\n";
    }
    
    //#ifdef COMPILE_WITH_HIP && UseHIP
    //#if defined(COMPILE_WITH_HIP) && defined(UseHIP)

    #ifdef UseHIP       
        hipGraphCreate(&hip_graph, 0);
        hip_nodeParams = {0};
        memset(&hip_nodeParams, 0, sizeof(hip_nodeParams));
    #endif

    #ifdef UseCUDA       
        cudaGraphCreate(&cuda_graph, 0);
        cuda_nodeParams = {0};
        memset(&cuda_nodeParams, 0, sizeof(cuda_nodeParams));
    #endif
    
}


template<typename Kernel, typename Input, typename Output>
    void Task::add_hip(const Kernel& kernel_function,
                                int numElems,
                                int iBegin,int iEnd,
                                Input* buffer,
                                Output* hostbuffer,
                                std::vector<int> links)
{            
    #ifdef UseHIP  
        bool qFlag=false;
        //BEGIN::Init new node
        hipGraphNode_t newKernelNode; M_hipGraphNode_t.push_back(newKernelNode);
        memset(&hip_nodeParams, 0, sizeof(hip_nodeParams));

        if (M_qViewInfo) { std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<"\n"; }
        //if (M_qViewInfo) { std::cout<<"[INFO]: M_hipGraphNode_t="<<M_hipGraphNode_t.size()<<" : "<<M_hipGraphNode_t[M_hipGraphNode_t.size()-1]<<"\n"; }
        if (M_qViewInfo) { std::cout<<"[INFO]: Num Graph Node = "<<M_hipGraphNode_t.size()<<"\n"; }

        //CRTL range
        if (iEnd>numElems) { iEnd=numElems; }
        if (iBegin<0)      { iBegin=0; }

        //printf("[%x]\n",M_hipGraphNode_t[M_hipGraphNode_t.size()-1]);

        hip_nodeParams.func   =  (void *)OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<Kernel,Input>;
        hip_nodeParams.gridDim        = dim3(M_numBlocksGPU, 1, 1);
        hip_nodeParams.blockDim       = dim3(M_nThPerBckGPU, 1, 1);
        hip_nodeParams.sharedMemBytes = 0;
        void *inputs[4];
        inputs[0]                     = (void *)&kernel_function;
        inputs[1]                     = (void *)&buffer;
        inputs[2]                     = (void *)&iBegin;
        inputs[3]                     = (void *)&iEnd;
        hip_nodeParams.kernelParams   = inputs;
        
        if (M_q_graph) { hip_nodeParams.extra          = NULL; }
        else           { hip_nodeParams.extra          = nullptr; }
        //END::Init new node

        //BEGIN::Dependencies part
        unsigned int nbElemLinks               = links.size();
        unsigned int nbElemKernelNode          = M_hipGraphNode_t.size();
        std::vector<hipGraphNode_t> dependencies;
        M_ListGraphDependencies.push_back(M_hipGraphNode_t.size()-1);
        for (int i = 0; i < nbElemLinks; i++) { 
            if (links[i]==-1) { qFlag=true; }
            if (links[i]!=-1) {
                dependencies.push_back(M_hipGraphNode_t[links[i]]); 
                M_ListGraphDependencies.push_back(links[i]);
            }
        }

        if (M_qViewInfo) {
            std::cout<<"[INFO]: Nb Elem Links  = "<<nbElemLinks<<"\n";
            std::cout<<"[INFO]: Link dependencies with >>> [";
                for (auto v: dependencies) { std::cout << v << " "; } std::cout<<"] <<<\n";
        }
        //END::Dependencies part

        //BEGIN::Add Node to kernel GPU
        if (M_q_graph) { hipGraphAddKernelNode(&M_hipGraphNode_t[M_hipGraphNode_t.size()-1],hip_graph,dependencies.data(),nbElemLinks, &hip_nodeParams); }
        else           { hipGraphAddKernelNode(&M_hipGraphNode_t[M_hipGraphNode_t.size()-1],hip_graph,nullptr,0, &hip_nodeParams); }
        //END::Add Node to kernel GPU

        M_q_graph=true;

        //BEGIN::Final node kernel GPU
        if (qFlag)
        {
            if (M_qViewInfo) {
                std::cout<<"[INFO]:"<<"\n";
                std::cout<<"[INFO]: List >>> [";
                for (auto v: M_hipGraphNode_t) { std::cout << v << " "; } std::cout<<"] <<<\n";
            }
            hipGraphNode_t copyBuffer;
            if (M_qViewInfo) { std::cout<<"[INFO]: Last M_hipGraphNode_t="<<M_hipGraphNode_t.size()<<" : "<<M_hipGraphNode_t[M_hipGraphNode_t.size()-1]<<"\n"; }
            std::vector<hipGraphNode_t> finalDependencies = { M_hipGraphNode_t[ M_hipGraphNode_t.size()-1] };
            hipGraphAddMemcpyNode1D(&copyBuffer,
                                    hip_graph,
                                    dependencies.data(),
                                    dependencies.size(),
                                    hostbuffer,
                                    buffer,
                                    numElems  * sizeof(typename std::decay<decltype(hostbuffer)>::type),
                                    hipMemcpyDeviceToHost);
        }
        //END::Final node kernel GPU

        dependencies.clear();
     #endif
}


template<typename Kernel, typename Input, typename Output>
    void Task::add_cuda(const Kernel& kernel_function,
                                int numElems,
                                int iBegin,int iEnd,
                                Input* buffer,
                                Output* hostbuffer,
                                std::vector<int> links)
{            
    #ifdef UseCUDA  
        bool qFlag=false;
        //BEGIN::Init new node
        cudaGraphNode_t newKernelNode; M_cudaGraphNode_t.push_back(newKernelNode);
        memset(&cuda_nodeParams, 0, sizeof(cuda_nodeParams));

        if (M_qViewInfo) { std::cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<"\n"; }
        //if (M_qViewInfo) { std::cout<<"[INFO]: M_cudaGraphNode_t="<<M_cudaGraphNode_t.size()<<" : "<<M_cudaGraphNode_t[M_cudaGraphNode_t.size()-1]<<"\n"; }
        if (M_qViewInfo) { std::cout<<"[INFO]: Num Graph Node = "<<M_cudaGraphNode_t.size()<<"\n"; }

        //CRTL range
        if (iEnd>numElems) { iEnd=numElems; }
        if (iBegin<0)      { iBegin=0; }

        //printf("[%x]\n",M_cudaGraphNode_t[M_cudaGraphNode_t.size()-1]);

        cuda_nodeParams.func   =  (void *)OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<Kernel,Input>;
        cuda_nodeParams.gridDim        = dim3(M_numBlocksGPU, 1, 1);
        cuda_nodeParams.blockDim       = dim3(M_nThPerBckGPU, 1, 1);
        cuda_nodeParams.sharedMemBytes = 0;
        void *inputs[4];
        inputs[0]                     = (void *)&kernel_function;
        inputs[1]                     = (void *)&buffer;
        inputs[2]                     = (void *)&iBegin;
        inputs[3]                     = (void *)&iEnd;
        cuda_nodeParams.kernelParams   = inputs;
        
        if (M_q_graph) { cuda_nodeParams.extra          = NULL; }
        else           { cuda_nodeParams.extra          = nullptr; }
        //END::Init new node

        //BEGIN::Dependencies part
        unsigned int nbElemLinks               = links.size();
        unsigned int nbElemKernelNode          = M_cudaGraphNode_t.size();
        std::vector<cudaGraphNode_t> dependencies;
        M_ListGraphDependencies.push_back(M_cudaGraphNode_t.size()-1);
        for (int i = 0; i < nbElemLinks; i++) { 
            if (links[i]==-1) { qFlag=true; }
            if (links[i]!=-1) {
                dependencies.push_back(M_cudaGraphNode_t[links[i]]); 
                M_ListGraphDependencies.push_back(links[i]);
            }
        }

        if (M_qViewInfo) {
            std::cout<<"[INFO]: Nb Elem Links  = "<<nbElemLinks<<"\n";
            std::cout<<"[INFO]: Link dependencies with >>> [";
                for (auto v: dependencies) { std::cout << v << " "; } std::cout<<"] <<<\n";
        }
        //END::Dependencies part

        //BEGIN::Add Node to kernel GPU
        if (M_q_graph) { cudaGraphAddKernelNode(&M_cudaGraphNode_t[M_cudaGraphNode_t.size()-1],cuda_graph,dependencies.data(),nbElemLinks, &cuda_nodeParams); }
        else           { cudaGraphAddKernelNode(&M_cudaGraphNode_t[M_cudaGraphNode_t.size()-1],cuda_graph,nullptr,0, &cuda_nodeParams); }
        //END::Add Node to kernel GPU

        M_q_graph=true;

        //BEGIN::Final node kernel GPU
        if (qFlag)
        {
            if (M_qViewInfo) {
                std::cout<<"[INFO]:"<<"\n";
                std::cout<<"[INFO]: List >>> [";
                for (auto v: M_cudaGraphNode_t) { std::cout << v << " "; } std::cout<<"] <<<\n";
            }
            cudaGraphNode_t copyBuffer;
            if (M_qViewInfo) { std::cout<<"[INFO]: Last M_cudaGraphNode_t="<<M_cudaGraphNode_t.size()<<" : "<<M_cudaGraphNode_t[M_cudaGraphNode_t.size()-1]<<"\n"; }
            std::vector<cudaGraphNode_t> finalDependencies = { M_cudaGraphNode_t[ M_cudaGraphNode_t.size()-1] };
            cudaGraphAddMemcpyNode1D(&copyBuffer,
                                    cuda_graph,
                                    dependencies.data(),
                                    dependencies.size(),
                                    hostbuffer,
                                    buffer,
                                    numElems  * sizeof(typename std::decay<decltype(hostbuffer)>::type),
                                    cudaMemcpyDeviceToHost);
        }
        //END::Final node kernel GPU

        dependencies.clear();
     #endif
}



void Task::run()
{
    if (M_qViewInfo) { std::cout<<"<=====================================================================>"<<"\n"; }
    if (M_qViewInfo) { std::cout<<"[INFO]: Run Graph"<<"\n"; }    
    M_t_begin = std::chrono::steady_clock::now();
    //Run HIP Graph
    if (M_q_graph) {   
        //#ifdef COMPILE_WITH_HIP && UseHIP
        //#if defined(COMPILE_WITH_HIP) && defined(UseHIP)
        #ifdef UseHIP
            //hipEventRecord            (hip_start); <<<===not using causes memory crashes
            hipGraphInstantiate       (&hip_graphExec, hip_graph, nullptr, nullptr, 0);
            hipStreamCreateWithFlags  (&hip_graphStream, hipStreamNonBlocking);
            hipGraphLaunch            (hip_graphExec, hip_graphStream);
            hipStreamSynchronize      (hip_graphStream);
            //hipEventRecord            (hip_stop);  <<<===not using causes memory crashes
            //hipEventElapsedTime       (&hip_milliseconds, hip_start, hip_stop);
        #endif
        #ifdef UseCUDA
            //cudaEventRecord           (cuda_start);  <<<===not using causes memory crashes
            cudaGraphInstantiate      (&cuda_graphExec, cuda_graph, nullptr, nullptr, 0);
            cudaStreamCreateWithFlags (&cuda_graphStream, cudaStreamNonBlocking);
            cudaGraphLaunch           (cuda_graphExec, cuda_graphStream);
            cudaStreamSynchronize     (cuda_graphStream);
            //cudaEventRecord           (cuda_stop);  <<<===not using causes memory crashes
            //cudaEventElapsedTime      (&cuda_milliseconds, cuda_start, cuda_stop);
        #endif
    }
    M_t_end = std::chrono::steady_clock::now();
    M_time_laps= std::chrono::duration_cast<std::chrono::microseconds>(M_t_end - M_t_begin).count();
    if (M_qViewInfo) { std::cout<<"<=====================================================================>"<<"\n"; }
}

void Task::close()
{
    if (M_q_graph) {   //HIP Graph
        //#ifdef COMPILE_WITH_HIP && UseHIP
        //#if defined(COMPILE_WITH_HIP) && defined(UseHIP)
        #ifdef UseHIP
            hipGraphExecDestroy     (hip_graphExec);
            hipGraphDestroy         (hip_graph);
            hipStreamDestroy        (hip_graphStream);
            if (M_qDeviceReset)     { hipDeviceReset(); }   // Not be used if Spex Hip AMD activated
            // Explicitly destroys and cleans up all resources associated with the current device in the current process. Any subsequent API call to this device will reinitialize the device.
        #endif
        #ifdef UseCUDA
            cudaGraphExecDestroy     (cuda_graphExec);
            cudaGraphDestroy         (cuda_graph);
            cudaStreamDestroy        (cuda_graphStream);
            if (M_qDeviceReset)      { cudaDeviceReset(); }  // Not be used if Spex Cuda NVidia activated
            // Explicitly destroys and cleans up all resources associated with the current device in the current process. Any subsequent API call to this device will reinitialize the device.
        #endif
        if (M_qViewInfo) { std::cout<<"[INFO]: Close Graph Hip"<<"\n"; }
    }
}



template<typename T>
class PtrTask:public Task
{
    private:
        #ifdef UseHIP
        bufferGraphHIP<T> BUFFER_HIP;   
        #endif 
        #ifdef UseCUDA
        bufferGraphCUDA<T> BUFFER_CUDA;   
        #endif 
        unsigned int bufferSizeBytes;  
    public:
        PtrTask();
        ~PtrTask();

        void set(int nb,T *v);
        void get(T *v);

        template<typename Kernel>
            void add(const Kernel& kernel_function,int iBegin,int iEnd,std::vector<int> links);

};


template<typename T>
PtrTask<T>::PtrTask()
{
    //...
}

template<typename T>
PtrTask<T>::~PtrTask()
{
    //...
}

template<typename T>
void PtrTask<T>::set(int nb,T *v)
{
    #ifdef UseHIP
    if (getNumType()==1) { 
        BUFFER_HIP.memoryInit(nb);  bufferSizeBytes=nb*sizeof(T); std::memcpy(&BUFFER_HIP.data, &v,bufferSizeBytes);  BUFFER_HIP.memmovHostToDevice();
    }
    #endif

    #ifdef UseCUDA
    if (getNumType()==1) { 
        BUFFER_CUDA.memoryInit(nb);  bufferSizeBytes=nb*sizeof(T); std::memcpy(&BUFFER_CUDA.data, &v,bufferSizeBytes); BUFFER_CUDA.memmovHostToDevice();
    }
    #endif
}

template<typename T>
void PtrTask<T>::get(T *v)
{
    #ifdef UseHIP
    if (getNumType()==1) {  BUFFER_HIP.memmovDeviceToHost(); std::memcpy(&v,&BUFFER_HIP.data,bufferSizeBytes); }
    #endif
    #ifdef UseCUDA
    if (getNumType()==2) {  BUFFER_CUDA.memmovDeviceToHost(); std::memcpy(&v,&BUFFER_CUDA.data,bufferSizeBytes); }
    #endif
}

template<typename T>
template<typename Kernel>
void PtrTask<T>::add(const Kernel& kernel_function,int iBegin,int iEnd,std::vector<int> links)
{
    #ifdef UseHIP
    if (getNumType()==1) { Task::add_hip(kernel_function,BUFFER_HIP.size,iBegin,iEnd,BUFFER_HIP.deviceBuffer,BUFFER_HIP.data,links); }
    #endif
    #ifdef UseCUDA
    if (getNumType()==2) { Task::add_cuda(kernel_function,BUFFER_CUDA.size,iBegin,iEnd,BUFFER_CUDA.deviceBuffer,BUFFER_CUDA.data,links); }
    #endif
}



} //End namespace Taskgpu
//--------------------------------------------------------------------------------------------------------------------------------




//================================================================================================================================
// GRAPH IMPLICIT

namespace Taskgpui {

#ifdef UseHIP

struct Task {
    enum class state_t      { capture, update };
    void add_kernel_node    (size_t key, hipKernelNodeParams params, hipStream_t s);
    void update_kernel_node (size_t key, hipKernelNodeParams params);
    state_t state()         { return M_state; }
    ~Task();

private:
    std::unordered_map<size_t, hipGraphNode_t> _node_map;
    state_t        M_state;
    hipGraph_t     M_graph;
    hipGraphExec_t M_graph_exec;
    bool M_qInstantiated = false;
    static void begin_capture  (hipStream_t stream);
    void end_capture           (hipStream_t stream);
    void launch_graph          (hipStream_t stream);

public:
    bool _always_recapture = false;
    template<class Obj>
        void wrap(Obj &o, hipStream_t stream);
};


Task::~Task() {
    if (M_qInstantiated) {
        hipGraphDestroy(M_graph);
        hipGraphExecDestroy(M_graph_exec);
        M_qInstantiated = false;
    }
}

void Task::begin_capture(hipStream_t stream) { hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal); }

void Task::end_capture(hipStream_t stream) {
    if (M_qInstantiated) { hipGraphDestroy(M_graph); }
    hipStreamEndCapture(stream, &M_graph);
    bool need_instantiation;

    if (M_qInstantiated) {
        hipGraphExecUpdateResult updateResult;
        hipGraphNode_t errorNode;
        hipGraphExecUpdate(M_graph_exec, M_graph, &errorNode, &updateResult);
        if (M_graph_exec == nullptr || updateResult != hipGraphExecUpdateSuccess) {
            hipGetLastError();
            if (M_graph_exec != nullptr) { hipGraphExecDestroy(M_graph_exec); }
            need_instantiation = true;
        } else {
            need_instantiation = false;
        }
    } else {
        need_instantiation = true;
    }

    if (need_instantiation) {
        hipGraphInstantiate(&M_graph_exec, M_graph, nullptr, nullptr, 0);
    }
    M_qInstantiated = true;
}

template<class Obj>
void Task::wrap(Obj &o, hipStream_t stream) 
{
    if (!_always_recapture && M_qInstantiated) {
        M_state = state_t::update;
        o(*this, stream);
    } 
    else
    {
        M_state = state_t::capture;
        begin_capture(stream);
        o(*this, stream);
        end_capture(stream);
    }
    launch_graph(stream);
}

void Task::launch_graph(hipStream_t stream) {
    if (M_qInstantiated) { hipGraphLaunch(M_graph_exec, stream);}
}

void Task::add_kernel_node(size_t key, hipKernelNodeParams params, hipStream_t stream)
{
    hipStreamCaptureStatus capture_status;
    hipGraph_t graph;
    const hipGraphNode_t *deps;
    size_t dep_count;
    hipStreamGetCaptureInfo_v2(stream, &capture_status, nullptr, &graph, &deps, &dep_count);
    hipGraphNode_t new_node;
    hipGraphAddKernelNode(&new_node, graph, deps, dep_count, &params);
    _node_map[key] = new_node;
    hipStreamUpdateCaptureDependencies(stream, &new_node, 1, 1);
}

void Task::update_kernel_node(size_t key, hipKernelNodeParams params)
{
    hipGraphExecKernelNodeSetParams(M_graph_exec, _node_map[key], &params);
}

#endif

} //End namespace Taskgpui
//--------------------------------------------------------------------------------------------------------------------------------



//================================================================================================================================
// CLASS Task: Unified Memory
//================================================================================================================================

// NOTA: With Unified Memory: GPU accesses data directly from the "host" may be used without a separate "host" allocation and no copy routine is required,
// greatly simplifying and reducing the size of the program. With:
// - System Allocated: no other changes required.
// - Managed Memory: data allocation changed to use (cuda or hip)-MallocManaged(),which returns a pointer valid from both host and device code.


namespace Taskgpuu {


//#if defined(useHIP) || defined(useCUDA)
#ifdef UseHIP

template<typename T>
class Task
{
    private:
        int         M_numBlocksGPU;
        int         M_nbThGPU;
        bool        M_isFree;
        bool        M_qViewInfo;
        bool        M_qDeviceReset;
        long int    M_time_laps;
        float       hip_milliseconds;
        float       cuda_milliseconds;
        int         M_nbStreams;
        int         M_numModel;
        int         M_numType;
        std::chrono::steady_clock::time_point M_t_begin,M_t_end;

        template<typename Kernel>
            void serial(const Kernel& kernel_function);

        template<typename Kernel>
            void stream(const Kernel& kernel_function);

    public:
        #ifdef UseHIP
        bufferGraphUnifiedHIP<T> BUFFER_HIP;   
        #endif 
        #ifdef UseCUDA
        bufferGraphUnifiedHIP<T> BUFFER_CUDA;   
        #endif 

        Task();
        ~Task();

        void open(int nbBlock,int NbTh);
        void close();
        void setViewInfo     (bool b)        {  M_qViewInfo = b; }
        void setDeviceReset  (bool b)        {  M_qDeviceReset = b; }
        void setNbStreams    (int v)         {  M_nbStreams    = v; }
        void setDeviceHIP    (int v);
        void setDeviceCUDA   (int v);
        void setNumType      (int v)         {  M_numType      = v; }
        int  getNbStreams    () const        {  return(M_nbStreams);    }
        bool isDeviceReset   () const        {  return(M_qDeviceReset); }
        int  getNumType      () const        {  return(M_numType);   }

        long int  getTimeLaps     () const        {  return(M_time_laps);   }

        template<typename Kernel>
            void run(const Kernel& kernel_function);

        void debriefing();
};

template<typename T>
Task<T>::~Task()
{   
    if (!M_isFree) { close(); }
}

template<typename T>
Task<T>::Task()
{
    M_numBlocksGPU    = 512;
    M_isFree          = false;
    M_qViewInfo       = true;
    M_qDeviceReset    = false;
    M_nbStreams       = 3;
    M_numType         = 1; // 1:HIP 2:CUDA
    M_numModel        = 1; // 1:Serial 2:Stream
    M_time_laps       = 0;
}

template<typename T>
void Task<T>::open(int nbBlock,int NbTh)
{
    M_nbThGPU= NbTh;
    M_numBlocksGPU=nbBlock;
    #ifdef UseHIP
    BUFFER_HIP.memoryInit(NbTh);
    #endif 
    #ifdef UseCUDA
    BUFFER_CUDA.memoryInit(NbTh);
    #endif 
}


template<typename T>
void Task<T>::close()
{
    #ifdef UseHIP
    hipFree(BUFFER_HIP.data); 
    if (M_qDeviceReset)      { hipDeviceReset(); }
    #endif   
    #ifdef UseCUDA
    cudaFree(BUFFER_CUDA.data); 
    if (M_qDeviceReset)      { cudaDeviceReset(); }
    #endif 
    M_isFree=true;
}

template<typename T>
void Task<T>::setDeviceHIP(int v)
{
    #ifdef UseHIP
    int numDevices=0; hipGetDeviceCount(&numDevices);
    if ((v>=0) && (v<=numDevices)) { hipSetDevice(v); }
    #endif
}

template<typename T>
void Task<T>::setDeviceCUDA(int v)
{
    #ifdef UseCUDA
    int numDevices=0; cudaGetDeviceCount(&numDevices);
    if ((v>=0) && (v<=numDevices)) { cudaSetDevice(v); }
    #endif
}


template<typename T>
void Task<T>::debriefing()
{
    if (M_qViewInfo) { 
        std::cout<<"<=====================================================================>"<<"\n"; 
        std::cout<<"[INFO]: Debriefing"<<"\n";  
        std::cout<<"[INFO]: Elapsed microseconds : "<<M_time_laps<< " us\n";
    }
}

template<typename T>
template<typename Kernel>
void Task<T>::serial(const Kernel& kernel_function)
{
    M_t_begin = std::chrono::steady_clock::now();
    int num_blocks = (BUFFER_HIP.size + M_numBlocksGPU - 1) / M_numBlocksGPU;
    dim3 thread_block(M_numBlocksGPU, 1, 1);
    dim3 grid(num_blocks, 1);
    #ifdef UseHIP
    hipLaunchKernelGGL(OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D,grid,thread_block,0,0,kernel_function,BUFFER_HIP.data,0,BUFFER_HIP.size);
    hipDeviceSynchronize();
    #endif 
    #ifdef UseCUDA
    OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<<<grid,thread_block,0,0>>>(kernel_function,BUFFER_CUDA.data,0,BUFFER_CUDA.size);
    cudaDeviceSynchronize();
    #endif 
    M_t_end = std::chrono::steady_clock::now();
    M_time_laps= std::chrono::duration_cast<std::chrono::microseconds>(M_t_end - M_t_begin).count();
}

template<typename T>
template<typename Kernel>
void Task<T>::stream(const Kernel& kernel_function)
{
    #ifdef UseHIP
        int   numVersion=1;
        float ms;
        int   i;

        //const int blockSize = 20;  //const int blockSize M_block_size;
        const int blockSize = M_numBlocksGPU;
        const int numElems=BUFFER_HIP.size;
        const int streamSize = numElems / M_nbStreams;
        const int streamBytes = streamSize * sizeof(decltype(BUFFER_HIP.data));
        const int sz = numElems * sizeof(decltype(BUFFER_HIP.data));
        int grid=streamSize/blockSize;
        if (M_qViewInfo) {
            std::cout<<"[INFO]: Block size :"<<M_numBlocksGPU<<"\n";
            std::cout<<"[INFO]: nb Streams :"<<M_nbStreams<<"\n";
        }
        if ((grid==0) && (M_qViewInfo)) {  std::cout<<"[INFO]: Grid size error : ==> auto correction"<<"\n"; }
        grid=max(streamSize/blockSize,1);
        hipEvent_t startEvent, stopEvent, dummyEvent;
        hipStream_t stream[M_nbStreams];
        hipEventCreate(&startEvent) ;
        hipEventCreate(&stopEvent) ;
        hipEventCreate(&dummyEvent) ;
        for (i = 0; i < M_nbStreams; ++i) { hipStreamCreate(&stream[i]) ; }

       //hipMallocManaged((void **) &data,sz,hipMemAttachHost);
        hipEventRecord(startEvent,0) ;

        if (numVersion==1) {
            for (i = 0; i < M_nbStreams; ++i) {
                int offset = i * streamSize;
                hipLaunchKernelGGL(OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D,grid,blockSize,0,stream[i],kernel_function,BUFFER_HIP.data,offset);
            }
            hipEventRecord(stopEvent, 0);
            hipEventSynchronize(stopEvent);
            //hipEventElapsedTime(&ms, startEvent, stopEvent);
            //printf("Time (ms): %f\n", ms);
        }

       
       // cleanup
        hipEventDestroy(startEvent) ;
        hipEventDestroy(stopEvent) ;
        hipEventDestroy(dummyEvent);
        for (int i = 0; i < M_nbStreams; ++i) { hipStreamDestroy(stream[i]); }
    #endif 
}


template<typename T>
template<typename Kernel>
void Task<T>::run(const Kernel& kernel_function)
{
    M_t_begin = std::chrono::steady_clock::now();
    if (M_numModel==1) { //SERIAL MODEL
        if (M_numType==1) {
            serial(kernel_function);
        }
    }

    if (M_numModel==2) { //STREAM MODEL
        if (M_numType==1) {
            stream(kernel_function);
        }
    }


    M_t_end = std::chrono::steady_clock::now();
    M_time_laps= std::chrono::duration_cast<std::chrono::microseconds>(M_t_end - M_t_begin).count();
}


#endif

}//End namespace Taskgpuu
//--------------------------------------------------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------------------------------------------

void resetAndDestroyAllGPU()
{
    //Reset and explicitly destroy all resources associated with the current device
    int numDevices=0;
    #ifdef UseHIP
        HIP_CHECK(hipGetDeviceCount(&numDevices));
        for (int i = 0; i < numDevices; i++) { HIP_CHECK(hipSetDevice(i)); HIP_CHECK(hipDeviceReset()); }
    #endif
    #ifdef UseCUDA
        CUDA_CHECK(cudaGetDeviceCount(&numDevices));
        for (int i = 0; i < numDevices; i++) { CUDA_CHECK(cudaSetDevice(i)); CUDA_CHECK(cudaDeviceReset()); }
    #endif
}

//--------------------------------------------------------------------------------------------------------------------------------



template <typename T>
bool isSpecxGPUFunctionBeta(T& fcv)
{
    std::string s1=typeid(fcv).name(); int l=s1.length();
    if (l>1) {
        if (s1.find("SpCallableType0EE") != std::string::npos) { return (true);  } //"SpCpu"
        else if (s1.find("SpCallableType2EE") != std::string::npos) { return (true);  } //"SpHip"
        else if (s1.find("SpCallableType1EE") != std::string::npos) { return (true);  } //"SpCuda
        
    }
    return (false);
}

template <typename T>
int numSpecxFunctionBeta(T& fcv)
{
    std::string s1=typeid(fcv).name(); int l=s1.length();
    //std::cout<<"s1="<<s1<<"\n";
    if (l>1) {
        if (s1.find("SpCallableType0EE") != std::string::npos) { return (1);  }           //"SpCpu"
        else if (s1.find("SpCallableType2EE") != std::string::npos) { return (2);  }      //"SpHip"
        else if (s1.find("SpCallableType1EE") != std::string::npos) { return (3);  }      //"SpCuda
        else if (s1.find("SpArrayView") != std::string::npos) { return (20);  }           //SpArrayView
        else if (s1.find("SpArrayAccessorIS1_EE") != std::string::npos) { return (21);  } //SpReadArray      
        else if (s1.find("SpDataAccessMode0") != std::string::npos) { return (10);  }     //"SpRead
        else if (s1.find("SpDataAccessMode1") != std::string::npos) { return (11);  }     //"SpWrite
        else if (s1.find("SpDataAccessMode3") != std::string::npos) { return (12);  }     //SpCommutativeWrite
    }
    return (0);
}


//================================================================================================================================
// THE END.
//================================================================================================================================

}