

//================================================================================================================================
// Get Information System CPU/GPU Hardware
//================================================================================================================================

namespace HARD {

void readFileViewInformation(char *filename) 
{
	FILE* FICH = NULL;
    int c = 0;
	FICH = fopen(filename, "r");
    if (FICH != NULL) { do { c = fgetc(FICH); printf("%c",c); } while (c != EOF); fclose(FICH); }
}


void scanInformationSystem()
{
	int Value;
	std::cout <<"\n";
	std::cout << "[INFO]: Scan Information System..."<<"\n";
	Value=std::system("lscpu>InfoSystemCPU.txt");
	Value=std::system("lshw -C display>InfoSystemGPU.txt");
	std::cout <<"\n";
    std::cout <<"\n";
}

void getInformationCPU()
{
	std::cout <<"\n";
	std::cout << "[INFO]: Information CPU"<<"\n";
    std::cout <<"\n";
	readFileViewInformation("InfoSystemCPU.txt");
	std::cout <<"\n";
    std::cout <<"\n";
}

void getInformationGPU()
{
	std::cout <<"\n";
	std::cout << "[INFO]: Information GPU"<<"\n";
    std::cout <<"\n";
	readFileViewInformation("InfoSystemGPU.txt");
	std::cout <<"\n";
    std::cout <<"\n";
}

#ifdef USE_MPI
void getMpiInformation(int argc, char *argv[])
{
	//BEGIN::INFO MPI
	bool qFullInfoSystem=false;
    MPI_Init(NULL, NULL);
    int world_rank,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name,&name_len);

    if (world_rank == 0) { 
	  	std::cout <<"\n";
      	int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
      	std::cout << "[INFO]: MPI Name worlds processor: "<<processor_name<<"\n";
      	std::cout << "[INFO]: MPI Nb CPU available: "<<numCPU<< "\n";
      	std::cout <<"\n";
      	std::cout << "[INFO]: MPI Scan..."<<"\n";
    }
    std::cout << "[INFO]: MPI Rank: "<<world_rank<<" out of "<<world_size<<"\n";
    MPI_Finalize();
	//END::INFO MPI
}
#endif

#ifdef USE_OpenMP
void getOpenMPInformation()
{
    std::cout << "[INFO]: OpenMP Nb num procs: "<<omp_get_num_procs( )<< "\n";
    std::cout << "[INFO]: OpenMP Nb max threads: "<<omp_get_max_threads()<< "\n";
}
#endif

void getShortInformationGPU()
{
	int deviceCount=0;
	std::cout <<"\n";
	std::cout << "[INFO]: Information GPU"<<"\n";

	#ifdef COMPILE_WITH_HIP && UseHIP
		hipGetDeviceCount(&deviceCount);
		if (deviceCount>0) {
			std::cout << "[INFO]: Number of available GPUs AMD: " << deviceCount << "\n";
			for (int deviceId = 0; deviceId < deviceCount; ++deviceId) {
				hipSetDevice(deviceId);
				std::cout << "[INFO]: GPU " << deviceId << " initialized and resources allocated." << "\n";
			}
		}
	#endif

    #if defined(COMPILE_WITH_CUDA) && defined(UseCUDA)
		cudaGetDeviceCount(&deviceCount);
		if (deviceCount>0) {
			std::cout << "[INFO]: Number of available GPUs NVIDIA: " << deviceCount << "\n";
			for (int deviceId = 0; deviceId < deviceCount; ++deviceId) {
				cudaSetDevice(deviceId);
				std::cout << "[INFO]: GPU " << deviceId << " initialized and resources allocated." << "\n";
			}
		}
	#endif
	std::cout <<"\n";
	if (deviceCount == 0) { std::cerr << "[INFO]: No GPUs found. Exiting." << "\n"; }
}


void getHipInformation()
{
  //BEGIN::INFO HIP AMD
     #if defined(COMPILE_WITH_HIP) && defined(UseHIP)
    std::cout<<std::endl;
    int numDevices=0;
    HIP_CHECK(hipGetDeviceCount(&numDevices));
    std::cout<<"[INFO]: Get numDevice                = "<<numDevices<<"\n";
    int deviceID=0;
    HIP_CHECK(hipGetDevice(&deviceID));
    std::cout<<"[INFO]: Get deviceID activated       = "<<deviceID<<"\n";
    deviceID=0;
    hipSetDevice(deviceID);

    hipDeviceProp_t devProp;
    for (int i = 0; i < numDevices; i++)
    {
                HIP_CHECK(hipSetDevice(i));
                HIP_CHECK(hipGetDeviceProperties(&devProp,i));
                std::cout<<"[INFO]:"<<std::endl;
                std::cout<<"[INFO]: DeviceID                     = "<<i<<std::endl;
                std::cout<<"[INFO]: Agent prop name              = "<< devProp.name<<std::endl;
                std::cout<<"[INFO]: System minor                 = "<< devProp.minor<<std::endl;
                std::cout<<"[INFO]: System major                 = "<< devProp.major<<std::endl;
                std::cout<<"[INFO]: Memory Clock Rate (KHz)      = "<< devProp.memoryClockRate<<std::endl;
                std::cout<<"[INFO]: Memory Bus Width (bits)      = "<< devProp.memoryBusWidth<<std::endl;
                std::cout<<"[INFO]: Peak Memory Bandwidth (GB/s) = "<< 2.0*devProp.memoryClockRate*(devProp.memoryBusWidth/8)/1.0e6<<std::endl;
                std::cout<<"[INFO]: max ThreadsPerBlock          = "<< devProp.maxThreadsPerBlock<<std::endl;
                std::cout<<"[INFO]: max ThreadsPerMultiProcessor = "<< devProp.maxThreadsPerMultiProcessor<<std::endl;
                std::cout<<"[INFO]: max ThreadsDim 3D            = "<< devProp.maxThreadsDim[0]<<" "<<devProp.maxThreadsDim[1]<<" "<<devProp.maxThreadsDim[2]<<std::endl;
                std::cout<<"[INFO]: max Grid Size 3D             = "<< devProp.maxGridSize[0]<<" "<<devProp.maxGridSize[1]<<" "<<devProp.maxGridSize[2]<<std::endl;
                std::cout<<"[INFO]: warpSize:                    = "<< devProp.warpSize << "\n";
                std::cout<<"[INFO]: regsPerBlock:                = "<< devProp.regsPerBlock << "\n";
                std::cout<<"[INFO]: concurrentKernels:           = "<< devProp.concurrentKernels << "\n";
                std::cout<<"[INFO]: total Global Mem             = "<< devProp.totalGlobalMem<<std::endl;
                std::cout<<"[INFO]: shared Mem Per Block         = "<< devProp.sharedMemPerBlock<<std::endl;
    }
    //(...)
    HIP_CHECK(hipSetDevice(0));
    std::cout<<std::endl;
    //END::INFO HIP AMD
    #endif
}

void getCudaInformation()
{
    //Nota: no code fusion because if hybrid CUDA and HIP system used
    //BEGIN::INFO CUDA NVIDIA
    #if defined(COMPILE_WITH_CUDA) && defined(UseCUDA)
        std::cout<<std::endl;
        int numDevices=0;
        CUDA_CHECK(cudaGetDeviceCount(&numDevices));
        std::cout<<"[INFO]: Get numDevice                = "<<numDevices<<"\n";
        int deviceID=0;
        CUDA_CHECK(cudaGetDevice(&deviceID));
        std::cout<<"[INFO]: Get deviceID activated       = "<<deviceID<<"\n";
        deviceID=0;
        cudaSetDevice(deviceID);

        hipDeviceProp_t devProp;
        for (int i = 0; i < numDevices; i++)
        {
                    CUDA_CHECK(cudaSetDevice(i));
                    CUDA_CHECK(cudaGetDeviceProperties(&devProp,i));
                    std::cout<<"[INFO]:"<<std::endl;
                    std::cout<<"[INFO]: DeviceID                     = "<<i<<std::endl;
                    std::cout<<"[INFO]: Agent prop name              = "<< devProp.name<<std::endl;
                    std::cout<<"[INFO]: System minor                 = "<< devProp.minor<<std::endl;
                    std::cout<<"[INFO]: System major                 = "<< devProp.major<<std::endl;
                    std::cout<<"[INFO]: Memory Clock Rate (KHz)      = "<< devProp.memoryClockRate<<std::endl;
                    std::cout<<"[INFO]: Memory Bus Width (bits)      = "<< devProp.memoryBusWidth<<std::endl;
                    std::cout<<"[INFO]: Peak Memory Bandwidth (GB/s) = "<< 2.0*devProp.memoryClockRate*(devProp.memoryBusWidth/8)/1.0e6<<std::endl;
                    std::cout<<"[INFO]: max ThreadsPerBlock          = "<< devProp.maxThreadsPerBlock<<std::endl;
                    std::cout<<"[INFO]: max ThreadsPerMultiProcessor = "<< devProp.maxThreadsPerMultiProcessor<<std::endl;
                    std::cout<<"[INFO]: max ThreadsDim 3D            = "<< devProp.maxThreadsDim[0]<<" "<<devProp.maxThreadsDim[1]<<" "<<devProp.maxThreadsDim[2]<<std::endl;
                    std::cout<<"[INFO]: max Grid Size 3D             = "<< devProp.maxGridSize[0]<<" "<<devProp.maxGridSize[1]<<" "<<devProp.maxGridSize[2]<<std::endl;
                    std::cout<<"[INFO]: warpSize:                    = "<< devProp.warpSize << "\n";
                    std::cout<<"[INFO]: regsPerBlock:                = "<< devProp.regsPerBlock << "\n";
                    std::cout<<"[INFO]: concurrentKernels:           = "<< devProp.concurrentKernels << "\n";
                    std::cout<<"[INFO]: total Global Mem             = "<< devProp.totalGlobalMem<<std::endl;
                    std::cout<<"[INFO]: shared Mem Per Block         = "<< devProp.sharedMemPerBlock<<std::endl;
        }
        //(...)
        CUDA_CHECK(cudaSetDevice(0));
        std::cout<<std::endl;
    //END::INFO CUDA NVIDIA
    #endif
}


void getInformationSystem()
{
    std::cout<<"[INFO]: ======================================================================================== "<<"\n";
    std::cout<<"[INFO]: Get Information System "<<"\n";
    getInformationCPU();
    scanInformationSystem();
    getInformationGPU();
    #if defined(COMPILE_WITH_HIP) && defined(UseHIP)
        getHipInformation();
    #endif
    #if defined(COMPILE_WITH_CUDA) && defined(UseCUDA)
        getCudaInformation();
    #endif
    std::cout<<"[INFO]: ======================================================================================== "<<"\n";
    std::cout<<"[INFO]: "<<"\n";
}

}//END::namespace