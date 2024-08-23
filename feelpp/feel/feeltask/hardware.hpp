/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

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

  @file
  @author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  @Date: 2024-06-04
  @copyright 2019 Feel++ Consortium
*/

namespace Feel
{

namespace SystemInformation
{

class Hardware
{
  private:
    void readFileViewInformation( std::string filename );
    void scanInformationSystem();
    void getInformationCPU();
    void getInformationGPU();
    void getShortInformationGPU();
    void getHipInformation();
    void getCudaInformation();
    void getOpenMPInformation();

  public:
    Hardware();
    ~Hardware();
    void getInformationSystem();

    // void getMpiInformation(int argc, char *argv[]);
};

Hardware::Hardware() {}
Hardware::~Hardware() {}

void Hardware::readFileViewInformation( std::string filename )
{
    std::ifstream inputFile( filename );
    std::string line;
    
    if ( inputFile )
    {
        while ( getline( inputFile, line ) )
        {
            VLOG( 1 ) << line << std::endl;
        }
        inputFile.close();
    }

}

void Hardware::scanInformationSystem()
{
    // Scan Information System..."
    int cmd1 = std::system( "lscpu>InfoSystemCPU.txt" );
    if ( cmd1 != 0 )
    {
        VLOG( 1 ) << "Error command system get information CPU" << std::endl;
    }
    int cmd2 = std::system( "lshw -C display>InfoSystemGPU.txt" );
    if ( cmd2 != 0 )
    {
        VLOG( 1 ) << "Error command system  get information GPU" << std::endl;
    }
}

void Hardware::getInformationCPU()
{
    // Get Information CPU
    VLOG( 1 ) << "Information CPU" << std::endl;
    readFileViewInformation( "InfoSystemCPU.txt" );
}

void Hardware::getInformationGPU()
{
    // Get Information GPU
    VLOG( 1 ) << "Information GPU" << std::endl;
    readFileViewInformation( "InfoSystemGPU.txt" );
}

/*
void Hardware::getMpiInformation(int argc, char *argv[])
{
    VLOG(1) <<"Information MPI"<< std::endl;
#ifdef USE_MPI
    bool qFullInfoSystem=false;
    MPI_Init(NULL, NULL);
    int world_rank,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name,&name_len);

    if (world_rank == 0) {
        int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
        VLOG(1) << "MPI Name worlds processor="<<processor_name<<"\n";
        VLOG(1) << "MPI Nb CPU available="<<numCPU<< "\n";
        VLOG(1) << "MPI Scan..."<<"\n";
    }
    VLOG(1) << "MPI Rank= "<<world_rank<<"out of "<<world_size<<"\n";
    MPI_Finalize();
#endif
}
*/

void Hardware::getOpenMPInformation()
{
#ifdef USE_OpenMP
    VLOG( 1 ) << "OpenMP Nb num procs=" << omp_get_num_procs() << "\n";
    VLOG( 1 ) << "OpenMP Nb max threads=" < < < omp_get_max_threads() << "\n";
#endif
}

void Hardware::getShortInformationGPU()
{
    int deviceCount = 0;
    VLOG( 1 ) << "Information GPU"
              << "\n";
#if defined( COMPILE_WITH_HIP ) && defined( UseHIP )
    hipGetDeviceCount( &deviceCount );
    if ( deviceCount > 0 )
    {
        VLOG( 1 ) << "Number of available GPUs AMD=" << deviceCount << "\n";
        for ( int deviceId = 0; deviceId < deviceCount; ++deviceId )
        {
            hipSetDevice( deviceId );
            VLOG( 1 ) << "GPU " << deviceId << "initialized and resources allocated."
                      << "\n";
        }
    }
#endif

#if defined( COMPILE_WITH_CUDA ) && defined( UseCUDA )
    cudaGetDeviceCount( &deviceCount );
    if ( deviceCount > 0 )
    {
        VLOG( 1 ) << "Number of available GPUs NVIDIA=" << deviceCount << "\n";
        for ( int deviceId = 0; deviceId < deviceCount; ++deviceId )
        {
            cudaSetDevice( deviceId );
            VLOG( 1 ) << "GPU " << deviceId << "initialized and resources allocated."
                      << "\n";
        }
    }
#endif
    if ( deviceCount == 0 )
    {
        VLOG( 1 ) << "No GPUs found. Exiting."
                  << "\n";
    }
}

void Hardware::getHipInformation()
{
#if defined( COMPILE_WITH_HIP ) && defined( UseHIP )
    int numDevices = 0;
    HIP_CHECK( hipGetDeviceCount( &numDevices ) );
    VLOG( 1 ) << "Get numDevice=" << numDevices << "\n";
    int deviceID = 0;
    HIP_CHECK( hipGetDevice( &deviceID ) );
    VLOG( 1 ) << "Get deviceID activated=" << deviceID << "\n";
    deviceID = 0;
    hipSetDevice( deviceID );

    hipDeviceProp_t devProp;
    for ( int i = 0; i < numDevices; i++ )
    {
        HIP_CHECK( hipSetDevice( i ) );
        HIP_CHECK( hipGetDeviceProperties( &devProp, i ) );
        VLOG( 1 ) << "DeviceID                     = " << i << std::endl;
        VLOG( 1 ) << "Agent prop name              = " << devProp.name << std::endl;
        VLOG( 1 ) << "System minor                 = " << devProp.minor << std::endl;
        VLOG( 1 ) << "System major                 = " << devProp.major << std::endl;
        VLOG( 1 ) << "Memory Clock Rate (KHz)      = " << devProp.memoryClockRate << std::endl;
        VLOG( 1 ) << "Memory Bus Width (bits)      = " << devProp.memoryBusWidth << std::endl;
        VLOG( 1 ) << "Peak Memory Bandwidth (GB/s) = " << 2.0 * devProp.memoryClockRate * ( devProp.memoryBusWidth / 8 ) / 1.0e6 << std::endl;
        VLOG( 1 ) << "max ThreadsPerBlock          = " << devProp.maxThreadsPerBlock << std::endl;
        VLOG( 1 ) << "max ThreadsPerMultiProcessor = " << devProp.maxThreadsPerMultiProcessor << std::endl;
        VLOG( 1 ) << "max ThreadsDim 3D            = " << devProp.maxThreadsDim[0] << "" << devProp.maxThreadsDim[1] << "" << devProp.maxThreadsDim[2] << std::endl;
        VLOG( 1 ) << "max Grid Size 3D             = " << devProp.maxGridSize[0] << "" << devProp.maxGridSize[1] << "" << devProp.maxGridSize[2] << std::endl;
        VLOG( 1 ) << "warpSize:                    = " << devProp.warpSize << "\n";
        VLOG( 1 ) << "regsPerBlock:                = " << devProp.regsPerBlock << "\n";
        VLOG( 1 ) << "concurrentKernels:           = " << devProp.concurrentKernels << "\n";
        VLOG( 1 ) << "total Global Mem             = " << devProp.totalGlobalMem << std::endl;
        VLOG( 1 ) << "shared Mem Per Block         = " << devProp.sharedMemPerBlock << std::endl;
    }
    HIP_CHECK( hipSetDevice( 0 ) );
#endif
}

void Hardware::getCudaInformation()
{
    // Nota: no code fusion because if hybrid CUDA and HIP system used
#if defined( COMPILE_WITH_CUDA ) && defined( UseCUDA )
    int numDevices = 0;
    CUDA_CHECK( cudaGetDeviceCount( &numDevices ) );
    VLOG( 1 ) << "Get numDevice                = " << numDevices << "\n";
    int deviceID = 0;
    CUDA_CHECK( cudaGetDevice( &deviceID ) );
    VLOG( 1 ) << "Get deviceID activated       = " << deviceID << "\n";
    deviceID = 0;
    cudaSetDevice( deviceID );
    cudaDeviceProp_t devProp;
    for ( int i = 0; i < numDevices; i++ )
    {
        CUDA_CHECK( cudaSetDevice( i ) );
        CUDA_CHECK( cudaGetDeviceProperties( &devProp, i ) );
        VLOG( 1 ) << "DeviceID                     = " << i << std::endl;
        VLOG( 1 ) << "Agent prop name              = " << devProp.name << std::endl;
        VLOG( 1 ) << "System minor                 = " << devProp.minor << std::endl;
        VLOG( 1 ) << "System major                 = " << devProp.major << std::endl;
        VLOG( 1 ) << "Memory Clock Rate (KHz)      = " << devProp.memoryClockRate << std::endl;
        VLOG( 1 ) << "Memory Bus Width (bits)      = " << devProp.memoryBusWidth << std::endl;
        VLOG( 1 ) << "Peak Memory Bandwidth (GB/s) = " << 2.0 * devProp.memoryClockRate * ( devProp.memoryBusWidth / 8 ) / 1.0e6 << std::endl;
        VLOG( 1 ) << "max ThreadsPerBlock          = " << devProp.maxThreadsPerBlock << std::endl;
        VLOG( 1 ) << "max ThreadsPerMultiProcessor = " << devProp.maxThreadsPerMultiProcessor << std::endl;
        VLOG( 1 ) << "max ThreadsDim 3D            = " << devProp.maxThreadsDim[0] << "" << devProp.maxThreadsDim[1] << "" << devProp.maxThreadsDim[2] << std::endl;
        VLOG( 1 ) << "max Grid Size 3D             = " << devProp.maxGridSize[0] << "" << devProp.maxGridSize[1] << "" << devProp.maxGridSize[2] << std::endl;
        VLOG( 1 ) << "warpSize:                    = " << devProp.warpSize << "\n";
        VLOG( 1 ) << "regsPerBlock:                = " << devProp.regsPerBlock << "\n";
        VLOG( 1 ) << "concurrentKernels:           = " << devProp.concurrentKernels << "\n";
        VLOG( 1 ) << "total Global Mem             = " << devProp.totalGlobalMem << std::endl;
        VLOG( 1 ) << "shared Mem Per Block         = " << devProp.sharedMemPerBlock << std::endl;
    }
    CUDA_CHECK( cudaSetDevice( 0 ) );
#endif
}

void Hardware::getInformationSystem()
{
    VLOG( 1 ) << "Get Information System"
              << "\n";
    getInformationCPU();
    scanInformationSystem();
    getInformationGPU();
    getShortInformationGPU();
#if defined( COMPILE_WITH_HIP ) && defined( UseHIP )
    getHipInformation();
#endif
#if defined( COMPILE_WITH_CUDA ) && defined( UseCUDA )
    getCudaInformation();
#endif
    getOpenMPInformation();
}

} // namespace SystemInformation

} // namespace Feel
