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

//================================================================================================================================
// Tools to manage memories between CPU and GPU
//================================================================================================================================

#ifdef UseHIP
struct VectorBuffer
{
    unsigned int dimension;
    unsigned int dimensionSizeof;
    unsigned int pitch;
    double* data;

    /////////////////////////////////////////////////////////////

    class DataDescr
    {
        std::size_t size;

      public:
        explicit DataDescr( const std::size_t inSize = 0 )
            : size( inSize ) {}

        auto getSize() const
        {
            return size;
        }
    };

    using DataDescriptor = DataDescr;

    std::size_t memmovNeededSize() const
    {
        return dimensionSizeof;
    }

    template <class DeviceMemmov>
    auto memmovHostToDevice( DeviceMemmov& mover, void* devicePtr, [[maybe_unused]] std::size_t size )
    {
        assert( size == dimensionSizeof );
        double* doubleDevicePtr = reinterpret_cast<double*>( devicePtr );
        mover.copyHostToDevice( doubleDevicePtr, data, dimensionSizeof );
        return DataDescr( dimension );
    }

    template <class DeviceMemmov>
    void memmovDeviceToHost( DeviceMemmov& mover, void* devicePtr, [[maybe_unused]] std::size_t size, const DataDescr& /*inDataDescr*/ )
    {
        assert( size == dimensionSizeof );
        double* doubleDevicePtr = reinterpret_cast<double*>( devicePtr );
        mover.copyDeviceToHost( data, doubleDevicePtr, dimensionSizeof );
    }
};

#endif

//================================================================================================================================
//
//================================================================================================================================

#ifdef UseHIP
// Structure to manage a buffer for HIP (Heterogeneous-compute Interface for Portability)
	template <typename ptrtype>
	struct bufferGraphHIP
	{
		unsigned int size;          // Size of the buffer
		ptrtype* data;             // Pointer to host memory
		ptrtype* deviceBuffer;     // Pointer to device memory (GPU)

		// Initializes the buffer size and allocates host memory
		void memoryInit(int dim)
		{
			size = dim;  // Set the size of the buffer
			data = (ptrtype*)malloc(sizeof(ptrtype) * size); // Allocate host memory
		}

		// Transfers data from host to device memory
		void memmovHostToDevice()
		{
			hipMalloc((void**)&deviceBuffer, sizeof(ptrtype) * size); // Allocate device memory
			hipMemcpy(deviceBuffer, data, sizeof(ptrtype) * size, hipMemcpyHostToDevice); // Copy data from host to device
		}

		// Transfers data from device to host memory and frees device memory
		void memmovDeviceToHost()
		{
			hipMemcpy(data, deviceBuffer, sizeof(ptrtype) * size, hipMemcpyDeviceToHost); // Copy data from device to host
			hipFree(deviceBuffer); // Free the allocated device memory
		}
	};
#endif

#ifdef UseCUDA
	// Structure to manage a buffer for CUDA (Compute Unified Device Architecture)
	template <typename ptrtype>
	struct bufferGraphCUDA
	{
		unsigned int size;          // Size of the buffer
		ptrtype* data;             // Pointer to host memory
		ptrtype* deviceBuffer;     // Pointer to device memory (GPU)

		// Initializes the buffer size and allocates host memory
		void memoryInit(int dim)
		{
			size = dim;  // Set the size of the buffer
			data = (ptrtype*)malloc(sizeof(ptrtype) * size); // Allocate host memory
		}

		// Transfers data from host to device memory
		void memmovHostToDevice()
		{
			cudaMalloc((void**)&deviceBuffer, sizeof(ptrtype) * size); // Allocate device memory
			cudaMemcpy(deviceBuffer, data, sizeof(ptrtype) * size, cudaMemcpyHostToDevice); // Copy data from host to device
		}

		// Transfers data from device to host memory and frees device memory
		void memmovDeviceToHost()
		{
			cudaMemcpy(data, deviceBuffer, sizeof(ptrtype) * size, cudaMemcpyDeviceToHost); // Copy data from device to host
			cudaFree(deviceBuffer); // Free the allocated device memory
		}
	};
#endif

#ifdef UseHIP
	// Structure to manage a unified buffer for HIP using managed memory (accessible from both host and device)
	template <typename ptrtype>
	struct bufferGraphUnifiedHIP
	{
		unsigned int size;          // Size of the unified buffer
		ptrtype* data;             // Pointer to managed memory

		// Initializes the unified buffer with specified dimensions using managed memory allocation
		void memoryInit(int dim)
		{
			size = dim;  // Set the size of the buffer
			hipMallocManaged(&data, sizeof(ptrtype) * size); // Allocate unified managed memory
		}

		// Initializes the unified buffer with specified dimensions and fills it with a given value using managed memory allocation
		void memoryInit(int dim, ptrtype v)
		{
			size = dim;  // Set the size of the buffer
			hipMallocManaged(&data, sizeof(ptrtype) * size); // Allocate unified managed memory

			// Fill the allocated buffer with the specified value 
			for (long int i = 0; i < dim; i++)
			{
				data[i] = v;
			}
		}
	};
#endif

#ifdef UseCUDA
	// Structure to manage a unified buffer for CUDA using managed memory (accessible from both host and device)
	template <typename ptrtype>
	struct bufferGraphUnifiedCUDA
	{
		unsigned int size;          // Size of the unified buffer
		ptrtype* data;             // Pointer to managed memory

		// Initializes the unified buffer with specified dimensions using managed memory allocation 
		void memoryInit(int dim)
		{
			size = dim;  // Set the size of the buffer 
			cudaMallocManaged(&data, sizeof(ptrtype) * size); // Allocate unified managed memory 
		}

		// Initializes the unified buffer with specified dimensions and fills it with a given value using managed memory allocation 
		void memoryInit(int dim, ptrtype v)
		{
			size = dim;  // Set the size of the buffer 
			cudaMallocManaged(&data, sizeof(ptrtype) * size); // Allocate unified managed memory 

			// Fill the allocated buffer with the specified value 
			for (long int i = 0; i < dim; i++)
			{
				data[i] = v;
			}
		}
	};
#endif



// GPU KERNEL FUNCTIONS

#if defined( UseHIP ) || defined( UseCUDA )

/**
 * @brief A generic 1D GPU kernel function that applies a kernel operation on input data.
 *
 * @tparam Kernel The type of the kernel function to be executed.
 * @tparam Input The type of input data.
 * @tparam Output The type of output data.
 * @param kernel_function The function to apply to each element.
 * @param n The total number of elements to process.
 * @param in Pointer to the input data array.
 * @param out Pointer to the output data array.
 */
	template <typename Kernel, typename Input, typename Output>
	__global__ void OP_IN_KERNEL_GPU_1D(const Kernel kernel_function, int n, Input* in, Output* out)
	{
		// Calculate the global thread index
		int i = threadIdx.x + blockIdx.x * blockDim.x;

		// Check if the index is within bounds
		if (i < n)
		{
			// Execute the kernel function for the current index
			kernel_function(i, in, out);
		}
	}

	/**
	 * @brief A 1D GPU kernel that applies a lambda operation on an input array.
	 *
	 * @tparam Kernel The type of the lambda operation.
	 * @tparam Input The type of input data.
	 * @param op The lambda operation to apply.
	 * @param A Pointer to the input data array.
	 * @param nb The number of elements to process.
	 */
	template <typename Kernel, typename Input>
	__global__ void OP_IN_KERNEL_LAMBDA_GPU_1D(Kernel op, Input* A, int nb)
	{
		// Calculate the global thread index
		int idx = blockIdx.x * blockDim.x + threadIdx.x;

		// Check if the index is within bounds
		if (idx < nb)
			op(idx, A); // Apply the lambda operation
		//__syncthreads(); // Uncomment for synchronization if needed
	}

	/**
	 * @brief A 1D GPU kernel that applies a lambda operation on three input arrays.
	 *
	 * @tparam Kernel The type of the lambda operation.
	 * @tparam Input The type of input data.
	 * @param op The lambda operation to apply.
	 * @param R Pointer to the first input data array (result).
	 * @param A Pointer to the second input data array.
	 * @param B Pointer to the third input data array.
	 * @param nb The number of elements to process.
	 */
	template <typename Kernel, typename Input>
	__global__ void OP_IN_KERNEL_LAMBDA_GPU_1D_3I(Kernel op, Input* R, Input* A, Input* B, int nb)
	{
		// Calculate the global thread index
		int idx = blockIdx.x * blockDim.x + threadIdx.x;

		// Check if the index is within bounds
		if (idx < nb)
			op(idx, R, A, B); // Apply the lambda operation with three inputs
		//__syncthreads(); // Uncomment for synchronization if needed
	}

	/**
	 * @brief A 1D GPU kernel that applies a lambda operation on an input array within specified bounds.
	 *
	 * @tparam Kernel The type of the lambda operation.
	 * @tparam Input The type of input data.
	 * @param op The lambda operation to apply.
	 * @param A Pointer to the input data array.
	 * @param iBegin Starting index for processing.
	 * @param iEnd Ending index for processing (exclusive).
	 */
	template <typename Kernel, typename Input>
	__global__ void OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D(Kernel op, Input* A, int iBegin, int iEnd)
	{
		// Calculate the global thread index
		int idx = blockIdx.x * blockDim.x + threadIdx.x;

		// Check if the index is within specified bounds
		if ((idx >= iBegin) && (idx < iEnd))
			op(idx, A); // Apply the lambda operation within bounds
		//__syncthreads(); // Uncomment for synchronization if needed
	}

	/**
	 * @brief A 1D GPU kernel that applies a lambda operation on an input array with an offset for indexing.
	 *
	 * @tparam Kernel The type of the lambda operation.
	 * @tparam Input The type of input data.
	 * @param op The lambda operation to apply.
	 * @param A Pointer to the input data array.
	 * @param offset Offset added to calculate global indices for processing elements.
	 */
	template <typename Kernel, typename Input>
	__global__ void OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D(Kernel op, Input* A, int offset)
	{
		// Calculate global thread index with an offset
		int idx = offset + blockIdx.x * blockDim.x + threadIdx.x;

		op(idx, A); // Apply the lambda operation at calculated index
		//__syncthreads(); // Uncomment for synchronization if needed
	}

	/**
	 * @brief A 1D GPU kernel that applies a lambda operation on three input arrays with an offset for indexing.
	 *
	 * @tparam Kernel The type of the lambda operation.
	 * @tparam Input The type of input data.
	 * @param op The lambda operation to apply.
	 * @param R Pointer to the first input data array (result).
	 * @param A Pointer to the second input data array.
	 * @param B Pointer to the third input data array.
	 * @param offset Offset added to calculate global indices for processing elements.
	 */
	template <typename Kernel, typename Input>
	__global__ void OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I(Kernel op, Input* R, Input* A, Input* B, int offset)
	{
		// Calculate global thread index with an offset
		int idx = offset + blockIdx.x * blockDim.x + threadIdx.x;

		op(idx, R, A, B); // Apply the lambda operation at calculated index with three inputs
	}

#endif // End of HIP or CUDA check

// CLASS Task: Provide to manipulate graph Hip/Cuda hybrid system Explicit and Implicit
// GRAPH EXPLICIT

namespace Taskgpu
{

class SingleTask
{
    private:
			// File name for saving task-related data
			std::string M_FileName;

			// Flag indicating whether to save task results
			bool M_qSave;

			// Type of processing: 1 for HIP, 2 for CUDA
			int M_numType;

			// Number of streams to be used for concurrent execution
			int M_nbStreams;

			// Flag for viewing information about the task (for debugging)
			bool M_qViewInfo;

			// Size of each block in GPU execution
			int M_block_size;

			// Model type for execution: 1 for Serial, 2 for Stream
			int M_numModel;

			// Flag indicating whether unified memory is being used
			bool M_qUsedUnifiedMemory;

			// Level of solving complexity (to be used in future implementations)
			int M_solveLevel;

			// Variable to track elapsed time in long format
			long int M_time_laps;

			// Number of GPUs used for execution
			int M_nbGPUused;

			// Time points for measuring execution duration
			std::chrono::steady_clock::time_point M_t_begin, M_t_end;

			// Private method to execute a kernel function serially using HIP
			template <typename Kernel, typename Input>
			void serial_hip(const Kernel& kernel_function,
				int numElems, Input* buffer);

			// Private method to execute a kernel function serially using CUDA
			template <typename Kernel, typename Input>
			void serial_cuda(const Kernel& kernel_function,
				int numElems, Input* buffer);

			// Private method to execute a kernel function using streams with HIP
			template <typename Kernel, typename Input>
			void stream_hip(const Kernel& kernel_function,
				int numElems, Input* buffer);

			// Private method to execute a kernel function using streams with CUDA
			template <typename Kernel, typename Input>
			void stream_cuda(const Kernel& kernel_function,
				int numElems, Input* buffer);

			// Private method for executing a kernel with three inputs serially using HIP
			template <typename Kernel, typename Input>
			void serial_hip_3I(const Kernel& kernel_function,
				int numElems,
				Input* bufferR, Input* bufferA, Input* bufferB);

			// Private method for executing a kernel with three inputs serially using CUDA
			template <typename Kernel, typename Input>
			void serial_cuda_3I(const Kernel& kernel_function,
				int numElems,
				Input* bufferR, Input* bufferA, Input* bufferB);

			// Private method for executing a kernel with three inputs using streams with HIP
			template <typename Kernel, typename Input>
			void stream_hip_3I(const Kernel& kernel_function,
				int numElems,
				Input* bufferR, Input* bufferA, Input* bufferB);

			// Private method for executing a kernel with three inputs using streams with CUDA
			template <typename Kernel, typename Input>
			void stream_cuda_3I(const Kernel& kernel_function,
				int numElems,
				Input* bufferR, Input* bufferA, Input* bufferB);

			// Private method for executing a kernel serially on multiple GPUs using HIP 
			template <typename Kernel, typename Input>
			void serial_hip_multi_GPU(const Kernel& kernel_function,
				int numElems,
				Input* buffer);

			// Private method for executing a kernel serially on multiple GPUs using CUDA 
			template <typename Kernel, typename Input>
			void serial_cuda_multi_GPU(const Kernel& kernel_function,
				int numElems,
				Input* buffer);

			// Private method for executing a kernel with three inputs on multiple GPUs using HIP 
			template <typename Kernel, typename Input>
			void serial_hip_3I_multi_GPU(const Kernel& kernel_function,
				int numElems,
				Input* bufferR, Input* bufferA, Input* bufferB);

			// Private method for executing a kernel with three inputs on multiple GPUs using CUDA 
			template <typename Kernel, typename Input>
			void serial_cuda_3I_multi_GPU(const Kernel& kernel_function,
				int numElems,
				Input* bufferR, Input* bufferA, Input* bufferB);

			// Private method for executing a multi-GPU version 2 kernel function using HIP 
			template <typename Kernel, typename Input>
			void serial_hip_multi_GPU_Vers2(const Kernel& kernel_function,
				int numElems,
				Input* buffer);

			// Private method for executing a multi-GPU version 2 kernel function using CUDA 
			template <typename Kernel, typename Input>
			void serial_cuda_multi_GPU_Vers2(const Kernel& kernel_function,
				int numElems,
				Input* buffer);

		public:

			// Sets the processing type (HIP or CUDA)
			void setNumType(int v) { M_numType = v; }

			// Sets the model type (Serial or Stream)
			void setNumModel(int v) { M_numModel = v; }

			// Sets the number of streams to be used 
			void setNbStreams(int v) { M_nbStreams = v; }

			// Sets the flag to view information about the task (for debugging)
			void setViewInfo(bool b) { M_qViewInfo = b; }

			// Sets the file name for saving task-related data 
			void setFileName(std::string s) { M_FileName = s; }

			// Sets the flag to save task results 
			void setSave(bool b) { M_qSave = b; }

			// Sets the block size for GPU execution 
			void setNbBlock(int v) { M_block_size = v; }

			// Sets the device identifier (HIP or CUDA)
			void setDevice(int v);

			// Sets the flag indicating whether unified memory is being used 
			void setUnifiedMemory(bool b) { M_qUsedUnifiedMemory = b; }

			// Sets the level of solving complexity (to be used in future implementations)
			void setSolveLevel(int v) { M_solveLevel = v; }

			// Sets the number of GPUs used for execution 
			void setNbGPUused(int v) { M_nbGPUused = v; }

			// Returns the processing type (HIP or CUDA)
			int getNumType() const { return (M_numType); }

			// Returns the model type (Serial or Stream)
			int getNumModel() const { return (M_numModel); }

			// Returns the number of streams currently set for execution 
			int getNbStreams() const { return (M_nbStreams); }

			// Returns the level of solving complexity 
			int getSolveLevel() const { return (M_solveLevel); }

			// Returns the number of GPUs currently being used 
			int getNbGPUused() const { return (M_nbGPUused); }

			// Checks if the save flag is set 
			bool isSave() const { return (M_qSave); }

			// Checks if unified memory is being used 
			bool isUnifiedMemory() const { return (M_qUsedUnifiedMemory); }

			// Returns elapsed time in long format since task started or last reset 
			long int getTimeLaps() const { return (M_time_laps); }

			// Constructor: Initializes a SingleTask object with default values
			SingleTask();

			// Destructor: Cleans up resources when a SingleTask object is destroyed
			~SingleTask();

			/**
			   * Runs a specified kernel function on data in a single input buffer.
			   * @tparam Kernel The type of the kernel function.
			   * @tparam Input The type of input data.
			   * @param kernel_function The function to run as a GPU kernel.
			   * @param numElems The number of elements to process.
			   * @param buffer Pointer to input data.
			   */
			template <typename Kernel, typename Input>
			void run(const Kernel& kernel_function,
				int numElems,
				Input* buffer);

			/**
			   * Runs a specified kernel function on data in three input buffers.
			   * @tparam Kernel The type of the kernel function.
			   * @tparam Input The type of input data.
			   * @param kernel_function The function to run as a GPU kernel.
			   * @param numElems The number of elements to process.
			   * @param bufferR Pointer to first input data.
			   * @param bufferA Pointer to second input data.
			   * @param bufferB Pointer to third input data.
			   */
			template <typename Kernel, typename Input>
			void run(const Kernel& kernel_function,
				int numElems,
				Input* bufferR,
				Input* bufferA,
				Input* bufferB);

			// Outputs debugging information about task execution 
			void debriefing();

			// Closes the task and releases resources 
			void close();
};

SingleTask::SingleTask()
{
    M_nbStreams = 3;
    M_numType = 1; // 1:HIP 2:CUDA
    M_qViewInfo = true;
    M_time_laps = 0;
    M_FileName = "NoName";
    M_qSave = false;
    M_block_size = 512;
    M_numModel = 1; // 1:Serial 2:Stream
    M_qUsedUnifiedMemory = false;
    M_solveLevel = 1;
    M_nbGPUused = 1;
}

SingleTask::~SingleTask()
{
}

void SingleTask::debriefing()
{
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Debriefing"
                  << "\n";
        VLOG( 1 ) << "Elapsed microseconds = " << M_time_laps << "us\n";
    }
    if ( M_qSave )
    {
        if ( M_qViewInfo )
        {
            VLOG( 1 ) << "Save Informations"
                      << "\n";
        }
        std::ofstream myfile;
        myfile.open( M_FileName + ".csv" );
        myfile << "Elapsed microseconds," << M_time_laps << "\n";
        myfile << "\n";
        myfile.close();
    }
}

void SingleTask::close()
{
}

void SingleTask::setDevice( int v )
{
    if ( M_numType == 1 )
    {
#ifdef UseHIP
        int numDevices = 0;
        hipGetDeviceCount( &numDevices );
        if ( ( v >= 0 ) && ( v <= numDevices ) )
        {
            hipSetDevice( v );
        }
#endif
    }

    if ( M_numType == 2 )
    {
#ifdef UseCUDA
        int numDevices = 0;
        cudaGetDeviceCount( &numDevices );
        if ( ( v >= 0 ) && ( v <= numDevices ) )
        {
            cudaSetDevice( v );
        }
#endif
    }
}

template <typename Kernel, typename Input>
void SingleTask::run( const Kernel& kernel_function,
                      int numElems,
                      Input* buffer )
{
    M_t_begin = std::chrono::steady_clock::now();
    if ( M_nbGPUused == 1 )
    {
        if ( M_numModel == 1 )
        { // SERIAL MODEL
            if ( M_numType == 1 )
            {
#ifdef UseHIP
                serial_hip( kernel_function, numElems, buffer );
#endif
            }

            if ( M_numType == 2 )
            {
#ifdef UseCUDA
                serial_cuda( kernel_function, numElems, buffer );
#endif
            }
        }

        if ( M_numModel == 2 )
        { // STREAM MODEL
            if ( M_numType == 1 )
            {
#ifdef UseHIP
                stream_hip( kernel_function, numElems, buffer );
#endif
            }

            if ( M_numType == 2 )
            {
#ifdef UseCUDA
                stream_cuda( kernel_function, numElems, buffer );
#endif
            }
        }
    }
    else
    {
        // SERIAL MODEL MULTI GPU
        if ( M_numModel == 1 )
        { // SERIAL MODEL
            if ( M_numType == 1 )
            {
#ifdef UseHIP
                serial_hip_multi_GPU( kernel_function, numElems, buffer );
#endif
            }

            if ( M_numType == 2 )
            {
#ifdef UseCUDA
                serial_cuda_multi_GPU( kernel_function, numElems, buffer );
#endif
            }
        }
    }

    M_t_end = std::chrono::steady_clock::now();
    M_time_laps = std::chrono::duration_cast<std::chrono::microseconds>( M_t_end - M_t_begin ).count();
}

template <typename Kernel, typename Input>
void SingleTask::serial_hip( const Kernel& kernel_function,
                             int numElems,
                             Input* buffer )
{
#ifdef UseHIP
    int sz = sizeof( decltype( buffer ) ) * numElems;
    int iEnd = numElems;
    int iBegin = 0;
    int num_blocks = ( numElems + M_block_size - 1 ) / M_block_size;
    dim3 thread_block( M_block_size, 1, 1 );
    dim3 grid( num_blocks, 1 );
    if ( M_qUsedUnifiedMemory )
    {
        hipLaunchKernelGGL( OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D, grid, thread_block, 0, 0, kernel_function, buffer, iBegin, iEnd );
        hipDeviceSynchronize();
    }
    else
    {
        Input* deviceBuffer;
        hipMalloc( (void**)&deviceBuffer, sz );
        hipMemcpy( deviceBuffer, buffer, sz, hipMemcpyHostToDevice );
        hipLaunchKernelGGL( OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D, grid, thread_block, 0, 0, kernel_function, deviceBuffer, iBegin, iEnd );
        hipMemcpy( buffer, deviceBuffer, sz, hipMemcpyDeviceToHost );
        hipFree( deviceBuffer );
    }
#endif
}

template <typename Kernel, typename Input>
void SingleTask::serial_cuda( const Kernel& kernel_function,
                              int numElems,
                              Input* buffer )
{
#ifdef UseCUDA
    int sz = sizeof( decltype( buffer ) ) * numElems;
    int iEnd = numElems;
    int iBegin = 0;
    int num_blocks = ( numElems + M_block_size - 1 ) / M_block_size;
    dim3 thread_block( M_block_size, 1, 1 );
    dim3 grid( num_blocks, 1 );
    if ( M_qUsedUnifiedMemory )
    {
        OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<<<gridthread_block>>>( kernel_function, buffer, iBegin, iEnd );
        cudaDeviceSynchronize();
    }
    else
    {
        Input* deviceBuffer;
        cudaMalloc( (void**)&deviceBuffer, sz );
        cudaMemcpy( deviceBuffer, buffer, sz, hipMemcpyHostToDevice );
        OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<<<gridthread_block>>>( kernel_function, deviceBuffer, iBegin, iEnd );
        cudaMemcpy( buffer, deviceBuffer, sz, hipMemcpyDeviceToHost );
        cudaFree( deviceBuffer );
    }
#endif
}

template <typename Kernel, typename Input>
void SingleTask::run( const Kernel& kernel_function,
                      int numElems,
                      Input* bufferR, Input* bufferA, Input* bufferB )
{
    M_t_begin = std::chrono::steady_clock::now();
    if ( M_nbGPUused == 1 )
    {
        if ( M_numModel == 1 )
        { // SERIAL MODEL
            if ( M_numType == 1 )
            {
#ifdef UseHIP
                serial_hip_3I( kernel_function, numElems, bufferR, bufferA, bufferB );
#endif
            }

            if ( M_numType == 2 )
            {
#ifdef UseCUDA
                serial_cuda_3I( kernel_function, numElems, bufferR, bufferA, bufferB );
#endif
            }
        }

        if ( M_numModel == 2 )
        { // STREAM MODEL
            if ( M_numType == 1 )
            {
#ifdef UseHIP
                stream_hip_3I( kernel_function, numElems, bufferR, bufferA, bufferB );
#endif
            }

            if ( M_numType == 2 )
            {
#ifdef UseCUDA
                stream_cuda_3I( kernel_function, numElems, bufferR, bufferA, bufferB );
#endif
            }
        }
    }
    else
    {
        if ( M_numModel == 1 )
        { // SERIAL MODEL MULTI GPU
            if ( M_numType == 1 )
            {
#ifdef UseHIP
                serial_hip_3I_multi_GPU( kernel_function, numElems, bufferR, bufferA, bufferB );
#endif
            }

            if ( M_numType == 2 )
            {
#ifdef UseCUDA
                serial_cuda_3I_multi_GPU( kernel_function, numElems, bufferR, bufferA, bufferB );
#endif
            }
        }
    }

    M_t_end = std::chrono::steady_clock::now();
    M_time_laps = std::chrono::duration_cast<std::chrono::microseconds>( M_t_end - M_t_begin ).count();
}

template <typename Kernel, typename Input>
void SingleTask::serial_hip_multi_GPU( const Kernel& kernel_function,
                                       int numElems,
                                       Input* buffer )
{
#ifdef UseHIP
    int numDevices = 0;
    hipGetDeviceCount( &numDevices );
    if ( M_nbGPUused > numDevices )
    {
        M_nbGPUused = numDevices;
    }
    int numElemsPerBlock = int( numElems / M_nbGPUused );
    int numElemsPerBlockRest = numElems - M_nbGPUused * numElemsPerBlock;

    for ( int i = 0; i < M_nbGPUused; i++ )
    {
        hipSetDevice( i );
        Input* d_B;
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        hipMalloc( (void**)&d_B, nbElemsPartiel );
        hipMemcpy( d_B, &buffer[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( buffer ) ), hipMemcpyHostToDevice );
        int num_blocks = ( nbElemsPartiel + M_block_size - 1 ) / M_block_size;
        dim3 thread_block( M_block_size, 1, 1 );
        dim3 grid( num_blocks, 1 );
        hipLaunchKernelGGL( OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D, grid, thread_block, 0, 0, kernel_function, d_B, 0, nbElemsPartiel );
        hipMemcpy( &buffer[i * numElemsPerBlock], d_B, nbElemsPartiel * sizeof( decltype( buffer ) ), hipMemcpyDeviceToHost );
        hipDeviceSynchronize();
        hipFree( d_B );
    }

#endif
}

template <typename Kernel, typename Input>
void SingleTask::serial_cuda_multi_GPU( const Kernel& kernel_function,
                                        int numElems,
                                        Input* buffer )
{
#ifdef UseCUDA
    int numDevices = 0;
    cudaGetDeviceCount( &numDevices );
    if ( M_nbGPUused > numDevices )
    {
        M_nbGPUused = numDevices;
    }
                int numElemsPerBlock=int(numElems/(M_nbGPUused);
		int numElemsPerBlockRest=numElems-M_nbGPUused*numElemsPerBlock;
		
        for (int i = 0; i < M_nbGPUused; i++) 
        {
        cudaSetDevice( i );
        Input* d_B;
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        cudaMalloc( (void**)&d_B, nbElemsPartiel );
        cudaMemcpy( d_B, &buffer[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( buffer ) ), cudaMemcpyHostToDevice );
        int num_blocks = ( nbElemsPartiel + M_block_size - 1 ) / M_block_size;
        dim3 thread_block( M_block_size, 1, 1 );
        dim3 grid( num_blocks, 1 );
        OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<<<grid, thread_block, 0, 0>>>( kernel_function, d_B, 0, nbElemsPartiel );
        cudaMemcpy( &buffer[i * numElemsPerBlock], d_B, nbElemsPartiel * sizeof( decltype( buffer ) ), cudaMemcpyDeviceToHost );
        cudaDeviceSynchronize();
        cudaFree( d_B );
		}
#endif
}

template <typename Kernel, typename Input>
void SingleTask::serial_hip_multi_GPU_Vers2( const Kernel& kernel_function,
                                             int numElems,
                                             Input* buffer )
{
#ifdef UseHIP
    int numDevices = 0;
    hipGetDeviceCount( &numDevices );
    if ( M_nbGPUused > numDevices )
    {
        M_nbGPUused = numDevices;
    }
    int numElemsPerBlock = int( numElems / ( M_nbGPUused ) );
    int numElemsPerBlockRest = numElems - M_nbGPUused * numElemsPerBlock;

    std::vector<Input> vdb;
    for ( int i = 0; i < M_nbGPUused; i++ )
    {
        hipSetDevice( i );
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        Input* d_B;
        vdb.push_back( std::move( d_B ) );
        hipMalloc( (void**)vdb[i], nbElemsPartiel );
        hipMemcpy( vdb[i], &buffer[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( buffer ) ), hipMemcpyHostToDevice );
    }

    for ( int i = 0; i < M_nbGPUused; i++ )
    {
        hipSetDevice( i );
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        int num_blocks = ( nbElemsPartiel + M_block_size - 1 ) / M_block_size;
        dim3 thread_block( M_block_size, 1, 1 );
        dim3 grid( num_blocks, 1 );
        hipLaunchKernelGGL( OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D, grid, thread_block, 0, 0, kernel_function, vdb[i], 0, nbElemsPartiel );
    }

    hipDeviceSynchronize();

    for ( int i = 0; i < M_nbGPUused; i++ )
    {
        hipSetDevice( i );
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        hipMemcpy( &buffer[i * numElemsPerBlock], vdb[i], nbElemsPartiel * sizeof( decltype( buffer ) ), hipMemcpyDeviceToHost );
        hipFree( vdb[i] );
    }
    vdb.clear();
#endif
}

template <typename Kernel, typename Input>
void SingleTask::serial_cuda_multi_GPU_Vers2( const Kernel& kernel_function,
                                              int numElems,
                                              Input* buffer )
{
#ifdef UseCUDA
    int numDevices = 0;
    cudaGetDeviceCount( &numDevices );
    if ( M_nbGPUused > numDevices )
    {
        M_nbGPUused = numDevices;
    }
    int numElemsPerBlock = int( numElems / ( M_nbGPUused ) );
    int numElemsPerBlockRest = numElems - M_nbGPUused * numElemsPerBlock;

    std::vector<Input> vdb;
    for ( int i = 0; i < M_nbGPUused; i++ )
    {
        cudaSetDevice( i );
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        Input* d_B;
        vdb.push_back( std::move( d_B ) );
        cudaMalloc( (void**)vdb[i], nbElemsPartiel );
        cudaMemcpy( vdb[i], &buffer[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( buffer ) ), cudaMemcpyHostToDevice );
    }

    for ( int i = 0; i < M_nbGPUused; i++ )
    {
        cudaSetDevice( i );
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        int num_blocks = ( nbElemsPartiel + M_block_size - 1 ) / M_block_size;
        dim3 thread_block( M_block_size, 1, 1 );
        dim3 grid( num_blocks, 1 );
        OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<<<grid, thread_block, 0, 0>>>( kernel_function, vdb[i], 0, nbElemsPartiel );
    }

    cudaDeviceSynchronize();

    for ( int i = 0; i < M_nbGPUused; i++ )
    {
        cudaSetDevice( i );
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        cudaMemcpy( &buffer[i * numElemsPerBlock], vdb[i], nbElemsPartiel * sizeof( decltype( buffer ) ), cudaMemcpyDeviceToHost );
        cudaFree( vdb[i] );
    }
    vdb.clear();
#endif
}

template <typename Kernel, typename Input>
void SingleTask::serial_hip_3I( const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR, Input* bufferA, Input* bufferB )
{
#ifdef UseHIP
    int sz = sizeof( decltype( bufferR ) ) * numElems;
    int iEnd = numElems;
    int iBegin = 0;
    int num_blocks = ( numElems + M_block_size - 1 ) / M_block_size;
    dim3 thread_block( M_block_size, 1, 1 );
    dim3 grid( num_blocks, 1 );
    if ( M_qUsedUnifiedMemory )
    {
        hipLaunchKernelGGL( OP_IN_KERNEL_LAMBDA_GPU_1D_3I, grid, thread_block, 0, 0, kernel_function,
                            bufferR, bufferA, bufferB, iBegin, iEnd );
        hipDeviceSynchronize();
    }
    else
    {
        Input *deviceBufferR, *deviceBufferA, *deviceBufferB;
        hipMalloc( (void**)&deviceBufferR, sz );
        hipMalloc( (void**)&deviceBufferA, sz );
        hipMalloc( (void**)&deviceBufferB, sz );
        hipMemcpy( deviceBufferR, bufferR, sz, hipMemcpyHostToDevice );
        hipMemcpy( deviceBufferA, bufferA, sz, hipMemcpyHostToDevice );
        hipMemcpy( deviceBufferB, bufferB, sz, hipMemcpyHostToDevice );
        hipLaunchKernelGGL( OP_IN_KERNEL_LAMBDA_GPU_1D_3I, grid, thread_block, 0, 0, kernel_function,
                            deviceBufferR, deviceBufferA, deviceBufferB, iBegin, iEnd );
        hipMemcpy( bufferR, deviceBufferR, sz, hipMemcpyDeviceToHost );
        hipFree( deviceBufferR );
        hipFree( deviceBufferA );
        hipFree( deviceBufferB );
    }
#endif
}

template <typename Kernel, typename Input>
void SingleTask::serial_hip_3I_multi_GPU( const Kernel& kernel_function,
                                          int numElems,
                                          Input* bufferR, Input* bufferA, Input* bufferB )
{
#ifdef UseHIP
    int numDevices = 0;
    hipGetDeviceCount( &numDevices );
    if ( M_nbGPUused > numDevices )
    {
        M_nbGPUused = numDevices;
    }
    int numElemsPerBlock = int( numElems / M_nbGPUused );
    int numElemsPerBlockRest = numElems - M_nbGPUused * numElemsPerBlock;

    for ( int i = 0; i < M_nbGPUused; i++ )
    {
        hipSetDevice( i );
        Input *d_B_R, *d_B_A, *d_B_B;
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        hipMalloc( (void**)&d_B_R, nbElemsPartiel );
        hipMalloc( (void**)&d_B_A, nbElemsPartiel );
        hipMalloc( (void**)&d_B_B, nbElemsPartiel );
        hipMemcpy( d_B_R, &bufferR[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( bufferR ) ), hipMemcpyHostToDevice );
        hipMemcpy( d_B_A, &bufferA[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( bufferR ) ), hipMemcpyHostToDevice );
        hipMemcpy( d_B_B, &bufferB[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( bufferR ) ), hipMemcpyHostToDevice );

        int num_blocks = ( nbElemsPartiel + M_block_size - 1 ) / M_block_size;
        dim3 thread_block( M_block_size, 1, 1 );
        dim3 grid( num_blocks, 1 );
        hipLaunchKernelGGL( OP_IN_KERNEL_LAMBDA_GPU_1D_3I, grid, thread_block, 0, 0, kernel_function,
                            d_B_R, d_B_A, d_B_B, 0, nbElemsPartiel );

        hipMemcpy( &bufferR[i * numElemsPerBlock], d_B_R, nbElemsPartiel * sizeof( decltype( bufferR ) ), hipMemcpyDeviceToHost );
        hipDeviceSynchronize();
        hipFree( d_B_R );
        hipFree( d_B_A );
        hipFree( d_B_B );
    }

#endif
}

template <typename Kernel, typename Input>
void SingleTask::serial_cuda_3I( const Kernel& kernel_function,
                                 int numElems,
                                 Input* bufferR, Input* bufferA, Input* bufferB )
{
#ifdef UseCUDA
    int sz = sizeof( decltype( bufferR ) ) * numElems;
    int iEnd = numElems;
    int iBegin = 0;
    int num_blocks = ( numElems + M_block_size - 1 ) / M_block_size;
    dim3 thread_block( M_block_size, 1, 1 );
    dim3 grid( num_blocks, 1 );
    if ( M_qUsedUnifiedMemory )
    {
        OP_IN_KERNEL_LAMBDA_GPU_1D_3I<<<gridthread_block>>>( kernel_function,
                                                             bufferR, bufferA, bufferB, iBegin, iEnd );
        cudaDeviceSynchronize();
    }
    else
    {
        Input *deviceBufferR, *deviceBufferA, *deviceBufferB;
        cudaMalloc( (void**)&deviceBufferR, sz );
        cudaMalloc( (void**)&deviceBufferA, sz );
        cudaMalloc( (void**)&deviceBufferB, sz );
        cudaMemcpy( deviceBufferR, bufferR, sz, cudaMemcpyHostToDevice );
        cudaMemcpy( deviceBufferA, bufferA, sz, cudaMemcpyHostToDevice );
        cudaMemcpy( deviceBufferB, bufferB, sz, cudaMemcpyHostToDevice );
        OP_IN_KERNEL_LAMBDA_GPU_1D_3I<<<gridthread_block>>>( kernel_function,
                                                             deviceBufferR, deviceBufferA, deviceBufferB, iBegin, iEnd );
        cudaMemcpy( bufferR, deviceBufferR, sz, cudaMemcpyDeviceToHost );
        cudaFree( deviceBufferR );
        cudaFree( deviceBufferA );
        cudaFree( deviceBufferB );
    }
#endif
}

template <typename Kernel, typename Input>
void SingleTask::serial_cuda_3I_multi_GPU( const Kernel& kernel_function,
                                           int numElems,
                                           Input* bufferR, Input* bufferA, Input* bufferB )
{
#ifdef UseCUDA
    int numDevices = 0;
    cudaGetDeviceCount( &numDevices );
    if ( M_nbGPUused > numDevices )
    {
        M_nbGPUused = numDevices;
    }
    int numElemsPerBlock = int( numElems / M_nbGPUused );
    int numElemsPerBlockRest = numElems - M_nbGPUused * numElemsPerBlock;

    for ( int i = 0; i < M_nbGPUused; i++ )
    {
        cudaSetDevice( i );
        Input *d_B_R, *d_B_A, *d_B_B;
        int nbElemsPartiel = numElemsPerBlock;
        if ( i == M_nbGPUused - 1 )
        {
            nbElemsPartiel = numElemsPerBlock + numElemsPerBlockRest;
        }
        cudaMalloc( (void**)&d_B_R, nbElemsPartiel );
        cudaMalloc( (void**)&d_B_A, nbElemsPartiel );
        cudaMalloc( (void**)&d_B_B, nbElemsPartiel );
        cudaMemcpy( d_B_R, &bufferR[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( bufferR ) ), cudaMemcpyHostToDevice );
        cudaMemcpy( d_B_A, &bufferA[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( bufferR ) ), cudaMemcpyHostToDevice );
        cudaMemcpy( d_B_B, &bufferB[i * numElemsPerBlock], nbElemsPartiel * sizeof( decltype( bufferR ) ), cudaMemcpyHostToDevice );

        int num_blocks = ( nbElemsPartiel + M_block_size - 1 ) / M_block_size;
        dim3 thread_block( M_block_size, 1, 1 );
        dim3 grid( num_blocks, 1 );
        OP_IN_KERNEL_LAMBDA_GPU_1D_3I<<<grid, thread_block, 0, 0>>>( kernel_function,
                                                                     d_B_R, d_B_A, d_B_B, 0, nbElemsPartiel );

        cudaMemcpy( &bufferR[i * numElemsPerBlock], d_B_R, nbElemsPartiel * sizeof( decltype( bufferR ) ), cudaMemcpyDeviceToHost );
        cudaDeviceSynchronize();
        cudaFree( d_B_R );
        cudaFree( d_B_A );
        cudaFree( d_B_B );
    }

#endif
}

template <typename Kernel, typename Input>
void SingleTask::stream_hip( const Kernel& kernel_function,
                             int numElems,
                             Input* buffer )
{
#ifdef UseHIP
    int numVersion = 1;
    float ms;
    int i;

    // const int blockSize = 20;  //const int blockSize M_block_size;
    const int blockSize = M_block_size;

    const int streamSize = numElems / M_nbStreams;
    const int streamBytes = streamSize * sizeof( decltype( buffer ) );
    const int sz = numElems * sizeof( decltype( buffer ) );
    int grid = streamSize / blockSize;
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Block size :" << M_block_size << "\n";
        VLOG( 1 ) << "nb Streams :" << M_nbStreams << "\n";
    }
    if ( ( grid == 0 ) && ( M_qViewInfo ) )
    {
        VLOG( 1 ) << "Grid size error : ==> auto correction"
                  << "\n";
    }
    grid = max( streamSize / blockSize, 1 );
    Input* d_buffer;
    hipMalloc( (void**)&d_buffer, sz );
    hipEvent_t startEvent, stopEvent, dummyEvent;
    hipStream_t stream[M_nbStreams];
    hipEventCreate( &startEvent );
    hipEventCreate( &stopEvent );
    hipEventCreate( &dummyEvent );
    for ( i = 0; i < M_nbStreams; ++i )
    {
        hipStreamCreate( &stream[i] );
    }
    hipEventRecord( startEvent, 0 );

    if ( numVersion == 1 )
    {
        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            hipMemcpyAsync( &d_buffer[offset], &buffer[offset], streamBytes, hipMemcpyHostToDevice, stream[i] );
            hipLaunchKernelGGL( OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D, grid, blockSize, 0, stream[i], kernel_function, d_buffer, offset );
            hipMemcpyAsync( &buffer[offset], &d_buffer[offset], streamBytes, hipMemcpyDeviceToHost, stream[i] );
        }
        hipEventRecord( stopEvent, 0 );
        hipEventSynchronize( stopEvent );
        // hipEventElapsedTime(&ms, startEvent, stopEvent);
        // printf("Time (ms): %f\n", ms);
    }

    if ( numVersion == 2 )
    {
        hipEventRecord( startEvent, 0 );
        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            hipMemcpyAsync( &d_buffer[offset], &buffer[offset], streamBytes, hipMemcpyHostToDevice, stream[i] );
        }

        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            hipLaunchKernelGGL( OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D, grid, blockSize, 0, stream[i], kernel_function, d_buffer, offset );
        }

        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            hipMemcpyAsync( &buffer[offset], &d_buffer[offset], streamBytes, hipMemcpyDeviceToHost, stream[i] );
        }
        hipEventRecord( stopEvent, 0 );
        hipEventSynchronize( stopEvent );
        // hipEventElapsedTime(&ms, startEvent, stopEvent);
        // printf("Time (ms): %f\n", ms);
    }

    // cleanup
    hipEventDestroy( startEvent );
    hipEventDestroy( stopEvent );
    hipEventDestroy( dummyEvent );
    for ( int i = 0; i < M_nbStreams; ++i )
    {
        hipStreamDestroy( stream[i] );
    }
    hipFree( d_buffer );
#endif
}

template <typename Kernel, typename Input>
void SingleTask::stream_cuda( const Kernel& kernel_function,
                              int numElems,
                              Input* buffer )
{
#ifdef UseCUDA
    int numVersion = 1;
    float ms;
    int i;

    // const int blockSize = 20;  //const int blockSize M_block_size;
    const int blockSize = M_block_size;

    const int streamSize = numElems / M_nbStreams;
    const int streamBytes = streamSize * sizeof( decltype( buffer ) );
    const int sz = numElems * sizeof( decltype( buffer ) );
    int grid = streamSize / blockSize;
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Block size :" << M_block_size << "\n";
        VLOG( 1 ) << "nb Streams :" << M_nbStreams << "\n";
    }
    if ( ( grid == 0 ) && ( M_qViewInfo ) )
    {
        VLOG( 1 ) << "Grid size error : ==> auto correction"
                  << "\n";
    }
    grid = max( streamSize / blockSize, 1 );
    Input* d_buffer;
    cudaMalloc( (void**)&d_buffer, sz );
    cudaEvent_t startEvent, stopEvent, dummyEvent;
    cudaStream_t stream[M_nbStreams];
    cudaEventCreate( &startEvent );
    cudaEventCreate( &stopEvent );
    cudaEventCreate( &dummyEvent );
    for ( i = 0; i < M_nbStreams; ++i )
    {
        cudaStreamCreate( &stream[i] );
    }
    cudaEventRecord( startEvent, 0 );

    if ( numVersion == 1 )
    {
        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            cudaMemcpyAsync( &d_buffer[offset], &buffer[offset], streamBytes, cudaMemcpyHostToDevice, stream[i] );
            OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D<<<grid, blockSize, 0, stream[i]>>>( kernel_function, d_buffer, offset );
            cudaMemcpyAsync( &buffer[offset], &d_buffer[offset], streamBytes, cudaMemcpyDeviceToHost, stream[i] );
        }
        cudaEventRecord( stopEvent, 0 );
        cudaEventSynchronize( stopEvent );
        // cudaEventElapsedTime(&ms, startEvent, stopEvent);
        // printf("Time (ms): %f\n", ms);
    }

    if ( numVersion == 2 )
    {
        cudaEventRecord( startEvent, 0 );
        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            cudaMemcpyAsync( &d_buffer[offset], &buffer[offset], streamBytes, cudaMemcpyHostToDevice, stream[i] );
        }

        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D<<<grid, blockSize, 0, stream[i]>>>( kernel_function, d_buffer, offset );
        }

        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            cudaMemcpyAsync( &buffer[offset], &d_buffer[offset], streamBytes, cudaMemcpyDeviceToHost, stream[i] );
        }
        cudaEventRecord( stopEvent, 0 );
        cudaEventSynchronize( stopEvent );
        // cudaEventElapsedTime(&ms, startEvent, stopEvent);
        // printf("Time (ms): %f\n", ms);
    }

    // cleanup
    cudaEventDestroy( startEvent );
    cudaEventDestroy( stopEvent );
    cudaEventDestroy( dummyEvent );
    for ( int i = 0; i < M_nbStreams; ++i )
    {
        cudaStreamDestroy( stream[i] );
    }
    cudaFree( d_buffer );

#endif
}

template <typename Kernel, typename Input>
void SingleTask::stream_hip_3I( const Kernel& kernel_function,
                                int numElems,
                                Input* bufferR, Input* bufferA, Input* bufferB )
{
#ifdef UseHIP
    int numVersion = 1;
    float ms;
    int i;

    // const int blockSize = 20;  //const int blockSize M_block_size;
    const int blockSize = M_block_size;

    const int streamSize = numElems / M_nbStreams;
    const int streamBytes = streamSize * sizeof( decltype( bufferR ) );
    const int sz = numElems * sizeof( decltype( bufferR ) );
    int grid = streamSize / blockSize;
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Block size :" << M_block_size << "\n";
        VLOG( 1 ) << "nb Streams :" << M_nbStreams << "\n";
    }
    if ( ( grid == 0 ) && ( M_qViewInfo ) )
    {
        VLOG( 1 ) << "Grid size error : ==> auto correction"
                  << "\n";
    }
    grid = max( streamSize / blockSize, 1 );
    Input* d_bufferR;
    hipMalloc( (void**)&d_bufferR, sz );
    Input* d_bufferA;
    hipMalloc( (void**)&d_bufferA, sz );
    Input* d_bufferB;
    hipMalloc( (void**)&d_bufferB, sz );
    hipEvent_t startEvent, stopEvent, dummyEvent;
    hipStream_t stream[M_nbStreams];
    hipEventCreate( &startEvent );
    hipEventCreate( &stopEvent );
    hipEventCreate( &dummyEvent );
    for ( i = 0; i < M_nbStreams; ++i )
    {
        hipStreamCreate( &stream[i] );
    }
    hipEventRecord( startEvent, 0 );

    if ( numVersion == 1 )
    {
        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            hipMemcpyAsync( &d_bufferR[offset], &bufferR[offset], streamBytes, hipMemcpyHostToDevice, stream[i] );
            hipMemcpyAsync( &d_bufferA[offset], &bufferA[offset], streamBytes, hipMemcpyHostToDevice, stream[i] );
            hipMemcpyAsync( &d_bufferB[offset], &bufferB[offset], streamBytes, hipMemcpyHostToDevice, stream[i] );
            hipLaunchKernelGGL( OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I, grid, blockSize, 0, stream[i], kernel_function,
                                d_bufferR, d_bufferA, d_bufferB, offset );
            hipMemcpyAsync( &bufferR[offset], &d_bufferR[offset], streamBytes, hipMemcpyDeviceToHost, stream[i] );
        }
        hipEventRecord( stopEvent, 0 );
        hipEventSynchronize( stopEvent );
    }

    if ( numVersion == 2 )
    {
        hipEventRecord( startEvent, 0 );
        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            hipMemcpyAsync( &d_bufferR[offset], &bufferR[offset], streamBytes, hipMemcpyHostToDevice, stream[i] );
            hipMemcpyAsync( &d_bufferA[offset], &bufferA[offset], streamBytes, hipMemcpyHostToDevice, stream[i] );
            hipMemcpyAsync( &d_bufferB[offset], &bufferB[offset], streamBytes, hipMemcpyHostToDevice, stream[i] );
        }

        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            hipLaunchKernelGGL( OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I, grid, blockSize, 0, stream[i], kernel_function,
                                d_bufferR, d_bufferA, d_bufferB, offset );
        }

        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            hipMemcpyAsync( &bufferR[offset], &d_bufferR[offset], streamBytes, hipMemcpyDeviceToHost, stream[i] );
        }
        hipEventRecord( stopEvent, 0 );
        hipEventSynchronize( stopEvent );
        // hipEventElapsedTime(&ms, startEvent, stopEvent);
        // printf("Time (ms): %f\n", ms);
    }
    // cleanup
    hipEventDestroy( startEvent );
    hipEventDestroy( stopEvent );
    hipEventDestroy( dummyEvent );
    for ( int i = 0; i < M_nbStreams; ++i )
    {
        hipStreamDestroy( stream[i] );
    }
    hipFree( d_bufferR );
    hipFree( d_bufferA );
    hipFree( d_bufferB );
#endif
}

template <typename Kernel, typename Input>
void SingleTask::stream_cuda_3I( const Kernel& kernel_function,
                                 int numElems,
                                 Input* bufferR, Input* bufferA, Input* bufferB )
{
#ifdef UseCUDA
    int numVersion = 1;
    float ms;
    int i;

    // const int blockSize = 20;  //const int blockSize M_block_size;
    const int blockSize = M_block_size;

    const int streamSize = numElems / M_nbStreams;
    const int streamBytes = streamSize * sizeof( decltype( bufferR ) );
    const int sz = numElems * sizeof( decltype( bufferR ) );
    int grid = streamSize / blockSize;
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Block size :" << M_block_size << "\n";
        VLOG( 1 ) << "nb Streams :" << M_nbStreams << "\n";
    }
    if ( ( grid == 0 ) && ( M_qViewInfo ) )
    {
        VLOG( 1 ) << "Grid size error : ==> auto correction"
                  << "\n";
    }
    grid = max( streamSize / blockSize, 1 );
    Input* d_bufferR;
    cudaMalloc( (void**)&d_bufferR, sz );
    Input* d_bufferA;
    cudaMalloc( (void**)&d_bufferA, sz );
    Input* d_bufferB;
    cudaMalloc( (void**)&d_bufferB, sz );
    cudaEvent_t startEvent, stopEvent, dummyEvent;
    cudaStream_t stream[M_nbStreams];
    cudaEventCreate( &startEvent );
    cudaEventCreate( &stopEvent );
    cudaEventCreate( &dummyEvent );
    for ( i = 0; i < M_nbStreams; ++i )
    {
        cudaStreamCreate( &stream[i] );
    }
    cudaEventRecord( startEvent, 0 );

    if ( numVersion == 1 )
    {
        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            cudaMemcpyAsync( &d_bufferR[offset], &bufferR[offset], streamBytes, cudaMemcpyHostToDevice, stream[i] );
            cudaMemcpyAsync( &d_bufferA[offset], &bufferA[offset], streamBytes, cudaMemcpyHostToDevice, stream[i] );
            cudaMemcpyAsync( &d_bufferB[offset], &bufferB[offset], streamBytes, cudaMemcpyHostToDevice, stream[i] );
            OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I<<<grid, blockSize, 0, stream[i]>>>( kernel_function,
                                                                                     d_bufferR, d_bufferA, d_bufferB, offset );
            cudaMemcpyAsync( &bufferR[offset], &d_bufferR[offset], streamBytes, cudaMemcpyDeviceToHost, stream[i] );
        }
        cudaEventRecord( stopEvent, 0 );
        cudaEventSynchronize( stopEvent );
    }

    if ( numVersion == 2 )
    {
        cudaEventRecord( startEvent, 0 );
        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            cudaMemcpyAsync( &d_bufferR[offset], &bufferR[offset], streamBytes, cudaMemcpyHostToDevice, stream[i] );
            cudaMemcpyAsync( &d_bufferA[offset], &bufferA[offset], streamBytes, cudaMemcpyHostToDevice, stream[i] );
            cudaMemcpyAsync( &d_bufferB[offset], &bufferB[offset], streamBytes, cudaMemcpyHostToDevice, stream[i] );
        }

        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            OP_IN_KERNEL_LAMBDA_STREAM_GPU_1D_3I<<<grid, blockSize, 0, stream[i]>>>( kernel_function,
                                                                                     d_bufferR, d_bufferA, d_bufferB, offset );
        }

        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            cudaMemcpyAsync( &bufferR[offset], &d_bufferR[offset], streamBytes, cudaMemcpyDeviceToHost, stream[i] );
        }
        cudaEventRecord( stopEvent, 0 );
        cudaEventSynchronize( stopEvent );
        // cudaEventElapsedTime(&ms, startEvent, stopEvent);
        // printf("Time (ms): %f\n", ms);
    }
    // cleanup
    cudaEventDestroy( startEvent );
    cudaEventDestroy( stopEvent );
    cudaEventDestroy( dummyEvent );
    for ( int i = 0; i < M_nbStreams; ++i )
    {
        cudaStreamDestroy( stream[i] );
    }
    cudaFree( d_bufferR );
    cudaFree( d_bufferA );
    cudaFree( d_bufferB );
#endif
}

//--------------------------------------------------------------------------------------------------------------------------------

class Task
{
    private:
  		// File name for saving task-related data
		std::string M_FileName;

		// Number of threads to be used in the task
		int M_nbTh;

		// Number of blocks to be used on the GPU
		int M_numBlocksGPU;

		// Number of threads per block on the GPU
		int M_nThPerBckGPU;

		// Flag indicating whether a graph has been created
		bool M_q_graph;

		// Flag for viewing information about the task (for debugging)
		bool M_qViewInfo;

		// Flag indicating whether to save task results
		bool M_qSave;

		// Variable to track elapsed time in long format
		long int M_time_laps;

		// Flag indicating if the device needs to be reset
		bool M_qDeviceReset;

		// Type identifier for the task (e.g., different processing types)
		int M_numType;

		// List of graph dependencies for managing execution order
		std::vector<int> M_ListGraphDependencies;

		// Time points for measuring execution duration
		std::chrono::steady_clock::time_point M_t_begin, M_t_end;

#if defined(COMPILE_WITH_HIP) && defined(UseHIP)
		// HIP-specific graph and execution objects
		hipGraph_t hip_graph;                // Graph representation for HIP operations
		hipGraphExec_t hip_graphExec;        // Executable graph representation for HIP operations
		hipStream_t hip_graphStream;         // Stream for executing HIP graphs
		hipKernelNodeParams hip_nodeParams;   // Parameters for kernel nodes in HIP graphs
		std::vector<hipGraphNode_t> M_hipGraphNode_t; // Vector to store HIP graph nodes
#endif

#ifdef UseCUDA
		// CUDA-specific graph and execution objects
		cudaGraph_t cuda_graph;              // Graph representation for CUDA operations
		cudaGraphExec_t cuda_graphExec;      // Executable graph representation for CUDA operations
		cudaStream_t cuda_graphStream;       // Stream for executing CUDA graphs
		cudaKernelNodeParams cuda_nodeParams; // Parameters for kernel nodes in CUDA graphs
		std::vector<cudaGraphNode_t> M_cudaGraphNode_t; // Vector to store CUDA graph nodes
#endif

	public:
		// Constructor: Initializes a Task object with default values
		Task();

		// Destructor: Cleans up resources when a Task object is destroyed
		~Task();

		// Sets the flag to save task results 
		void setSave(bool b)
		{
			M_qSave = b;
		}

		// Sets the flag to view information about the task (for debugging)
		void setViewInfo(bool b)
		{
			M_qViewInfo = b;
		}

		// Sets the device reset flag (to control device state)
		void setDeviceReset(bool b)
		{
			M_qDeviceReset = b;
		}

		// Sets the type identifier for this task 
		void setNumType(int v)
		{
			M_numType = v;
		}

		// Sets the device identifier for HIP processing 
		void setDeviceHIP(int v);

		// Sets the device identifier for CUDA processing 
		void setDeviceCUDA(int v);

		// Sets the file name for saving task-related data 
		void setFileName(std::string s)
		{
			M_FileName = s;
		}

		// Returns the type identifier for this task 
		int getNumType() const
		{
			return (M_numType);
		}

		// Checks if the save flag is set 
		bool isSave() const
		{
			return (M_qSave);
		}

		// Checks if the device reset flag is set 
		bool isDeviceReset() const
		{
			return (M_qDeviceReset);
		}

		// Opens a task with specified number of blocks and threads 
		void open(int nbBlock, int NbTh);

		// Adds a HIP kernel function to the task with specified parameters 
		template <typename Kernel, typename Input, typename Output>
		void add_hip(const Kernel& kernel_function,
			int numElems,
			int iBegin, int iEnd,
			Input* buffer,
			Output* hostbuffer,
			std::vector<int> links);

		// Adds a CUDA kernel function to the task with specified parameters 
		template <typename Kernel, typename Input, typename Output>
		void add_cuda(const Kernel& kernel_function,
			int numElems,
			int iBegin, int iEnd,
			Input* buffer,
			Output* hostbuffer,
			std::vector<int> links);

		// Executes the added tasks or kernels 
		void run();

		// Closes the task and releases resources 
		void close();

		// Outputs debugging information about task execution 
		void debriefing();

};

Task::Task()
{
    M_numType = 1; // 1: HIP  2: CUDA
    M_nbTh = 1;
    M_numBlocksGPU = 1;
    M_q_graph = false;
    M_qViewInfo = true;
    M_time_laps = 0;
    M_FileName = "NoName";
    M_qSave = true;
    M_qDeviceReset = false;
    M_ListGraphDependencies.clear();
}

Task::~Task()
{
    M_ListGraphDependencies.clear();
}

void Task::debriefing()
{
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Debriefing"
                  << "\n";
        VLOG( 1 ) << "Elapsed microseconds = " << M_time_laps << "us\n";
        VLOG( 1 ) << "List Graph Dependencie  >>> [";
        for ( int i = 0; i < M_ListGraphDependencies.size(); i++ )
        {
            std::cout << M_ListGraphDependencies[i];
        }
        std::cout << "] <<<\n";
    }

    if ( M_qSave )
    {
        if ( M_qViewInfo )
        {
            VLOG( 1 ) << "Save Informations"
                      << "\n";
        }
        std::ofstream myfile;
        myfile.open( M_FileName + ".csv" );
        myfile << "Elapsed microseconds," << M_time_laps << "\n";
        myfile << "Nb Thread," << M_nbTh << "\n";
        myfile << "Nb Block used," << M_numBlocksGPU << "\n";
        myfile << "Nb Th/Block," << M_nThPerBckGPU << "\n";
        myfile << "List Graph Dependencie,";
        for ( int i = 0; i < M_ListGraphDependencies.size(); i++ )
        {
            myfile << M_ListGraphDependencies[i];
        }
        myfile << "\n";
        myfile.close();
    }
}

void Task::setDeviceHIP( int v )
{
#ifdef UseHIP
    int numDevices = 0;
    hipGetDeviceCount( &numDevices );
    if ( ( v >= 0 ) && ( v <= numDevices ) )
    {
        hipSetDevice( v );
    }
#endif
}

void Task::setDeviceCUDA( int v )
{
#ifdef UseCUDA
    int numDevices = 0;
    cudaGetDeviceCount( &numDevices );
    if ( ( v >= 0 ) && ( v <= numDevices ) )
    {
        cudaSetDevice( v );
    }
#endif
}

void Task::open( int nbBlock, int NbTh )
{
    M_q_graph = false;
    M_nbTh = NbTh;
    M_numBlocksGPU = nbBlock;
    M_nThPerBckGPU = M_nbTh / M_numBlocksGPU;
    if ( M_qViewInfo )
    {
        std::cout << "<=====================================================================>"
                  << "\n";
        VLOG( 1 ) << "Open Graph"
                  << "\n";
        VLOG( 1 ) << "nb Thread        = " << M_nbTh << "\n";
        VLOG( 1 ) << "nb Block         = " << M_numBlocksGPU << "\n";
        VLOG( 1 ) << "Thread Per Block = " << M_nThPerBckGPU << "\n";
        std::cout << "<=====================================================================>"
                  << "\n";
    }

    //#if defined(COMPILE_WITH_HIP) && defined(UseHIP)

#ifdef UseHIP
    hipGraphCreate( &hip_graph, 0 );
    hip_nodeParams = { 0 };
    memset( &hip_nodeParams, 0, sizeof( hip_nodeParams ) );
#endif

#ifdef UseCUDA
    cudaGraphCreate( &cuda_graph, 0 );
    cuda_nodeParams = { 0 };
    memset( &cuda_nodeParams, 0, sizeof( cuda_nodeParams ) );
#endif
}

template <typename Kernel, typename Input, typename Output>
void Task::add_hip( const Kernel& kernel_function,
                    int numElems,
                    int iBegin, int iEnd,
                    Input* buffer,
                    Output* hostbuffer,
                    std::vector<int> links )
{
#ifdef UseHIP
    bool qFlag = false;
    // BEGIN::Init new node
    hipGraphNode_t newKernelNode;
    M_hipGraphNode_t.push_back( newKernelNode );
    memset( &hip_nodeParams, 0, sizeof( hip_nodeParams ) );

    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Num Graph Node = " << M_hipGraphNode_t.size() << "\n";
    }

    // CRTL range
    if ( iEnd > numElems )
    {
        iEnd = numElems;
    }
    if ( iBegin < 0 )
    {
        iBegin = 0;
    }

    // printf("[%x]\n",M_hipGraphNode_t[M_hipGraphNode_t.size()-1]);

    hip_nodeParams.func = (void*)OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<Kernel, Input>;
    hip_nodeParams.gridDim = dim3( M_numBlocksGPU, 1, 1 );
    hip_nodeParams.blockDim = dim3( M_nThPerBckGPU, 1, 1 );
    hip_nodeParams.sharedMemBytes = 0;
    void* inputs[4];
    inputs[0] = (void*)&kernel_function;
    inputs[1] = (void*)&buffer;
    inputs[2] = (void*)&iBegin;
    inputs[3] = (void*)&iEnd;
    hip_nodeParams.kernelParams = inputs;

    if ( M_q_graph )
    {
        hip_nodeParams.extra = NULL;
    }
    else
    {
        hip_nodeParams.extra = nullptr;
    }
    // END::Init new node

    // BEGIN::Dependencies part
    unsigned int nbElemLinks = links.size();
    unsigned int nbElemKernelNode = M_hipGraphNode_t.size();
    std::vector<hipGraphNode_t> dependencies;
    M_ListGraphDependencies.push_back( M_hipGraphNode_t.size() - 1 );
    for ( int i = 0; i < nbElemLinks; i++ )
    {
        if ( links[i] == -1 )
        {
            qFlag = true;
        }
        if ( links[i] != -1 )
        {
            dependencies.push_back( M_hipGraphNode_t[links[i]] );
            M_ListGraphDependencies.push_back( links[i] );
        }
    }

    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Nb Elem Links  = " << nbElemLinks << "\n";
        VLOG( 1 ) << "Link dependencies with >>> [";
        for ( auto v : dependencies )
        {
            VLOG( 1 ) << v << "";
        }
        VLOG( 1 ) << "] <<<\n";
    }
    // END::Dependencies part

    // BEGIN::Add Node to kernel GPU
    if ( M_q_graph )
    {
        hipGraphAddKernelNode( &M_hipGraphNode_t[M_hipGraphNode_t.size() - 1], hip_graph, dependencies.data(), nbElemLinks, &hip_nodeParams );
    }
    else
    {
        hipGraphAddKernelNode( &M_hipGraphNode_t[M_hipGraphNode_t.size() - 1], hip_graph, nullptr, 0, &hip_nodeParams );
    }
    // END::Add Node to kernel GPU

    M_q_graph = true;

    // BEGIN::Final node kernel GPU
    if ( qFlag )
    {
        if ( M_qViewInfo )
        {
            VLOG( 1 ) << ""
                      << "\n";
            VLOG( 1 ) << "List >>> [";
            for ( auto v : M_hipGraphNode_t )
            {
                std::cout << v << "";
            }
            std::cout << "] <<<\n";
        }
        hipGraphNode_t copyBuffer;
        if ( M_qViewInfo )
        {
            VLOG( 1 ) << "Last M_hipGraphNode_t=" << M_hipGraphNode_t.size() << "= " << M_hipGraphNode_t[M_hipGraphNode_t.size() - 1] << "\n";
        }
        std::vector<hipGraphNode_t> finalDependencies = { M_hipGraphNode_t[M_hipGraphNode_t.size() - 1] };
        hipGraphAddMemcpyNode1D( &copyBuffer,
                                 hip_graph,
                                 dependencies.data(),
                                 dependencies.size(),
                                 hostbuffer,
                                 buffer,
                                 numElems * sizeof( typename std::decay<decltype( hostbuffer )>::type ),
                                 hipMemcpyDeviceToHost );
    }
    // END::Final node kernel GPU

    dependencies.clear();
#endif
}

template <typename Kernel, typename Input, typename Output>
void Task::add_cuda( const Kernel& kernel_function,
                     int numElems,
                     int iBegin, int iEnd,
                     Input* buffer,
                     Output* hostbuffer,
                     std::vector<int> links )
{
#ifdef UseCUDA
    bool qFlag = false;
    // BEGIN::Init new node
    cudaGraphNode_t newKernelNode;
    M_cudaGraphNode_t.push_back( newKernelNode );
    memset( &cuda_nodeParams, 0, sizeof( cuda_nodeParams ) );

    // if (M_qViewInfo) { VLOG(1) << "M_cudaGraphNode_t="<<M_cudaGraphNode_t.size()<<"= "<<M_cudaGraphNode_t[M_cudaGraphNode_t.size()-1]<<"\n"; }
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Num Graph Node = " << M_cudaGraphNode_t.size() << "\n";
    }

    // CRTL range
    if ( iEnd > numElems )
    {
        iEnd = numElems;
    }
    if ( iBegin < 0 )
    {
        iBegin = 0;
    }

    // printf("[%x]\n",M_cudaGraphNode_t[M_cudaGraphNode_t.size()-1]);

    cuda_nodeParams.func = (void*)OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<Kernel, Input>;
    cuda_nodeParams.gridDim = dim3( M_numBlocksGPU, 1, 1 );
    cuda_nodeParams.blockDim = dim3( M_nThPerBckGPU, 1, 1 );
    cuda_nodeParams.sharedMemBytes = 0;
    void* inputs[4];
    inputs[0] = (void*)&kernel_function;
    inputs[1] = (void*)&buffer;
    inputs[2] = (void*)&iBegin;
    inputs[3] = (void*)&iEnd;
    cuda_nodeParams.kernelParams = inputs;

    if ( M_q_graph )
    {
        cuda_nodeParams.extra = NULL;
    }
    else
    {
        cuda_nodeParams.extra = nullptr;
    }
    // END::Init new node

    // BEGIN::Dependencies part
    unsigned int nbElemLinks = links.size();
    unsigned int nbElemKernelNode = M_cudaGraphNode_t.size();
    std::vector<cudaGraphNode_t> dependencies;
    M_ListGraphDependencies.push_back( M_cudaGraphNode_t.size() - 1 );
    for ( int i = 0; i < nbElemLinks; i++ )
    {
        if ( links[i] == -1 )
        {
            qFlag = true;
        }
        if ( links[i] != -1 )
        {
            dependencies.push_back( M_cudaGraphNode_t[links[i]] );
            M_ListGraphDependencies.push_back( links[i] );
        }
    }

    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Nb Elem Links  = " << nbElemLinks << "\n";
        VLOG( 1 ) << "Link dependencies with >>> [";
        for ( auto v : dependencies )
        {
            VLOG( 1 ) << v << "";
        }
        VLOG( 1 ) << "] <<<\n";
    }
    // END::Dependencies part

    // BEGIN::Add Node to kernel GPU
    if ( M_q_graph )
    {
        cudaGraphAddKernelNode( &M_cudaGraphNode_t[M_cudaGraphNode_t.size() - 1], cuda_graph, dependencies.data(), nbElemLinks, &cuda_nodeParams );
    }
    else
    {
        cudaGraphAddKernelNode( &M_cudaGraphNode_t[M_cudaGraphNode_t.size() - 1], cuda_graph, nullptr, 0, &cuda_nodeParams );
    }
    // END::Add Node to kernel GPU

    M_q_graph = true;

    // BEGIN::Final node kernel GPU
    if ( qFlag )
    {
        if ( M_qViewInfo )
        {
            VLOG( 1 ) << ""
                      << "\n";
            VLOG( 1 ) << "List >>> [";
            for ( auto v : M_cudaGraphNode_t )
            {
                std::cout << v << "";
            }
            std::cout << "] <<<\n";
        }
        cudaGraphNode_t copyBuffer;
        if ( M_qViewInfo )
        {
            VLOG( 1 ) << "Last M_cudaGraphNode_t=" << M_cudaGraphNode_t.size() << "= " << M_cudaGraphNode_t[M_cudaGraphNode_t.size() - 1] << "\n";
        }
        std::vector<cudaGraphNode_t> finalDependencies = { M_cudaGraphNode_t[M_cudaGraphNode_t.size() - 1] };
        cudaGraphAddMemcpyNode1D( &copyBuffer,
                                  cuda_graph,
                                  dependencies.data(),
                                  dependencies.size(),
                                  hostbuffer,
                                  buffer,
                                  numElems * sizeof( typename std::decay<decltype( hostbuffer )>::type ),
                                  cudaMemcpyDeviceToHost );
    }
    // END::Final node kernel GPU

    dependencies.clear();
#endif
}

void Task::run()
{
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Run Graph"
                  << "\n";
    }
    M_t_begin = std::chrono::steady_clock::now();
    // Run HIP Graph
    if ( M_q_graph )
    {

//#if defined(COMPILE_WITH_HIP) && defined(UseHIP)
#ifdef UseHIP
        // hipEventRecord            (hip_start); <<<===not using causes memory crashes
        hipGraphInstantiate( &hip_graphExec, hip_graph, nullptr, nullptr, 0 );
        hipStreamCreateWithFlags( &hip_graphStream, hipStreamNonBlocking );
        hipGraphLaunch( hip_graphExec, hip_graphStream );
        hipStreamSynchronize( hip_graphStream );
        // hipEventRecord            (hip_stop);  <<<===not using causes memory crashes
        // hipEventElapsedTime       (&hip_milliseconds, hip_start, hip_stop);
#endif

#ifdef UseCUDA
        // cudaEventRecord           (cuda_start);  <<<===not using causes memory crashes
        cudaGraphInstantiate( &cuda_graphExec, cuda_graph, nullptr, nullptr, 0 );
        cudaStreamCreateWithFlags( &cuda_graphStream, cudaStreamNonBlocking );
        cudaGraphLaunch( cuda_graphExec, cuda_graphStream );
        cudaStreamSynchronize( cuda_graphStream );
        // cudaEventRecord           (cuda_stop);  <<<===not using causes memory crashes
        // cudaEventElapsedTime      (&cuda_milliseconds, cuda_start, cuda_stop);
#endif
    }

    M_t_end = std::chrono::steady_clock::now();
    M_time_laps = std::chrono::duration_cast<std::chrono::microseconds>( M_t_end - M_t_begin ).count();
}

void Task::close()
{
    if ( M_q_graph )
    {
        //#if defined(COMPILE_WITH_HIP) && defined(UseHIP)
#ifdef UseHIP
        hipGraphExecDestroy( hip_graphExec );
        hipGraphDestroy( hip_graph );
        hipStreamDestroy( hip_graphStream );
        if ( M_qDeviceReset )
        {
            hipDeviceReset();
        } // Not be used if Spex Hip AMD activated
          // Explicitly destroys and cleans up all resources associated with the current device in the current process. Any subsequent API call to this device will reinitialize the device.
#endif

#ifdef UseCUDA
        cudaGraphExecDestroy( cuda_graphExec );
        cudaGraphDestroy( cuda_graph );
        cudaStreamDestroy( cuda_graphStream );
        if ( M_qDeviceReset )
        {
            cudaDeviceReset();
        } // Not be used if Spex Cuda NVidia activated
          // Explicitly destroys and cleans up all resources associated with the current device in the current process. Any subsequent API call to this device will reinitialize the device.
#endif
        if ( M_qViewInfo )
        {
            VLOG( 1 ) << "Close Graph Hip"
                      << "\n";
        }
    }
}

template <typename T>
class PtrTask : public Task
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

    void set( int nb, T* v );
    void get( T* v );

    template <typename Kernel>
    void add( const Kernel& kernel_function, int iBegin, int iEnd, std::vector<int> links );
};

template <typename T>
PtrTask<T>::PtrTask()
{
    //...
}

template <typename T>
PtrTask<T>::~PtrTask()
{
    //...
}

template <typename T>
void PtrTask<T>::set( int nb, T* v )
{
#ifdef UseHIP
    if ( getNumType() == 1 )
    {
        BUFFER_HIP.memoryInit( nb );
        bufferSizeBytes = nb * sizeof( T );
        std::memcpy( &BUFFER_HIP.data, &v, bufferSizeBytes );
        BUFFER_HIP.memmovHostToDevice();
    }
#endif

#ifdef UseCUDA
    if ( getNumType() == 1 )
    {
        BUFFER_CUDA.memoryInit( nb );
        bufferSizeBytes = nb * sizeof( T );
        std::memcpy( &BUFFER_CUDA.data, &v, bufferSizeBytes );
        BUFFER_CUDA.memmovHostToDevice();
    }
#endif
}

template <typename T>
void PtrTask<T>::get( T* v )
{
#ifdef UseHIP
    if ( getNumType() == 1 )
    {
        BUFFER_HIP.memmovDeviceToHost();
        std::memcpy( &v, &BUFFER_HIP.data, bufferSizeBytes );
    }
#endif
#ifdef UseCUDA
    if ( getNumType() == 2 )
    {
        BUFFER_CUDA.memmovDeviceToHost();
        std::memcpy( &v, &BUFFER_CUDA.data, bufferSizeBytes );
    }
#endif
}

template <typename T>
template <typename Kernel>
void PtrTask<T>::add( const Kernel& kernel_function, int iBegin, int iEnd, std::vector<int> links )
{
#ifdef UseHIP
    if ( getNumType() == 1 )
    {
        Task::add_hip( kernel_function, BUFFER_HIP.size, iBegin, iEnd, BUFFER_HIP.deviceBuffer, BUFFER_HIP.data, links );
    }
#endif
#ifdef UseCUDA
    if ( getNumType() == 2 )
    {
        Task::add_cuda( kernel_function, BUFFER_CUDA.size, iBegin, iEnd, BUFFER_CUDA.deviceBuffer, BUFFER_CUDA.data, links );
    }
#endif
}

} // End namespace Taskgpu
//--------------------------------------------------------------------------------------------------------------------------------

//================================================================================================================================
// GRAPH IMPLICIT

namespace Taskgpui
{

#ifdef UseHIP

struct Task
{
		// Enumeration to represent the state of the task (capturing or updating)
		enum class state_t
		{
			capture, // State when capturing a new graph
			update   // State when updating an existing graph
		};

		// Adds a kernel node to the graph with specified parameters and stream
		void add_kernel_node(size_t key, hipKernelNodeParams params, hipStream_t s);

		// Updates an existing kernel node with new parameters
		void update_kernel_node(size_t key, hipKernelNodeParams params);

		// Returns the current state of the task
		state_t state() { return M_state; }

		// Destructor: Cleans up resources when a Task object is destroyed
		~Task();

	private:
		// Map to store kernel nodes with their associated keys
		std::unordered_map<size_t, hipGraphNode_t> _node_map;

		// Current state of the task (capture or update)
		state_t M_state;

		// Graph representation for HIP operations
		hipGraph_t M_graph;

		// Executable graph representation for HIP operations
		hipGraphExec_t M_graph_exec;

		// Flag indicating if the graph has been instantiated
		bool M_qInstantiated = false;

		// Begins capturing commands into a graph on the specified stream
		static void begin_capture(hipStream_t stream);

		// Ends capturing commands and finalizes the graph on the specified stream
		void end_capture(hipStream_t stream);

		// Launches the instantiated graph on the specified stream
		void launch_graph(hipStream_t stream);

	public:
		// Flag to control whether to always recapture the graph
		bool _always_recapture = false;

		// Wraps an object and captures or updates the task based on its state
		template <class Obj>
		void wrap(Obj& o, hipStream_t stream);
};

Task::~Task()
{
    if ( M_qInstantiated )
    {
        hipGraphDestroy( M_graph );
        hipGraphExecDestroy( M_graph_exec );
        M_qInstantiated = false;
    }
}

void Task::begin_capture( hipStream_t stream )
{
    hipStreamBeginCapture( stream, hipStreamCaptureModeGlobal );
}

void Task::end_capture( hipStream_t stream )
{
    if ( M_qInstantiated )
    {
        hipGraphDestroy( M_graph );
    }
    hipStreamEndCapture( stream, &M_graph );
    bool need_instantiation;

    if ( M_qInstantiated )
    {
        hipGraphExecUpdateResult updateResult;
        hipGraphNode_t errorNode;
        hipGraphExecUpdate( M_graph_exec, M_graph, &errorNode, &updateResult );
        if ( M_graph_exec == nullptr || updateResult != hipGraphExecUpdateSuccess )
        {
            hipGetLastError();
            if ( M_graph_exec != nullptr )
            {
                hipGraphExecDestroy( M_graph_exec );
            }
            need_instantiation = true;
        }
        else
        {
            need_instantiation = false;
        }
    }
    else
    {
        need_instantiation = true;
    }

    if ( need_instantiation )
    {
        hipGraphInstantiate( &M_graph_exec, M_graph, nullptr, nullptr, 0 );
    }
    M_qInstantiated = true;
}

template <class Obj>
void Task::wrap( Obj& o, hipStream_t stream )
{
    if ( !_always_recapture && M_qInstantiated )
    {
        M_state = state_t::update;
        o( *this, stream );
    }
    else
    {
        M_state = state_t::capture;
        begin_capture( stream );
        o( *this, stream );
        end_capture( stream );
    }
    launch_graph( stream );
}

void Task::launch_graph( hipStream_t stream )
{
    if ( M_qInstantiated )
    {
        hipGraphLaunch( M_graph_exec, stream );
    }
}

void Task::add_kernel_node( size_t key, hipKernelNodeParams params, hipStream_t stream )
{
    hipStreamCaptureStatus capture_status;
    hipGraph_t graph;
    const hipGraphNode_t* deps;
    size_t dep_count;
    hipStreamGetCaptureInfo_v2( stream, &capture_status, nullptr, &graph, &deps, &dep_count );
    hipGraphNode_t new_node;
    hipGraphAddKernelNode( &new_node, graph, deps, dep_count, &params );
    _node_map[key] = new_node;
    hipStreamUpdateCaptureDependencies( stream, &new_node, 1, 1 );
}

void Task::update_kernel_node( size_t key, hipKernelNodeParams params )
{
    hipGraphExecKernelNodeSetParams( M_graph_exec, _node_map[key], &params );
}

#endif

} // End namespace Taskgpui
//--------------------------------------------------------------------------------------------------------------------------------

// CLASS Task: Unified Memory
// NOTA: With Unified Memory: GPU accesses data directly from the "host"may be used without a separate "host"allocation and no copy routine is required,
// greatly simplifying and reducing the size of the program. With:
// - System Allocated: no other changes required.
// - Managed Memory: data allocation changed to use (cuda or hip)-MallocManaged(),which returns a pointer valid from both host and device code.

namespace Taskgpuu
{

//#if defined(useHIP) || defined(useCUDA)
#ifdef UseHIP

template <typename T>
class Task
	{
	private:
		// Number of blocks to be used on the GPU
		int M_numBlocksGPU;

		// Number of threads per block on the GPU
		int M_nbThGPU;

		// Flag indicating if resources are free for reuse
		bool M_isFree;

		// Flag for viewing information about the task
		bool M_qViewInfo;

		// Flag indicating if the device needs to be reset
		bool M_qDeviceReset;

		// Variable to track elapsed time in long format
		long int M_time_laps;

		// Variables to store elapsed time in milliseconds for HIP and CUDA
		float hip_milliseconds;
		float cuda_milliseconds;

		// Number of streams for concurrent execution
		int M_nbStreams;

		// Number of models being processed
		int M_numModel;

		// Type identifier for the task (e.g., different processing types)
		int M_numType;

		// Time points for measuring execution duration
		std::chrono::steady_clock::time_point M_t_begin, M_t_end;

		// Private method for executing kernels serially
		template <typename Kernel>
		void serial(const Kernel& kernel_function);

		// Private method for executing kernels using streams for concurrency
		template <typename Kernel>
		void stream(const Kernel& kernel_function);

	public:
#ifdef UseHIP
		// Buffer for Unified Memory management on HIP devices
		bufferGraphUnifiedHIP<T> BUFFER_HIP;
#endif

#ifdef UseCUDA
		// Buffer for Unified Memory management on CUDA devices
		bufferGraphUnifiedHIP<T> BUFFER_CUDA;
#endif

		// Constructor: Initializes a Task object
		Task();

		// Destructor: Cleans up resources when a Task object is destroyed
		~Task();

		// Opens a task with specified number of blocks and threads
		void open(int nbBlock, int NbTh);

		// Closes the task and releases resources
		void close();

		// Sets the view information flag (for debugging or logging)
		void setViewInfo(bool b)
		{
			M_qViewInfo = b;
		}

		// Sets the device reset flag (to control device state)
		void setDeviceReset(bool b)
		{
			M_qDeviceReset = b;
		}

		// Sets the number of streams for concurrent execution
		void setNbStreams(int v)
		{
			M_nbStreams = v;
		}

		// Sets the device identifier for HIP processing
		void setDeviceHIP(int v);

		// Sets the device identifier for CUDA processing
		void setDeviceCUDA(int v);

		// Sets the type identifier for this task
		void setNumType(int v)
		{
			M_numType = v;
		}

		// Returns the number of streams currently set for execution
		int getNbStreams() const
		{
			return (M_nbStreams);
		}

		// Checks if the device reset flag is set
		bool isDeviceReset() const
		{
			return (M_qDeviceReset);
		}

		// Returns the type identifier for this task
		int getNumType() const
		{
			return (M_numType);
		}

		// Returns elapsed time in long format since task started or last reset 
		long int getTimeLaps() const
		{
			return (M_time_laps);
		}

		// Runs a kernel function on the GPU using either serial or stream execution 
		template <typename Kernel>
		void run(const Kernel& kernel_function);

		// Outputs debugging information about task execution 
		void debriefing();
};



template <typename T>
Task<T>::~Task()
{
    if ( !M_isFree )
    {
        close();
    }
}

template <typename T>
Task<T>::Task()
{
    M_numBlocksGPU = 512;
    M_isFree = false;
    M_qViewInfo = true;
    M_qDeviceReset = false;
    M_nbStreams = 3;
    M_numType = 1;  // 1:HIP 2:CUDA
    M_numModel = 1; // 1:Serial 2:Stream
    M_time_laps = 0;
}

template <typename T>
void Task<T>::open( int nbBlock, int NbTh )
{
    M_nbThGPU = NbTh;
    M_numBlocksGPU = nbBlock;
#ifdef UseHIP
    BUFFER_HIP.memoryInit( NbTh );
#endif
#ifdef UseCUDA
    BUFFER_CUDA.memoryInit( NbTh );
#endif
}

template <typename T>
void Task<T>::close()
{
#ifdef UseHIP
    hipFree( BUFFER_HIP.data );
    if ( M_qDeviceReset )
    {
        hipDeviceReset();
    }
#endif

#ifdef UseCUDA
    cudaFree( BUFFER_CUDA.data );
    if ( M_qDeviceReset )
    {
        cudaDeviceReset();
    }
#endif
    M_isFree = true;
}

template <typename T>
void Task<T>::setDeviceHIP( int v )
{
#ifdef UseHIP
    int numDevices = 0;
    hipGetDeviceCount( &numDevices );
    if ( ( v >= 0 ) && ( v <= numDevices ) )
    {
        hipSetDevice( v );
    }
#endif
}

template <typename T>
void Task<T>::setDeviceCUDA( int v )
{
#ifdef UseCUDA
    int numDevices = 0;
    cudaGetDeviceCount( &numDevices );
    if ( ( v >= 0 ) && ( v <= numDevices ) )
    {
        cudaSetDevice( v );
    }
#endif
}

template <typename T>
void Task<T>::debriefing()
{
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Debriefing"
                  << "\n";
        VLOG( 1 ) << "Elapsed microseconds = " << M_time_laps << "us\n";
    }
}

template <typename T>
template <typename Kernel>
void Task<T>::serial( const Kernel& kernel_function )
{
    M_t_begin = std::chrono::steady_clock::now();
    int num_blocks = ( BUFFER_HIP.size + M_numBlocksGPU - 1 ) / M_numBlocksGPU;
    dim3 thread_block( M_numBlocksGPU, 1, 1 );
    dim3 grid( num_blocks, 1 );
#ifdef UseHIP
    hipLaunchKernelGGL( OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D, grid, thread_block, 0, 0, kernel_function, BUFFER_HIP.data, 0, BUFFER_HIP.size );
    hipDeviceSynchronize();
#endif
#ifdef UseCUDA
    OP_IN_KERNEL_GRAPH_LAMBDA_GPU_1D<<<grid, thread_block, 0, 0>>>( kernel_function, BUFFER_CUDA.data, 0, BUFFER_CUDA.size );
    cudaDeviceSynchronize();
#endif
    M_t_end = std::chrono::steady_clock::now();
    M_time_laps = std::chrono::duration_cast<std::chrono::microseconds>( M_t_end - M_t_begin ).count();
}

template <typename T>
template <typename Kernel>
void Task<T>::stream( const Kernel& kernel_function )
{
#ifdef UseHIP
    int numVersion = 1;
    float ms;
    int i;

    // const int blockSize = 20;  //const int blockSize M_block_size;
    const int blockSize = M_numBlocksGPU;
    const int numElems = BUFFER_HIP.size;
    const int streamSize = numElems / M_nbStreams;
    const int streamBytes = streamSize * sizeof( decltype( BUFFER_HIP.data ) );
    const int sz = numElems * sizeof( decltype( BUFFER_HIP.data ) );
    int grid = streamSize / blockSize;
    if ( M_qViewInfo )
    {
        VLOG( 1 ) << "Block size :" << M_numBlocksGPU << "\n";
        VLOG( 1 ) << "nb Streams :" << M_nbStreams << "\n";
    }
    if ( ( grid == 0 ) && ( M_qViewInfo ) )
    {
        VLOG( 1 ) << "Grid size error : ==> auto correction"
                  << "\n";
    }
    grid = max( streamSize / blockSize, 1 );
    hipEvent_t startEvent, stopEvent, dummyEvent;
    hipStream_t stream[M_nbStreams];
    hipEventCreate( &startEvent );
    hipEventCreate( &stopEvent );
    hipEventCreate( &dummyEvent );
    for ( i = 0; i < M_nbStreams; ++i )
    {
        hipStreamCreate( &stream[i] );
    }

    // hipMallocManaged((void **) &data,sz,hipMemAttachHost);
    hipEventRecord( startEvent, 0 );

    if ( numVersion == 1 )
    {
        for ( i = 0; i < M_nbStreams; ++i )
        {
            int offset = i * streamSize;
            hipLaunchKernelGGL( OP_IN_KERNEL_GRAPH_LAMBDA_STREAM_GPU_1D, grid, blockSize, 0, stream[i], kernel_function, BUFFER_HIP.data, offset );
        }
        hipEventRecord( stopEvent, 0 );
        hipEventSynchronize( stopEvent );
        // hipEventElapsedTime(&ms, startEvent, stopEvent);
        // printf("Time (ms): %f\n", ms);
    }

    // cleanup
    hipEventDestroy( startEvent );
    hipEventDestroy( stopEvent );
    hipEventDestroy( dummyEvent );
    for ( int i = 0; i < M_nbStreams; ++i )
    {
        hipStreamDestroy( stream[i] );
    }
#endif
}

template <typename T>
template <typename Kernel>
void Task<T>::run( const Kernel& kernel_function )
{
    M_t_begin = std::chrono::steady_clock::now();
    if ( M_numModel == 1 )
    {
        // SERIAL MODEL
        if ( M_numType == 1 )
        {
            serial( kernel_function );
        }
    }

    if ( M_numModel == 2 )
    { // STREAM MODEL
        if ( M_numType == 1 )
        {
            stream( kernel_function );
        }
    }

    M_t_end = std::chrono::steady_clock::now();
    M_time_laps = std::chrono::duration_cast<std::chrono::microseconds>( M_t_end - M_t_begin ).count();
}

#endif

} // End namespace Taskgpuu
//--------------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------------

void resetAndDestroyAllGPU()
{
    // Reset and explicitly destroy all resources associated with the current device
    int numDevices = 0;
#ifdef UseHIP
    HIP_CHECK( hipGetDeviceCount( &numDevices ) );
    for ( int i = 0; i < numDevices; i++ )
    {
        HIP_CHECK( hipSetDevice( i ) );
        HIP_CHECK( hipDeviceReset() );
    }
#endif
#ifdef UseCUDA
    CUDA_CHECK( cudaGetDeviceCount( &numDevices ) );
    for ( int i = 0; i < numDevices; i++ )
    {
        CUDA_CHECK( cudaSetDevice( i ) );
        CUDA_CHECK( cudaDeviceReset() );
    }
#endif
}

//--------------------------------------------------------------------------------------------------------------------------------

// Temporary functions intended to control the types of variables. Must be removed in the future

template <typename T>
bool isSpecxGPUFunctionBeta( T& fcv )
{
    std::string s1 = typeid( fcv ).name();
    int l = s1.length();
    if ( l > 1 )
    {
        if ( s1.find( "SpCallableType0EE" ) != std::string::npos )
        {
            return ( true );
        } //"SpCpu"
        else if ( s1.find( "SpCallableType2EE" ) != std::string::npos )
        {
            return ( true );
        } //"SpHip"
        else if ( s1.find( "SpCallableType1EE" ) != std::string::npos )
        {
            return ( true );
        } //"SpCuda
    }
    return ( false );
}

template <typename T>
int numSpecxFunctionBeta( T& fcv )
{
    std::string s1 = typeid( fcv ).name();
    int l = s1.length();
    // std::cout<<"s1="<<s1<<"\n";
    if ( l > 1 )
    {
        if ( s1.find( "SpCallableType0EE" ) != std::string::npos )
        {
            return ( 1 );
        } //"SpCpu"
        else if ( s1.find( "SpCallableType2EE" ) != std::string::npos )
        {
            return ( 2 );
        } //"SpHip"
        else if ( s1.find( "SpCallableType1EE" ) != std::string::npos )
        {
            return ( 3 );
        } //"SpCuda
        else if ( s1.find( "SpArrayView" ) != std::string::npos )
        {
            return ( 20 );
        } // SpArrayView
        else if ( s1.find( "SpArrayAccessorIS1_EE" ) != std::string::npos )
        {
            return ( 21 );
        } // SpReadArray
        else if ( s1.find( "SpDataAccessMode0" ) != std::string::npos )
        {
            return ( 10 );
        } //"SpRead
        else if ( s1.find( "SpDataAccessMode1" ) != std::string::npos )
        {
            return ( 11 );
        } //"SpWrite
        else if ( s1.find( "SpDataAccessMode3" ) != std::string::npos )
        {
            return ( 12 );
        } // SpCommutativeWrite
    }
    return ( 0 );
}

} // namespace Feel
