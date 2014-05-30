// Enable double support: not needed for 1.2
#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#else
    #error "Double precision floating point not supported by OpenCL implementation."
#endif

/* In this first version, we make as many blocks as there are value in BetaAq (Attention to block size limit) */
/* we apply this kernel to matrix and vector (same way stored) */
/* size is the total number of elements of the vector or the matrix (NxN) */
/* Size limit if we want to copy NxN matrix into local memory (32Ko on Tahiti => 64x64 matrix) */
/* Constraints: Block count limit (Tahiti: 256) */
__kernel void SVProd(__global double * S, __global double * V, const int Vsz)
{
    /* get beta value */
    __private double beta = S[get_group_id(0)];

    /* compute size and indexes */
    __private int i = 0;
    __private int idx = get_group_id(0) * Vsz;
    __private int nIt = Vsz / get_local_size(0);
    __private int rIt = Vsz % get_local_size(0);
    __private int lidx = 0;

    /* Check if the matrix size is a multiple of the number of work items */
    if(rIt > 0)
    {
        nIt++;
    }

    /* iterate over the matrix */
    for(i = 0; i < nIt; i++)
    {
        /* get current element index */
        lidx = i * get_local_size(0) + get_local_id(0);
        
        /* check if we are still inside valid memory */
        if(lidx < Vsz)
        {
            V[idx + lidx] = beta * V[idx + lidx];
        }
    }
}

/* Computes the sum of vectors or matrix (not optimized) */
/* Dedicates one block to it (one vector must at least fit in shared memory: 32Kb for Tahiti => 4096 elements max */
/* V: storage for vector/matrix */
/* Vsz: Nb elements in vector/matrix */
/* nV: Number of vectors/matrices */
__kernel void VSum(__global double * out, __global double * V, const int Vsz, const int nV) 
{
    // forced to alloc statically
    __local double lV[2048];

    /* compute size and indexes */
    __private int i = 0;
    __private int j = 0;
    __private int nIt = Vsz / get_local_size(0);
    __private int rIt = Vsz % get_local_size(0);
    __private int lidx = 0;

    /* Check if the matrix size is a multiple of the number of work items */
    if(rIt > 0)
    {
        nIt++;
    }

    /* iterate over all vectors */
    for(j = 0; j < nV; j++)
    {
        /* iterate over the components of the vector */
        for(i = 0; i < nIt; i++)
        {
            /* get current element index */
            lidx = i * get_local_size(0) + get_local_id(0);

            //printf("%d/%d %d/%d %d\n", get_group_id(0), get_num_groups(0), get_local_id(0), get_local_size(0), lidx);

            /* check if we are still inside valid memory */
            if(lidx < Vsz)
            {
                /* first time, we only fetch the data */
                if(j == 0)
                {
                    lV[lidx] = V[lidx];
                }
                /* otherwise we sum with previous results */
                else
                {
                    lV[lidx] = lV[lidx] + V[Vsz * j + lidx];
                }
            }
        }
    }

    /* send back data to global memory */
    for(i = 0; i < nIt; i++)
    {
        /* get current element index */
        lidx = i * get_local_size(0) + get_local_id(0);

        /* check if we are still inside valid memory */
        if(lidx < Vsz)
        {
            //printf("%f\n", lV[lidx]);
            out[lidx] = lV[lidx];
        }
    }
}

