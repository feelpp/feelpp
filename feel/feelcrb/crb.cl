
__kernel void VMProd(__global double * V, __global double * M, int nComp, int nProd)
{
    int i = get_global_id(0);
}

__kernel void MVProd(__global double * M, __global double * V, int nComp, int numElements)
{
    int i = get_global_id(0);
}
