/*
 * Kernel to compute watershed seeded from one pixel
 *
 * David Pfau, 2014
 */

void __global__ watershedKernel(const float * A, const float * B, const int * seedIdx, const int N, const int ndims, const int * dims)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        B[i] = 2.0 * A[i];
    }
}