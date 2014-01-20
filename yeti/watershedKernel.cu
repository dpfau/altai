/*
 * Kernel to compute watershed seeded from one pixel
 *
 * David Pfau, 2014
 */

__global__ void watershedKernel(float * A, float * B, const int * seedIdx, const int N, const int ndims, const int * dims)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        B[i] = 2.0f * A[i];
    }
}