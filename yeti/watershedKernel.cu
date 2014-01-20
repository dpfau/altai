/*
 * Kernel to compute watershed seeded from one pixel
 *
 * David Pfau, 2014
 */

#define MEM_BLOCK 512

__global__ void watershedKernel(const float * A, int * B, const int * seedIdx, const int N, const int ndims, const int * dims)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N) {
        int * idx = (int *)malloc(sizeof(int) * MEM_BLOCK); // This array will grow dynamically as we add new indices to the list
        int npix = 1; // number of pixels in the watershed
        int dpix = 1; // change in number of pixels in this loop
        idx[0] = seedIdx[i];
        while (dpix > 0) {
            dpix = 0;
            for (int k=0; k<npix; k++) {
               for (int ii=-1; ii<=1; ii++) {
                  for (int jj=-1; jj<=1; jj++) {
                      switch(ndims) {
                          case 2:
                              if ( (idx[k] % dims[0]) + ii >= 0 && (idx[k] % dims[0]) + ii < dims[0]
                                && (idx[k] / dims[0]) + jj >= 0 && (idx[k] / dims[0]) + jj < dims[1] ) { // check that index is within bounds
                                  int new_idx = (idx[k] % dims[0]) + ii + dims[0] * ((idx[k] / dims[0]) + jj);
                                  for (int kk=0; kk<npix; kk++) { 
                                  // Check that index is not in the list. 
                                  // Linear search is stupid and slow, but I'm not sure how to put a hash table inside CUDA device code
                                      if (idx[kk]==new_idx) {
                                          new_idx = -1;
                                          break;
                                      }
                                  }
                                  if (new_idx != -1) {
                                      if (A[new_idx] > 0.0f) {
                                          int neighbor_idx = new_idx;
                                          for (int iii=-1; iii<=1; iii++) {
                                              for (int jjj=-1; jjj<=1; jjj++) {
                                                 if ( (new_idx % dims[0]) + iii >= 0 && (new_idx % dims[0]) + iii < dims[0]
                                                   && (new_idx / dims[0]) + jjj >= 0 && (new_idx / dims[0]) + jjj < dims[1] ) {
                                                    int this_idx = (new_idx % dims[0]) + iii + dims[0] * ((new_idx / dims[0]) + jjj);
                                                    if (A[this_idx] > A[neighbor_idx]) {
                                                        neighbor_idx = this_idx;
                                                    }
                                                 } 
                                              }
                                          }
                                          for (int kk=0; kk < npix; kk++) {
                                              if (idx[kk]==neighbor_idx) {
                                                  if ( (npix + dpix) % MEM_BLOCK == 0 ) {
                                                      int * bigger_idx = (int *)realloc(idx, npix + dpix + MEM_BLOCK);
                                                      if (bigger_idx != NULL) {
                                                          idx = bigger_idx;
                                                      } else {
                                                          free(idx);
                                                          puts("Error reallocating memory");
                                                          exit(1);
                                                      }
                                                  }
                                                  idx[npix + dpix++] = new_idx;
                                                  break;
                                              }
                                          }
                                      }
                                  }
                              }
                          case 3:
                      }
                  }
               }
            }
            npix += dpix;
        }

        for (int k=0; k < npix; k++) {
            B[idx[k]] = i;
        }

        free(idx);
    }
}
