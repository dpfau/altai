/*
 * Kernel to compute watershed seeded from one pixel
 *
 * David Pfau, 2014
 */

#define MEM_BLOCK 512

__device__ int * myrealloc(int * old, int oldsize, int newsize)
{
    int * newT = (int *) malloc (newsize*sizeof(int));

    for(int i=0; i<oldsize; i++)
    {
        newT[i] = old[i];
    }

    free(old);
    return newT;
}

__device__ int checkIndex(const int idx, const int * offset, const int * dims, const int ndims) {
    switch (ndims) {
        case 2:
            return (idx % dims[0]) + offset[0] >= 0 && (idx % dims[0]) + offset[0] < dims[0]
                && (idx / dims[0]) + offset[1] >= 0 && (idx / dims[0]) + offset[1] < dims[1];
        case 3:
            return (idx % dims[0])             + offset[0] >= 0 && (idx % dims[0])             + offset[0] < dims[0]
                && ((idx / dims[0]) % dims[1]) + offset[1] >= 0 && ((idx / dims[0]) % dims[1]) + offset[1] < dims[1] 
                && (idx / (dims[0] * dims[1])) + offset[2] >= 0 && (idx / (dims[0] * dims[1])) + offset[2] < dims[2];
    }
    return 0;
}

__device__ int offsetIndex(const int idx, const int * offset, const int * dims, const int ndims) {
    switch (ndims) {
        case 2:
            return (idx % dims[0]) + offset[0] + dims[0] * ((idx / dims[0]) + offset[1]);
        case 3:
            return                  (idx % dims[0])              + offset[0]  + 
                dims[0] *           (((idx / dims[0]) % dims[1]) + offset[1]) + 
                dims[0] * dims[1] * ((idx / (dims[0] * dims[1])) + offset[2]);
    }
    return 0;
}

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
            for (int q=0; q<npix; q++) {
               for (int ii=-1; ii<=1; ii++) {
                  for (int jj=-1; jj<=1; jj++) {
                      switch(ndims) {
                          case 2:
                              int offset[] = {ii, jj};
                              if ( checkIndex(idx[q], offset, dims, ndims) ) { // check that index is within bounds
                                  int new_idx = offsetIndex(idx[q], offset, dims, ndims);
                                  for (int qq=0; qq<npix; qq++) { 
                                  // Check that index is not in the list. 
                                  // Linear search is stupid and slow, but I'm not sure how to put a hash table inside CUDA device code
                                      if (idx[qq]==new_idx) {
                                          new_idx = -1;
                                          break;
                                      }
                                  }
                                  if (new_idx != -1) {
                                      if (A[new_idx] > 0.0f) {
                                          int neighbor_idx = new_idx;
                                          for (int iii=-1; iii<=1; iii++) {
                                              for (int jjj=-1; jjj<=1; jjj++) {
                                                 if ( checkIndex(new_idx, offset, dims, ndims) ) {
                                                    offset[0] = iii; offset[1] = jjj;
                                                    int this_idx = offsetIndex(new_idx, offset, dims, ndims);
                                                    if (A[this_idx] > A[neighbor_idx]) {
                                                        neighbor_idx = this_idx;
                                                    }
                                                 } 
                                              }
                                          }
                                          for (int qq=0; qq < npix; qq++) {
                                              if (idx[qq]==neighbor_idx) {
                                                  if ( (npix + dpix) % MEM_BLOCK == 0 ) {
                                                      idx = myrealloc(idx, npix + dpix, npix + dpix + MEM_BLOCK);
                                                  }
                                                  idx[npix + dpix++] = new_idx;
                                                  break;
                                              }
                                          }
                                      }
                                  }
                              }
                          case 3:
                              for (int kk=-1; kk<=1; kk++) {
                                  int offset[] = {ii, jj, kk};
                                  if ( checkIndex(idx[q], offset, dims, ndims) ) { // check that index is within bounds
                                      int new_idx = offsetIndex(idx[q], offset, dims, ndims);
                                      for (int qq=0; qq<npix; qq++) { 
                                      // Check that index is not in the list. 
                                      // Linear search is stupid and slow, but I'm not sure how to put a hash table inside CUDA device code
                                          if (idx[qq]==new_idx) {
                                              new_idx = -1;
                                              break;
                                          }
                                      }
                                      if (new_idx != -1) {
                                          if (A[new_idx] > 0.0f) {
                                              int neighbor_idx = new_idx;
                                              for (int iii=-1; iii<=1; iii++) {
                                                  for (int jjj=-1; jjj<=1; jjj++) {
                                                      for (int kkk=-1; kkk<=1; kkk++) {
                                                          if ( checkIndex(new_idx, offset, dims, ndims)) {
                                                              offset[0] = iii; offset[1] = jjj; offset[2] = kkk;
                                                              int this_idx = offsetIndex(new_idx, offset, dims, ndims);
                                                              if (A[this_idx] > A[neighbor_idx]) {
                                                                  neighbor_idx = this_idx;
                                                              }
                                                          } 
                                                      }
                                                  }
                                              }
                                              for (int qq=0; qq < npix; qq++) {
                                                  if (idx[qq]==neighbor_idx) {
                                                      if ( (npix + dpix) % MEM_BLOCK == 0 ) {
                                                          idx = myrealloc(idx, npix + dpix, npix + dpix + MEM_BLOCK);
                                                      }
                                                      idx[npix + dpix++] = new_idx;
                                                      break;
                                                  }
                                              }
                                          }
                                      }
                                  }
                              }
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
