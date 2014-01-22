/*
* Compute the watershed regions greater than zero. Useful when the watershed regions are
* quite a small part of a large image, as it avoids doing any computation outside a
* relevant region. Could be accelerated by using Boost's unordered_set insted of C++ STL
* sets (constant vs. logarithmic insertion and checking), but works well for small watersheds.
*
*
* David Pfau, 2014
*/

#include <set>
#include <utility>
#include "mex.h"

using namespace std;

int checkIndex(const int idx, const int * offset, const int * dims, const int ndims) {
    switch (ndims) {
        case 2:
            return ((idx % dims[0]) + offset[0] >= 0) && ((idx % dims[0]) + offset[0] < dims[0])
                && ((idx / dims[0]) + offset[1] >= 0) && ((idx / dims[0]) + offset[1] < dims[1]);
        case 3:
            return (idx % dims[0])                 + offset[0] >= 0 && (idx % dims[0])             + offset[0] < dims[0]
                    && ((idx / dims[0]) % dims[1]) + offset[1] >= 0 && ((idx / dims[0]) % dims[1]) + offset[1] < dims[1]
                    && (idx / (dims[0] * dims[1])) + offset[2] >= 0 && (idx / (dims[0] * dims[1])) + offset[2] < dims[2];
    }
    return 0;
}

int offsetIndex(const int idx, const int * offset, const int * dims, const int ndims) {
    switch (ndims) {
        case 2:
            return (idx % dims[0]) + offset[0] + dims[0] * ((idx / dims[0]) + offset[1]);
        case 3:
            return                      (idx % dims[0])              + offset[0]  +
                    dims[0] *           (((idx / dims[0]) % dims[1]) + offset[1]) +
                    dims[0] * dims[1] * ((idx / (dims[0] * dims[1])) + offset[2]);
    }
    return 0;
}

int checkOneNeighbor(const float * data, int idx, int old_idx, const int * offset, const int * dims, const int ndims) {
    if ( checkIndex(idx, offset, dims, ndims) ) {
        int new_idx = offsetIndex(idx, offset, dims, ndims);
        if (data[new_idx] > data[old_idx]) {
            return new_idx;
        } else {
            return old_idx;
        }
    }
    return old_idx;
}

int checkOneIndex(const float * data, set<int> * watershed, int idx, int * offset, const int * dims, const int ndims) {
    if ( checkIndex(idx, offset, dims, ndims) ) { // check that index is within bounds
        int new_idx = offsetIndex(idx, offset, dims, ndims);
        if (watershed->find(new_idx) == watershed->end()) { // logarithmic rather than linear. Could do this in constant time if I had Boost or C++11
            if (data[new_idx] > 0.0f) {
                int neighbor_idx = new_idx;
                for (int iii=-1; iii<=1; iii++) {
                    for (int jjj=-1; jjj<=1; jjj++) {
                        switch (ndims) {
                            case 2:
                            {
                                if (iii != 0 || jjj != 0) {
                                    int offset2[] = {iii,jjj};
                                    neighbor_idx = checkOneNeighbor(data, new_idx, neighbor_idx, offset2, dims, ndims);
                                }
                                break;
                            }
                            case 3:
                            {
                                for (int kkk=-1; kkk<=1; kkk++) {
                                    if (iii != 0 || jjj != 0 || kkk != 0) {
                                        int offset2[] = {iii,jjj,kkk};
                                        neighbor_idx = checkOneNeighbor(data, new_idx, neighbor_idx, offset2, dims, ndims);
                                    }
                                }
                            }
                        }
                    }
                }
                if (watershed->find(neighbor_idx) != watershed->end())
                {
                    watershed->insert(new_idx);
                    return 1;
                } else {
                    return 0;
                }
            }
        }
    }
    return 0;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    float *data = (float*) mxGetData( prhs[0] );
    int *output = (int*) mxGetData( prhs[1] );
    int *seeds = (int*) mxGetData( prhs[2] );
    
    int N = mxGetNumberOfElements( prhs[2] );
    int ndims = mxGetNumberOfDimensions( prhs[0] );
    const int * dims = mxGetDimensions( prhs[0] );
    
    set<int> * watersheds = new set<int>[N];
    set<int>::iterator iter;
    
    for (int i=0; i < N; i++) {
        int dpix = 1; // change in number of pixels in this loop
        watersheds[i].insert(seeds[i]-1); // because we pass Matlab indices in here, subtract one
        while (dpix > 0) {
            dpix = 0;
            for (iter = watersheds[i].begin(); iter != watersheds[i].end(); ++iter) {
                for (int ii=-1; ii<=1; ii++) {
                    for (int jj=-1; jj<=1; jj++) {
                        switch(ndims) {
                            case 2:
                            {
                                if (ii != 0 || jj != 0) {
                                    int offset[] = {ii, jj};
                                    dpix += checkOneIndex(data, watersheds + i, *iter, offset, dims, ndims);
                                }
                                break;
                            }
                            case 3:
                            {
                                for (int kk=-1; kk<=1; kk++) {
                                    if (ii != 0 || jj != 0 || kk != 0) {
                                        int offset[] = {ii, jj, kk};
                                        dpix += checkOneIndex(data, watersheds + i, *iter, offset, dims, ndims);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    for (int i=0; i < N; i++) {
        for (iter = watersheds[i].begin(); iter != watersheds[i].end(); ++iter) {
            output[*iter] = i+1;
        }
    }
}
