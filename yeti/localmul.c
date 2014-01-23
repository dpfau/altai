#include "localmul.h"
#include <stdlib.h>

void indToSub(int idx, int * sub, const int * dims, const int ndims) {
    sub[0] = idx % dims[0];
    int i;
    for (i = 1; i < ndims; i++) {
        idx = idx / dims[i-1];
        sub[i] = idx % dims[i];
    }
}

int localToGlobal(int * sub, const int * dims, int * offset, const int ndims) {
    int new_idx = 0;
    int i;
    for (i = ndims-1; i > 0; i--) {
        new_idx += sub[i] + offset[i];
        new_idx *= dims[i-1];
    }
    new_idx += sub[0] + offset[0];
    return new_idx;
}

void localMultiply(int mode, int m, int n, double x[], double y[], void *UsrWrk) {
    // Wrapper around matrix multiplication in the format required by lsqr
    param * p = (param *) UsrWrk;
    int q = 1;
    int i, j;
    for (i = 0; i < p->ndims; i++) {
        q *= p->localDims[i];
    }

    int * foo = (int *) calloc(p->ndims,sizeof(int)); // the zero offset
    int * sub = (int *) calloc(p->ndims,sizeof(int)); // subscripts at each step of the loop. better to allocate it once.
    double A;
    int yIdx;
    for (i = 0; i < q; i++) { // If we wanted to really knock things out of the park, maybe could parallelize this loop on GPU? But only in R2013b...
        for (j = 0; j < p->nROI; j++) {
            A = p->ROIs[i + j * q];
            if (A != 0) {
                indToSub(i, sub, p->localDims, p->ndims);
                yIdx = localToGlobal(sub, p->globalDims, p->loc + (p->ndims * j), p->ndims);
                if (mode == 1) {
                    y[yIdx] += A * x[j];
                } else if (mode == 2) {
                    x[j] += A * y[yIdx];
                }
            }
        }
    }
    free(foo);
    free(sub);
}