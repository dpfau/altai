#ifndef _LOCALMUL_H_
#define _LOCALMUL_H_

typedef struct {
    int ndims;
    const int * localDims;
    const int * globalDims;
    int * loc;
    int nROI;
    double * ROIs;
} param;

void indToSub(int idx, int * sub, const int * dims, const int ndims);
int localToGlobal(int * sub, const int * dims, int * offset, const int ndims);
void localMultiply(int mode, int m, int n, double x[], double y[], void *UsrWrk);

#endif