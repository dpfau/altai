#include "mex.h"
#include "lsqr.h"

typedef struct {
    int ndims;
    const int * localDims;
    const int * globalDims;
    int * loc;
    int nROI;
    double * ROIs;
} param;

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

void aprod(int mode, int m, int n, double x[], double y[], void *UsrWrk) {
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

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    param p;
    p.ndims = mxGetNumberOfDimensions( prhs[0] );
    p.localDims  = mxGetDimensions( prhs[1] );
    p.globalDims = mxGetDimensions( prhs[0] );
    p.loc = (int *)mxGetData( prhs[2] );
    p.nROI = *(int *)( mxGetData( prhs[3] ) );
    p.ROIs = mxGetPr( prhs[1] );

    double * frame = mxGetPr( prhs[0] );
    double * ROIs = mxGetPr( prhs[1] );

    plhs[0] = mxCreateDoubleMatrix( (mwSize)p.nROI, 1, mxREAL );
    double * rates = mxGetPr( plhs[0] );

    // create workspace variables
    double * v = (double *)malloc(sizeof(double) * p.nROI);
    double * w = (double *)malloc(sizeof(double) * p.nROI);

    // Output-only variables requires by lsqr
    int istop_out;
    int itn_out;
    double anorm_out;
    double acond_out;
    double rnorm_out;
    double arnorm_out;
    double xnorm_out;

    lsqr( mxGetNumberOfElements( prhs[0] ), p.nROI, &aprod, 0.0, (void *) &p, frame, v, w, rates, NULL, 1e-9, 1e-9, 1e8, 10, NULL,
        &istop_out, &itn_out, &anorm_out, &acond_out, &rnorm_out, &arnorm_out, &xnorm_out );

    free(v);
    free(w);
}
