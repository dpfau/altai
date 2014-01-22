#include "mex.h"
#include "lsqr.h"



void aprod(int mode, int m, int n, double x[], double y[], void *UsrWrk) {
    // Wrapper around matrix multiplication in the format required by lsqr
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    int numPixels = mxGetNumberOfElements( prhs[0] );
    double * frame = mxGetPr( prhs[0] );
    int numROI = *(int *)( mxGetData( prhs[2] ) );

    plhs[1] = mxCreateDoubleArray( (mwSize)numROI, 1, mxREAL);
    double * rates = mxGetPr( plhs[1] );

    // Output-only variables requires by lsqr
    int istop_out;
    int itn_out;
    double anorm_out;
    double acond_out;
    double rnorm_out;
    double arnorm_out;
    double xnorm_out;

    lsqr( numPixels, numROI, &aprod, 0.0, (void *) p, y, v, w, x, NULL, 1e-9, 1e-9, 1e8, 10, NULL,
        &istop_out, &itn_out, &anorm_out, &acond_out, &rnorm_out, &arnorm_out, &xnorm_out );
}
