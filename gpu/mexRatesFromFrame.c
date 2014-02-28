#include "mex.h"
#include "lsqr.h"
#include "localmul.h"

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

    double * v = (double *)malloc(sizeof(double) * p.nROI);
    double * w = (double *)malloc(sizeof(double) * p.nROI);

    int istop_out;
    int itn_out;
    double anorm_out;
    double acond_out;
    double rnorm_out;
    double arnorm_out;
    double xnorm_out;

    lsqr( mxGetNumberOfElements( prhs[0] ), p.nROI, &localMultiply, 0.0, (void *) &p, frame, v, w, rates, NULL, 1e-9, 1e-9, 1e8, 10, NULL,
        &istop_out, &itn_out, &anorm_out, &acond_out, &rnorm_out, &arnorm_out, &xnorm_out );

    free(v);
    free(w);
}