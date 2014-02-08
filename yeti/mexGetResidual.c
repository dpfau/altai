#include "mex.h"
#include "localmul.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
	/* Because each frame is stored to single precision to fit in GPU memory, but LSQR works on double precision,
	 * we have to make lots of copies of the data, which is slow. To quickly compute the residual after running
	 * LSQR, we use already-allocated memory, even though that leads to an odd mixing of float and double types
     */
	param p;
    p.ndims = mxGetNumberOfDimensions( prhs[0] );
    p.localDims  = mxGetDimensions( prhs[2] );
    p.globalDims = mxGetDimensions( prhs[0] );
    p.loc = (int *)mxGetData( prhs[3] );
    p.nROI = mxGetNumberOfElements( prhs[4] );
    p.ROIs = mxGetPr( prhs[2] );

    float * frame = (float *)mxGetData( prhs[0] );
    double * residual = mxGetPr( prhs[1] );
    double * rates = mxGetPr( prhs[4] );

    int i;
    for (i = 0; i < mxGetNumberOfElements( prhs[0] ); i++) {
    	residual[i] = (double)frame[i];
    }

    double * negRates = (double *) malloc( p.nROI * sizeof(double) );

    for (i = 0; i < p.nROI; i++) {
    	negRates[i] = -1.0 * rates[i];
    }

    localMultiply(1, 0, 0, negRates, residual, (void *) &p);

    free(negRates);
}