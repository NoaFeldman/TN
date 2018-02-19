/* ==================================================================
 * tstmex.cc
 *
 *    C++ test routine.
 *
 * This is a MEX-file for MATLAB.
 * ================================================================== */

char USAGE[]=
"   Usage: [E,U] = wbeig(H);                                     \n\
                                                                 \n\
       eigenvalue decomposition based on BLAS dsyev              \n\
       for (genearlized) symmetric rank-2 object.                \n\
                                                                 \n\
   AWb C Sep 2009                                                \n";

#include "wblib.h"


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }

    if (nargin!=1 || nargout>2) usage(FL,
    "ERR invalid number of I/O arguments");

    wbarray<double> H(argin[0]),U;
    wbvector<double> E;

    wbEigenS(H,U,E);

    if (nargout>0) argout[0]=E.toMx();
    if (nargout>1) argout[1]=U.toMx();
}

