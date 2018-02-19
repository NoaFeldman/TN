char USAGE[] =
/* ==================================================================
 * vecdata2Data.cc */

"   Usage: A.data = vecdata2Data(vecdata, sizeA)                     \n\
                                                                     \n\
   Inverse operation to vectorizeData()                              \n\
   i.e. put the vectorized data back into its (cell) block structure \n\
                                                                     \n\
   AWb © Feb 2006                                                    \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include <mex.h>
#include <math.h>
#include <string.h>
#include "wblib.h"

#include "qspace_lib.cc"

void mexFunction(
    int nargout, mxArray *plhs[],
    int nargin, const mxArray *prhs[]
){
    wbMatrix<double> QA;
    wbvector<unsigned> Is,dqa,iv;
    wbvector<double> Avec;
    wbMatrix<int> sda;

    if (nargin) if (isHelpIndicator(prhs[0])) { usage(); return; }


    if (nargin!=2)  usage(FLINE, "Invalid number of input arguments.");
    if (nargout>1)  usage(FLINE, "Invalid number of output arguments.");
    if (nargout==0) usage(FLINE, "No ouput argument.");

    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
    usage(FLINE, "Input structures must be double arrays.");

    sda.init0(prhs[1]);


    vecdata2Data(prhs[0], sda, plhs[0]);

    return;
}


