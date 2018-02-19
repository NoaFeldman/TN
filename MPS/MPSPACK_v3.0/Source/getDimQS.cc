char USAGE[] =
/* ==================================================================
 * mpsGetDim.cc */

"   Usage: size = mpsGetDim(A)                                    \n\
                                                                  \n\
       get full dimension of given QSpace object                  \n\
                                                                  \n\
   (C) AW : May 2006 ; Oct 2014                                   \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsGetDim"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    wbvector<INDEX_T> D,D2;
    mxArray *a;


    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    if (nargin!=1 || nargout>1) usage(FL,
       "ERR Invalid number of I/O arguments.");

    str[0]=0;

    try { mxIsQSpace(FL,argin[0],'c'); }
    catch (...) {
       if (str[0]) printf("\n%s\n\n",str);
       wblog(FL,"ERR input not a valid QSpace object.");
    }

    if (mxIsQSpace(argin[0])) {
       const QSpace<gTQ,gTD> A(argin[0],'r');
       A.getDim(D,&D2);
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0]);
       A.getDim(D,&D2);
    }

    if (D!=D2) {
       wbvector< wbvector<INDEX_T> const* > dd(2); dd[0]=&D; dd[1]=&D2;
       a=wbMatrix<INDEX_T>().CAT(1,dd).toMx();
    }
    else a=D.toMx();

    argout[0]=a;
}

