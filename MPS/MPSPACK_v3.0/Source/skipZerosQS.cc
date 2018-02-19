/* ==================================================================
 * tstmex.cc
 *
 *    C++ test routine.
 *
 * This is a MEX-file for MATLAB.
 * ================================================================== */

char USAGE[]="\n\
Usage: [A,n] = mpsSkipZeros(A [,eps]);                        \n\
                                                              \n\
   Skip data blocks that have all (close to) zero entries     \n\
   and return resulting object. The number of blocks skipped  \n\
   is returned as 2nd argument. Note that elements with       \n\
   norm<eps are considered zero (default: 1E-14).             \n\
                                                              \n\
Wb,May14,10                                                   \n\
";
/* ================================================================== */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsSkipZeros"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    unsigned k,r=-1; char isr; double eps=1E-14;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }
    if (!nargin || nargin>2 || nargout>2) usage(FL,
       "ERR Invalid number of I/O arguments.");

    if (nargin==2)  if (mxGetNumber(argin[1],eps))
    wblog(FL,"ERR reading 2nd argument eps\n%s", str);

    isr=mxIsQSpace(argin[0]); str[0]=0;
    if (!isr) {
       mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
       "invalid QSpace as 1st argument");
    }

    if (isr) {
       QSpace<gTQ,double> A;
       A.init(FL,argin[0],'r'); k=A.SkipZeroData(eps);
       argout[0]=A.toMx();
    }
    else {
       QSpace<gTQ,wbcomplex> A;
       A.init(FL,argin[0]); k=A.SkipZeroData(eps);
       argout[0]=A.toMx();
    }

    if (nargout>1) argout[1]=numtoMx(k);
}

