/* ==================================================================
 * tstmex.cc
 *
 *    C++ test routine.
 *
 * This is a MEX-file for MATLAB.
 * ================================================================== */

char USAGE[]="\n\
  function S=sparse2wb(S)                                    \n\
                                                             \n\
     convert matlab sparse array to wbsparse format.         \n\
                                                             \n\
  (C) Wb,Dec09,11                                            \n\
";
/* ================================================================== */

#include "wblib.h"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "sparse2wb"
#endif

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   if (nargin && isHelpIndicator(argin[0])) { usage(); return; }
   if (nargin!=1 || nargout>1 || !mxIsDouble(argin[0])) wblog(FL,
      "ERR %s() invalid usage",PROG);

   wbsparray<unsigned,double> A(FL,argin[0]);
   argout[0]=A.toMx();
};

