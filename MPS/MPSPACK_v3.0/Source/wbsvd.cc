/* ==================================================================
 * tstmex.cc
 *
 *    C++ test routine.
 *
 * This is a MEX-file for MATLAB.
 * ================================================================== */

char USAGE[]=
"   Usage: [U,S,V] = wbsvd(H [,I]);                                \n\
                                                                 \n\
       singular value decomposition based on BLAS dgesvd         \n\
       for genearlized rank-r object.                            \n\
       index set I specifies which dimensions are combined       \n\
       into first fused dimension. For rank-2 objects            \n\
       and I omitted, this corresponds to the regular SVD        \n\
       on matrizes.                                              \n\
                                                                 \n\
   AWb C Sep 2009                                                \n";

#include "wblib.h"


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }

    if (nargin<1 ||nargin>2 || nargout>3) usage(FL,
    "ERR wbsvd() - invalid number of I/O arguments");

    wbarray<double> A(argin[0]),U,V;
    wbvector<double> S;
    wbvector<unsigned> I;
    wbperm P;

    unsigned r=A.SIZE.len;

    if (nargin>1) {
       I.init(FL,argin[1]);
       if (I.len<1 || I.len>r || I.anyST(1) || I.anyGT(r)) wblog(FL,
       "ERR %s - invalid dim-specs (%s)",__FUNCTION__,I.toStr().data);
       I-=1;
    }

    wbSVD(A,U,S,V,I);

    if (nargout<=1) argout[0]=S.toMx();
    else {
       argout[0]=U.toMx();
       argout[1]=S.toMx(); if (nargout>2) {
       argout[2]=V.toMx(); }
    }
}

