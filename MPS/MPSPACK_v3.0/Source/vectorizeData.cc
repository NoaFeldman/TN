char USAGE[] =
/* ==================================================================
 * vectorizeData.cc */

"   Usage: [Q, iA, vecA, sizeA] = vectorizeData(A)                \n\
                                                                  \n\
   vectorize A.data for iterative routines such as bicg();        \n\
                                                                  \n\
   NB! since this routine may be used in different contexts,      \n\
   A.Q may show different ordering; thus the convetion taken      \n\
   here is, that A.data is always resorted with respect to        \n\
   ordered A.Q.                                                   \n\
                                                                  \n\
   Input: A structure with fields {Q, data [, ...]}               \n\
   Output:                                                        \n\
                                                                  \n\
      Q       sorted A.Q data                                     \n\
      iA      index of re-ordering A.Q data to obtain Q           \n\
      vecA    vectorized A.data(iA)                               \n\
      sizeA   original size of A.data arrays                      \n\
                                                                  \n\
   AWb © Feb 2006                                                 \n";

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
    unsigned ma;
    int iD, iQ;

    wbMatrix<double> QA;
    wbvector<unsigned> dqa,iv;
    wbvector<double> Avec;
    wbMatrix<int> sda;
    wbperm Is;

    mxArray *AQ, *ADATA;

    if (nargin) if (isHelpIndicator(prhs[0])) { usage(); return; }


    if (nargin!=1)  usage(FLINE, "Invalid number of input arguments.");
    if (nargout>4)  usage(FLINE, "Invalid number of output arguments.");
    if (nargout==0) usage(FLINE, "No ouput argument.");

    iQ=mxGetFieldNumber(prhs[0], "Q");
    iD=mxGetFieldNumber(prhs[0], "data");

    if (iQ<0 || iD<0)
    usage(FLINE, "Input argument must be structure {Q,data,...}.");

    AQ   =mxGetFieldByNumber(prhs[0],0, iQ);
    ADATA=mxGetFieldByNumber(prhs[0],0, iD);

    check_qdim(AQ, ma, dqa);

    getDataSize(ADATA, dqa.len, ma, sda);
    Data2vec(ADATA, sda, Avec);
 
    getQall(AQ, ma, dqa, QA);
    QA.sortRecs(Is);


    switch (nargout) {
       case 4: sda.mat2mx0(plhs[3]);
       case 3: Avec.mat2mx(plhs[2], 1);
       case 2: iv=Is+1; iv.mat2mx(plhs[1], 1);
       case 1: QtoCell(QA, dqa, plhs[0], -1 );
    }

    return;
}


