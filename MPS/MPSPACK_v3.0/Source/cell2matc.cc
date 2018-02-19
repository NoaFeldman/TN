char USAGE[] = 
/* =============================================================
 * cell2matc.cc */

"   Usage: M = cell2mat(<cellarray>)                           \n\
                                                               \n\
      C++ program that mimics the MatLab function cell2mat.m   \n\
      gain in speed: 2..10                                     \n\
                                                               \n\
   AWb © Jan 2006                                              \n";

/* This is a MEX-file for MATLAB.
 * ============================================================= */

#include <mex.h>
#include <math.h>
#include <string.h>
#include <string.h>

#include "wblib.h"

void check_dim(
    const mxArray *C,
    unsigned *i1, unsigned n1,
    unsigned *i2, unsigned n2,
    unsigned &M, unsigned &N
){

    unsigned i,j,d1,d2;
    mxArray *a;

    M=N=0;

    for (i=0; i<n1; i++)
    for (j=0; j<n2; j++) {
        a=mxGetCell(C,n1*j+i);
        d1=mxGetM(a);
        d2=mxGetN(a);

        if (j>0) {
           if (i1[i]!=d1) {
              char istr[128]; sprintf(istr,
              "Dimension mismatch on input cell array (%d,%d): %dx%d <> %dx%d",
              i+1,j+1, i1[i], j>0 ? i2[j] : d2, d1, d2);
              usage(FLINE, istr);
           }
        }
        else { i1[i]=d1; M+=d1; }

        if (i>0) {
           if (i2[j]!=d2) {
              char istr[128]; sprintf(istr,
              "Dimension mismatch on input cell array (%d,%d): %dx%d <> %dx%d",
              i+1,j+1, i1[i], i2[j], d1, d2);
              usage(FLINE, istr);
           }
        }
        else { i2[j]=d2; N+=d2; }
    }
}

void mexFunction(
    int nargout, mxArray *plhs[],
    int nargin, const mxArray *prhs[]
){
    unsigned int i,j,k,s1,M,N,n1,n2,d1,d2,ir,ic;

    double *dd, *DD, *DI;

    if (nargin!=1)      usage(FLINE, "One cell array as input required.");
    else if (nargout>1) usage(FLINE, "Too many output arguments.");
    else if (!mxIsCell(prhs[0]))
    usage(FLINE, "One cell array as input required.");

    if (nargout==0) return;

    n1=mxGetM(prhs[0]);
    n2=mxGetN(prhs[0]);

    for (k=n1*n2,i=0; i<k; i++)
    if (mxGetClassID(mxGetCell(prhs[0],i)) != mxDOUBLE_CLASS)
    usage(FLINE, "Cell array of type double as input required.");

    unsigned i1[n1], i2[n2];
    check_dim(prhs[0],i1,n1,i2,n2,M,N);

    plhs[0]=mxCreateDoubleMatrix(M,N,mxREAL);

    if (M==0 || N==0) return;

    DD=mxGetPr(plhs[0]);
    if (DD==NULL) {
        wblog(FLINE, "ERR mxCreateDoubleMatrix failed (%dx%d).", M,N);
        wberror(FLINE,"");
    }

    for (ir=i=0; i<n1; i++) {
         d1=i1[i];
         for (ic=j=0; j<n2; j++) {
              d2=i2[j];

              dd=mxGetPr(mxGetCell(prhs[0], j*n1+i));

              DI = DD+M*ic+ir;
              s1 = d1*sizeof(double);

              for (k=0; k<d2; k++)
              memcpy(DI+k*M, dd+k*d1, s1);

              ic+=d2;
         }
         ir+=d1;
    }

    return;
}

