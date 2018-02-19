char USAGE[] = 
/* =============================================================
 * mat2cellc.cc */

"   Usage: A = mat2cellc(cell_array)                           \n\
      C++ program that mimics the MatLab function mat2cell.m   \n\
      gain in speed: ~2                                        \n\
   AWb © Jan 2006                                              \n";

/* This is a MEX-file for MATLAB.
 * ============================================================= */

#include "wblib.h"

inline void check_rdim(double *r, unsigned n, unsigned ntot) {

    unsigned i,k,m=0;
    for (i=0; i<n; i++) {
       k=(unsigned)round(r[i]);
       if (r[i]!=k) usage(FLINE,
         "Input argument 2 and 3 must be positive integer vectors.");
       m+=k;
    }
    if (m!=ntot)
    usage(FLINE, "Input argument 2 and 3 do not match index range.");
}

/* The gateway routine.  */
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
){
    unsigned int i,j,k,s1,M,N,m,n,n1,n2,d1,d2,ir,ic;

    double *i1, *i2, *dd, *DD, *DI, iaux[]={-1,-1};
    int classID;

 /* Check proper input and output */
    if (nrhs!=3)     usage(FLINE, "Three input required.");
    else if (nlhs>1) usage(FLINE, "Too many output arguments.");
    else if (!mxIsDblMat(prhs[0]) ||
             !mxIsDblMat(prhs[1]) ||
             !mxIsDblMat(prhs[2]))
    usage(FLINE, "Input argument 1 must be of type double.");

    M=mxGetM(prhs[0]);
    N=mxGetN(prhs[0]);

    if (mxIsEmpty(prhs[1])) {
       i1=iaux; i1[0]=M; n1=1;
    }
    else {
       m=mxGetM(prhs[1]); n=mxGetN(prhs[1]);
       if (m!=1 && n!=1)
       usage(FLINE, "Input argument 2 and 3 must be (int) vectors.");
       i1=mxGetPr(prhs[1]);
       n1=MAX(m,n);
    }

    if (mxIsEmpty(prhs[2])) {
       i2=iaux+1; i2[0]=N; n2=1;
    }
    else {
       m=mxGetM(prhs[2]); n=mxGetN(prhs[2]);
       if (m!=1 && n!=1)
       usage(FLINE, "Input argument 2 and 3 must be (int) vectors.");
       i2=mxGetPr(prhs[2]);
       n2=MAX(m,n);

       for (i=0; i<n2; i++)
       if (i2[i]!=round(i2[i]))
       usage(FLINE, "Input argument 2 and 3 must be (int) vectors.");
    }

    check_rdim(i1,n1,M);
    check_rdim(i2,n2,N);

 /* input check is ok. */

 /* split data into cell array
  * see also /opt/DIR/matlab-14.3/extern/examples/mx/mxcreatecellmatrix.c */

    plhs[0] = mxCreateCellMatrix(n1,n2);

    DD=mxGetPr(prhs[0]);
    classID=mxGetClassID(prhs[0]);

    for (ir=i=0; i<n1; i++) {
         d1=(unsigned)i1[i];
         for (ic=j=0; j<n2; j++) {
              d2=(unsigned)i2[j];

              k=j*n1+i;
              mxSetCell(plhs[0], k, mxCreateDoubleMatrix(d1,d2,mxREAL));
              dd=mxGetPr(mxGetCell(plhs[0],k));

              if (dd==NULL)
              usage(FLINE, "MatLab could not allocate memory for cell array.");

              DI = DD+M*ic+ir; /* matlab is column major!! */
              s1 = d1*sizeof(double);

           /* mexPrintf("... %d,%d [%d,%d]: %dx%d  %8.4f %8.4f ...\n",
              i,j, ir,ic, d1, d2, DI[0], DI[1]);
            */

              for (k=0; k<d2; k++)
              memcpy(dd+k*d1, DI+k*M, s1);

              ic+=d2;
         }
         ir+=d1;
    }

 /* mexPrintf("\ndim(A) = %dx%d\n", M, N);

    mexPrintf(  "dim(1) = %d:", n1);
    for (i=0; i<n1; i++) mexPrintf(" %2g", i1[i]);

    mexPrintf("\ndim(2) = %d:", n2);
    for (i=0; i<n2; i++) mexPrintf(" %2g", i2[i]);
    mexPrintf("\n\n");
  */

    return;
}

