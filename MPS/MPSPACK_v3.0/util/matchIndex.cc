char USAGE[] =
/* ==================================================================
 * matchIndex.cc
 *
 *    C++ program that gets index sets for generic contractQ() */

"   Usage: [Ia,Ib,I] = matchIndex(QA, QB [,OPTS])                   \n\
                                                                    \n\
   QA or QB are matrizes whos rows are matched; consequently        \n\
   their row dimension must be equal. If Q[AB] are QSpaces,         \n\
   all of their Q{} data is merged into a single matrix.            \n\
                                                                    \n\
   Options                                                          \n\
                                                                    \n\
      '-s'    sort indizes with respect to A                        \n\
      '-q'    quiet mode, i.e. do not issue warning when NaN's      \n\
              are encountered.                                      \n\
      'eps',..specify epsilon within which to accept match (0)      \n\
      '-f'    use float instead of double (gets rid of numerical    \n\
              noise in case of double precision input data)         \n\
                                                                    \n\
   Output: CIDX is a structure with elements                        \n\
                                                                    \n\
      Ia      index into QA (QA.Q)                                  \n\
      Ib      index into QB (QB.Q)                                  \n\
      m       first m of Ia and Ib do overlap                       \n\
                                                                    \n\
   AWb * Feb 2006                                                   \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include <mex.h>
#include <math.h>
#include <string.h>
#include "wblib.h"

#include "qspace_lib.cc"

/* ----------------------------------------------------------------------- *
 * THE GATEWAY ROUTINE
 * ----------------------------------------------------------------------- */
void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    wbMatrix<double> QA,QB;
    wbMatrix<wbcomplex> ZA,ZB; // leave here! QA,QB initialized to ref if used!
    wbindex In1,In2,ia,ib;
    unsigned nan1,nan2;
    double eps=0;

    char sflag=0, fflag=0, qflag=0;

    if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }

    if (nargin<2 || nargout>3)
    usage(FL,"Invalid number of I/O arguments.");

    if (nargin>2) {
       OPTS opts; opts.init(argin+2, nargin-2);

       fflag=opts.getOpt("-f");
       qflag=opts.getOpt("-q");
       sflag=opts.getOpt("-s");
       opts.getOpt("eps", eps);

       opts.checkAnyLeft(); // also closes parameter file if open
    }

 /* if (mxIsQSpace(argin[0])) { // real
       if (mxIsStruct(argin[0])) getQall(argin[0],QA); else QA.init(argin[0]);
       if (mxIsStruct(argin[1])) getQall(argin[1],QB); else QB.init(argin[1]);

       if (QA.dim1==1 && QB.dim1==1 && QA.dim2 != QB.dim2) {
           SWAP(QA.dim1,QA.dim2);
           SWAP(QB.dim1,QB.dim2);
       }
    } else */
    if (!mxIsComplex(argin[0]) && !mxIsComplex(argin[1])) {
       QA.init(argin[0]);
       QB.init(argin[1]);

       if (QA.dim1==1 && QB.dim1==1 && QA.dim2 != QB.dim2) {
           SWAP(QA.dim1,QA.dim2);
           SWAP(QB.dim1,QB.dim2);
       }
    }
    else {
    // hard to compare complex numbers -> map complex to real array
       ZA.init(argin[0]);
       ZB.init(argin[1]);

       if (ZA.dim1==1 && ZB.dim1==1 && ZA.dim2 != ZB.dim2) {
           SWAP(ZA.dim1,ZA.dim2);
           SWAP(ZB.dim1,ZB.dim2);
       }

       QA.init2ref(ZA.dim1,2*ZA.dim2,(double*)ZA.data);
       QB.init2ref(ZB.dim1,2*ZB.dim2,(double*)ZB.data);
    }

    if (QA.dim2!=QB.dim2) wblog(FL,
    "ERR column size mismatch (%d/%d)", QA.dim2, QB.dim2);

    if (fflag) {
       QA.SkipTiny_float(); // 0 = plain double->float->double conversion
       QB.SkipTiny_float();
    }

    nan1=QA.skipNanRecs(In1);
    nan2=QB.skipNanRecs(In2);

    if (!qflag) if (nan1 || nan2) wblog(FL,
       "WRN skipping %d record(s) containing NaN",nan1+nan2);

    matchIndex(QA,QB,ia,ib,1,NULL,NULL,eps);

    if (sflag) {
       wbperm p; Wb::hpsort(ia,p); ib.Select(p);

    // in case of degenerate index ia, sub-sort ib
       for (unsigned i=1,j=0; i<ia.len; i++) {
          if (ia[i-1]!=ia[i]) {
             if (j+1<i) {
                wbindex jj; jj.init2ref(i-j, ib.data+j);
                Wb::hpsort(jj,p);
             }
             j=i; continue;
          }
       }
    }

    if (nan1) { In1.Select(ia); In1.save2(ia); }
    if (nan2) { In2.Select(ib); In2.save2(ib); }

    if (nargout>=3) {
       mxArray *S=mxCreateStructMatrix(1,1,0,NULL);
       wbindex iax,ibx; wbvector<unsigned> sz(2);
       
    // get index sets for extra space (0: i[ab] possibly not unique)
       invertIndex(ia, QA.dim1, iax, 0);
       invertIndex(ib, QB.dim1, ibx, 0);

    // NB! toMx() below already makes index 1-based! // Wb,Feb10,11
    // iax+=1; ibx+=1; // make 1-based

       sz[0]=QA.dim1;
       sz[1]=QB.dim1;

       mxAddField2Scalar(FL,S,"n",   numtoMx(ia.len));
       mxAddField2Scalar(FL,S,"dim", sz.toMx());
       mxAddField2Scalar(FL,S,"ix1", iax.toMx());
       mxAddField2Scalar(FL,S,"ix2", ibx.toMx());

       argout[2]=S;
    }

 // NB! toMx() below already makes index 1-based! // Wb,Feb10,11
 // ia+=1; ib+=1; // 1-based index for MatLab

    argout[0]=ia.toMx(); // MX_USING_ANS
    if (nargout>1) argout[1]=ib.toMx();

} // EOM


/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

