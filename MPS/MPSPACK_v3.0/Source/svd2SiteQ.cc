char USAGE[] =
/* ==================================================================
 * svd2SiteQ - similar to setup2SiteQ() */

"   Usage: CIDX = svd2SiteQ(A4)                                   \n\
                                                                  \n\
     A4 must be a Matrix Product State object containing          \n\
     two sites, i.e. of the structure (L,R,sigma,tau)!            \n\
                                                                  \n\
     The input must obey symmetries for this routine to make      \n\
     sense, e.g. it must hold L+sigma+R+tau = Qtot = const        \n\
     (this is double checked within the routine).                 \n\
                                                                  \n\
     OUTPUT: CIDX is a structure containing the fields            \n\
                                                                  \n\
       Q1    sum of combined index (L+sigma)    dim: ng x 1       \n\
       Q2    sum of combined index (R+tau)      dim: ng x 1       \n\
       IC    cell array of indizes              dim: ng x 1       \n\
                                                                  \n\
     Here I{i} contains the index array                           \n\
                                                                  \n\
         [ i11 i12 ...                                            \n\
           i21 i22 ...                                            \n\
           ... ... imn ]                                          \n\
                                                                  \n\
     such that the rows label different (L,sigma) to the same     \n\
     Q1=L+sigma, and the rows labe different (R,tau) to the same  \n\
     Q2=R+tau. Q1+Q2 adds up to Qtot.                             \n\
                                                                  \n\
     See also: setup2SiteQ.                                       \n\
                                                                  \n\
   AWb © Feb 2006                                                 \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include <mex.h>

#define __WBDEBUG__

#include "wblib.h"
#include "qspace_lib.cc"


void checkQ4SubIdx(
   const wbMatrix<double> &ql,
   const wbMatrix<double> &qr,
   const wbMatrix<double> &qs,
   const wbMatrix<double> &qt,
   const wbvector<unsigned> &ik,
   const unsigned DQ,
   wbMatrix<int> &IM,
   wbMatrix<double> &Q1b,
   wbMatrix<double> &Q2b
){
    unsigned i,k,l,d,d1,d2,i0,s=DQ*sizeof(double);
    double *q;

    wbMatrix<double> QS(qs), QT(qt), Qst, QTk;
    wbvector<unsigned> dg,ii,il,ix,tg;
    wbperm is,it;

    Qst.Cat(2,qs,qt);
    hpsort(Qst,is);

    QS.recPermute(is);
    QS.groupSortedRecs(dg); d1=QS.dim1;

    QT.groupRecs(it,tg); d2=QT.dim1; 

    IM.init(d1,d2); IM.set(-1);

    for (i0=k=0; k<d1; k++,i0+=d) {
         d=dg[k];

         is.get(i0, i0+d-1, il);
         qt.getRecs(il, QTk);

         matchUSortedIdx(QT,QTk,ix,ii);

         if (ii.len!=il.len)
         wberror(FLINE,"Index mismatch !?");

         for (i=0; i<il.len; i++)
         IM(k,ix[i])=ik[il[i]];
    }

    Q1b.init(d1,3*DQ);
    for (i=k=0; k<d1; k++, i+=d) { d=dg[k]; l=is[i]; q=Q1b.rec(k);
        memcpy(q,    ql.rec(l), s);
        memcpy(q+DQ, qs.rec(l), s);
        addRange(q,q+DQ,q+2*DQ,DQ);
    }
   
    Q2b.init(d2,3*DQ);
    for (i=k=0; k<d2; k++, i+=d) { d=tg[k]; l=it[i]; q=Q2b.rec(k);
        memcpy(q,    qr.rec(l), s);
        memcpy(q+DQ, qt.rec(l), s);
        addRange(q,q+DQ,q+2*DQ,DQ);
    }

    return;
}

void mexFunction(
    int nargout, mxArray *plhs[],
    int nargin, const mxArray *prhs[]
){

    unsigned i,l,k,d,ma,n,DQ;
    int e=0;

    wbvector<unsigned> dq4,dg1,ik;
    wbMatrix<int> sd4;
    wbMatrix<double> QL, Qs, QR, Qt, Q1, Q2, QX;
    wbvector< wbMatrix<int> > IM;
    wbvector< wbMatrix<double> > Q1b, Q2b;
    wbperm is1;

    mxArray *A4Q, *A4DATA;


    if (nargin==1) if (isHelpIndicator(prhs[0])) { usage(); return; }


    if (nargin!=1) usage(FLINE, "Invalid number of input arguments.");
    if (nargout>1) usage(FLINE, "Invalid number of output arguments.");
    if (!nargout){ usage(FLINE, "No output argument - done."); return; }

    n=mxGetNumberOfElements(prhs[0]);
    if (!mxIsStruct(prhs[0]) || n!=1) {
        sprintf(str,"Input argument must be%s structure {Q,data,...}.",
        n==1 ? "" : " SINGLE");
        usage(FLINE,str);
    }

    A4Q   =mxGetField(prhs[0],0,"Q");
    A4DATA=mxGetField(prhs[0],0,"data");

    if (A4Q==NULL || A4DATA==NULL) usage(FLINE,
    "Input argument must be structure with fields {Q,data,...}.");

    if (!mxIsCell(A4Q) || !mxIsCell(A4DATA)) usage(FLINE,
    "Input argument must be structure of cells {Q{},data{},...}.");

    check_qdim(A4Q, ma, dq4); DQ=dq4[0];

    getDataSize(A4DATA, dq4.len, ma, sd4);

    for (e=i=0; i<sd4.dim1; i++) 
    if (sd4(i,2)!=1 || sd4(i,3)!=1) e++;

sd4.put("sd4");
dq4.put("dq4"); wblog(FL,"TST e=%d", e);

    if (e || dq4.len!=4)
    usage(FLINE, "Input MPS data must have (L,R,sigma,tau) Q space order.");

    for (i=0; i<4; i++) if (dq4[i]!=DQ)
    usage(FLINE, "Input MPS must have equal Q-dim for every index.");

    getQBlock(A4Q,0,QL);
    getQBlock(A4Q,1,QR);
    getQBlock(A4Q,2,Qs); Q1=QL; Q1.add(Qs);
    getQBlock(A4Q,3,Qt); Q2=QR; Q2.add(Qt);

    QX=Q1; QX.add(Q2);
    for (i=1; i<QX.dim1; i++)
    if (QX.recCompare(i,0)) {
        QX.recPrint(0, "Qtot"); QX.recPrint(i, "Qtot"); usage(FLINE,
        "Input MPS must preserve symmetry (Qtot is not unique).");
    }

    Q1.groupRecs(is1, dg1); n=dg1.len;

    IM.init(n); Q1b.init(n); Q2b.init(n);

    for (l=k=0; k<n; k++,l+=d) {
         d=dg1[k];
         ik.init(d, is1.data+l);

         checkQ4SubIdx(
           QL.getRecs(ik), QR.getRecs(ik),
           Qs.getRecs(ik), Qt.getRecs(ik), ik, DQ, IM[k], Q1b[k], Q2b[k]
         );
    }


  { const char *fnames[]={"Q1","Q2","IC"};
    plhs[0]=mxCreateStructMatrix(1,1,3,fnames);
  }

    matvec2cell(Q1b, plhs[0], 0);
    matvec2cell(Q2b, plhs[0], 1);
    matvec2cell(IM,  plhs[0], 2);

    return;
}


