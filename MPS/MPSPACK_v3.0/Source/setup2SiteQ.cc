char USAGE[] =
/* ==================================================================
 * setup2SiteQ
 *
 *    Similar to contractIndex() but with complete local basis
 *    taking into account global symmetry Qtot */

"USAGE                                                          \n\
                                                                \n\
   [QAB,CIDX,permA,permB] = setup2SiteQ(A,B,Qsigma,Qtau,Qtot)   \n\
                                                                \n\
   A and B must be Matrix Product State objects                 \n\
   i.e. A.Q={L1, R1, sigma}                                     \n\
        B.Q={L2, R2, tau  }                                     \n\
   that is contracted over R1:L2                                \n\
   resulting in the object {L(1),R(2),sigma,tau}                \n\
                                                                \n\
   Qsigma and Qtau are unique arrays providing full local       \n\
   state space.                                                 \n\
                                                                \n\
   The index set {L,R,sigma,tau} will be expanded/reduced       \n\
   in {sigma,tau} space to get all combinations in {sigma,tau}  \n\
   that provide given Qtot.                                     \n\
                                                                \n\
OUTPUT                                                          \n\
                                                                \n\
   QAB    contracted and completed combined Q space of A and B  \n\
    .Q    actual Q space with order (L,R,sigma, tau)            \n\
    .data size information                                      \n\
    .Sab  size information                                      \n\
                                                                \n\
   CIDX   index set IG with the following data in its rows      \n\
    .IG  [ number of already existing entries in Q4             \n\
           index 1 into matched (degenerate) Q4                 \n\
           index 2 into matched (degenerate) Q4                 \n\
           ... ]                                                \n\
    .Ia   in debug mode, index match on A (prior to Q check)    \n\
    .Ib   in debug mode, index match on B (prior to Q check)    \n\
    .Iu   in debug mode, group sizes in Ia and Ib               \n\
                                                                \n\
   permA  permuation that needs to be done on A data            \n\
   permB  permuation that needs to be done on A data            \n\
                                                                \n\
AWb © Jan 2006                                                  \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include <mex.h>

#define __WBDEBUG__

#include "wblib.h"
#include "qspace_lib.cc"

void initSigmaSpace(
    const mxArray *a,
    const unsigned &DQ,
    wbMatrix<double> &QSIG
){
    unsigned i,j,k=0,d1=mxGetM(a), d2=mxGetN(a);
    double* d;

    if (!mxIsDouble(a) || d1<1 || d2!=DQ)
    usage(FLINE,"Dimension mismatch on local state spaces.");

    QSIG.init(d1,d2); d=mxGetPr(a);

    for (j=0; j<d2; j++)
    for (i=0; i<d1; i++) QSIG(i,j)=d[k++]; 

    return;
}


void checkQtot(const mxArray *a, wbvector<double> &Qtot) {

    if (!mxIsDouble(a) || mxGetM(a)!=1 && mxGetN(a)!=1)
    usage(FLINE,"Third argument must be double vector");

    Qtot.init(mxGetNumberOfElements(a), mxGetPr(a));

    for (unsigned i=0; i<Qtot.len; i++) {
       if (Qtot[i]<0 || floor(Qtot[i])!=Qtot[i])
       usage(FLINE,"Qtot must be index vector.");
    }
}


int complete2SiteIndex(
    wbMatrix<double> Q4,          
    const wbMatrix<double> &QSIG,
    const wbMatrix<double> &QTAU,
    const wbvector<double> &Qtot, 
    wbvector<unsigned> &g4,
    wbMatrix<double> &QNEW,
    wbMatrix<unsigned> &IG
){
    unsigned i,j,k,kk,l,n,ng,K, d,d2, s,s2,s4, iQ, *iq,*c, e=0;
    unsigned ns=QSIG.dim1, nt=QTAU.dim1, D4=Q4.dim2, NQ, DQ=Qtot.len;

    wbvector<unsigned> mark,dLR,dST,dd,ig4,Ist,uu;
    wbvector<int> u1,u2;
    wbvector<double> Q1;
    wbperm is4, IST;

    wbvector<unsigned> miLR;
    wbMatrix<double> QQ,QLR,QST, QS, QLRs;

    double *q;

    if (D4!=4*DQ || QTAU.dim2!=DQ || QSIG.dim2!=DQ)
    wberror(FLINE,"Dimensions mismatch in Q space.");

    if (!isUnique(QSIG) || !isUnique(QTAU))
    wberror(FLINE,"Local state space must have unique Q labels.");


    Q4.groupRecs(is4,g4); K=g4.max();

    Q4.getBlocks(0,1,DQ,QLR);
    QLR.groupSortedRecs(dLR);
    QLR.blockSum(0,1,DQ,QLRs);

    NQ=Q4.dim1; Ist.init(NQ);
    ig4.init(g4.len);
    for (i=1; i<ig4.len; i++) ig4[i]=ig4[i-1]+g4[i-1];

    Q4.getBlock(2,DQ, QQ); e+=QQ.findRecsInSet(QSIG, u1);
    Q4.getBlock(3,DQ, QQ); e+=QQ.findRecsInSet(QTAU, u2);

    if (e) wberror(FLINE, "sigma or tau not in local state space.");

    for (k=0; k<NQ; k++) {
        i=(unsigned)u1[k]; j=(unsigned)u2[k];
        Ist[k]=i*nt+j;
    }

    QST.init(ns*nt, 2*DQ); s=DQ*sizeof(double);

    for (k=i=0; i<ns; i++)
    for (  j=0; j<nt; j++) { q=QST.rec(k++);
        memcpy(q,    QSIG.rec(i), s);
        memcpy(q+DQ, QTAU.rec(j), s);
    }

    QST.blockSum(0,1,DQ, QS);
    QS.groupRecs(IST,dST); 

    QST.recPermute(IST);

    getIPerm(IST,uu);
    for (i=0; i<NQ; i++) Ist[i]=uu[Ist[i]];


    mark.init(IST.len); miLR.init(IST.len);

    IG.init(K*MIN(ns,nt)*NQ, 2+K); 

    iQ=0; ng=dLR.len;

    for (kk=k=0; k<ng; k++) { 
        mark.reset(); miLR.reset();

        Q1=Qtot; Q1.subtract(QLRs.rec(k));

        for (l=i=0; i<QS.dim1; i++, l+=d) {
            d=dST[i];

            if (Q1.isEqual(QS.rec(i)))
            for (c=mark.data+l, j=0; j<d; j++) c[j]+=1;
        }

        d=dLR[k];
        for (i=0; i<d; i++) {
             j=Ist[kk+i];  
             if (mark[j]==0) continue;

             mark[j]++;
             miLR[j]=kk+i;

             if (mark[j]>2) wblog(FLINE, "ERR Index mismatch.");
        }

        for (i=0; i<mark.len; i++) {
            if (mark[i]==0) continue;


            iq=IG.rec(iQ++); iq[0]=i;

            if (mark[i]==1) {
                iq[2]=kk; 
            }
            else {
                j=miLR[i]; iq[1]=g4[j]; iq[2]=j;
            }

            if (iQ>=IG.dim1)  
            wblog(FLINE, "ERR Running out of allocated space (%dx%dx%d=%d)",
            K, MIN(ns,nt), NQ, K*MIN(ns,nt)*NQ);
        }

        kk+=d;
    }

    IG.Resize(iQ, IG.dim2);

    n=IG.dim1; d2=2*DQ; s2=d2*sizeof(double); s4=4*DQ*sizeof(double);

    QNEW.init(n,D4);

    for (k=0; k<n; k++) {
        q=QNEW.rec(k); iq=IG.rec(k);

        if (iq[1]) {
            memcpy(q, Q4.rec(iq[2]), s4);
        }
        else {
            memcpy(q,    Q4 .rec(iq[2]), s2);
            memcpy(q+d2, QST.rec(iq[0]), s2);
        }

        j=iq[2]; l=ig4[j]; d=iq[1]; if (d==0) d=1;

        for (iq+=2,i=0; i<d; i++) iq[i]=is4[l+i];
    }

    IG.set2Cols(1,IG.dim2-1); 

    return 0;
}

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
){
    unsigned i,k,e,m,ma,mb, DQ;
    int ifield_Q, ifield_D;

    wbvector<unsigned> ica,icb,ika,ikb,ikP,isa,isb,isk,dqa,dqb,dab;
    wbvector<unsigned> Isa, Isb, ia, ib, Ia, Ib, Iu, uu;
    wbvector<unsigned> iL, iR, is, it, isnew;
    wbvector<double> Qtot;

    wbMatrix<double> QL,Q1,Qs, Q2,QR,Qt, QAB, Q4, QX, QSIG, QTAU;
    wbMatrix<unsigned> IG,UU;
    wbMatrix<int> sda,sdb;

    mxArray *AQ, *ADATA, *BQ, *BDATA;


    if (nrhs==1)  if (isHelpIndicator(prhs[0])) { usage(); return; }


    if (nrhs!=5) usage(FLINE, "Three input arguments required.");
    if (nlhs> 2) usage(FLINE, "Max two output arguments available.");

    if (!mxIsStruct(prhs[0]) || !mxIsStruct(prhs[1])) usage(FLINE,
    "First two input arguments must be structures {Q,data,...}.");

    if (mxGetNumberOfElements(prhs[0])!=1 ||
        mxGetNumberOfElements(prhs[1])!=1) usage(FLINE,
    "First two input arguments must be SINGLE structures {Q,data,...}.");

    ifield_Q=mxGetFieldNumber(prhs[0], "Q");
    ifield_D=mxGetFieldNumber(prhs[0], "data");

    if (ifield_Q<0 || ifield_D<0) usage(FLINE,
    "First two input arguments must be structures {Q,data,...}.");

    AQ   =mxGetFieldByNumber(prhs[0], 0, ifield_Q);
    ADATA=mxGetFieldByNumber(prhs[0], 0, ifield_D);

    ifield_Q=mxGetFieldNumber(prhs[1], "Q");
    ifield_D=mxGetFieldNumber(prhs[1], "data");

    if (ifield_Q<0 || ifield_D<0) usage(FLINE,
    "First two input arguments must be structures {Q,data,...}.");

    BQ   =mxGetFieldByNumber(prhs[1], 0, ifield_Q);
    BDATA=mxGetFieldByNumber(prhs[1], 0, ifield_D);

    check_qdim(AQ, ma, dqa); DQ=dqa[0];
    check_qdim(BQ, mb, dqb);

    if (!mxIsCell(ADATA) || !mxIsCell(BDATA)) usage(FLINE,
    "First two input arguments must be structures {Q{},data{},...}.");

    getDataSize(ADATA,dqa.len,ma,sda);
    getDataSize(BDATA,dqb.len,mb,sdb);

    e=0; k=2;
    for (i=0; i<sda.dim1; i++) if (sda(i,k)!=1) e++;
    for (i=0; i<sdb.dim1; i++) if (sdb(i,k)!=1) e++;
    if (e) {
        sda.put("sda"); sdb.put("sdb");
        wberror(FLINE, "Local state space has dimension > 1 !?");
    }

    initSigmaSpace(prhs[2],DQ,QSIG);
    initSigmaSpace(prhs[3],DQ,QTAU);

    checkQtot(prhs[4],Qtot);

    if (Qtot.len!=DQ) usage(FLINE,"Length mismatch on third argument.");

    if (dqa.len!=3 || dqb.len!=3)
    usage(FLINE, "A and B must be 3D MPS block structures.");

    for (i=1; i<3; i++)
    if (dqa[i]!=DQ || dqb[i]!=DQ)
    usage(FLINE, "Q-dim mismatch on A or B.");

    ica.init(1); ica[0]=1;
    icb.init(1); icb[0]=0;

    invertIndex(ica,dqa.len,ika);
    invertIndex(icb,dqb.len,ikb);

  { unsigned pp[4]={0,2,1,3}; 
    ikP.init(4,pp); }


    getQBlock(AQ,0,QL);
    getQBlock(AQ,1,Q1);
    getQBlock(AQ,2,Qs);

    getQBlock(BQ,0,Q2);
    getQBlock(BQ,1,QR);
    getQBlock(BQ,2,Qt);

    e=0; QX.Cat(2,QL,Q1,Qs); if (!isUnique(QX)) { QX.put("Q4a"); e+=1; }
         QX.Cat(2,Q2,QR,Qt); if (!isUnique(QX)) { QX.put("Q4b"); e+=2; }
    if (e) wberror(FLINE, "State spaces (LRs) must be unique.");

    matchIndex(Q1,Q2,Ia,Ib);

    Q4.Cat(2,
    QL.getRecs(Ia), QR.getRecs(Ib), Qs.getRecs(Ia), Qt.getRecs(Ib));

    complete2SiteIndex(
       Q4, QSIG, QTAU, Qtot,
       Iu, QAB, IG    
    );


    dab.init(ika.len+ikb.len); m=ika.len;
    for (i=0; i<ika.len; i++) dab[i  ]=dqa[ika[i]];
    for (i=0; i<ikb.len; i++) dab[i+m]=dqb[ikb[i]];

    if (nlhs>=1) {
        const char *fnames[]={"Q", "Sab", "SAB"};

        plhs[0]=mxCreateStructMatrix(1,1,3,fnames);

        QtoCell(QAB,dab, plhs[0], 0 );

        Sab2Site2Cell(
           sda, sdb, Ia, Ib, IG, ika, ikb, ikP,
           plhs[0],1,2
        );

    if (nlhs>=2) {  

      #ifdef __WBDEBUG__
        const char *fnames[]={"IG","Ia","Ib","Iu"};
        plhs[1]=mxCreateStructMatrix(1,1,4,fnames);

        uu=Ia+1; uu.mat2mxs(plhs[1],1);
        uu=Ib+1; uu.mat2mxs(plhs[1],2);
                 Iu.mat2mxs(plhs[1],3);
      #else
        const char *fnames[]={"IG"};
        plhs[1]=mxCreateStructMatrix(1,1,1,fnames);
      #endif

        for (UU=IG,k=0; k<UU.dim1; k++)
        for (i=1; i==1 || i<=UU(k,0); i++) UU(k,i)+=1;
        UU.mat2mxs(plhs[1], 0 );

    if (nlhs>=3) {
        uu=ika; uu.cat(ica); uu+=1; uu.mat2mx(plhs[2]);

    if (nlhs>=4) {
        uu=icb; uu.cat(ikb); uu+=1; uu.mat2mx(plhs[3]);
    }}}}


    return;
}


