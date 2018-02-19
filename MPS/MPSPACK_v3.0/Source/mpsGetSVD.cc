
char USAGE[] = 
/* =================================================================
 * mpsGetSVD.cc */

"   Usage 1:                                                            \n\
   [U,S,V',[,info]] = mpsGetSVD(Psi1, Psi2, ic1, ic2 [,OPTS]);          \n\
                                                                        \n\
     C++ program to obtain SVD decomposition of QSpaceRM object           \n\
     for the `left' (Psi1, LR) combined with the `right' (Psi2, RL)     \n\
     QSpaceRM.                                                            \n\
                                                                        \n\
   Usage 2:                                                             \n\
   [U,S,V',[,info]] = mpsGetSVD(PSI,idx [,OPTS]);                       \n\
                                                                        \n\
     SVD wrt. indizes given in idx and separate into Psi1.              \n\
                                                                        \n\
   VARIABLES                                                            \n\
                                                                        \n\
     Psi[12]  nearest neighbor MPS blocks    (rank-3, real)             \n\
     ic[12]   indizes that connect Psi1 to Psi2                         \n\
                                                                        \n\
   OPTIONS                                                              \n\
                                                                        \n\
    'disp'     display data to stdout                                   \n\
                                                                        \n\
   OUTPUT                                                               \n\
                                                                        \n\
     Psi1(2)   updated and properly orthonormalized Psi1(2)             \n\
     info      info structure containing                                \n\
      .SVD     singular value data blocked within QS symmetry space     \n\
      .svd     vectorized singular value data                           \n\
                                                                        \n\
   Wb,Sep04,08                                                          \n";

/* This is a MEX-file for MATLAB ================================ //
 *
 *   As the C-code is row-major indizes with dimension 1 (sigma1(2))
 *   will be moved to the first position on the composite object PSI
 *   as this keeps the (L',R') data blocks together
 *
 * AWb © May 2006
 *
 * ============================================================== */

#define _TQ int

#include "wblib.h"
#include "mpsortho.cc"

template<class T1, class T2, class T3>
void MPS_ORTHO(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   const QSpaceRM<_TQ,T1> &Psi1,
   const QSpaceRM<_TQ,T2> &Psi2,
   QSpaceRM<_TQ,T3> &PSI
);

template<class TD>
void MPS_ORTHO(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   QSpaceRM<_TQ,TD> &PSI
);


void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   unsigned two; int r=-1;

   if (nargin) if (isHelpIndicator(argin[0])) { usage(); return; }

   if (nargin<2 || nargout<3 || nargout>4)
   wberror(FL,"Invalid number of I/O arguments.");

   if (!mxIsQSpace(argin[0],r,'C',0,'i')) {
      wbstring istr(str); sprintf(str,
      "%s:%d ERR argument #1 must be a valid QSpaceRM\n%s",
      FL, istr.data); wberror(FL,str);
   }

   r=-1;
   two = mxIsQSpace(argin[1],r,'C',0,'i');

   if (two && nargin<4)
   wberror(FL,"Invalid number of I/O arguments.");


   if (two) {
      const unsigned
      isra=mxIsQSpace(argin[0]),
      isrb=mxIsQSpace(argin[1]);

      if (isra) {
         const QSpaceRM<_TQ,double> Psi1(argin[0],'r',"ref");
         if (isrb) {
            const QSpaceRM<_TQ,double> Psi2(argin[1],'r',"ref");
            QSpaceRM<_TQ,double> PSI;
            MPS_ORTHO(nargout,argout,nargin,argin,Psi1,Psi2,PSI);
         }
         else {
            const QSpaceRM<_TQ,wbcomplex> Psi2(argin[1],'r');
            QSpaceRM<_TQ,wbcomplex> PSI;
            MPS_ORTHO(nargout,argout,nargin,argin,Psi1,Psi2,PSI);
         }
      }
      else {
         const QSpaceRM<_TQ,wbcomplex> Psi1(argin[0],'r');
         if (isrb) {
            const QSpaceRM<_TQ,double> Psi2(argin[1],'r',"ref");
            QSpaceRM<_TQ,wbcomplex> PSI;
            MPS_ORTHO(nargout,argout,nargin,argin,Psi1,Psi2,PSI);
         }
         else {
            const QSpaceRM<_TQ,wbcomplex> Psi2(argin[1],'r');
            QSpaceRM<_TQ,wbcomplex> PSI;
            MPS_ORTHO(nargout,argout,nargin,argin,Psi1,Psi2,PSI);
         }
      }
   }
   else {
      const unsigned
      isra=mxIsQSpace(argin[0]);

      if (isra) {
         QSpaceRM<_TQ,double> PSI(argin[0],'r');
         MPS_ORTHO(nargout,argout,nargin,argin,PSI);
      }
      else {
         QSpaceRM<_TQ,wbcomplex> PSI(argin[0],'r');
         MPS_ORTHO(nargout,argout,nargin,argin,PSI);
      }
   }
}


template<class T1, class T2, class T3>
void MPS_ORTHO(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   const QSpaceRM<_TQ,T1> &Psi1,
   const QSpaceRM<_TQ,T2> &Psi2,
   QSpaceRM<_TQ,T3> &PSI
){
   unsigned i,ic1,ic2,r1,r2, K, disp=0;
   mxArray *S;

   QSpaceRM<_TQ,T3> A1, A2;
   QSpaceRM<_TQ,double> SV;
   OPTS opts;

   if (Psi1.isEmpty()) wblog(FL,"ERR Argument #1 (Psi1) is empty!");
   if (Psi2.isEmpty()) wblog(FL,"ERR Argument #2 (Psi2) is empty!");

   r1=Psi1.rank(); r2=Psi2.rank();

   for (i=2; i<4; i++) if (!mxIsNumber(argin[i])) wblog(FL,
   "ERR argument #%d must be a valid index.", i+1);

   if (mxGetNumber(argin[2], ic1)) wblog(FL,"ERR %s",str);
   if (mxGetNumber(argin[3], ic2)) wblog(FL,"ERR %s",str);

   if (ic1 && ic1<=r1 && ic2 && ic2<=r2) {
       ic1=r1-ic1; ic2=r2-ic2;
   }
   else wblog(FL,
   "ERR arguments #3,#4 - Index out of bounds (%d,%d)",ic1,ic2);

   opts.init(argin+4,nargin-4);

   disp    = opts.getOpt("disp");

   opts.checkAnyLeft();


   twoSiteInit(Psi2, Psi1, PSI, K, ic2, ic1, 0,-1);

   S=mpsGetSVD(
     PSI, A2, SV, A1, K, disp ? "dS" : "s"
   );

   twoSiteFinal(Psi2, Psi1, A2, A1, ic2, ic1, 0,-1);

   argout[0]=A1.toMx('r');
   argout[1]=SV.toMx('r'); if (nargout>2)
   argout[2]=A2.toMx('r'); if (nargout>3)
   argout[3]=S; else if (S) mxDestroyArray(S);
}


template<class TD>
void MPS_ORTHO(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   QSpaceRM<_TQ,TD> &PSI
){
   unsigned i,r,K, disp;
   wbindex I; wbperm P,pa,pb;
   wbvector<unsigned> D;

   QSpaceRM<_TQ,TD> A1, A2;
   QSpaceRM<_TQ,double> SV;
   const mxArray *a; mxArray *S;
   OPTS opts;

   if (PSI.isEmpty())
   wblog(FL,"ERR Argument #1 (PSI) is empty!");

   r=PSI.rank(FL);

   a=argin[1]; {
      unsigned m=mxGetM(a), n=mxGetN(a);

      if (!mxIsDblMat(a) || m>1 && n>1)
      wblog(FL,"ERR (%d,%d)\n%s",m,n,str);

      I.init(FL,argin[1],"invalid index arg #2",'T');
      if (I.len>r) wblog(FL,"ERR Invalid index set arg #2 (%d/%d)",I.len,r);

      for (i=0; i<I.len; i++) if (I[i]<1 || I[i]>r)
      wblog(FL,"ERR Index out of bounds (%d/%d)",I[i],r);
      I-=1;
   }

   opts.init(argin+2,nargin-2);

   disp=opts.getOpt("disp");

   opts.getOpt("permA",S); if (S) pa.init(FL,S,1);
   opts.getOpt("permB",S); if (S) pb.init(FL,S,1);

   if (pa.len>r || pb.len>r || pa.len && pb.len && pa.len+pb.len!=r+2)
   wblog(FL,"ERR Invalid permutations perm[AB] !??");

   opts.checkAnyLeft();

   if (I.len==0 || I.len==r)
   wblog(FL,"ERR invalid index set");


   if (I.extend2Perm(r,P))
   wblog(FL,"ERR Invalid index set%N%N%s%N",str);

   P.FlipIdx(); PSI.Permute(P);
   K=I.len;

   S=mpsGetSVD(PSI, A1, SV, A2, K, disp ? "dSr":"sr");

   if (!pa.isEmpty()) A1.Permute(pa);
   if (!pb.isEmpty()) A2.Permute(pb);

   argout[0]=A1.toMx(); if (nargout>1)
   argout[1]=SV.toMx(); if (nargout>2)
   argout[2]=A2.toMx(); if (nargout>3)
   argout[3]=S; else if (S) mxDestroyArray(S);
}


