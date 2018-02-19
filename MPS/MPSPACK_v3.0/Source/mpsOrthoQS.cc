
char USAGE[] = 
/* =================================================================
 * mpsOrthoQS.cc */

"   Usage 1:                                                            \n\
   [Psi1,Psi2,[,info]] = mpsOrthoQS(                                    \n\
       Psi1, Psi2, ic1, ic2, lrdir [,OPTS] );                           \n\
                                                                        \n\
     C++ program that orthonormalizes MPS state using Abelian           \n\
     symmetries for the `left' (Psi1, LR) or the `right' (Psi2, RL)     \n\
     QSpace.                                                            \n\
                                                                        \n\
   Usage 2:                                                             \n\
   [Psi1,Psi2,[,info]] =mpsOrthoQS(PSI,idx, lrdir [,OPTS]);             \n\
                                                                        \n\
     Orthonormalize indizes given in idx and separate into Psi1.        \n\
                                                                        \n\
   VARIABLES                                                            \n\
                                                                        \n\
     Psi[12]  nearest neighbor MPS blocks    (rank-3, real)             \n\
     ic[12]   indizes that connects Psi1 to Psi2                        \n\
                                                                        \n\
     lrdir = 'LR' or 'RL', 12 or 21, +1 or -1, '>>' or '<<'             \n\
                                                                        \n\
        direction (left to right, right to left) that specifies order   \n\
        of SVD on output.                                               \n\
                                                                        \n\
   OPTIONS                                                              \n\
                                                                        \n\
    'Nkeep',.. maximum dimension connecting Psi1 to Psi2 (-1=all)       \n\
    'stol',... tolerance on SVD (i.e. on sqrt(rho_i); 1E-8)             \n\
    'disp'     display data to stdout                                   \n\
                                                                        \n\
   OUTPUT                                                               \n\
                                                                        \n\
     Psi1(2)   updated and properly orthonormalized Psi1(2) using lrflag\n\
     info      info structure containing                                \n\
      .SVD     singular value data blocked within QS symmetry space     \n\
      .svd     vectorized singular value data                           \n\
                                                                        \n\
   AWb (C) Jun 2006 - 2009                                              \n";

/* This is a MEX-file for MATLAB ================================ //
 *
 *   As the C-code is row-major indizes with dimension 1 (sigma1(2))
 *   will be moved to the first position on the composite object PSI
 *   as this keeps the (L',R') data blocks together
 *
 * AW (C) May 2006
 *
 * ============================================================== */

#define QTYPE double
#define _TQ QTYPE

#include "wblib.h"
#include "mpsortho.cc"

template<class TQ, class TD>
void MPS_ORTHO(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   const QSpace<TQ,TD> &Psi1,
   const QSpace<TQ,TD> &Psi2,
   QSpace<TQ,TD> &PSI
);

template<class TQ, class TD>
void MPS_ORTHO(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   QSpace<TQ,TD> &PSI
);


void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   unsigned got12,r=-1;

   if (!nargin || isHelpIndicator(argin[0])) { usage(); return; }

   if (nargin<3 || nargout<2 || nargout>3) wblog(FL,
      "ERR invalid number of I/O arguments (%d/%d)",nargin,nargout);

   mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
      "argument #1 must be a valid QSpace");

   r=-1;
   got12 = mxIsQSpace(0,0,argin[1],r,'c');

   if ((got12 && nargin<5) || nargout>3) wblog(FL, 
      "ERR invalid number of I/O arguments (%d/%d)",nargin,nargout);


   if (got12) {
      const unsigned
      isra=mxIsQSpace(argin[0]),
      isrb=mxIsQSpace(argin[1]);

      if (isra && isrb) {
         const QSpace<_TQ,double> Psi1(argin[0],'r');
         const QSpace<_TQ,double> Psi2(argin[1],'r');
         QSpace<_TQ,double> PSI;
         MPS_ORTHO(nargout,argout,nargin-2,argin+2,Psi1,Psi2,PSI);
      }
      else {
         const QSpace<_TQ,wbcomplex> Psi1(argin[0]);
         const QSpace<_TQ,wbcomplex> Psi2(argin[1]);
         QSpace<_TQ,wbcomplex> PSI;
         MPS_ORTHO(nargout,argout,nargin-2,argin+2,Psi1,Psi2,PSI);
      }
   }
   else {
      const unsigned
      isra=mxIsQSpace(argin[0]);

      if (isra) {
         QSpace<_TQ,double> PSI(argin[0],'r');
         MPS_ORTHO(nargout,argout,nargin-1,argin+1,PSI);
      }
      else {
         QSpace<_TQ,wbcomplex> PSI(argin[0]);
         MPS_ORTHO(nargout,argout,nargin-1,argin+1,PSI);
      }
   }
}

template<class TQ, class TD>
void MPS_ORTHO(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   const QSpace<TQ,TD> &Psi1,
   const QSpace<TQ,TD> &Psi2,
   QSpace<TQ,TD> &PSI
){
   unsigned i,ic1,ic2,r1,r2, K, disp=0;
   double stol=1E-8; int Nkeep=-1;
   char ldir;
   mxArray *S;

   QSpace<TQ,TD> A1, A2;
   OPTS opts;

   if (Psi1.isEmpty()) wblog(FL,"ERR input arg #1 is empty (Psi1) !??");
   if (Psi2.isEmpty()) wblog(FL,"ERR input arg #2 is empty (Psi2) !??");

   r1=Psi1.rank(); r2=Psi2.rank();

   for (i=0; i<2; i++) if (!mxIsNumber(argin[i])) wblog(FL,
   "ERR input arg #%d must be a valid index", i+1);

   if (mxGetNumber(argin[0], ic1)) wblog(FL,"ERR %s",str);
   if (mxGetNumber(argin[1], ic2)) wblog(FL,"ERR %s",str);

   if (ic1 && ic1<=r1 && ic2 && ic2<=r2) {
       ic1--; ic2--;
   }
   else wblog(FL,
   "ERR args #3,#4: index out of bounds (%d/%d,%d/%d)",ic1,r1,ic2,r2);

   ldir=mxGetLRDir(FL,argin[2],4);

   opts.init(argin+3,nargin-3);

   opts.getOpt("stol",  stol );
   opts.getOpt("Nkeep", Nkeep);

   disp    = opts.getOpt("disp");

   opts.checkAnyLeft();

   twoSiteInit(Psi1,Psi2,PSI, K,ic1,ic2, 0,ldir);

   if (Nkeep<0) {
       wbvector<unsigned> D;
       PSI.getDim(D);
       Nkeep=prodRange(D.data,K);
   }

   S=mpsOrthoQS(
      PSI, A1, A2,
      K, Nkeep, stol, disp ? "dS" : "s"
   );

   twoSiteFinal(Psi1,Psi2, A1,A2, ic1,ic2, 0,ldir);

   argout[0]=A1.toMx();
   argout[1]=A2.toMx();

   if (nargout>2) argout[2]=S; else
   if (S) mxDestroyArray(S);
}


template<class TQ, class TD>
void MPS_ORTHO(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   QSpace<TQ,TD> &PSI
){
   unsigned i,r,K, disp; char ldir; int Nkeep=-1;
   wbindex I; wbperm P,pa,pb;
   wbvector<unsigned> D;
   double stol=1E-8;

   QSpace<TQ,TD> A1, A2;
   const mxArray *a; mxArray *S;
   OPTS opts;

   if (PSI.isEmpty())
   wblog(FL,"ERR input arg #1 (PSI) is empty!");

   r=PSI.rank(FL);

   ldir=mxGetLRDir(FL,argin[1],2);

   a=argin[0]; {
      unsigned m=mxGetM(a), n=mxGetN(a);

      if (!mxIsDblMat(a) || (m>1 && n>1))
      wblog(FL,"ERR (%d,%d)\n%s",m,n,str);

      I.init(FL,argin[0],"invalid index arg #2",'T');
      if (I.len>r) wblog(FL,"ERR invalid index set arg #2 (%d/%d)",I.len,r);

      for (i=0; i<I.len; i++) if (I[i]<1 || I[i]>r)
      wblog(FL,"ERR index out of bounds (%d/%d)",I[i],r);
      I-=1;
   }

   opts.init(argin+2,nargin-2);

   opts.getOpt("stol",  stol );
   opts.getOpt("Nkeep", Nkeep);

   disp=opts.getOpt("disp");

   opts.getOpt("permA",S); if (S) pa.init(FL,S,1);
   opts.getOpt("permB",S); if (S) pb.init(FL,S,1);

   if (pa.len>r || pb.len>r || (pa.len && pb.len && pa.len+pb.len!=r+2))
   wblog(FL,"ERR invalid permutations perm[AB] !??");

   opts.checkAnyLeft();

   if (I.len==0 || I.len==r) {
      argout[0]=PSI.toMx();
      argout[1]= A2.toMx();
      
      if (I.len==0) SWAP(argout[0],argout[1]);

      if (nargout>2) {
         S=mxCreateStructMatrix(1,1,0,NULL); sprintf(str,
         "no orthonormalization required (n=%d/%d)",I.len,r);
         mxAddField2Scalar(FL, S, "istr",wbstring(str).toMx());
         argout[2]=S;
      }
      return;
   }

   if (ldir>=0) {
      if (I.extend2Perm(r,P))
      wblog(FL,"ERR invalid index set%N%N%s%N",str);

      PSI.Permute(P); K=I.len;

      if (Nkeep<0) {
         PSI.getDim(D);
         Nkeep=prodRange(D.data,K);
      }

      S=mpsOrthoQS(PSI,A1,A2, K, Nkeep, stol, disp ? "dS":"s");

      if (!pa.isEmpty()) A1.Permute(pa);
      if (!pb.isEmpty()) A2.Permute(pb);

      argout[0]=A1.toMx();
      argout[1]=A2.toMx();
   }
   else {
      if (I.extend2Perm(r,P,"end"))
      wblog(FL,"ERR invalid index set%N%N%s%N",str);

      PSI.Permute(P); K=r-I.len;

      if (Nkeep<0) {
         PSI.getDim(D);
         Nkeep=prodRange(D.data,K);
      }

      S=mpsOrthoQS(PSI,A2,A1, K, Nkeep, stol, disp ? "dS":"s");

      P.initFirstTo(I.len,I.len+1);
      if (pa.isEmpty()) pa=P; else P.times(pa,pa);
      A1.Permute(pa);

      P.initLastTo(0,K+1);
      if (pb.isEmpty()) pb=P; else P.times(pb,pb);
      A2.Permute(pb);

      argout[0]=A1.toMx();
      argout[1]=A2.toMx();
   }

   if (nargout>2) argout[2]=S; else
   if (S) mxDestroyArray(S);
}


