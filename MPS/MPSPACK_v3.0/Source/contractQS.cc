
// => update usage

char USAGE[] =
/* ================================================================== */
"                                                                    \n\
   Usage #1: S=contractQS(A, ica, B, icb [, perm, OPTS ]);             \n\
                                                                     \n\
       Plain contraction of two QSpaces, A and B, with respect       \n\
       to given index sets while also respecting Clebsch Gordan      \n\
       coefficient spaces if present.                                \n\
                                                                     \n\
       Options specific to usage #1:                                 \n\
                                                                     \n\
         'conjA'  use complex conjugate of A for contraction         \n\
         'conjB'  use complex conjugate of B for contraction         \n\
                                                                     \n\
       Contraction indices ica and icb can also be specified           \n\
       as strings which allows to simultaneously also specify        \n\
       conj flag [using semantics of ctrIdx() class]:                \n\
                                                                     \n\
         contractQS(A,[1 3],B,[1 3],'conjA')  isequivalent to        \n\
         contractQS(A,'1,3;*',B,[1 3])                               \n\
                                                                     \n\
   Usage #2: S=contractQS({A,{B,C}},... [, perm, OPTS ]);            \n\
                                                                     \n\
       Generalized cell-contraction of tensors based on their        \n\
       itags (index labels; see info.itags). These must be unique    \n\
       such that matching index spaces define the contraction.       \n\
       The trailing options affect the entire contraction.           \n\
       In particular, perm stands for the permutation to be used     \n\
       with the final object.                                        \n\
                                                                     \n\
       The order of contractions can be controlled by grouping       \n\
       contractions pairwise using cell structures (thus the         \n\
       name cell-contraction).                                       \n\
                                                                     \n\
       For every operator, additional optional strings can be        \n\
       specified, appearing right after the affected QSpace (e.g. A):\n\
                                                                     \n\
         A,'!kk'  do not contract indices specified by kk            \n\
                  despite they share common matching itags.          \n\
         A,'*'    apply overall (complex) conjugation of input       \n\
                  tensor A (note that this also includes reverting   \n\
                  the directions of all indizes!)                    \n\
         A,'!kk*' both of the above (with * always trailing).        \n\
                                                                     \n\
       Further options for individual QSpaces:                       \n\
                                                                     \n\
         A,'-op:<tag>'                                               \n\
            this sets operator itags '<tag>;<tag>*[;op]' for QSpace A\n\
            and issues a warning, if existing itags are overwritten. \n\
                                                                     \n\
       Usage #2 also allows set of sequential contractions           \n\
       such as in S=contractQS(A,B,C,... [, perm, OPTS ]); this is   \n\
       interpreted as S=contractQS({A,{B, {C,...}}}, [,perm,OPTS ]), \n\
       i.e. the sequential contractions are started from the end     \n\
       onwards to the beginning of the set.                          \n\
                                                                     \n\
       Options specific to usage #2:                                 \n\
                                                                     \n\
         '-v'  verbose mode that shows level of cell contraction     \n\
               together with actual contractions performed.          \n\
                                                                     \n\
   Note that mixed usage of #2 and #1 is not possible.               \n\
   General options                                                   \n\
                                                                     \n\
       perm    permutation of dimensions after contraction           \n\
                                                                     \n\
   AW (C) May 2010 ; Aug 2012 ; Apr 2013 ; Dec 2014                  \n";





#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "contractQS"
#endif


#define LOAD_CGC_QSPACE
#include "wblib.h"


template<class TA, class TB, class TC>
QSpace<gTQ,TC>& CONTRACT_CG(
    const QSpace<gTQ,TA> &A, const ctrIdx &ica,
    const QSpace<gTQ,TB> &B, const ctrIdx &icb,
    QSpace<gTQ,TC> &C,
    wbperm P
);

char isComplexQS(unsigned na, const mxArray *aa[], unsigned level=0);


class icFlags {
 public:

   icFlags(const mxArray* ain=NULL) : a(ain), conj(0) {
      checkQSpace(a); };

  ~icFlags() {};

   icFlags& checkQSpace(const mxArray *a);

   char check_arg(const mxArray *ain) { a=ain;
      if (a) {
         if (mxIsCell(a)) { return 3; }
         if (mxIsChar(a)) { return (set(a,0) ? 0:1); }
         if (mxIsQSpace(0,0,a,'c')) return 2;
      }
      return 0;
   };

   icFlags& init(const mxArray *ain) {
      if (!ain) wblog(FL,"ERR %s() got NULL QSpace or cell !??",FCT);
      a=ain; checkQSpace(a);
      ktags.init(); otags.init(); conj=0;
      return *this;
   };

   char set(const mxArray *a, const char *istr="");
   
   template <class TQ, class TD>
   void apply(QSpace<TQ,TD> &X);

   const mxArray *a;
   wbstring ktags;
   wbstring otags;"<otag>;<otag>*;op"
   char conj;


 protected:
 private:

};


icFlags& icFlags::checkQSpace(const mxArray *a) {

   if (a) {
      if (mxIsCell(a)) {
         if (mxGetNumberOfElements(a)<2) wblog(FL,
           "ERR %s() invalid cell-contraction\n"
           "(got cell with fewer than two entries (%d)!)",
            FCT, mxGetNumberOfElements(a)
         );
      }
      else {
         if (!mxIsQSpace(FL,a,'c')) wblog(FL,
           "ERR %s() invalid QSpace in cell-contraction (%s)",
            FCT, mxGetClassName(a));
         if (mxGetNumberOfElements(a)!=1) wblog(FL,
           "ERR %s() invalid QSpace(%d) in cell-contraction",
            FCT, mxGetNumberOfElements(a)
         );
      }
   }
   return *this;
};


char icFlags::set(const mxArray *a, const char *istr) {

   wbstring S(FL,a); const char *s=S.data;

   if (!s || !s[0]) {
      if (istr) wblog(FL,"ERR %s() %s\n"
        "invalid string option for cell-contraction (%s)",
         FCT, istr, s? "empty":"null");
      else return 1;
   }

   if (!strncasecmp(s,"-op:",4)) { otags=s+1; }
   else {
      unsigned l=0;
      if (s[0]=='!') {
         for (++s; s[l]; ++l) { if (!isdigit(s[l])) break; }
      }
      if (s[l]) {
         if (s[l]==CC_ITAG && !s[l+1]) { conj=1; }
         else {
            if (istr) wblog(FL,
               "ERR %s() %s\ninvalid cell-contraction option\n"
               "(expecting '[!##][*]': got %s)",FCT,istr,S.data);
            else return 1;
         }
      }
      ktags=s;
      if (conj) ktags.data[l]=0;
   }
      
   return 0;
};


template <class TQ, class TD>
void icFlags::apply(QSpace<TQ,TD> &X) {


   if (otags.len) {
      unsigned r=X.rank(FL);
      if (r!=2 && !(r==3 && X.otype==QS_OPERATOR)) wblog(FL,
         "ERR %s() cannot set operator itags\n"
         "for rank-%d QSpace (%s)", FCT,r,X.otype2Str().data
      );
      X.init_itags(FL,otags.data,-1);
   }

   if (ktags.len) {
      unsigned i=0, k, r=X.rank(FL);
      if (r!=X.itags.len) wblog(FL,
         "ERR %s() valid set of info.itags required (%s; %d)",
         FCT,IT2STR(X),r
      );

      for (; i<ktags.len; ++i) { k=ktags[i]; if (!k) break;
         k-='0';
         if (!k || k>r) wblog(FL,"ERR %s() index in cc-string "
            "out of bounds (%d: '%s'; %d/%d)",FCT,i,ktags.data,k,r);
         X.SetFlag(FL,k,1);
      }
   }
};


mxArray* contractQS_plain(int nargin, const mxArray *argin[]);

template<class TD>
unsigned contractQS_itags(
   int nargin, const mxArray *argin[],
   unsigned level, unsigned vflag, QSpace<gTQ,TD> &C
);


void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){

    if (nargin && isHelpIndicator(argin[0])) { usage(); return; }

#ifdef __WBDEBUG__
    Wb::MemCheck(FL,"start");
#endif

    if (nargin<1) wblog(FL,
       "ERR invalid number of input arguments (%d)",nargin);
    if (nargout>1) wblog(FL,"ERR invalid number of output arguments");

    mxArray *a=NULL;

    if (nargin>1 && isCtrIdx(argin[1])) {
       a=contractQS_plain(nargin,argin);
    }
    else {
       unsigned i,
          nc=nargin,
          level=0;"base")
       char vflag=0; wbperm P;

       icFlags Q;
       for (i=0; i<(unsigned)nargin; ++i) {
          if (!Q.check_arg(argin[i])) { nc=i; break; }
       }

       if (nc<(unsigned)nargin) {
          OPTS opts; unsigned l=nc;

          if (nargin>(int)l && !mxIsChar(argin[l])) {
             try {
                P.init(FL,argin[l++],1);
             }
             catch (...) { wblog(FL,"ERR %s() expecting permutation\n"
               "(got %s at input arg #%d/%d)", FCT,
                mxGetClassName(argin[l]),l+1,nargin);
             }
          }
          opts.init(argin+l,nargin-l);

          vflag = opts.getOpt("-v");


          opts.checkAnyLeft();
       }

       if (isComplexQS(nc,argin)) {
          QSpace<gTQ,wbcomplex> C;
          contractQS_itags(nc,argin,level,vflag,C);
          C.Permute(P); a=C.toMx();
       }
       else {
          QSpace<gTQ,double> C;
          contractQS_itags(nc,argin,level,vflag,C);
          C.Permute(P); a=C.toMx();
       }
    }

    argout[0]=a;

#ifdef __WBDEBUG__
    Wb::MemCheck(FL,"stop");
#endif
};


mxArray* contractQS_plain(int nargin, const mxArray *argin[]) {

   char isra=1, isrb=1;
   unsigned i, l=4, k=2;

   ctrIdx ica, icb;
   mxArray *a=NULL;

   OPTS opts;
   wbperm P;

   if (nargin<int(l)) wblog(FL,
      "ERR invalid number of input arguments (%d/%d; usage #1)",nargin,l);
   ica.init(FL,argin[1]);
   icb.init(FL,argin[3]);

    if (nargin>(int)l && !mxIsChar(argin[l])) {
       P.init(FL,argin[l++],1);
    }

   opts.init(argin+l,nargin-l);

   if (opts.getOpt("conjA")) { ica.Conj(); }
   if (opts.getOpt("conjB")) { icb.Conj(); }


   opts.checkAnyLeft();

   for (i=0; i<2; ++i) { unsigned j=(i==1 ? k:i);
      try { mxIsQSpace(FL,argin[j],'c'); }
      catch (...) {
        wblog(FL,"ERR %s() invalid QSpace at arg #%d (%s)",
        FCT,i+1,mxGetClassName(argin[j]));
      }
   }

   isra=mxIsQSpace(argin[0]); if (!isra) mxIsQSpace(FL,argin[0],'c');
   isrb=mxIsQSpace(argin[k]); if (!isrb) mxIsQSpace(FL,argin[k],'c');

   if (isra) {
      QSpace<gTQ,double> A(argin[0],'r');

      if (isrb) {
         QSpace<gTQ,double> B(argin[k],'r'); QSpace<gTQ,double> C;
         CONTRACT_CG(A,ica,B,icb,C,P); a=C.toMx();
      }
      else {
         QSpace<gTQ,wbcomplex> B(argin[k]), C;
         CONTRACT_CG(A,ica,B,icb,C,P); a=C.toMx();
      }
   }
   else {
      QSpace<gTQ,wbcomplex> A(argin[0]);

      if (isrb) {
         QSpace<gTQ,double> B(argin[k],'r'); QSpace<gTQ,wbcomplex> C;
         CONTRACT_CG(A,ica,B,icb,C,P); a=C.toMx();
      }
      else {
         QSpace<gTQ,wbcomplex> B(argin[k]), C;
         CONTRACT_CG(A,ica,B,icb,C,P); a=C.toMx();
      }
   }

   return a;
};


template<class TA, class TB, class TC>
QSpace<gTQ,TC>& CONTRACT_CG(
    const QSpace<gTQ,TA> &A, const ctrIdx &ica,
    const QSpace<gTQ,TB> &B, const ctrIdx &icb, QSpace<gTQ,TC> &C,
    wbperm P
){
    if (A.isEmpty() || B.isEmpty()) { C.init(); return C; }

    unsigned i, ra=A.rank(FL), rb=B.rank(FL);

    for (i=0; i<ica.len; ++i) if (ica[i]>=ra) wblog(FL,
       "ERR index out of bounds (%d; %d)",ica[i],ra);
    for (i=0; i<icb.len; ++i) if (icb[i]>=rb) wblog(FL,
       "ERR index out of bounds (%d; %d)",icb[i],rb);

    A.contract(FL,ica,B,icb, C, P);

    C.SkipZeroData();

    return C;
};


template<class TD>
unsigned contractQS_itags(
   int nargin, const mxArray *argin[],
   unsigned level, unsigned vflag,
   QSpace<gTQ,TD> &C
){
   
   unsigned i, l=-1, len=0; char mark[nargin];
   icFlags q;


   for (i=0; i<(unsigned)nargin; ++i) { mark[i]=q.check_arg(argin[i]);
      if (mark[i]>1) ++len; else
      if (!mark[i] || (mark[i]==1 && !i)) wblog(FL,
         "ERR %s() invalid QSpace cell-structure (%d,%d)\n"
         "hint usage: ({QSpace} [,string opts])", FCT,level,i
      );
   }

   wbvector<icFlags> Q(len);
   wbvector<QSpace<gTQ,TD> > X(len);

   for (i=0; i<(unsigned)nargin; ++i) {
      if (mark[i]>1) {
         icFlags &q=Q[++l]; q.init(argin[i]);
         if (mark[i]==2) { X[l].init(FL,argin[i],'r'); }
         else {
            if (!mxIsCell(q.a)) wblog(FL,"ERR %s() got non-cell !??",FCT);
            unsigned j=0, nc=mxGetNumberOfElements(q.a);
               const mxArray* ac[nc];
               for (; j<nc; ++j) ac[j]=mxGetCell(q.a,j);
            contractQS_itags(nc,ac,level+1,vflag,X[l]);
         }
      }
      else if (mark[i]==1) { Q[l].set(argin[i]); }
      else if (mark[i]) {
         wblog(FL,"ERR %s() mark[%d]=%d !??",FCT,i,mark[i]);
      }
   }

   if (len<2) wblog(FL,"ERR %s() invalid usage #2: at least\n"
      "two QSpaces required for each contraction (%d)",FCT,len);

   unsigned k=len-1, ra,rb;
   ctrIdx ica,icb; wbperm P;

   QSpace<gTQ,TD> &B=X[k];

   Q[k].apply(B); if (Q[k].conj) icb.conj=1;

   for (--k; k<len; --k) {
      QSpace<gTQ,TD> &A=X[k];
      Q[k].apply(A); if (Q[k].conj) ica.conj=1;

      if (A.isEmpty() || B.isEmpty()) { 
         if (vflag) {
            char i[2]={ A.isEmpty(), B.isEmpty() };
            sprintf(str,"%s(#%d) %s empty", PROG, level,
            i[0] ? (i[1] ? "both QSpaces are": "QSpace 1 is")
                 : (i[1] ? "QSpace 2 is" : "NEITHER (!?) QSpace is"));
            wblog(FL,"%s",str);
         }
         C.init(); return nargin;
      }

      ra=A.rank(FL); rb=B.rank(FL);
      A.matchITags(FL,B,ica,icb);

      for (i=0; i<ica.len; ++i) if (ica[i]>=ra) wblog(FL,
         "ERR index out of bounds (%d; %d)",ica[i],ra);
      for (i=0; i<icb.len; ++i) if (icb[i]>=rb) wblog(FL,
         "ERR contract() index out of bounds (%d; %d)",icb[i],rb);

      if (vflag) wblog(FL,
         " *  %s(#%d)\n       %-15s -> %s\n    <> %-15s -> %s", PROG, level,
         IT2STR(A), ica.toStr().data,
         IT2STR(B), icb.toStr().data
      );

      try {
         A.contract(FL,ica,B,icb, C, P);
      }
      catch (...) { wblog(FL,"ERR %s(#%d) "
        "invalid cell-contraction (k=%d/%d)", PROG,level,k,len);
      }

      C.UnsetFlags();
      C.SkipZeroData();

      if (k) { C.save2(B); icb.conj=0; }
   }


   return nargin;
};


char isComplexQS(
   unsigned na, const mxArray *aa[],
   unsigned level
){
   icFlags Q;

   for (unsigned i=0; i<na; ++i) { const mxArray *a=aa[i];
      if (mxIsCell(a)) {
         unsigned i=0, nc=mxGetNumberOfElements(a);
            const mxArray* ac[nc];
            for (; i<nc; ++i) ac[i]=mxGetCell(a,i);
         if (isComplexQS(nc,ac,level+1)) return 1; 
      }
      else if (mxIsChar(a)) {
         if (Q.set(a,0)!=2) { wbstring s(a); wblog(FL,
            "ERR %s() invalid QSpace cell-structure ('%s'; l=%d)",
            FCT,s.data,level);
         }
      }
      else {
         try {
            if (!mxIsQSpace(FL,a) || mxGetNumberOfElements(a)!=1)
               throw "Wb: invalid QSpace";
            if (mxIsQSpace(FL,a,'c')) return 1;
         }
         catch (...) {
            wblog(FL,"ERR %s() invalid QSpace cell-structure (%s; l=%d)",
            FCT,mxGetClassName(a),level);
         }
      }
   }

   return 0;
};


