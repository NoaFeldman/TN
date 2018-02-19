#ifndef __WB_SPARRAY_CC__
#define __WB_SPARRAY_CC__

/*-------------------------------------------------------------------//
// A.sameSize(B [,strict])
// return > 0 if size is the same (up to trailing singletons)
       1: exactly the same size (use strict to use just this)
      11: same size, but only A is in DIAG_REP
      12: same size, but only B is in DIAG_REP
      13: same size, where both, A and B, are in DIAG_REP
      21: same size, but A.SIZE contains trailing singletons
      22: same size, but B.SIZE contains trailing singletons
      23: same size, but A or B contain intermediate singletons
// Wb,Dec20,11
//-------------------------------------------------------------------*/

template <class TD>
char wbsparray<TD>::sameSize(const wbsparray<TD> &B, char strict) const {

    if (SIZE.len!=B.SIZE.len && !strict) {
       if (  isDiag(FL)) return (B.isSMatrix(FL,  D.len) ? 11 : 0);
       if (B.isDiag(FL)) return (  isSMatrix(FL,B.D.len) ? 12 : 0);

       char s=sameSizeUp2Singletons(B.SIZE);
       if (s) {
          if (s>2)
             return 23;
          else {
             return (SIZE.len>B.SIZE.len ? 21 : 22);
          }
       }
       else return 0;
    }

    if (SIZE==B.SIZE) {
       if (!SIZE.len && D.len) 
            return (D.len==B.D.len ? 13 : 0);
       else return 1;
    }

    return 0;
};


template <class TD> inline
char wbsparray<TD>::sameSizeUp2Singletons(
  const wbvector<SPIDX_T> &S, const char strict
) const {
  
   if (SIZE==S) return 1;
   if (numel()!=S.prod(0)) return 0;

   if (isDiag()) {
      if (S.len<2 || S[0]!=S[1] || S[0]!=D.len) return 0;
      if (S.len==2) return 1;
      if (!strict) {
         for (unsigned i=2; i<S.len; ++i) if (S[i]!=1) return 0;
         return 2;
      }
      return 0;
   }

   {  char ok=1;
      unsigned n,N, i=0; const SPIDX_T *s;

      if (SIZE.len<S.len)
           { n=SIZE.len; N=   S.len; s=   S.data; }
      else { n=   S.len; N=SIZE.len; s=SIZE.data; }

      for (; i<n; ++i) { if (SIZE[i]!=S[i]) { ok=0; break; }}
      if (ok) {
      for (; i<N; ++i) if (s[i]!=1) { ok=0; break; }}

      if (ok) return 2;
   }

   if (!strict) {
      unsigned i=0, j=0, ra=SIZE.len, rb=S.len;
      const SPIDX_T *sa=SIZE.data, *sb=S.data;

      for (; i<ra; ++i) { if (sa[i]!=1) {
         for (; j<rb; ++j) { if (sb[j]!=1) break; }
         if (j<rb && sa[i]==sb[j]) { ++j; }
         else return 0;
      }}
      for (; j<rb; ++j) { if (sb[j]!=1) return 0; }
      return 3;
   }

   return 0;
};


template <class TD> inline
char wbsparray<TD>::matchWithSingletons(
   const char *F, int L, const wbvector<SPIDX_T> &S, wbvector<int> &M
 ) const {

   unsigned i=0, j=0, e=0, ra=SIZE.len, rb=S.len;
   const SPIDX_T *sa=SIZE.data, *sb=S.data;

   if (!S.len) {
      if (SIZE.len) wblog(FL,"ERR %s() got %s <> ()",FCT,sizeStr().data);
      M.init(); return 0;
   }

   M.init(S.len).set(-1);

   for (; i<ra; ++i) { if (sa[i]!=1) {
      for (; j<rb; ++j) { if (sb[j]!=1) break; }
      if (j<rb && sa[i]==sb[j]) { M[j++]=i; }
      else { ++e; break; }
   }}

   if (!e)
   for (; j<rb; ++j) { if (sb[j]!=1) { ++e; break; }}

   if (e) {
      if (F) wblog(F,L,"ERR %s() got size mismatch (%s <> %s)",
      FCT, SIZE.toStrD().data, S.toStrD().data);
      return e;
   }

   return 0;
};


template <class TD> inline
bool wbsparray<TD>::hasSize(SPIDX_T d1) const {
   return (SIZE.len==1 && SIZE[0]==d1);
};

template <class TD> inline
bool wbsparray<TD>::hasSize(SPIDX_T d1, SPIDX_T d2) const {
   return (SIZE.len==2 && SIZE[0]==d1 && SIZE[1]==d2);
};

template <class TD> inline
bool wbsparray<TD>::hasSize(SPIDX_T d1, SPIDX_T d2, SPIDX_T d3) const {
   return (SIZE.len==3 && SIZE[0]==d1 && SIZE[1]==d2 && SIZE[2]==d3);
};


template <class TD> inline
wbsparray<TD>& wbsparray<TD>::setRand(double x, char pnflag){

    static int firstcall=1;
    if (firstcall) {
       ::srand((unsigned)time((time_t*)NULL)); firstcall=0;
    }

    SPIDX_T n=numel(); if (!n) return *this;
    n=MIN(n,SPIDX_T(x*pow(double(n),1./rank(FL))));

    init_nnz(n); if (!n) return *this;

    SPIDX_T i,j, m=SIZE.len, *I=IDX.data;
    wbvector<double> Sfac(m);
    double fac=1./(double)RAND_MAX, dfac=fac, *sfac=Sfac.data;
    wbindex I1;

    for (j=0; j<m; ++j) { sfac[j]=fac*SIZE[j]; }

    if (!WbUtil<TD>::isFloat()) dfac*=100;

    for (i=0; i<n; ++i, I+=m) {
       for (j=0; j<m; ++j) {
          I[j]= SPIDX_T(sfac[j]*::rand());
          if (I[j]>=SIZE[j]) wblog(FL,"ERR %s() index out of bounds "
             "(%d,%d: %d; %s)",FCT,i+1,j+1,I[j],sizeStr().data
          );
       }
       D.data[i] = (TD)(dfac*::rand());
    }

    IDX.makeUnique1(I1); D.Select(I1); n=D.len;

    if (pnflag) {
       if (pnflag=='d' || pnflag=='D') {
          double x,a=2./(dfac*double(RAND_MAX));
          for (i=0; i<n; ++i) {
             x=a*D.data[i]-1;
             if (fabs(x)<1) { D.data[i]=atanh(x); }
             else wblog(FL,"ERR %s() got x=%.4g (%.4g) !??",FCT,x,D.data[i]);
          }
       }
       else {
          dfac *= (0.5*double(RAND_MAX));
          for (i=0; i<n; ++i) D.data[i] -= dfac;
       }
    }

    return Compress();
};


template <class TD> inline
SPIDX_T wbsparray<TD>::findRecSortedP(
   SPIDX_T *idx, SPIDX_T n, char lex
){
   unsigned i=0;

   if (long(n)<0) n=IDX.dim2;
   else if (n!=SIZE.len) {
      if (lex>0) {
         if (n<SIZE.len) {
            for (i=n; i<SIZE.len; ++i) { if (SIZE[i]!=1) { i=-1; break; }}
         }
         else if (n>SIZE.len) {
            for (i=SIZE.len; i<n; ++i) { if (idx[i]!=0) { i=-1; break; }}
            n=SIZE.len;
         }
      }
      else {
         wblog(FL,"ERR %s() got n=%d/%d (%d)",FCT,n,SIZE.len,lex);
      }
   }
      
   if (int(i)<0 || !n || !SIZE.len) wblog(FL,
      "ERR %s() invalid index [%s] (%s)",
      FCT, wbvector<SPIDX_T>(n,idx,'r'), sizeStr().data
   );

   return IDX.findRecSorted(idx,n,lex);
};


template<class TD> inline
char wbsparray<TD>::sameUptoFac(
  const wbsparray<TD> &B, TD *fac_, char lex, TD eps
) const {

   INDEX_T i=-1, n=D.len; TD a0,b0,fac=1;
   const TD *a=D.data, *b=B.D.data;
   wbindex Ia,Ib;
   char iA=0, iB=0, sab=sameSize(B);

   if (!sab) return 1;
   if (!n) { 
      if (fac_) { (*fac_)=(B.D.len ? 0. : 1.); }
      return (B.D.len ? 0 : 2);
   }

   if (sab>=20) wblog(FL,"ERR %s() got singletons (%s <> %s)",
      FCT,sizeStr().data, B.sizeStr().data); 

   if (sab>=10) {
      iA=((sab-10) & 1);
      iB=((sab-10) & 2);
   }

   a0=D.aMax(&i);
   
   if (sINDEX_T(i)<0) wblog(FL,
      "ERR %s() got zero data (%g;%d)",FCT,double(a0),i
   );

   if (iB) {
      if (iA) b0=B.D[i];
      else {
         const SPIDX_T *I=IDX.ref(i);
         if (IDX.dim2!=2) wblog(FL,
            "ERR %s() %dx%d !??",FCT,IDX.dim1,IDX.dim2);
         b0=(I[0]==I[1] ? B.D[I[0]] : TD(0.));
      }

   }
   else {
      b0=B.value(IDX.rec(i),IDX.dim2);
   }

   if (a0<=eps) {
      if (fac_) { (*fac_)=((a0==0 && b0) || ABS(b0)>eps ? 0. : 1.); }
      return (ABS(b0)>eps ? 0 : 3);
   }
   if (ABS(b0)<=eps) { return 4; }

   fac=a[i]/b0; if (fac_) (*fac_)=fac;

   if (sab>=10) {
      if (iA && iB) {
         for (i=0; i<D.len; ++i) { if (ABS(a[i]-fac*b[i])>eps) return 11; }
         return 0;
      }
      else if (iA) {
         const SPIDX_T *I=B.IDX.data; wbvector<char> mark(D.len);
         if (B.IDX.dim2!=2) wblog(FL,
            "ERR %s() %dx%d !??",FCT,B.IDX.dim1,B.IDX.dim2);
         for (i=0; i<B.D.len; ++i, I+=2) {
            if (I[0]==I[1]) { mark.el(I[0])=1;
                 if (ABS(a[I[0]]-fac*b[i])>eps) return 12; }
            else if (ABS(b[i])>eps) return 13;
         }
         for (i=0; i<D.len; ++i) {
            if (!mark[i] && ABS(a[i])>eps) return 14;
         }
         return 0;
      }
      else {
         const SPIDX_T *I=IDX.data; wbvector<char> mark(B.D.len);
         if (IDX.dim2!=2) wblog(FL,
            "ERR %s() %dx%d !??",FCT,IDX.dim1,IDX.dim2);
         for (i=0; i<D.len; ++i, I+=2) {
            if (I[0]==I[1]) { mark.el(I[0])=1;
                 if (ABS(a[i]-fac*b[I[0]])>eps) return 21; }
            else if (ABS(a[i])>eps) return 22;
         }
         for (i=0; i<B.D.len; ++i) {
            if (!mark[i] && ABS(b[i])>eps) return 23;
         }
         return 0;
      }
   }

   matchSortedIdxU(FL,IDX,B.IDX,Ia,Ib,-1,lex);
   for (n=Ia.len, i=0; i<n; ++i) {
      if (ABS(a[Ia[i]]-fac*b[Ib[i]])>eps) return 31;
   }

   Ia.Invert(  IDX.dim1,'u');
   for (i=0; i<Ia.len; ++i) { if (ABS(a[Ia[i]])>eps) return 32; }

   Ib.Invert(B.IDX.dim1,'u');
   for (i=0; i<Ib.len; ++i) { if (ABS(b[Ib[i]])>eps) return 33; }

   return 0;
};


template<class TD> inline
bool wbsparray<TD>::isZero(TD eps, char bflag) const {

   if (!D.len) return 1;
   Wb::scale_eps(eps,D.data,D.len);

   if (eps==0) {
      for (SPIDX_T i=0; i<D.len; ++i) if (D[i]!=0) return 0;
   }
   else {
      if (bflag) { return (D.norm()<eps); }
      else {
         for (SPIDX_T i=0; i<D.len; ++i) {
         if (ABS(D[i])>eps) return 0; }
      }
   }

   return 1;
};


template <class TD>
bool wbsparray<TD>::isDiagMatrix(TD eps) const {

   if (!isDiag()) { unsigned r=rank(FL); 

      if (r!=2) return 0;
      if (IDX.dim2==1) return 1;

      SPIDX_T i=0; const SPIDX_T *id=IDX.data;

      if (IDX.dim2!=2) sperror_this(FLF);
      Wb::scale_eps(eps,D.data,D.len);

      if (eps==0) {
         for (; i<IDX.dim1; ++i, id+=2) { if (id[0]!=id[1]) return 0; }
      }
      else {
         for (; i<IDX.dim1; ++i, id+=2) if (id[0]!=id[1]) {
            if (ABS(D[i])>eps) return 0;
         }
      }
   }

   return 1;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::getDiag(
   const char *F, int L, wbsparray<TD> &A
){
   if (isDiag(F_L)) { A=*this; }
   else {
      SPIDX_T i=0, n=-1; const SPIDX_T *I=IDX.data;

      if (!isSMatrix(F_L,&n) || IDX.dim2!=2) wblog(F_L,
         "ERR %s() requires rank-2 object (%s)",FCT,sizeStr().data);

      A.initDiag(n);

      for (; i<IDX.dim1; ++i, I+=2) if (I[0]==I[1]) {
      A.D.el(I[0]) += D[i]; }
   }

   return A;
};

template <class TD>
wbvector<TD>& wbsparray<TD>::getDiag(
   const char *F, int L, wbvector<TD> &dd
){
   SPIDX_T i=0, n=-1; const SPIDX_T *I=IDX.data;

   if (!isSMatrix(F_L,&n) || IDX.dim2!=2) wblog(F_L,
      "ERR %s() requires rank-2 object (%s)",FCT,sizeStr().data);

   dd.init(n);

   for (; i<IDX.dim1; ++i, I+=2) if (I[0]==I[1]) {
   dd.el(*I) += D[i]; }

   return dd;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::diag2reg(
  const char *F, int L, wbsparray<TD> &A, unsigned r
) const {

   if (int(r)<0 || r==2) {
      if (!D.len) { checkSize(F_L); A=(*this); return A; }
      if (!isDiag(FL)) sperror_this(F_LF);

      SPIDX_T i=0, n=D.len, *I; const SPIDX_T s[]={n,n};
      wbvector<SPIDX_T> S(2,s);

      A.init(S,n); A.D=D; I=A.IDX.data;
      for (; i<n; ++i, I+=2) { I[0]=I[1]=i; }
   }
   else {
      if (!r) wblog(FL,"ERR %s() invalid rank (%d)",FCT,r);
      if (!isScalar()) sperror_this(F_LF);
      A.SIZE.init(r).set(1);
      A.IDX.init(D.len,r); A.D=D; A.isref=0;
   }

   return A;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::diag2reg(
  const char *F, int L, unsigned r
){
   if (isref) wblog(FL,"ERR %s() got reference (%d)",FCT,isref);

   if (int(r)<0 || r==2) {
      if (!D.len) { checkSize(F_L); return *this; }
      if (!isDiag(FL)) sperror_this(F_LF);

      SPIDX_T i=0, n=D.len, *I;
      SIZE.init(2); SIZE[0]=SIZE[1]=n;
      IDX.init(n,SIZE.len);
      
      for (I=IDX.data; i<n; ++i, I+=2) { I[0]=I[1]=i; }
   }
   else {
      if (!r) wblog(FL,"ERR %s() invalid rank (%d)",FCT,r);
      if (!isScalar()) sperror_this(F_LF);
      SIZE.init(r).set(1);
      IDX.init(D.len,r);
   }

   return *this;
};


template <class TD>
bool wbsparray<TD>::isProptoId(TD &q, TD eps) const {

   if (!D.len) return 0;
   if (isDiag(FL)) {
      if (q==0) q=D[0];
      for (SPIDX_T i=1; i<D.len; ++i) { if (ABS(D[i]-q)>eps) return 0; }
      return 1;
   }

   unsigned j, r=rank(FL), r2=r/2;
   SPIDX_T i=0,d=1; const SPIDX_T *id=IDX.data, *s=SIZE.data;

   if (IDX.dim2!=r) sperror_this(FLF);
   Wb::scale_eps(eps,D.data,D.len);

   if (!r || r%2) return 0;
   for (j=0; j<r2; ++j) { d*=s[j];
      if (!s[j] || s[j]!=s[j+r2]) return 0;
   }

   if (double(eps)<=0) {
      if (D.len!=d) return 0;
      for (; i<IDX.dim1; ++i, id+=r) {
         for (j=0; j<r2; ++j) { if (id[j]!=id[j+r2]) return 0; }
         if (i || q!=0) { if (ABS(D[i]-q)>eps) return 0; }
         else q=D[i];
      }
   }
   else {
      SPIDX_T n=0;
      for (; i<IDX.dim1; ++i, id+=r) {
         for (j=0; j<r2; ++j) { if (id[j]!=id[j+r2]) {
            if (ABS(D[i])>eps) return 0; break;
         }}
         if (j==r2) {
            if ((++n)!=1 || q!=0) { if (ABS(D[i]-q)>eps) return 0; }
            else q=D[i];
         }
      }
      if (!n || n!=d) return 0;
   }

   return 1;
};


template <class TD>
char wbsparray<TD>::isIdentity(TD eps) const {

   if (!D.len) return 1;
   if (isDiag(FL)) {
      for (SPIDX_T i=1; i<D.len; ++i) { if (ABS(D[i]-1)>eps) return 2; }
      return 0;
   }

   unsigned r=rank(FL), r2=r/2;
   SPIDX_T j,i=0,d=1; const SPIDX_T *id=IDX.data, *s=SIZE.data;

   if (IDX.dim2!=r) sperror_this(FLF);
   Wb::scale_eps(eps,D.data,D.len);

   if (!r || r%2) return 3;
   for (j=0; j<r2; ++j) { d*=s[j];
      if (!s[j] || s[j]!=s[j+r2]) return 4;
   }

   if (eps<=0) {
      if (D.len!=d) return 5;
      for (; i<IDX.dim1; ++i, id+=r) {
         for (j=0; j<r2; ++j) { if (id[j]!=id[j+r2]) return 6; }
         if (ABS(D[i]-1)>eps) return 7;
      }
   }
   else {
      SPIDX_T n=0;
      for (; i<IDX.dim1; ++i, id+=r) {
         for (j=0; j<r2; ++j) { if (id[j]!=id[j+r2]) {
            if (ABS(D[i])>eps) return 8; break;
         }}
         if (j==r2) { ++n;
            if (ABS(D[i]-1)>eps) return 9;
         }
      }
      if (!n || n!=d) return 10;
   }

   return 0;
};


template <class TD>
char wbsparray<TD>::isIdentity(
   const wbperm &P, TD *dval_, TD eps) const {

   TD dval=1;

   if (!D.len) return 1;
   if (isDiag(FL)) {
      if (dval_) { (*dval_)=dval=D.avg(); eps*=dval; }
      for (SPIDX_T i=1; i<D.len; ++i) { if (ABS(D[i]-dval)>eps) return 2; }
      return 0;
   }

   unsigned j1,j2, r=rank(FL);
   const SPIDX_T *id=IDX.data;

   if (IDX.dim2!=r) sperror_this(FLF);

   if (P.len) {
      if (P.len!=r || r<2) wblog(FL,"ERR %s() "
         "invalid permutation P=[%s] (r=%d)",FCT,P.toStr().data,r);
      j1=P[0]; j2=P[1];
   }
   else {
      unsigned i=0; j1=j2=SIZE.len;
      for (   ; i<SIZE.len; ++i) { if (SIZE[i]>1) { j1=i; break; }}
      for (++i; i<SIZE.len; ++i) { if (SIZE[i]>1) { j2=i; break; }}
      for (++i; i<SIZE.len; ++i) { if (SIZE[i]>1) {       break; }}
      if (j1==SIZE.len) {
         if (SIZE.allEqual(1)) {
            if (D.len!=1) return 3;
            dval=D[0]; if (dval_) { (*dval_)=dval; }
            return 0;
         }
         else return 4;
      }
      if (j2==SIZE.len || i<SIZE.len) return 5;
   }

   if (j1>=SIZE.len || j2>=SIZE.len || j1==j2) wblog(FL,
      "ERR %s() got (%d,%d/%d)",FCT,j1,j2,SIZE.len);
   if (SIZE[j1]!=SIZE[j2]) return 6;

   Wb::scale_eps(eps,D.data,D.len);

   if (dval_) { 
      SPIDX_T i=0, l=0; dval=0;
      for (; i<IDX.dim1; ++i, id+=r) {
         if (id[j1]==id[j2]) { dval+=D[i]; ++l; }
      }
      if (l!=SIZE[j1]) return 7;
      if (l) { dval/=l; eps*=dval; }; (*dval_)=dval;
      id=IDX.data;
   }

   if (eps<=TD(0)) { if (D.len!=SIZE[j1]) return 8;
      for (SPIDX_T i=0; i<IDX.dim1; ++i, id+=r) {
         if (id[j1]!=id[j2])     return 12;"off-diag"
         if (ABS(D[i]-dval)>eps) return 11;"diag"
      }
   }
   else {
      SPIDX_T i=0, l=0;
      for (; i<IDX.dim1; ++i, id+=r) {
         if (id[j1]==id[j2]) { ++l;
                if (ABS(D[i]-dval)>eps) return 11;  }
         else { if (ABS(D[i])     >eps) return 12; }
      }
      if (!l || l!=SIZE[j1]) return 7;
   }

   return 0;
};


template <class TD>
bool wbsparray<TD>::isSym_aux(
  const char *F, int L, const char *fct,
  const wbsparray<TD> &B0, TD eps, TD* xref,
  const char symflag,
  const char lflag
) const {

   wbsparray<TD> B; B0.transpose(F_L,B);
   SPIDX_T i,r,s; INDEX_T m,ma,mb;
   wbindex Ia,Ib;
   wbvector<SPIDX_T> Ja(IDX.dim1), Jb(B.IDX.dim1);
   TD x;

   checkSize(FL,"A: "); if (&B!=this) checkSize(FL,"B: ");

   if (IDX.dim2!=B.IDX.dim2 || IDX.dim2%2) {
      if (lflag) sprintf(str,
         "%s %s() only applies to even-rank objects (%ld;%ld).",
          shortFL(F,L),fct, SIZE.len, B.SIZE.len);
      return 0;
   }
   if (isEmpty() && B.isEmpty()) return 1;

   m=matchIndex(IDX,B.IDX,Ia,Ib,1,&ma,&mb);
   if (m) wblog(FL,"ERR %s() input sparse matrices "
      "not compressed (%d,%d,%d) !??",FCT,m,ma,mb);

   if (xref) (*xref)=0;

   if (symflag=='s') {
      for (i=0; i<Ia.len; ++i) { r=Ia[i]; s=Ib[i];
         if ((++Ja[r])>1 || (++Jb[s])>1) {
            wblog(FL,"ERR %s()",FCT);
         }
         x=ABS(D[r]-CONJ(B.D[s])); if (xref) *xref=MAX(*xref,x);
         if (x>eps) return 0;
      }
   }
   else if (symflag=='a') {
      for (i=0; i<Ia.len; ++i) { r=Ia[i]; s=Ib[i];
         if ((++Ja[r])>1 || (++Jb[s])>1) {
            wblog(FL,"ERR %s()",FCT);
         }
         x=ABS(D[r]+CONJ(B.D[s])); if (xref) *xref=MAX(*xref,x);
         if (x>eps) return 0;
      }
   }
   else wblog(FL,"ERR invalid symflag=%c<%d>",symflag,symflag);

   for (i=0; i<Ja.len; ++i) { if (!Ja[i] && ABS(  D[i])>eps) return 0; }

   if (&B0!=this)
   for (i=0; i<Jb.len; ++i) { if (!Jb[i] && ABS(B.D[i])>eps) return 0; }

   return 1;
};


template <class TD>
TD wbsparray<TD>::dotProd(
  const char *F, int L, const wbsparray<TD> &B, char tnorm
) const {

   char cflag = ((tnorm || WbUtil<TD>::hasNoConj()) ? 0 : 1);
   char iA=0, iB=0, sab=sameSize(B);
   TD x=TD();

   if (!sab || sab>=20) wblog(F_L,"ERR %s() "
      "size mismatch (%s <> %s)",FCT,sizeStr().data,B.sizeStr().data);

   if (sab>=10) {
      iA=((sab-10) & 1);
      iB=((sab-10) & 2);

      if (iA && iB) {
         const TD *a=D.data, *b=B.D.data; SPIDX_T i=0;
         if (cflag)
              { for (; i<D.len; ++i) x+=CONJ(a[i])*b[i]; }
         else { for (; i<D.len; ++i) x+=     a[i] *b[i]; }
      }
      else if (iA)
           { x=B.dotProd_full_diag(*this, cflag ? 2:0); }
      else { x=  dotProd_full_diag(B,     cflag ? 1:0); }
   }
   else {
      wbindex Ia,Ib; SPIDX_T i=0; INDEX_T m,ma=0,mb=0;
      const TD *a=D.data, *b=B.D.data;

      m=matchIndex(IDX,B.IDX,Ia,Ib,1,&ma,&mb);
      if (m) wblog(F_L,"ERR %s() input sparse "
        "matrices not compressed (%d,%d,%d) !??",FCT,m,ma,mb);

      if (cflag)
           { for (; i<Ia.len; ++i) x+=CONJ(a[Ia[i]])*b[Ib[i]]; }
      else { for (; i<Ia.len; ++i) x+=     a[Ia[i]] *b[Ib[i]]; }
   }

   return x;
};


template <class TD> inline
TD wbsparray<TD>::dotProd_full_diag(const wbsparray<TD> &B, char cflag
) const {

   const TD *a=D.data, *b=B.D.data;
   const SPIDX_T *I=IDX.data;
   TD x=TD();

   if (IDX.dim2!=2 || B.IDX.dim1 || B.D.len!=SIZE[0]) wblog(FL,
      "ERR %s() %s <> %s !??",FCT, sizeStr().data, B.sizeStr().data
   );

   if (cflag==0) {
      for (SPIDX_T i=0; i<D.len; ++i, I+=2) {
         if (I[0]==I[1]) { x += a[i] * b[I[0]]; }
      }
   }
   else if (cflag==1) {
      for (SPIDX_T i=0; i<D.len; ++i, I+=2) {
         if (I[0]==I[1]) { x += CONJ(a[i]) * b[I[0]]; }
      }
   }
   else if (cflag==2) {
      for (SPIDX_T i=0; i<D.len; ++i, I+=2) {
         if (I[0]==I[1]) { x += a[i] * CONJ(b[I[0]]); }
      }
   }
   else wblog(FL,"ERR %s() invalid cflag=%d",FCT,cflag);

   return x;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::TimesEl(
   const wbsparray<TD> &B, char tnorm
){
   char cflag = ((tnorm || WbUtil<TD>::hasNoConj()) ? 0 : 1);
   char iA=0, iB=0, sab=sameSize(B);

   if (!sab || sab>=20) wblog(FL,"ERR %s() "
      "size mismatch (%s <> %s)",FCT,sizeStr().data,B.sizeStr().data);

   if (sab>=10) { wbsparray<TD> X;

      iA=((sab-10) & 1);
      iB=((sab-10) & 2);

      if (iA && iB) {
         SPIDX_T i=0; const TD *a=D.data; TD *x;
         X=B; x=X.D.data;
         if (cflag)
              { for (; i<D.len; ++i) x[i]*=CONJ(a[i]); }
         else { for (; i<D.len; ++i) x[i]*=     a[i] ; }
      }
      else if (iA)
           { B.timesEl_full_diag(*this,X, cflag ? 2:0); }
      else {   timesEl_full_diag(B,    X, cflag ? 1:0); }
   }
   else {
      wbindex Ia,Ib; SPIDX_T i=0, m,ma=0,mb=0;
      const TD *a=D.data, *b=B.D.data;

      m=matchIndex(IDX,B.IDX,Ia,Ib,1,&ma,&mb);
      if (m) wblog(FL,"ERR %s() input sparse "
        "objects not compressed (%d,%d,%d) !??",FCT,m,ma,mb);

      wbsparray<TD> X(SIZE,Ia.len);
      TD *x=X.D.data;

      for (i=0; i<Ia.len; ++i) { X.IDX.recSetP(i, IDX.rec(Ia[i])); }
      if (cflag)
           { for (i=0; i<Ia.len; ++i) { x[i]=CONJ(a[Ia[i]])*b[Ib[i]]; }}
      else { for (i=0; i<Ia.len; ++i) { x[i]=     a[Ia[i]] *b[Ib[i]]; }}

      X.save2(*this).Compress();
   }

   return *this;
};


template <class TD> inline
wbsparray<TD>& wbsparray<TD>::timesEl_full_diag(
  const wbsparray<TD> &B, wbsparray<TD> &X, char cflag
) const {

   const TD *a=D.data, *b=B.D.data;
   const SPIDX_T *I=IDX.data;

   X.initDiag(B.D.len); TD *x=X.D.data;

   if (IDX.dim2!=2 || B.IDX.dim1 || B.D.len!=IDX.dim1) wblog(FL,
      "ERR %s() (%d,%d)x%d <> (%d,%d)x%d !??",FCT,
      D.len, IDX.dim1,IDX.dim2, B.D.len, B.IDX.dim1,B.IDX.dim2
   );

   if (cflag==0) {
      for (SPIDX_T i=0; i<IDX.dim1; ++i, I+=2) { if (I[0]==I[1]) {
         x[I[0]]=a[i]*b[I[0]];
      }}
   }
   else if (cflag==1) {
      for (SPIDX_T i=0; i<IDX.dim1; ++i, I+=2) { if (I[0]==I[1]) {
         x[I[0]]=CONJ(a[i])*b[I[0]];
      }}
   }
   else if (cflag==2) {
      for (SPIDX_T i=0; i<IDX.dim1; ++i, I+=2) { if (I[0]==I[1]) {
         x[I[0]]=a[i]*CONJ(b[I[0]]);
      }}
   }
   else wblog(FL,"ERR %s() invalid cflag=%d",FCT,cflag);
   return X;
};


template <class TD> inline
wbsparray<TD>& wbsparray<TD>::setCol(
   const char *F, int L, SPIDX_T k, const wbsparray<TD> &v
){
   unsigned rv=v.rank(F_L);

   if (rank(F_L)!=2 || rv<1 || rv>2 || SIZE[0]!=v.SIZE[0] ||
      (rv==2 && v.SIZE[1]!=1)
    ) wblog(F_L,"ERR %s() invalid input (%s <> %s)",
      FCT, sizeStr().data, v.sizeStr().data
   );

   if (k>=SIZE[1]) wblog(F_L,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,SIZE[1]);
   if ((IDX.dim1 && IDX.dim2!=2) || v.IDX.dim2!=rv) wblog(F_L,
      "ERR %s() %dx%d <> %dx%d/%d !??",FCT,IDX.dim1,IDX.dim2,
      v.IDX.dim2,v.IDX.dim2,rv);

   SPIDX_T i=0, l=0, l2, n=IDX.dim1, *I2=IDX.data+1;

   for (; i<n; ++i, I2+=2) { if ((*I2)!=k) {
       if (l<i) { setRec(l,i); }; ++l;
   }}

   if (l<i && F) wblog(F_L,
      "WRN %s() *this already contains data at col=%d/%d (%d/%d)", 
      FCT,k,SIZE[1],l,SIZE[0]
   );
   if (!l) {
      if (v.IDX.dim2==1)
           { v.IDX.resize(v.IDX.dim1, 2, IDX); }
      else { IDX=v.IDX; }
      IDX.setLastCol(k); 

      D=v.D; return *this;
   }

   n=v.D.len; l2=l+n;
   if (l2) {
      if (l2<=IDX.dim1)
           { IDX.dim1=l2;             D.len=l2; }
      else { IDX.Resize(l2,IDX.dim2); D.Resize(l2); }
   }
   else { IDX.init(0,IDX.dim2); D.init(); }

   if (n) {
      const SPIDX_T *Iv=v.IDX.data;
      cpyRange(v.D.data, D.data+l, n); I2=IDX.rec(l);
      for (i=0; i<n; ++i, I2+=2, Iv+=rv) {
         I2[0]=(*Iv);"rank-2")
         I2[1]=k;
      }
   }
   return *this;
};


template <class TD> inline
TD wbsparray<TD>::Normalize(char tnorm, char qflag){
   TD x = SQRT(overlap(D.data,D.data,D.len,1,tnorm));

   if (ABS(x)>TD(1E-15)) timesRange(D.data,TD(1)/x,D.len);
   else if (!qflag) {
      wblog(FL,"ERR %s() got vector with norm %.3g !??",FCT,double(x));
   }

   return x;
};

template <class TD>
TD wbsparray<TD>::NormalizeCol(SPIDX_T k, char tnorm, char qflag){

   if (rank()!=2 || IDX.dim2!=2) wblog(FL,
      "ERR %s() requires rank-2 (got %s)",FCT,sizeStr().data);
   if (k>=SIZE[1]) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,SIZE[1]);
   checkSize(FLF);

   SPIDX_T i=0, n=IDX.dim1, *I2=IDX.data+1;
   opFlags<TD> xflags(tnorm ? 'N' : 'C');
   TD x=0;

   if (xflags.conj())
        { for (; i<n; ++i, I2+=2) if ((*I2)==k) { x+=CONJ(D[i])*D[i]; }}
   else { for (; i<n; ++i, I2+=2) if ((*I2)==k) { x+=D[i]*D[i]; }}

   x=SQRT(x);
   if (ABS(x)>1E-15) { x=1/x;
      for (i=0; i<n; ++i, I2+=2) if ((*I2)==k) { D[i]*=x; }
   }
   else if (!qflag) {
      wblog(FL,"ERR %s() got vector with norm %.3g !??",FCT,x);
   }

   return x;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::OrthoNormalizeCols(const char *F, int L,
   char tnorm,
   char qflag,
   TD eps,
   unsigned np
){
   if (rank()!=2) wblog(FL,
      "ERR %s() requires rank-2 array (%s)",FCT,sizeStr().data);
  
   SPIDX_T i,j,ip;
   char xflag=(qflag=='x' || qflag=='X');
   wbvector< wbsparray<TD> > X;
   wbvector<SPIDX_T> I0;
   TD a,z;

   splitSparseCM(F_L,X, xflag ? NULL : &I0);

   WbUtil<TD>().adjust_tnorm(tnorm);

   for (i=0; i<X.len; ++i) {
      for (ip=0; ip<np; ++ip) {
         for (j=0; j<i; ++j) { if (X[j].D.len) {
            z=X[i].dotProd(FL,X[j],tnorm); if (z!=0) {
            X[i].Plus(FL,X[j],-z); }
         }}
      }

      z=X[i].norm2(tnorm); a=ABS(z);

      if (a>eps) { X[i]*=(TD(1)/SQRT(a)); }
      else {
         if (qflag) {
            X[i].initz();
         }
         else if (!xflag) {
            MXPut(FL,"x").add(*this,"V").add(X,"X")
              .add(qflag,"qflag").add(xflag,"xflag")
              .add(i+1,"i").add(a,"a").add(eps,"eps");
            wblog(F_L,"ERR |U(:,%d/%d)| = %.3g (%g) %s !??", i+1, X.len,
               double(a), double(eps), tnorm ? " (using t-norm)" : ""
            );
         }
      }
   }

   wbvector< wbsparray<TD>* > xp(X.len);
   if (xflag) {
      for (j=i=0; i<X.len; ++i) { if (X[i].D.len) {
         xp[j++]=(&X[i]);
      }}
      if (j<i) {
         if (!j) wblog(F_L,"ERR %s() got all null-vectors",FCT);
         xp.len=j;
      }
   }
   else {
      for (i=0; i<X.len; ++i) { xp[i]=(&X[i]); }
   }

   initCAT(FL,xp,&I0);
   if (!IDX.dim1) wblog(F_L,"ERR %s() got all null-vectors",FCT);

   return *this;
};


template <class TD>
template <class T2>
wbvector<T2>& wbsparray<TD>::norm2vec(unsigned k,
   wbvector<T2> &A, char tnorm
 ) const {

   if (k>=SIZE.len) wblog(FL,
      "ERR %s() index out of bounds (k=%d/%d)",FCT,k,SIZE.len);
   A.init(SIZE[k]);

   unsigned m=IDX.dim2;
   SPIDX_T i=0, *ip=IDX.data+k;
   T2 *a=A.data;

   if (tnorm) {
      for (; i<D.len; ++i, ip+=m) { a[*ip] += D[i]*D[i]; }
   }
   else {
      for (; i<D.len; ++i, ip+=m) { a[*ip] += CONJ(D[i])*D[i]; }
   }

   return A;
};


template <class TD>
TD wbsparray<TD>::trace() const {

   TD x=0;
   
   unsigned r=rank(FL); if (!r) return x;

   if (r==2 && !IDX.dim1) return D.sum();
   if (r%2 || r!=SIZE.len || r!=IDX.dim2) sperror_this(FLF);

   SPIDX_T i=0; unsigned r2=r/2, n=r2*sizeof(SPIDX_T);
   const SPIDX_T *idx=IDX.data;

   for (; i<r2; ++i) {
      if (SIZE[i]!=SIZE[i+r2]) wblog(FL,
      "ERR %s() got non-symmetric tensor (%s)",FCT,sizeStr().data);
   }

   for (i=0; i<IDX.dim1; ++i, idx+=r) {
      if (memcmp(idx,idx+r2,n)==0) x+=D[i];
   }

   return x;
};


template <class TD>
template <class T2>
wbvector<T2>& wbsparray<TD>::trace(
   unsigned k, wbvector<T2> &A
 ) const {

   SPIDX_T i=0;
   unsigned r=IDX.dim2, r2=(r-1)/2, n=r2*sizeof(SPIDX_T);
   const SPIDX_T *idx=IDX.data+(k?0:1), *sz=SIZE.data+(k?0:1);
   int l=(k ? r-1 : -1);

   if (k && k+1!=SIZE.len) wblog(FL,
      "ERR %s() only accepts k=1 or k=rank (%d/%d)",FCT,k+1,r);
   if (SIZE.len%2!=1) wblog(FL,"ERR %s() "
      "requires odd-rank tensor (r=%d; k=%d)",FCT,SIZE.len,k);
   if (SIZE.len!=IDX.dim2 || IDX.dim1!=D.len) sperror_this(FLF);

   A.init(SIZE[k]);
   T2 *a=A.data;

   for (; i<r2; ++i) {
      if (sz[i]!=sz[i+r2]) wblog(FL,"ERR %s() "
      "got non-symmetric tensor (%s; k=%d)",FCT,sizeStr().data,k+1);
   }

   for (i=0; i<IDX.dim1; ++i, idx+=r) {
      if (memcmp(idx,idx+r2,n)==0) {
         if (idx[l]>=A.len) wblog(FL,"ERR %s() "
            "index out of bounds (%d,%d: %d/%d)",FCT,i,l,idx[l],A.len);
         a[idx[l]]+=D[i];
      }
   }

   return A;
};


bool mxIsWbsparray(
   const char *F, int L, const mxArray *a, unsigned k
){

   if (!a || mxIsEmpty(a)) { return 0; }

   SPIDX_T n=mxGetNumberOfElements(a);
   if (k>=n) { if (F) wblog(F,L,
      "ERR sparse::%s() index exceeds dimension (%d/%d)",FCT,k,n);
      return 0;
   }

   int e[3] { 
      mxGetFieldNumber(a,"S"),
      mxGetFieldNumber(a,"idx"), mxGetFieldNumber(a,"data")
   };

   if (e[0]<0 || e[1]<0 || e[2]<0) { if (F) wblog(FL,
      "ERR %s() invalid sparse array (missing fields%s%s%s)",FCT,
       e[0]?" S":"", e[1]?" idx":"", e[2]?" data":"");
      return 0;
   }

   return 1;
};


bool mxIsWbsparray(
   const char *F, int L, const mxArray *a, unsigned k,
   const mxArray **as_, const mxArray **ad_, const mxArray **ai_
){

   if (!a || mxIsEmpty(a)) { return 0; }

   SPIDX_T n=mxGetNumberOfElements(a);
   if (k>=n) { if (F) wblog(F,L,
      "ERR sparse::%s() index exceeds dimension (%d/%d)",FCT,k,n);
      return 0;
   }

   const mxArray 
      *as=mxGetField(a,k,"S"),
      *ai=mxGetField(a,k,"idx"),
      *ad=mxGetField(a,k,"data");

   if (!as || !ad) { if (F) {
      if (as && mxGetNumberOfElements(as)) wblog(F_L,
         "ERR %s() invalid sparse array (0x%lX, 0x%lX: S (%d)",
          FCT,as,ad, as ? mxGetNumberOfElements(as):-1);
      if (ad && mxGetNumberOfElements(ad)) wblog(F_L,
         "ERR %s() invalid sparse array (0x%lX, 0x%lX: data (%d)",
          FCT,as,ad, ad ? mxGetNumberOfElements(ad):-1);
      }
      return 0;
   }

   if (as_) (*as_)=as;
   if (ad_) (*ad_)=ad;
   if (ai_) (*ai_)=ai;

   return 1;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::toScalar(const char *F, int L, unsigned r) {

   if (D.len>1 || (SIZE.len && !SIZE.allEqual(1))) wblog(F_L,
      "ERR %s() invalid scalar %s",FCT,info2Str().data);

   if (!D.len) { D.init(1); D[0]=0; }

   if (int(r)<=0) {
      if (SIZE.data) SIZE.init();
      if (IDX.dim1 || IDX.dim2) IDX.init();
   }
   else {
      if (SIZE.len) { if (SIZE.len!=r) wblog(FL,
         "ERR %s() got rank mismatch (%d/%d)",FCT,SIZE.len,r); }
      else {
         SIZE.init(r).set(1);
         IDX.init(1,r);
      }
   }

   return *this;
};


template <class TD>
TD wbsparray<TD>::getScalar(const char *F, int L) const {

   if (D.len>1 || (SIZE.len && !SIZE.allEqual(1))) wblog(F_L,
      "ERR %s() invalid scalar %s",FCT,info2Str().data);

   return (D.len ? D[0] : TD(0));
};


template <class TD>
int wbsparray<TD>::checkTrailingSingletons(
   const char *F, int L, unsigned *r
 ) const {

   if (SIZE.len!=IDX.dim2) wblog(F_L,
      "ERR %s() size mismatch (%d/%d)",FCT,SIZE.len,IDX.dim2);
   if (!SIZE.len) {
      if (r && int(*r)>=0) wblog(FL,
         "ERR %s() out of bounds (%d/%d)",FCT,*r,SIZE.len);
      return 0;
   }

   unsigned l=SIZE.len-1;
   for (; l<SIZE.len; --l) { if (SIZE[l]!=1) break; }; ++l;

   if (r) {
      if (int(*r)>=0) {
         if ((*r)>SIZE.len) wblog(FL,
            "ERR %s() out of bounds (%d/%d)",FCT,*r,SIZE.len);
         if (l>(*r)) wblog(FL,"ERR %s() "
            "got non-singletons for dim>=%d (%d,%d)",FCT,*r,l,SIZE.len);
         l=(*r);
      }
      else (*r)=l;
   }

   return (SIZE.len-l);
};


template <class TD>
int wbsparray<TD>::skipTrailingSingletons(
   const char *F, int L, wbsparray<TD> &C, unsigned r
 ) const {

   C.init();

   if (SIZE.len && checkTrailingSingletons(F_L,&r)) {
      SPIDX_T i=0, n=IDX.dim1, zero=0;
      unsigned m=SIZE.len-r;
         
      #ifndef NSAFEGUARDS
         for (; i<n; ++i) if (Wb::anyUnequal(IDX.ref(i,r),m,zero)) {
            wblog(F_L,"ERR %s() IDX out of bounds (%d/%dx%d: %d/%d) !?",
            FCT,i+1,IDX.dim1,IDX.dim2,r,SIZE.len); 
         }
      #endif

      C.SIZE.init(r,SIZE.data);
      IDX.resize(IDX.dim1,r, C.IDX);
      C.D=D;

      return m;
   }

   if (int(r)>=0) wblog(FL,
      "ERR %s() out of bounds (%d/%d)",FCT,r,SIZE.len);

   return 0;
};


template <class TD>
int wbsparray<TD>::skipTrailingSingletons(
   const char *F, int L, unsigned r
){

   if (SIZE.len && checkTrailingSingletons(F_L,&r)) {
      SPIDX_T i=0, n=IDX.dim1, zero=0; unsigned m=SIZE.len-r;
         
      #ifndef NSAFEGUARDS
         for (; i<n; ++i) if (Wb::anyUnequal(IDX.ref(i,r),m,zero)) {
            wblog(F_L,"ERR %s() IDX out of bounds (%d/%dx%d: %d/%d) !?",
            FCT,i+1,IDX.dim1,IDX.dim2,r,SIZE.len); 
         }
      #endif

      if (r) SIZE.len=r; else SIZE.init();
      IDX.Resize(IDX.dim1,r);

      return m;
   }

   if (int(r)>=0) wblog(FL,
      "ERR %s() out of bounds (%d/%d)",FCT,r,SIZE.len);

   return 0;
};


template <class TD> inline
bool wbsparray<TD>::hasSingletons(const wbindex &I) const {

   if (!I.len) wblog(FL,"ERR %s() got emtpy index set",FCT);
   for (unsigned i=0; i<I.len; ++i) {
      if (I[i]>=SIZE.len || SIZE[I[i]]!=1) return 0;
   }


   return 1;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::skipSingletons(
   const char *F, int L, const wbindex &I, wbsparray<TD> &A
 ) const {

   if (!I.len) { A=*this; return A; }
   if (&A==this) {
      wbsparray<TD> X; A.save2(X);
      return X.skipSingletons(F,L,I,A);
   }

   if (!hasSingletons(I)) wblog(FL,"ERR %s() got S=(%s) "
      "given I=[%s] !??", FCT,sizeStr().data, I.toStr().data);

   wbindex I2; I.invert(IDX.dim2,I2,'u');
   
   SIZE.select(I2,A.SIZE); A.isref=0;
   IDX.getCols(I2,A.IDX); A.D=D;

   return A;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::addTrailingSingletons(
   const char *F, int L, unsigned r
){
   unsigned rk=rank(FL);
   if (int(r)<0 || r<rk) wblog(F_L,
      "ERR %s() invalid rank (%d/%d)",FCT,r,rk);

   if (r==rk) return *this;
   if (r<rk) wblog(F_L,"ERR %s() cannot reduce rank (%d/%d)",FCT,r,rk);

   if (isDiag()) {
      if (IDX.data) wblog(F_L,
         "ERR %s() got IDX data (%s) !??",FCT,IDX.dim1,IDX.dim2);
      SIZE.init(r); SIZE[0]=SIZE[1]=D.len;
         for (unsigned i=2; i<r; ++i) SIZE[i]=1;
      IDX.init(D.len,r);
         for (SPIDX_T *I=IDX.data, i=0; i<D.len; ++i, I+=r) {
         I[0]=I[1]=i;
      }
   }
   else {
      if (SIZE.len) {
         unsigned j=SIZE.len; SIZE.Resize(r);
         for (; j<r; ++j) SIZE[j]=1;
      }

      if (IDX.dim1) {
         if (!SIZE.len) wblog(F_L,
            "ERR %s() got empty size (%d/%d)",FCT,SIZE.len,r);
         IDX.Resize(IDX.dim1,r);
      }
   }

   return *this;
};


template <class TD>
char wbsparray<TD>::checkSize2(
  const char *F, int L, const wbsparray<TD> &B
) const {

    char sab=sameSize(B);
    this->checkSize(F_L,"A:"); B.checkSize(F_L,"B:");

    if (sab==1) {
       if (IDX.dim2!=B.IDX.dim2){ wblog(F_L,
             "ERR sparse::%s() IDX size mismatch (%dx%d <> %dx%d)",
             FCT,IDX.dim1,IDX.dim2,B.IDX.dim1,B.IDX.dim2);
          return -1;
       }
    }
    else if (sab<=0) {
       wblog(F_L,"ERR %s() mismatch %s <> %s (%d)",
          FCT,sizeStr().data,B.sizeStr().data,sab);
    }

    return sab;
};


template <class TD> inline
char wbsparray<TD>::checkSize(
  const char *F, int L, const char *fct, const char *istr
) const {

   if (!SIZE.len) {
      if (IDX.dim1 || IDX.dim2) {
         if (F) sperror_this(F_L_F,istr); else return 1;
      }
   }
   else {
      if (IDX.dim1!=D.len || 
         (IDX.dim2!=SIZE.len && (IDX.dim1 || IDX.dim2))
      ){
         if (F) sperror_this(F_L_F,istr); else return 2;
      }
   }
   return 0;
};

template <class TD> inline
void wbsparray<TD>::check_IDX_range(const char *F, int L) const {

   if (!SIZE.len) {
      if (IDX.dim1 || IDX.dim2) sperror_this(F_LF);
      return;
   }

   if (IDX.dim1!=D.len || 
      (IDX.dim2!=SIZE.len && (IDX.dim1 || IDX.dim2))
   ) sperror_this(F_LF);

   if (IDX.dim1) {
      SPIDX_T i=0; unsigned j=0, m=IDX.dim2;
      const SPIDX_T *s=SIZE.data, *I=IDX.data;

      for (; i<IDX.dim1; ++i, I+=m) {
         for (j=0; j<m; ++j) {
            if (I[j]>=s[j]) wblog(FL,
               "ERR %s() index out of bounds (%d: %d <> %s)",
               FCT,j+1,I[j]+1, SIZE.toStrD().data
            );
         }
      }
   }
};


template <class TD>
mxArray* wbsparray<TD>::toMxSp() const {

   if (!isBaseType(typeid(TD))) {
      double x;
      for (SPIDX_T i=0; i<D.len; ++i) { x=double(D[i]);
         if (fabs(double(TD(x)-D[i])/x)>1E-20) {
             MXPut(FL,"ans").add(*this,"data").add(i+1,"i").add(x,"x");
             wblog(FL,"ERR %s() got data type `%s' (d[%d]=%g @ %.3g)",
                FCT,getName(typeid(TD)).data, i+1,x,double(TD(x)-D[i])
             );
         }
      }
   }

   if (isScalar()) {
      if (D.len) {
          if (D.len>1) sperror_this(FLF);
          return numtoMx(D[0]);
      }
      return numtoMx(0);
   }

   if (isDiag(FL)) {
      wbMatrix<SPIDX_T> IJ; SPIDX_T d1=1,d2=1;
      get2DIndex(IJ,d1,d2);

      return Wb::mxCreateSparse(FL,d1,d2,IJ, D);
   }

   if (!SIZE.len) {
      if (IDX.isEmpty() && D.isEmpty())
         return Wb::mxCreateSparse(FL,0,0,IDX,D);
      else wblog(FL,
     "ERR %s() invalid wbsparray (%s;%d)",FCT,sizeStr().data,D.len);
   }

   if (SIZE.len==2) {
      return Wb::mxCreateSparse(FL,SIZE[0],SIZE[1],IDX,D);
   }

   if (SIZE.len==1) {
      wbMatrix<SPIDX_T> IX(IDX.dim1,2);
      SPIDX_T i=0, *I=IX.data;

      if (IDX.dim2!=1) sperror_this(FLF);
      for (; i<IDX.dim1; ++i, I+=2) { I[0]=0; I[1]=IDX.data[i]; }

      return Wb::mxCreateSparse(FL,1,SIZE[0],IX,D);
   }

   wblog(FL,"ERR %s() requires rank-2 wbsparray (%s)",
   FCT,sizeStr().data); return 0;
};


template <class TD>
mxArray* wbsparray<TD>::mxCreateStructX(
   unsigned m, unsigned n, int ma
) const {
   const char *fields[]={"S","idx","data","info","A"};
   return mxCreateStructMatrix(m,n,ma!=-99 ? 5:4,fields);
};

template <class TD>
void wbsparray<TD>::add2MxStructX(mxArray *a, unsigned i, int ma) const {

   char isId=isIdentityMatrix(), isc=isScalar();

   mxSetFieldByNumber(a,i,0, SIZE.toMx());

   if (isc) {
      mxSetFieldByNumber(a,i,1, IDX.toMx());
      mxSetFieldByNumber(a,i,2, D.toMx('t'));
      mxSetFieldByNumber(a,i,3, wbstring("scalar").toMx());
   }
   else if (!IDX.isEmpty() && IDX.allEqual(0) && D.allEqual(0)) {
      mxSetFieldByNumber(a,i,1,mxCreateSparse(IDX.dim1,IDX.dim2,0,mxREAL));
      mxSetFieldByNumber(a,i,2,mxCreateSparse(1,D.len,0,mxREAL));
      mxSetFieldByNumber(a,i,3, wbstring("(init)").toMx());
   }
   else if (isId) {
      mxSetFieldByNumber(a,i,2, numtoMx(D[0]));
      mxSetFieldByNumber(a,i,3, wbstring("identity").toMx());
   }
   else if (isDiag(FL)) {
      mxSetFieldByNumber(a,i,2, D.toMx('t'));
      mxSetFieldByNumber(a,i,3, wbstring("diag").toMx());
   }
   else if (IDX.dim2!=1 && isDiagMatrix()) {
      mxSetFieldByNumber(a,i,1, IDX.getCols(0,0).toMx());
      mxSetFieldByNumber(a,i,2, D.toMx('t'));
   }
   else {
      mxSetFieldByNumber(a,i,1, IDX.toMx());
      mxSetFieldByNumber(a,i,2, D.toMx('t'));
   };

   if (ma!=-99) {
      wbMatrix<SPIDX_T> IJ; SPIDX_T d1=1,d2=1;

      if (!SIZE.len) { d1=d2=0; }
      else ma=get2DIndex(IJ,d1,d2,ma);

      mxSetFieldByNumber(a,i,4,
         Wb::mxCreateSparse(FL,d1,d2,IJ,D)
      );
   }
};


template <class TD>
mxArray* wbsparray<TD>::mxCreateSTRUCT(
   unsigned m, unsigned n
) const {
   const char *fields[]={"S","idx","data","A"};
   return mxCreateStructMatrix(m,n,4,fields);
};

template <class TD>
mxArray* wbsparray<TD>::add2MxSTRUCT(mxArray *a, unsigned i) const {

   wbMatrix<SPIDX_T> IJ; SPIDX_T d1=1,d2=1; char tflag='t';
   get2DIndex(IJ,d1,d2);

   mxSetFieldByNumber(a,i,0, SIZE.toMx());
   mxSetFieldByNumber(a,i,1, IDX.toMx());

#ifdef USE_WB_MPFR
   if (typeid(TD)==typeid(Wb::quad)) {
      tflag=30;
   }
#endif

   mxSetFieldByNumber(a,i,2, D.toMx(tflag));

   mxSetFieldByNumber(a,i,3,Wb::mxCreateSparse(FL,d1,d2,IJ,D));

   return a;
};


template <class TD> inline
mxArray* wbsparray<TD>::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[]={"S","data"};
   return mxCreateStructMatrix(m,n,2,fields);
};


#ifdef USE_WB_MPFR

template <> inline mxArray*
wbsparray<Wb::quad>::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[]={"S","idx","data"};
   return mxCreateStructMatrix(m,n,3,fields);
};

#endif


template <class TD>
void wbsparray<TD>::add2MxStruct(mxArray *a, unsigned i, int ma) const {

   wbMatrix<SPIDX_T> IJ;
   SPIDX_T d1,d2; get2DIndex(IJ,d1,d2,ma);
   int e=0;

   mxSetFieldByNumber(a,i,0, SIZE.toMx());
   mxSetFieldByNumber(a,i,1, Wb::mxCreateSparse(FL,d1,d2,IJ,D,NULL,&e));

   if (e && !isDiag()) {
      static unsigned ic=0;
      sprintf(str,"i%02d",ic++);
      MXPut(FL,str).add(*this,"A").addP(toMX(),"A0").add(e,"e");
      wblog(FL,"ERR %s() got e=%d",FCT,e);
   }
};


#ifdef USE_WB_MPFR

template <>
void wbsparray<Wb::quad>::add2MxStruct(
   mxArray *a, unsigned i, int ma  __attribute__ ((unused))
 ) const {

   mxSetFieldByNumber(a,i,0, SIZE.toMx());
   mxSetFieldByNumber(a,i,1, IDX.toMxT());
   mxSetFieldByNumber(a,i,2, D.toMx());
};

#endif


template <class TD> inline
mxArray* wbsparray<TD>::mxCreateCell(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
};

template <class TD>
void wbsparray<TD>::add2MxCell(mxArray *a, unsigned i) const {

   if (rank()>2) wblog(FL,
      "ERR %s() cell array expects rank <=2 (%s)",FCT,sizeStr().data);
   mxSetCell(a,i,toMxSp());

};


template <class TD>
mxArray* wbsparray<TD>::IDtoMx() const {

   const char *fields[]={"S","nnz","norm","stat"};
   mxArray *S=mxCreateStructMatrix(1,1,4,fields);

   const char *field2[]={"istat","idx","d"};
   mxArray *as=mxCreateStructMatrix(1,1,3,field2);

   if (isDiag()) {
      SPIDX_T s[2]={D.len,D.len};
      mxSetFieldByNumber(S,0,0, wbvector<SPIDX_T>(2,s).toMx());
   } else 
   mxSetFieldByNumber(S,0,0, SIZE.toMx());

   mxSetFieldByNumber(S,0,1, numtoMx(nnz()));
   mxSetFieldByNumber(S,0,2, numtoMx(double(norm())));

   if (D.len>3) {
      unsigned i=0; const unsigned m=3;
      const TD *d=D.data; TD a;

      wbindex I(m); INDEX_T *k=I.data;
      wbvector<double> x_(m); double *x=x_.data;
         x[0]=d[0];
         x[1]=ABS(d[0]);
         x[2]=d[0];

      for (++i; i<D.len; ++i) { a=ABS(d[i]);
         if (x[0]>d[i]) { x[0]=d[i]; k[0]=i; }
         if (x[1]>a   ) { x[1]=a;    k[1]=i; }
         if (x[2]<d[i]) { x[2]=d[i]; k[2]=i; }
      }

      mxSetFieldByNumber(as,0,0, I  .toMx());
      mxSetFieldByNumber(as,0,1, IDX.getRecs(I).toMx());
      mxSetFieldByNumber(as,0,2, x_ .toMx());
   }
   else if (D.len) {
      mxSetFieldByNumber(as,0,1, IDX.toMx());
      mxSetFieldByNumber(as,0,2, D  .toMx());
   }

   mxSetFieldByNumber(S,0,3,as);

   return S;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::initCAT(
   const char *F, int L, wbvector< wbsparray<TD>* > &X,
   wbvector<SPIDX_T> *I0, char del
){
   SPIDX_T i,id,l, r=0, i0=-1, i2=0, n=0, nz=0;
   wbvector<SPIDX_T> S;

   for (i=0; i<X.len; ++i) {
      if (!X[i] || X[i]->isEmpty()) continue;
      if (X[i]==this) {
         wbsparray<TD> A; A.initCAT(F,L,X,I0,del);
         return A.save2(*this);
      }
      X[i]->checkSize(F_LF);
      nz+=X[i]->nnz();

      if ((++n)==1) { i0=i; r=X[i]->SIZE.len+1;
         S.init(r,X[i]->SIZE);
      }; i2=i;

      if (i!=i0 && !X[i0]->sameSize(*X[i])) wblog(F_L,
         "ERR %s() severe size mismatch (%d: %s <> %d: %s)",FCT,
         i0+1, X[i0]->sizeStr().data,
         i+1, X[i]->sizeStr().data
      );
   }
   if (long(i0)<0) { init(); return *this; }

   char gotI0=(I0 && I0->len);
   if (gotI0 && I0->len!=X.len+1) wblog(F_L,
      "ERR %s() size mismatch (%d/%d+1)",FCT,I0->len,X.len);

   S.end()=(gotI0 ? I0->end() : n);
   init(S,nz);

   for (id=-1, l=0, i=i0; i<=i2; ++i) {
      if (X[i] && !X[i]->isEmpty()) {
         if (gotI0) { id=I0->data[i];
            if (id>=S.end()) { wblog(F_L,
               "ERR %s() index out of bounds (%d/%d)",FCT,id,S.end());
            }
         }
         else ++id;

         if (X[i]->nnz()) l+=catRecs(l,*X[i], id);
      }
   }

   if (l!=D.len) wblog(F_L,"ERR %s() %d/%d !??",FCT,l,D.len);

   if (del) for (i=0; i<X.len; ++i) {
      if (X[i]) { WB_DELETE_1(X[i]); }
      else break;
   }

   return Compress();
};


template <class TD>
SPIDX_T wbsparray<TD>::catRecs(
   SPIDX_T i0, const wbsparray<TD> &a, SPIDX_T id
){
   SPIDX_T n=a.IDX.dim1;

   if (IDX.dim2!=a.IDX.dim2+1) wblog(FL,
      "ERR %s() size mismatch (%d+1/%d)",FCT,IDX.dim2,a.IDX.dim2);
   if (i0+n>IDX.dim1) wblog(FL,
      "ERR %s() index out of bounds (%d+%d = %d ?)",FCT,i0,n,IDX.dim1);

   if (!n && a.D.len==1) {
      if (IDX.dim2!=1 || IDX.dim1!=D.len) wblog(FL,"ERR %s() "
         "%s <> %dx%d (%d)",FCT,sizeStr().data,IDX.dim1,IDX.dim2,i0); 
      if (a.D[0])
           { IDX(i0,0)=id; D[i0]=a.D[0]; return 1; }
      else { return 0; }
   }

   if (a.D.len!=n) wblog(FL,"ERR %s() %s <> %dx%d (%d)",
      FCT, a.sizeStr().data, a.IDX.dim1, a.IDX.dim2, i0);

   cpyRange(a.D.data, D.data+i0, n);
   Wb::cpyStride(
      IDX.rec(i0), a.IDX.data, a.IDX.dim2, a.IDX.dim1,
      IDX.dim2
   );

   for (SPIDX_T l=IDX.dim2-1, i=0; i<n; ++i) {
      IDX(i0+i,l)=id;
   }

   return n;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::Cat(
   unsigned dim,
   const wbsparray<TD> &B, const char *F, int L
){
   if (!dim || dim>SIZE.len) wblog(F_L,
      "ERR %s() dim out of bounds (%d/%d)",FCT,dim,SIZE.len);
   --dim;

   if (SIZE.len!=B.SIZE.len) wblog(F_L,
      "ERR %s() rank mismatch (%s <> %s; %d)",
      FCT, sizeStr().data, B.sizeStr().data, dim+1);
   for (unsigned i=0; i<SIZE.len; ++i) {
      if (i!=dim && SIZE[i]!=B.SIZE[i]) wblog(F_L,
      "ERR %s() size mismatch (%s <> %s; %d/%d)",
      FCT, sizeStr().data, B.sizeStr().data, i+1,dim+1);
   }
   checkSize(FLF); B.checkSize(FLF);

   SPIDX_T D1=SIZE[dim], d1=D.len, d2=B.D.len, d12=d1+d2;

   SIZE[dim]+=B.SIZE[dim]; D.Resize(d12, B.D.data);

   if (!IDX.dim2) {
      if (D1) wblog(F_L,
         "ERR %s() got %dx%d (%s)",FCT,IDX.dim1,IDX.dim1,sizeStr().data);
      IDX=B.IDX;
   }
   else {
      if (IDX.dim2!=B.IDX.dim2) wblog(F_L,
         "ERR %s() %d/%d",FCT,IDX.dim2,B.IDX.dim2);
      IDX.Resize(d12,IDX.dim2, B.IDX.data);
      for (SPIDX_T i=d1; i<d12; ++i) IDX(i,dim)+=D1;
   }

   return Compress();
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::initCAT(
   const char *F, int L, const wbvector< const wbsparray<TD>* > &X,
   unsigned dim
){
   if (!X.len) { init(); return *this; }

   SPIDX_T i,j, i0=-1, i2=0, n=0;
   wbvector<SPIDX_T> S0, DD(X.len), dd(X.len);

   for (i=0; i<X.len; ++i) { if (X[i] && !X[i]->isEmpty()) {
      const wbvector<SPIDX_T> &S2=X[i]->SIZE;
      if (X[i]==this) {
         wbsparray<TD> A; A.initCAT(F,L,X,dim);
         return A.save2(*this);
      }

      if (S2.len!=X[i]->IDX.dim2) X[i]->sperror_this(F_L);
      X[i]->checkSize(F_LF);

      if ((++n)==1) { i0=i; S0=S2;
         if (!dim || dim>S2.len) wblog(FL,"ERR %s() dim out of bounds "
            "(%s; %d)",FCT,S2.toStrD().data,dim);
         --dim;
      }
      else {
         if (S0.len!=S2.len) wblog(FL,"ERR %s() rank mismatch "
            "(%d: %s <> %d/%d: %s; %d)",FCT,i0+1,S0.toStrD().data,
            i+1,X.len, S2.toStrD().data, dim+1);
         for (j=0; j<S0.len; ++j) {
            if (j!=dim && S0[j]!=S2[j]) wblog(FL,"ERR %s() size mismatch "
            "(%d: %s <> %d/%d: %s; %d/%d)",FCT,i0+1,S0.toStrD().data,
            i+1,X.len, S2.toStrD().data, j+1,dim+1);
         }
      }

      DD[i]=S2[dim];
      dd[i]=X[i]->D.len; i2=i;
   }}

   SPIDX_T dall=dd.sum(), Dall=DD.sum();

   if (!n) { init(); return *this; }
   S0[dim]=Dall; init(S0,dall); if (!dall) return *this;

   SPIDX_T *pi=IDX.data,N,s, d1=0, d2, D1=0, m=S0.len;
   TD *pd=D.data;

   for (i=i0; i<=i2; ++i, d1+=n, D1+=N) { n=dd[i]; N=DD[i];
      s=n*m; MEM_CPY<SPIDX_T>(pi,s,X[i]->IDX.data); pi+=s;
      s=n  ; MEM_CPY<TD>(pd,s,X[i]->D.data); pd+=s;
      for (d2=d1+n, j=d1; j<d2; ++j) IDX(j,dim)+=D1;
   }
   if (d1!=IDX.dim1) wblog(FL,"ERR %s() %d/%d",FCT,d1,IDX.dim1);

   return Compress();
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::init(
   const char *F, int L, const mxArray *a, unsigned k
){

   init(); if (!a) {
      if (k) wblog(FL,"ERR %s() got k=%d for null mxArray !??",FCT,k);
      return *this;
   }

   if (mxIsEmpty(a) && mxIsDouble(a)) {
      SPIDX_T i,r=mxGetNumberOfDimensions(a);
      if (r) {
         const mwSize *sp=mxGetDimensions(a);
         wbvector<SPIDX_T> S(r);

         for (i=0; i<r; ++i) S[i]=(SPIDX_T)sp[i];
         init(S,0);
      }
      return *this;
   }

   if (mxIsStruct(a)) {
      const mxArray *as,*ad,*ai;
      if (!mxIsWbsparray(F_L,a,k,&as,&ad,&ai))wblog(FL,
         "ERR %s() invalid sparse array",FCT);

      SIZE.init(F_L,as);
      IDX .init(F_L,ai);
      D   .init(F_L,ad);

      check_IDX_range(F_L);

      if (IDX.isEmpty()) {
         if (!SIZE.isEmpty() && !IDX.dim2) {
            if (IDX.dim1) wblog(FL,
               "ERR %s() got size %s !??",FCT,sizeStr().data);
            IDX.init(IDX.dim1,SIZE.len);
         }
      }
      if (!isDiag()) {
         if (IDX.dim1!=D.len || (IDX.dim1 && IDX.dim2!=SIZE.len))
         sperror_this(F_L);
      }

      return Compress();
   }

   if (mxIsCell(a)) {
      a=mxGetCell(a,k); if (!a) return *this;
   }
   else if (k) wblog(FL,"ERR %s() got k=%d !??",FCT,k);

   if (mxIsSparse(a)) {

      SPIDX_T j,k,l,m,n,d, nnz=mxGetNzmax(a);
      int cflag=mxIsComplex(a);
      double *dr,*di=0;
      mwIndex *Ir,*Jc;

      if (!mxIsDouble(a) || mxGetNumberOfDimensions(a)!=2) wblog(FL,
         "ERR %s() invalid input (got %s)",FCT,mxTypeSize2Str(a).data);
      if (cflag && typeid(TD)!=typeid(wbcomplex)) { wblog(FL,
         "WRN %s() got complex sparse array for %s",
          FCT,getName(typeid(TD)).data); cflag=-cflag;
      }

      m=mxGetM(a); n=mxGetN(a);
      initz(m,n,nnz);

      Ir=mxGetIr(a);  dr=mxGetPr(a);
      Jc=mxGetJc(a);  di=mxGetPi(a);

      for (l=j=0; j<n; ++j) {
         d=SPIDX_T(Jc[j+1]-Jc[j]); if (d==0) continue;

         for (k=0; k<d; ++k,++l) {
            IDX(l,0)=SPIDX_T(Ir[l]); IDX(l,1)=j;
            if (dr[l]!=0.) D[l]=dr[l]; if (cflag) {
            if (di[l]!=0.) DSET_IMAG(D[l],di[l]); }
         }
      }
      if (l!=nnz && l>0) wblog(FL,
         "ERR severe sparse inconsistency (mex: %d,%d)",l,nnz
      );
   }
   else {

      if (!mxIsDouble(a)) wblog(FL,
         "ERR %s() got `%s' data",FCT,mxGetClassName(a));
      if (typeid(TD)==typeid(wbcomplex)) wblog(FL,
         "WRN %s() ignores complex input data",FCT);

      SPIDX_T r=mxGetNumberOfDimensions(a);
      const mwSize *sp=mxGetDimensions(a);
      const double *dr=mxGetPr(a);

      wbvector<SPIDX_T> S(r), I(r);
      SPIDX_T i,k,N, nz=0, iz=0, l=r-1, *s=S.data, *ip=I.data;

      for (i=0; i<r; ++i) S[i]=(SPIDX_T)sp[i];

      N=S.prod(0);
      for (i=0; i<N; ++i) { if (dr[i]!=0) ++nz; }
      init(S,nz);

      for (i=0; i<N; ++i) {
          if (dr[i]!=0) {
             if (iz>=nz) wblog(FL,
                "ERR %s() index out of bounds (%d/%d) !??",FCT,iz+1,nz);
             D[iz]=dr[i]; IDX.recSetP(iz,ip); ++iz;
          }

          k=0; ++ip[0];
          while (ip[k]>=s[k] && k<l) { ip[k]=0; ++ip[++k]; }
      }
   }

   return *this;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::init(
   const char *F, int L, const wbarray<TD> &A, TD eps
){
   unsigned k, r=A.rank(), l=r-1;
   wbvector<SPIDX_T> I(r);
   SPIDX_T i,N, nz=0, iz=0, *ip=I.data; const size_t *s=A.SIZE.data;

   N=A.numel();
   for (i=0; i<N; ++i) { if (ABS(A.data[i])>eps) ++nz; }
   init(A.SIZE,nz);

   for (i=0; i<N; ++i) {
       if (ABS(A.data[i])>eps) {
          if (iz>=nz) wblog(F_L,
             "ERR %s() index out of bounds (%d/%d) !??",FCT,iz+1,nz);
          D[iz]=A.data[i]; IDX.recSetP(iz,ip); ++iz;
       }

       k=0; ++ip[0];
       while (ip[k]>=s[k] && k<l) { ip[k]=0; ++ip[++k]; }
   }

   return *this;
};


template <class TD>
wbarray<TD>& wbsparray<TD>::toFull(wbarray<TD> &A) const {

   if (isDiag()) { A.initDiag(D); }
   else {
      unsigned k, r=rank(FL), l=r-1;
      SPIDX_T i,j; const SPIDX_T *s=SIZE.data, *ip=IDX.data;

      A.init(SIZE);
      if (IDX.dim1 && A.data!=NULL) {
         for (i=0; i<IDX.dim1; ++i, ip+=IDX.dim2) {
             for (j=ip[l], k=l-1; k<r; --k) { j=j*s[k]+ip[k]; }
             A.data[j]=D.data[i];
         }
      }
   }
   return A;
};


template <class TD>
template <class T2>
wbvector<T2>& wbsparray<TD>::toFull(wbvector<T2> &A) const {

   if (SIZE.len) { 
      unsigned m=IDX.dim2, r=rank(FL);
      SPIDX_T i=0; const SPIDX_T *s=SIZE.data, *ip=0;

      for (; i<r; ++i) { if (s[i]!=1) {
         if (!ip) { ip=IDX.data+i; }
         else wblog(FL,
            "ERR %s() invalid usage (got %s tensor)",
            FCT,sizeStr().data
         );
      }}
      if (!ip) ip=IDX.data;

      A.init(numel());
      for (i=0; i<D.len; ++i, ip+=m) {
         if ((*ip)>=A.len) wblog(FL,
            "ERR %s() invalid sparse setting (%d/%d)",FCT,*ip,A.len);
         if (A[*ip]) wblog(FL,
            "ERR %s() invalid sparse setting (%d: %g)",FCT,*ip,A[*ip]);
         A[*ip]=D.data[i];
      }
   }
   else { A.init(1);
      if (D.len) {
         if (D.len>1) wblog(FL,
            "ERR %s() invalid scalar (len=%d)",FCT,D.len);
         A.data[0]=D.data[0];
      }
      else { A.data[0]=0; }
   }
   return A;
};


template <class TD> inline
wbsparray<TD>& wbsparray<TD>::reshape( 
   const wbvector<SPIDX_T> &S, wbsparray<TD> &A
 ) const {

   if (S==SIZE) { A=*this; return A; }
   if (sameSizeUp2Singletons(S)) {
      unsigned j; SPIDX_T i=0, *I; const SPIDX_T *I0=IDX.data; int *ms;
     
      wbvector<int> Ms; matchWithSingletons(FL,S,Ms);
      A.init(S); A.D=D; A.IDX.init(IDX.dim1,S.len);
      ms=Ms.data; I=A.IDX.data;
      
      for (; i<IDX.dim1; ++i, I+=S.len, I0+=IDX.dim2) {
         for (j=0; j<S.len; ++j) {
            if (ms[j]>=0) { I[j]=I0[ms[j]]; }
         }
      }
      return A;
   }

   if (numel()!=S.prod(0)) wblog(FL,
      "ERR %s() got size mismatch (%s => %s)",
      FCT,sizeStr().data, S.toStrD().data
   );

   A.init(S); A.D=D; A.IDX.init(IDX.dim1,S.len);

   if (!IDX.isEmpty()) {
      SPIDX_T i,k, *I=A.IDX.data; unsigned r0=SIZE.len, r=S.len;
      const SPIDX_T *I0=IDX.data, *S0=SIZE.data, *S=A.SIZE.data;

      for (i=0; i<IDX.dim1; ++i, I0+=r0, I+=r) {
         k=Wb::sub2ind(S0,I0,r0);
         Wb::ind2sub(k,S,I,r);   
      }
   }
      
   return A;
};


template <class TD>
template <class DB> inline
wbsparray<TD>& wbsparray<TD>::init(
   const char *F, int L,
   SPIDX_T d1, SPIDX_T d2, const wbsparray<DB> &B
){
   SPIDX_T s[]= {d1,d2};
   wbvector<SPIDX_T> S; S.init(2,s);

   if (B.numel()!=S.prod(0)) wblog(F_L,"ERR %s() size mismatch "
      "(%s <> %s)",FCT,S.toStrD().data,B.sizeStr().data);
   (*this)=B; return Reshape(S);
};


template <class TD>
template <class DB> inline
wbsparray<TD>& wbsparray<TD>::init(
   const char *F, int L,
   SPIDX_T d1, SPIDX_T d2, SPIDX_T d3, const wbsparray<DB> &B
){
   SPIDX_T s[]= {d1,d2,d3};
   wbvector<SPIDX_T> S; S.init(3,s);

   if (B.numel()!=S.prod(0)) wblog(F_L,"ERR %s() size mismatch "
      "(%s <> %s)",FCT,S.toStrD().data,B.sizeStr().data);
   (*this)=B; return Reshape(S);
};


template <class TD>
template <class DB> inline
wbsparray<TD>& wbsparray<TD>::init(
   const char *F, int L,
   SPIDX_T d1, SPIDX_T d2, SPIDX_T d3, SPIDX_T d4, const wbsparray<DB> &B
){
   SPIDX_T s[]= {d1,d2,d3,d4};
   wbvector<SPIDX_T> S; S.init(4,s);

   if (B.numel()!=S.prod(0)) wblog(F_L,"ERR %s() size mismatch "
      "(%s <> %s)",FCT,S.toStrD().data,B.sizeStr().data);
   (*this)=B; return Reshape(S);
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::select0(
   const WBINDEX &P, unsigned dim,
   wbsparray<TD> &X
) const {

   if (&X==this) {
      wbsparray<TD> Q; X.save2(Q);
      return Q.select0(P,dim,X);
   }
   if (dim>=SIZE.len) wblog(FL,
      "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len);

   SPIDX_T i,l,k, n=0, d0=SIZE[dim], *I, *j;

   wbvector<SPIDX_T> I_(d0); I_.set(-1); I=I_.data;

   for (i=0; i<P.len; ++i) {
      if (P[i]>=d0) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,P[i],d0);
      if (long(I[P[i]])>=0) wblog(FL,
         "ERR %s() input index not unique (%d/%d)",FCT,P[i],d0);
      I[P[i]]=i;
   }

   for (j=IDX.data+dim, i=0; i<IDX.dim1; ++i, j+=IDX.dim2) {
      if (long(I[*j])>=0) ++n;
   }

   X.init(SIZE,n); X.SIZE[dim]=P.len;
   
   for (j=IDX.data+dim, l=i=0; i<IDX.dim1; ++i, j+=IDX.dim2) { k=I[*j];
      if (long(k)>=0) {
         X.IDX.recSetP(l,IDX.rec(i)); X.IDX(l,dim)=k;
         X.D[l++]=D[i];
      }
   }

   return X;
};


template <class TD> inline
wbsparray<TD>& wbsparray<TD>::Permute(const wbperm &P, char iflag){

   if (!SIZE.len && !IDX.dim1 && D.len<=1) return *this;
   if (!P.isEmpty()) {
      SIZE.Permute(P,iflag); IDX.colPermute(P,iflag);
      Sort();
   }
   return *this;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::permute(
   const wbperm &P, wbsparray<TD>&B, char iflag
) const {

   if (this==&B) return B.Permute(P,iflag);

   if (P.isEmpty() || (!SIZE.len && !IDX.dim1 && D.len<=1)) { B=(*this); }
   else {
      SIZE.permute(B.SIZE,P,iflag); B.D=D; B.isref=0;
      IDX.colPermute(P,B.IDX,iflag);
      B.Sort();
   }
   return B;
};


template <class TD> inline
wbsparray<TD>& wbsparray<TD>::Permute(const char* s, char iflag){
   return Permute(Str2Idx(s,1),iflag);
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::permute(
   const char* s, wbsparray<TD>&Q, char iflag
) const { return permute(Q,Str2Idx(s,1),iflag); };


template <class TD> inline
wbsparray<TD>& wbsparray<TD>::ColPermute(const wbperm &P0, char iflag){

   if (isIdentityPerm(P0)) return *this;

   if (SIZE.len!=2 || IDX.dim2!=2) wblog(FL,
      "ERR %s() requires rank-2 (%s; %d)",FCT,sizeStr().data,IDX.dim2);
   if (P0.len!=SIZE[1]) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,P0.len,SIZE[1]);

   SPIDX_T k,i=0, l=IDX.dim2-1; wbperm P;

   if (!iflag)
        P0.getIPerm(P);
   else P.wbvector<PERM_T>::init2ref(P0);

   for (; i<IDX.dim1; ++i) { k=IDX(i,l);
      if (k>=P.len) wblog(FL,"ERR %s() %d/%d !??",FCT,k,P.len); 
      IDX(i,l)=P[k];
   }
   return Sort();
};


template <class TD> inline
wbsparray<TD>& wbsparray<TD>::MatPermute(const wbperm &P0, char iflag){

   if (isIdentityPerm(P0)) return *this;

   SPIDX_T d=-1;
   if (!isSMatrix(FL,&d)) wblog(FL,
      "ERR %s() requires rank-2 (%s; %d)",FCT,sizeStr().data,IDX.dim2);
   if (P0.len!=d) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,P0.len,d);

   if (isDiag()) { D.Permute(P0); }
   else {
      wbperm P;
      if (!iflag)
           P0.getIPerm(P);
      else P.wbvector<PERM_T>::init2ref(P0);

      SPIDX_T *I=IDX.data, n=IDX.numel(); PERM_T *p=P.data;
      for (SPIDX_T i=0; i<n; ++i) {
         if (I[i]>=P.len) wblog(FL,"ERR %s() %d/%d !??",FCT,I[i],P.len); 
         I[i]=p[I[i]];
      }
   }

   return Sort();
};


template <class TD>
void wbsparray<TD>::setSIZE_kron(
   const char *F, int L,
   const wbvector<SPIDX_T> &sa, const wbvector<SPIDX_T> &sb,
   char kflag
){
   wbvector<SPIDX_T> S;
   if (sa.len!=sb.len) wblog(F_L,
      "ERR %s() rank mismatch (%d/%d)",FCT,sa.len,sb.len);

   if (kflag) { S.init(sa.len); SPIDX_T *s=S.data;
      for (unsigned i=0; i<sa.len; ++i) {
         s[i]=sa[i]*sb[i];
      }
   }
   else { S.init(sa.len+sb.len); SPIDX_T *s=S.data;
      for (unsigned i=0, l=0; i<sa.len; ++i) {
         s[l++]=sa[i];
         s[l++]=sb[i];
      }
   }
   S.save2(SIZE);
};

template <class TD>
void wbsparray<TD>::setSIZE_kron(const char *F, int L, char kflag) {

   if (kflag) {
      unsigned l=0, i=0, r=IDX.dim2; SPIDX_T *s=SIZE.data;
      if (SIZE.len!=2*r) wblog(F_L,
         "ERR %s() unexpected size (%dx%d <> %d/(%d/2))",
         FCT,IDX.dim1,IDX.dim2,D.len,SIZE.len
      );
      if (r) {
         for (; i<r; ++i,l+=2) { s[i]=s[l]*s[l+1]; }
         SIZE.len=r;
      }
   }
   else if (SIZE.len!=IDX.dim2) sperror_this(F_L);
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::setRec_kron(
   SPIDX_T l, TD x, unsigned r, SPIDX_T *i1, SPIDX_T *i2,
   unsigned kflag
){
   if (l>IDX.dim1) wblog(FL,"ERR %s() index "
      "out of bounds (%d,%d <> %dx%d)",FCT,l,r,IDX.dim1,IDX.dim2);
   
   if (kflag) {
      if (l==0 && (SIZE.len!=2*r || IDX.dim2!=r))
         wblog(FL,"ERR %s() 2*%d/2*%d/%d",FCT,r,IDX.dim2,SIZE.len);
      SPIDX_T *q=IDX.rec(l);
      for (unsigned i=0; i<r; ++i) {
         q[i]=i1[i] + SIZE[2*i]*i2[i];
      }
   }
   else {
      if (l==0 && (SIZE.len!=2*r || SIZE.len!=IDX.dim2))
         wblog(FL,"ERR %s() index out of bounds (%d,%d <> %dx%d)",
         FCT,l,r,IDX.dim1,IDX.dim2
      );
   
      SPIDX_T *q=IDX.rec(l);
      for (unsigned j=0, i=0; i<r; ++i) {
          q[j++]=i1[i];
          q[j++]=i2[i];
      }
   }
   D[l]=x; return *this;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::transpose(
   const char *F, int L, wbsparray<TD> &A
 ) const {

   unsigned r=rank(); wbperm P;

   if (r%2) wblog(F_L,"ERR %s() "
      "applies to even-rank arrays only! (%s)",FCT,sizeStr().data);

   P.initTranspose(r);
   permute(P,A);

   return A.Sort();
};


template <class TD> inline
wbsparray<TD>& wbsparray<TD>::tensorProd(const char *F, int L,
   const wbsparray<TD>& B0, wbsparray<TD>& C,
   char aflag0, char bflag0,
   char kflag
 ) const {

   char iA=isDiag(FL), iB=B0.isDiag(FL);

   if (iA || iB) {
      if (iA && iB) {
          if (kflag) {
             C.initDiag(D.len*B0.D.len); D.tensorProd(B0.D,C.D);
             return C;
          }
          wbsparray<TD> A2;    diag2reg(FL,A2);
          wbsparray<TD> B2; B0.diag2reg(FL,B2);
          return A2.tensorProd(F,L,B2,C,aflag0,bflag0,kflag);
      }
      else if (iA) {
          wbsparray<TD> A2; diag2reg(FL,A2, B0.rank());
          return A2.tensorProd(F,L,B0,C,aflag0,bflag0,kflag);
      }
      else {
          wbsparray<TD> B2; B0.diag2reg(FL,B2, rank());
          return tensorProd(F,L,B2,C,aflag0,bflag0,kflag);
      }
   }

   if ((void*)this==(void*)(&C) || (void*)(&B0)==(void*)(&C)) {
       wbsparray<TD> X; tensorProd(F_L,B0,X,aflag0,bflag0,kflag);
       return X.save2(C);
   }

   unsigned ra=rank(), rb=B0.rank();

   if (!ra || ra%2 || ra!=rb) wblog(F_L,"ERR %s() invalid input to "
      "kron(%s,%s)",FCT,sizeStr().data,B0.sizeStr().data); 
   if (SIZE.len!=IDX.dim2 || B0.SIZE.len!=B0.IDX.dim2) wblog(FL,
      "ERR %s() data inconsistency (%s -> %dx%d, %s -> %dx%d)",FCT,
         sizeStr().data,   IDX.dim1,   IDX.dim2,
      B0.sizeStr().data,B0.IDX.dim1,B0.IDX.dim2);

   opFlags<TD> aflag(aflag0), bflag(bflag0);
   wbsparray<TD> A(*this),B(B0);

   A.Compress(); aflag.apply(F_L,A);
   B.Compress(); bflag.apply(F_L,B);

   SPIDX_T i,j, l=0, na=A.IDX.dim1, nb=B.IDX.dim1;
   C.init_kron(SIZE,B.SIZE,na*nb,kflag);
   
   if (C.IDX.dim1) {
      for (i=0; i<na; ++i)
      for (j=0; j<nb; ++j,++l) {
         C.setRec_kron(
           l, A.D[i]*B.D[j], ra, A.IDX.rec(i), B.IDX.rec(j), kflag
         );
      }
   }

   C.setSIZE_kron(F_L,kflag);


   return C;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::Plus(
   const char *F, int L,
   const wbsparray<TD> &B, TD bfac, TD afac, char strict
){
   if (isref) wblog(FL,"ERR %s() got isref=%d !??",FCT,isref);
   if ((void*)this==(void*)(&B)) {
      TD x=afac+bfac;
      if (x== 0) { D.init(); IDX.init(0,IDX.dim2); isref=0;    } else
      if (x==-1) { for (SPIDX_T i=0; i<D.len; ++i) D[i]=-D[i];} else
      if (x!=+1) { for (SPIDX_T i=0; i<D.len; ++i) D[i]*=x;   }
      return *this;
   }

   char sab=checkSize2(F_L,B);

   if (afac!=1) {
      if (afac== 0) { D.init(); IDX.init(0,IDX.dim2); isref=0;    } else
      if (afac==-1) { for (SPIDX_T i=0; i<D.len; ++i) D[i]=-D[i];} else
      if (afac!=+1) { for (SPIDX_T i=0; i<D.len; ++i) D[i]*=afac;}
      afac=1;
   }

   if (bfac==0) return *this;

   if (sab==1) {
      SPIDX_T n1=D.len, n2=B.D.len;
      IDX.cat(1,B.IDX); D.cat(B.D);

      if (bfac!=1) { TD* d=D.data+n1;
         if (bfac==-1)
              { for (SPIDX_T i=0; i<n2; ++i) d[i]=-d[i];}
         else { for (SPIDX_T i=0; i<n2; ++i) d[i]*=bfac;}
      }
      Compress(F_L);
   }
   else if (strict) wblog(F_L,"ERR %s() got size mismatch",FCT);
   else if (sab>=20) {
      wblog(FL,"ERR %s() got singletons %s <> %s",
         FCT,sizeStr().data,B.sizeStr().data);

   }
   else if (sab>=10) {
      char iA=((sab-10) & 1),
           iB=((sab-10) & 2);
      if (sab==10 || sab>13) wblog(FL,"ERR %s() got sab=%d",FCT,sab);

      if (iA && iB) D.Plus(B.D,bfac);
      else if (iA) {
          wbsparray<TD> A2; diag2reg(FL,A2);
          A2.Plus(F,L,B,bfac,afac,1);
          A2.save2(*this);
      }
      else {
          wbsparray<TD> B2; B.diag2reg(FL,B2);
          Plus(F,L,B2,bfac,afac,1);
      }
   }
   else wblog(FL,"ERR %s() got sab=%d !??",FCT,sab);

   return *this;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::Compress(
   const char *F, int L, TD eps, char lex, const wbperm *Pc
){
   wbperm P; WBINDEX dg; rank(F_L);

   checkSize(F_L);
   if (!IDX.dim1) {
      if (Pc && Pc->len!=2) wblog(FL,
         "ERR %s() invalid permutation (%s)",FCT,Pc->toStr().data); 
      if (double(eps)>0) {
         for (SPIDX_T k=0; k<D.len; ++k) {
            TD &x=D.data[k]; if (x && ABS(x)<=eps) x=0;
         }
      }
      return *this;
   }

   if (!Pc || Pc->isIdentityPerm()) {
      if (double(eps)<=0 && IDX.isUniqueSorted(+1,lex)) { return *this; }
      IDX.groupRecs(P,dg,-1,lex);
   }
   else {
      wbMatrix<SPIDX_T> X; IDX.colPermute(*Pc,X);
      X.groupRecs(P,dg,-1,lex);
   }

#ifdef WB_SPARSE_CLOCK
   Wb::UseClock spf(&wbc_sparse_cmpr);
#endif

   SPIDX_T i,j,k,d; PERM_T *p=P.data;
   wbvector<TD> D0; D.save2(D0); D.init(dg.len);
   const TD *d0=D0.data;

   for (j=k=0; j<dg.len; ++j,++k) {
      TD &x=D.data[k]; x=d0[p[0]];
      for (d=dg.data[j], i=1; i<d; ++i) { x+=d0[p[i]]; }

      p+=d;
      if (x==0) { --k; } else if (k<j) { IDX.recSet(k,j); }
   }

   if (k!=IDX.dim1) {
      if (k>IDX.dim1) wblog(FL,"ERR %s() !??",FCT);
      if (k) { IDX.dim1=k; D.len=k; }
      else { IDX.init(0,SIZE.len); D.init(); }
   }

   SkipTiny(eps); return *this;
};


template <class TD>
int wbsparray<TD>::splitSparseCM(
   const char *F, int L, wbvector< wbsparray<TD> > &X,
   wbvector<SPIDX_T> *I0
){
   SPIDX_T l=0, n;
   wbsparray<TD> a;

   iterSparseCMref<TD> ISP(F_L,*this); n=ISP.numIter();

   if (!n) { X.init(); return X.len; }

   X.init(n);
   while (ISP.next(a)) { X[l++]=a; }
   if (l!=X.len) wblog(F_L,"ERR %s() %d/%d !??",FCT,l,X.len);


   if (I0) { (*I0)=ISP.IX; }
   return X.len;
};



template <class T>
inline int checkContract(const char* F, int L,
   const wbvector<T> &SA, const wbindex &ica,
   const wbvector<T> &SB, const wbindex &icb,
   const wbperm *pfinal, wbvector<T> *SC
){
   unsigned i,e=0;

   if (ica.len!=icb.len) {
      if (F) wblog(F,L,
         "ERR %s() invalid index set (length mismatch, [%s], [%s])",
          FCT,ica.toStr().data,icb.toStr().data);
      return 1;
   }
   if (!ica.len) {
      if (F) wblog(F,L,"ERR contract() got empty index set");
      return 1;
   }

   if (!e)
   for (i=0; i<ica.len; ++i) {
      if ((ica[i]<SA.len ? SA[ica[i]] : 1) !=
          (icb[i]<SB.len ? SB[icb[i]] : 1) ){ ++e; break; }
   }

   if (e) { wblog(FL,
      "ERR invalid contraction %s (%s) * %s (%s) !??",
       SA.len ? SA.toStrD().data : "[]", (ica+1).toStr().data,
       SB.len ? SB.toStrD().data : "[]", (icb+1).toStr().data);
      return 2;
   }

   if (SC) {
      wbvector<T> sC; sC.initI(SA,ica,SB,icb);

      if (SC->isEmpty()) { sC.save2(*SC); }
      else {
         if (pfinal && pfinal->len) sC.Permute(*pfinal);
         if ((*SC)!=sC) { e=6; if (F) wblog(FL,
           "ERR %s() severe size mismatch (C: %s <> %s)",
            FCT, sC.toStrD().data, SC->toStrD().data);
         }
      }

      if (!SC->prod(0)) return -1;
   }

   return e;
};


template <class TD> inline
char wbsparray<TD>::contract_scalar(
   const char* F, int L,  const wbindex &ica,
   const wbsparray<TD> &B, const wbindex &icb, wbsparray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
 ) const {

   char isa=isScalar(), isb=B.isScalar(); unsigned r=-1;
   if (!isa || !isb) return 0;

   TD x=(D.len && B.D.len ? afac*D[0]*B.D[0] : TD(0));
   if (cfac && C.D.len) { x+=cfac*C.D[0]; }

   if (ica.len!=icb.len) wblog(F_L,"ERR %s() length mismatch: "
      "%s <> %s",FCT, ica.toStr().data, icb.toStr().data); 
   if (cfac && (!C.isScalar() || C.D.len>1)) wblog(F_L,
      "ERR %s() invalid scalar %s",FCT,C.info2Str("C").data);
   if (SIZE.len!=IDX.dim2 || B.SIZE.len!=B.IDX.dim2) wblog(F_L,
      "ERR %s() invalid scalars %s <> %s", FCT,
      info2Str("A").data, B.info2Str("B").data); 

   if (SIZE.len) {
      if (B.SIZE.len) { 
         checkContract(F_L, SIZE,ica,B.SIZE,icb,
            (cfac? &pfinal : NULL), (cfac? &C.SIZE : NULL));
         r=(SIZE.len + B.SIZE.len - 2*ica.len);
      }
      else {
         wbvector<SPIDX_T> Sb(2); Sb[0]=Sb[1]=B.D.len;
         checkContract(F_L, SIZE,ica,Sb,icb,
            (cfac? &pfinal : NULL), (cfac? &C.SIZE : NULL));
         r=(SIZE.len + 2 - 2*ica.len);
      }
   }
   else {
      if (B.SIZE.len) { 
         wbvector<SPIDX_T> Sa(2); Sa[0]=Sa[1]=D.len;
         checkContract(F_L, Sa,ica,B.SIZE,icb,
            (cfac? &pfinal : NULL), (cfac? &C.SIZE : NULL));
         r=(2 + B.SIZE.len - 2*ica.len);
      }
   }

   C.initScalar(x,r);

   return 1;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::contract_diag_diag(
   const char* F, int L,  const wbindex &ica,
   const wbsparray<TD> &B, const wbindex &icb, wbsparray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
 ) const {

   SPIDX_T i,e=0, ra=2, rb=2;

   if (SIZE.len || B.SIZE.len) wblog(F_L,
      "ERR %s() (%s) <> (%s)",FCT,sizeStr().data,B.sizeStr().data);

   if  (ica.len!=icb.len || ica.len>2) e=2; else
   for (i=0; i<ica.len; ++i) if (ica[i]>=ra || icb[i]>=rb) e=3;
   if (!e && ica.len==2 && (ica[0]==ica[1] || icb[0]==icb[1])) e=4;
   if (e) wblog(F_L,
      "ERR %s() size mismatch %s (%s) <> %s (%s)",
      FCT,sizeStr().data,ica.toStr().data,B.sizeStr().data,
      icb.toStr().data
   );

   wbvector<TD> x(D.len);
   if (afac) {
      for (i=0; i<D.len; ++i) { x[i]=D[i]*B.D[i]; }
   }

   if (ica.len==2) { TD q=afac*x.sum();
      if (cfac) {
         if (!C.isScalar()) wblog(F_L,
            "ERR %s() invalid scalar C (%s)",FCT,C.sizeStr().data);
         C.D*=cfac; C[0]+=q;
      }
      else C.initScalar(q);
   }
   else {
      if (cfac) {
         if (!C.isDiag() || C.D.len!=D.len) wblog(F_L,"ERR %s() "
            "invalid diag C (%s; %d)",FCT,C.sizeStr().data,D.len);
         C.D*=cfac; C.D.Plus(x,afac);
      }
      else {
         if (afac) x*=afac;
         C.initDiag(D.len); x.save2(C.D);
      }
   }

   if (!pfinal.isEmpty() && pfinal.isIdentityPerm())
   wblog(F_L,"ERR %s() got pfinal=[%s]",FCT,pfinal.toStr().data); 

   return C;
};


template <class TD>
char wbsparray<TD>::contract_check_2full(
   const char* F, int L,  const wbindex &ica,
   const wbsparray<TD> &B, const wbindex &icb, wbsparray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac,
   char flag, double eps, SPIDX_T *nnz_C
) const {

   if (flag=='s' || flag=='S'
       || typeid(TD)!=typeid(double)
   ) return 0;

   if (flag && flag!='f' && flag!='F')
   wblog(F_L,"ERR %s() invalid flag=%c<%d>",FCT,flag,flag);

   if (!flag || nnz_C!=NULL) {
      SPIDX_T NA=numel(), NB=B.numel(), M=SIZE.prod(ica.data,ica.len);
      SPIDX_T NC=(NA/M)*(NB/M);
      double pa=nnz()/double(NA), pb=B.nnz()/double(NB);

      double pc=M*pa*pb; if (pc>1) pc=1;

      if (!NA || !NB || !M) wblog(F_L,
         "ERR %s() got\n%s [%s] <> %s [%s] (%d)",FCT,sizeStr().data,
         ica.toStr().data, B.sizeStr().data, icb.toStr().data, M);

      if ((pc>0.50 && NC<(1<<24)) ||
          (pc>0.20 && NC<(1<<16)) ||
          (pc>0.10 && NC<(1<<12))
       ) { flag='F';
         wbvector<SPIDX_T> sC; sC.initI(SIZE,ica,B.SIZE,icb);
      }

      if (nnz_C) {
         double nc=pc*NC;
         if (NC>128) (nc)*=1.20; else
         if (NC> 32) (nc)*=1.50; else
         if (NC>  4) (nc)*=2; else nc=NC;
         (*nnz_C)=MIN(SPIDX_T(nc),NC);
      }

   }

   if (!flag) return 0;

#ifdef WB_SPARSE_CLOCK
   Wb::UseClock spf(&wbc_sparse_cntf);
#endif


   wbarray<TD> Af,Bf,ABf; this->toFull(Af); B.toFull(Bf);

   Af.contract(F_L,ica,Bf,icb,ABf,pfinal,afac);

   if (cfac!=0 && C.D.len) { C*=cfac;
      wbarray<TD> Cf; C.toFull(Cf);
      ABf+=Cf;
   }

   C.init(F_L,ABf,eps); return 1;
};




template <class TD>
wbsparray<TD>& wbsparray<TD>::contract(
   const char* F, int L,  const wbindex &ica,
   const wbsparray<TD> &B, const wbindex &icb, wbsparray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac,
   char sflag,
   double eps
) const {

#ifdef WB_SPARSE_CLOCK
   Wb::UseClock spc(&wbc_sparse_cont);
#endif

   if (D.len<=1 && B.D.len<=1) {
   if (contract_scalar(F_L,ica,B,icb,C,pfinal,afac,cfac)) {
      return C;
   }}

   if (!sflag && typeid(TD)!=typeid(double)) sflag='s';

   if (this==&C || &B==&C) { wbsparray<TD> X(C);
      contract(F,L,ica,B,icb,X,pfinal,afac,cfac,sflag,eps);
      return X.save2(C);
   }


   char iA=isDiag(FL), iB=B.isDiag(FL);
   if (iA || iB) {
      if (iA && iB) {
         return contract_diag_diag(F_L,ica,B,icb,C,pfinal,afac,cfac);
      }
      else if (iA) {
          wbsparray<TD> A2; diag2reg(FL,A2);
          return A2.contract(F,L,ica,B,icb,C,pfinal,afac,cfac,sflag,eps);
      }
      else if (iB) {
          wbsparray<TD> B2; B.diag2reg(FL,B2);
          return contract(F,L,ica,B2,icb,C,pfinal,afac,cfac,sflag,eps);
      }
   }

   wbvector<SPIDX_T> sC; if (cfac) sC=C.SIZE; int e=0;

   if ((e=checkContract(F_L,SIZE,ica,B.SIZE,icb, &pfinal, &sC))) {
      if (e<0) {
         if (sC.len) {
            C.init(sC);
            return C;
         }
      }
      else wblog(FL,"ERR %s()",FCT);
   }

   if (!D.len || !B.D.len) {
      if (cfac) { C.D*=cfac; } else { C.init(sC); }
      if (e<0) {
         if (C.SIZE.len || C.IDX.dim2) wblog(FL,"ERR %s() !??",FCT);
         if (!C.D.len) { C.D.init(1); C.D[0]=0; }
      }
      return C;
   }



   if (contract_check_2full(F_L,
         ica,B,icb,C,pfinal,afac,cfac,sflag,eps)) { return C; }

#ifdef WB_SPARSE_CLOCK
   Wb::UseClock sp2(&wbc_sparse_cnt2);
#endif

   SPIDX_T i,l,ja,jb,ic,is, nnzc=0;
   INDEX_T ka,la,kb,lb, kgb=-1, jga=-1, nc=ica.len;
   unsigned r1=SIZE.len-nc, r2=B.SIZE.len-nc;
   const PERM_T *pa,*pb;
   bool sortflag;
   const char lCM=-1;"col-major")
   TD bjk;

   groupIndex<INDEX_T> gIJ,gIA,gJK,gJB,gKC;
   wbMatrix<SPIDX_T> IJ,JK,IA,JA,JB;
   wbvector<SPIDX_T> u,ii;
   wbindex Ja,Jb,gia,gjb;
   wbvector<TD> w;
   wbsparray<TD> C0;

   if (cfac!=0) { C.D*=cfac; } 
   if (afac==0) {
      if (!cfac) C.init(sC.Permute(pfinal));
      return C; 
   }
   if (cfac!=0) C.save2(C0);





#ifdef WB_SPARSE_CLOCK
   sp2.done();
   Wb::UseClock sp3(&wbc_sparse_cnt3);
#endif

     IDX.cols2End  (ica,IJ).groupRecs(gIJ,nc,lCM,&pa);
   B.IDX.cols2Front(icb,JK).groupRecs(gJK,r2,lCM,&pb);

   IJ.getCols(0,r1-1,IA).groupRecs(gIA,-1,lCM,&gia);
   JK.getCols(0,nc-1,JB).groupRecs(gJB,-1,lCM,&gjb);"col-major" C=AB

   gIJ.getRecs0(IJ,JA,-nc);
   matchSortedIdxU(FL,JA,JB,Ja,Jb,-1,lCM);

   u.init(JB.dim1).set(-1);
   for (i=0; i<Jb.len; ++i) { u[Jb[i]]=Ja[i]; }
   for (i=0; i<gjb.len; ++i) { gjb[i]=u[gjb[i]]; }


#ifdef WB_SPARSE_CLOCK
   sp3.done();
   Wb::UseClock sp4(&wbc_sparse_cnt4);
#endif

   gKC.init(gJK.numGroups()); u.init(gIA.numGroups());

   while (gJK.groupIter(kgb,kb,lb)) {
      for (jb=kb; jb<lb; ++jb) { jga=gjb[jb];
          if (!gIJ.getGroup(jga,ka,la)) continue;
          for (ja=ka; ja<la; ++ja) { ic=gia[ja];
             if (u[ic]<kgb+1) {
                u[ic]=kgb+1;
                ++nnzc;
             }
          }
      }
      gKC.init_g(kgb,nnzc);
   }

   if (!nnzc) {
      if (C0.isEmpty()) C.init(sC); else C0.save2(C);
      return C;
   }

   w.init(u.len);
   ii.init(u.len);
   u.set(0); l=0;

   C.init(sC,nnzc);

   while (gJK.groupIter(kgb,kb,lb)) {

      sortflag=((lb-kb)>SPIDX_T(ceil(0.2*u.len))); is=0;

      for (jb=kb; jb<lb; ++jb) { jga=gjb[jb];
          if (!gIJ.getGroup(jga,ka,la)) continue;
          bjk=B.D[pb[jb]];

          for (ja=ka; ja<la; ++ja) { ic=gia[ja];
             if (u[ic]<kgb+1) { u[ic]=kgb+1;
                if (sortflag) ii[is++]=ic;
                w[ic]=bjk * D[pa[ja]];
             }
             else {
                w[ic] += bjk * D[pa[ja]];
             }
          }
      }


      if (sortflag) {
         if (is) {
            ii.len=is; ii.Sort();
            for (i=0; i<ii.len; ++i) { ic=ii[i];
               C.IDX.recSetP(l,IA.rec(ic),r1,JK.ref(kb,nc),r2);
               C.D[l++]=w[ic];
            }
            ii.len=u.len;
         }
      }
      else {
         for (ic=0; ic<u.len; ++ic) {
            if (u[ic]==kgb+1) {
               C.IDX.recSetP(l,IA.rec(ic),r1,JK.ref(kb,nc),r2);
               C.D[l++]=w[ic];
            }
         }
      }
   }

   if (l!=C.D.len) wblog(FL,"ERR %s() %d/%d/%d !??",FCT,l,C.D.len,nnzc);

#ifdef WB_SPARSE_CLOCK
   sp4.done();
#endif

   C.Permute(pfinal); if (afac) C*=afac;
   if (!C0.isEmpty()) { C+=C0; }



   return C.Compress();
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::contract(
   const char *F, int L,
   unsigned k, const wbvector<TD> &B, wbsparray<TD> &C
) const {

   if (k>=SIZE.len) wblog(F_L,
      "ERR %s() index out of bounds (%d/%d)",FCT,k+1,SIZE.len);
   if (SIZE[k]!=B.len) wblog(F_L,
      "ERR %s() size mismatch (%s ; %d @ %d)",FCT,sizeStr().data,B.len,k+1);
   if (IDX.dim1!=D.len || IDX.dim2!=SIZE.len) wblog(FL,
      "ERR %s() size inconsistency (%dx%d <> %dx%d) !?",
      FCT,IDX.dim1,IDX.dim2,D.len,SIZE.len); 

   if (IDX.dim2>1) {
      const SPIDX_T *idx=IDX.data+k;
      const TD *b=B.data; TD *c; wbindex Ik;

      Ik.Index_ex(IDX.dim2,k); C.init();
      IDX.getCols(Ik,C.IDX); C.D=D; c=C.D.data;
      SIZE.select(Ik,C.SIZE);

      for (SPIDX_T n=IDX.dim2, i=0; i<D.len; ++i, idx+=n) { 
         if ((*idx)>=B.len) wblog(F_L,
            "ERR %s() index out of bounds (%d/%d)",FCT,(*idx)+1,B.len);
         c[i]*=b[*idx];
      }

      C.Compress();
   }
   else if (IDX.dim2==1) {
      if (int(k)>=0 && k) wblog(FL,
         "ERR %s() ctrIndex out of bounds (%d/1)",FCT,k+1);
      C.initScalar();

      const SPIDX_T *idx=IDX.data;
      const TD *a=D.data, *b=B.data;
      TD &c=C.D[0]; c=0;

      for (SPIDX_T i=0; i<D.len; ++i) { 
         if (idx[i]>=B.len) wblog(F_L,
            "ERR %s() index out of bounds (%d/%d)",FCT,idx[i]+1,B.len);
         c+=a[i]*b[idx[i]];
      }
   }
   else {
      info(FL); wblog(FL,"ERR %s()",FCT);
   }

   return C;
};


template <class TD>
wbarray<TD>& wbsparray<TD>::contract(
   const char* F, int L,  const wbindex &ica,
   const wbsparray<TD> &B, const wbindex &icb, wbarray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
) const {

   wbarray<TD> Af,Bf,X; 

   if (cfac && !C.isEmpty()) { C.save2(X); X*=cfac; }

   this->toFull(Af); B.toFull(Bf);
   Af.contract(FL,ica,Bf,icb,C,pfinal,afac);

   if (!X.isEmpty()) {
      if (!X.hasSameSize(C)) wblog(FL,
         "ERR %s() severe size mismatch (%s <> %s; %g)",
         FCT, C.sizeStr().data, X.sizeStr().data, cfac);
      C+=X;
   }
   return C;
};


template <class TD>
wbarray<TD>& contract(const char* F, int L,  
   const wbarray<TD>   &A, const wbindex &ic1,
   const wbsparray<TD> &B, const wbindex &ic2, wbarray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
){
   unsigned m=ic1.len, r1=A.rank(), r2=B.rank(FL); int e=0;
   wbarray<TD> X;

   if ((void*)(&A)==(void*)(&C)) wblog(F,L,
      "ERR sparse::%s() must have distinct target space",FCT);

   wbvector<SPIDX_T> sC; if (cfac) sC=C.SIZE;
   if ((e=checkContract(F_L,A.SIZE,ic1,B.SIZE,ic2, &pfinal, &sC))) {
      if (e<0) { C.init(sC); if (sC.len) return C; }
      else wblog(FL,"ERR %s()",FCT);
   }

   if (afac==0 || !B.D.len) {
      if (cfac!=0) { C*=cfac; } else { C.init(sC.Permute(pfinal)); }
      return C;
   }

   if (cfac!=0) { C*=cfac; C.save2(X); }
   C.init(sC);

   wbvector<SPIDX_T> s1,D2;
   wbMatrix<SPIDX_T> I2;
   wbperm Pa,P2;

   A.SIZE.getI(ic1,s1); Pa.init2End(ic1,r1);
   B.IDX.cols2End(ic2,I2).groupRecs(P2,D2,r2-m);

   wbindex Ia(r1),Ic(sC.len);
   SPIDX_T i, ig,d,k,l, k1=r1-ic1.len, k2=r2-ic2.len;
   SPIDX_T *i2, *ia=Ia.data, *ic=Ic.data, *pa=Pa.data;
   wbIndex K1(s1);

   for (l=ig=0; ig<D2.len; ++ig) { d=D2[ig];
   for (k=0; k<d; ++k,++l) {
       TD b=B.D.data[P2[l]]; if (b==0) continue; if (afac!=1) b*=afac;

       i2=I2.ref(l); if (k==0) {
       for (i=0; i<k2; ++i) ic[k1+i]=i2[i]; }
       for (i=0; i<m; ++i) ia[pa[k1+i]]=i2[k2+i];

       for (K1.reset(); ++K1;) {
          for (i=0; i<k1; ++i) { ia[pa[i]]=ic[i]=K1.data[i]; }
          if (b== 1) { C(Ic)+=A(Ia); } else
          if (b==-1) { C(Ic)-=A(Ia); } else { C(Ic)+=b*A(Ia); }
       }
   }}

   C.Permute(pfinal); if (!X.isEmpty()) C+=X;
   return C;
};


template <class TD>
wbarray<TD>& contract(const char* F, int L,  
   const wbsparray<TD> &A, const wbindex &ic1,
   const wbarray<TD>   &B, const wbindex &ic2, wbarray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
){
   unsigned m=ic1.len, r1=A.rank(FL), r2=B.rank(); int e=0;
   wbarray<TD> X;

   if ((void*)(&B)==(void*)(&C)) wblog(F,L,
      "ERR sparse::%s() must have distinct target space",FCT);

   wbvector<SPIDX_T> sC; if (cfac) sC=C.SIZE;
   if ((e=checkContract(F_L,A.SIZE,ic1,B.SIZE,ic2, &pfinal, &sC))) {
      if (e<0) { C.init(sC); if (sC.len) return C; }
      else wblog(FL,"ERR %s()",FCT);
   }

   if (afac==0 || !A.D.len) {
      if (cfac!=0) { C*=cfac; } else { C.init(sC.Permute(pfinal)); }
      return C;
   }

   if (cfac!=0) { C*=cfac; C.save2(X); }
   C.init(sC);

   wbvector<SPIDX_T> s2,D1;
   wbMatrix<SPIDX_T> I1;
   wbperm Pb,P1;

   A.IDX.cols2End(ic1,I1).groupRecs(P1,D1,r1-m);
   B.SIZE.getI(ic2,s2); Pb.init2End(ic2,r2);

   wbindex Ib(r2),Ic(sC.len);
   SPIDX_T i, ig,d,k,l, k1=r1-ic1.len, k2=r2-ic2.len, *i1,
      *ib=Ib.data, *ic=Ic.data, *pb=Pb.data;
   wbIndex K2(s2);

   for (l=ig=0; ig<D1.len; ++ig) { d=D1[ig];
   for (k=0; k<d; ++k,++l) {
       TD a=A.D.data[P1[l]]; if (a==0) continue; if (afac!=1) a*=afac;

       i1=I1.ref(l); if (k==0) {
       for (i=0; i<k1; ++i) ic[i]=i1[i]; }
       for (i=0; i<m; ++i) ib[pb[k2+i]]=i1[k1+i];

       for (K2.reset(); ++K2;) {
          for (i=0; i<k2; ++i) { ib[pb[i]]=ic[k1+i]=K2.data[i]; }
          if (a== 1) { C(Ic)+=B(Ib); } else
          if (a==-1) { C(Ic)-=B(Ib); } else { C(Ic)+=a*B(Ib); }
       }
   }}

   C.Permute(pfinal); if (!X.isEmpty()) C+=X;
   return C;
};


template <class TD>
TD wbVMatVprod(
  const wbsparray<TD> &v1, const wbsparray<TD> &M, const wbsparray<TD> &v2
){
   wbsparray<TD> x,X; wbindex i1(1), i2(1);
   i1[0]=0; i2[0]=1;


   M.contract(FL,i2,v2,i1,X);
   v1.contract(FL,i1,X,i1,x); 

   x.checkSize(FL);
   if (x.D.len && !x.isScalar()) {
      MXPut(FL,"q").add(M,"M").add(v1,"v1").add(v2,"v2")
        .add(X,"X").add(x,"x");
      wblog(FL,"ERR %s() got %s",FCT,x.sizeStr().data);
   }

   return (x.D.len ? x[0] : 0.);
};

template <class TD>
wbsparray<TD>& wbMatProd(
  const wbsparray<TD> &A, const wbsparray<TD> &B, wbsparray<TD> &C,
  char aflag0, char bflag0, TD afac, TD cfac
){
   if (&C==&A || &C==&B) {
      wbsparray<TD> X; if (cfac) X=C;
      wbMatProd(A,B,X,aflag0,bflag0,afac,cfac);
      return X.save2(C);
   }
   if (A.rank()!=2 || B.rank()!=2) wblog(FL,
      "ERR %s() expects matrices (%s; %s)",
      FCT,A.sizeStr().data, B.sizeStr().data
   );

   opFlags<TD> aflag(aflag0), bflag(bflag0);
   const wbsparray<TD>
      &AX = aflag.applyConjOrRef(FL,A),
      &BX = bflag.applyConjOrRef(FL,B);
   wbindex ia(1),ib(1);

   ia[0]=(aflag.trans() ? 0:1);
   ib[0]=(bflag.trans() ? 1:0);


   AX.contract(FL,ia,BX,ib,C, wbperm(), afac, cfac);


   if (&AX!=&A) { delete &AX; }
   if (&BX!=&B) { delete &BX; }; return C;
};


template <class TD>
wbsparray<TD>& wbsparray<TD>::comm(const char *F, int L,
  const wbsparray<TD> &B, wbsparray<TD> &C, char aflag, char bflag
) const {

   if (&C==this || &C==&B) {
      wbsparray<TD> X; comm(F,L,B,X,aflag,bflag);
      return X.save2(C);
   }

   SPIDX_T d=-1;
   if (!isSMatrix(FL,&d) || !B.isSMatrix(FL,&d)) wblog(F_L,
      "ERR %s() invalid operators for [A,B] (%s; %s)", FCT,
      sizeStr().data, B.sizeStr().data
   );

   wbMatProd(*this,B,C,aflag,bflag);
   wbMatProd(B,*this,C,bflag,aflag,TD(-1.),TD(1.));

   return C;
};


#endif

