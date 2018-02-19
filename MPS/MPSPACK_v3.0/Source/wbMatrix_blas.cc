#ifndef __WB_BLAS_CC__
#define __WB_BLAS_CC__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //
// wrapper routine

template<class TA, class TB, class TC>
void MMDIAG(
    const wbMatrix<TA> &A,
    const wbMatrix<TB> &B, wbMatrix<TC> &C,
    char aflag='N', char bflag='N',
    const TA afac=1., const TC cfac=0.
);


template<class TA, class TB, class TC>
wbMatrix<TC>& MMDIAG(
    const wbvector<TA> &A,
    const wbMatrix<TB> &B, wbMatrix<TC> &C,
    char aflag='N', char bflag='N',
    const TA afac=1., const TC cfac=0.
);


template<class TA, class TB, class TC>
wbMatrix<TC>& MMDIAG(
    const wbMatrix<TA> &A,
    const wbvector<TB> &B, wbMatrix<TC> &C,
    char aflag='N', char bflag='N',
    const TA afac=1., const TC cfac=0.
);

template<class T>
void DZGEMM_VEC(
    const wbMatrix<T> &A,
    const wbvector<T> &B,
    wbvector<T> &C, int m, int n,
    char aflag='N', const T afac=1., const T cfac=0.
);


template<class TD>
void DZGEMM(
    const wbMatrix<TD> &A,
    const wbMatrix<TD> &B, wbMatrix<TD> &C,
    const unsigned k, 
    char aflag='N', char bflag='N',
    const TD afac=1., const TD cfac=0.
);



template<>
inline void DZGEMM<double>(
    const wbMatrix<double> &A,
    const wbMatrix<double> &B, wbMatrix<double> &C,
    const unsigned k, 
    char aflag, char bflag,
    const double afac, const double cfac
){

#ifdef WB_CLOCK
   wbc_dgemm.resume();

   if (aflag=='N')
        stat_dgemm.account(A.dim1*A.dim2*(A.dim2==B.dim1 ? B.dim2 : B.dim1));
   else stat_dgemm.account(A.dim1*A.dim2*(A.dim1==B.dim1 ? B.dim2 : B.dim1));
#endif


    dgemm (
      bflag, aflag, (int)C.dim2, (int)C.dim1, (int)k,
      afac, B.data, (int)B.dim2, A.data, (int)A.dim2,
      cfac, C.data, (int)C.dim2
    );

#ifdef WB_CLOCK
   wbc_dgemm.stop();
#endif
}

template<>
inline void DZGEMM<wbcomplex>(
    const wbMatrix<wbcomplex> &A,
    const wbMatrix<wbcomplex> &B, wbMatrix<wbcomplex> &C,
    const unsigned k, 
    char aflag, char bflag,
    const wbcomplex afac, const wbcomplex cfac
){
#ifdef WB_CLOCK
    wbc_zgemm.resume();

    if (bflag=='N')
         stat_dgemm.account(A.dim1*A.dim2*(A.dim2==B.dim1 ? B.dim2 : B.dim1));
    else stat_dgemm.account(A.dim1*A.dim2*(A.dim1==B.dim1 ? B.dim2 : B.dim1));
#endif

    zgemm (
      bflag, aflag, (int)C.dim2, (int)C.dim1, (int)k,
      afac, B.data, (int)B.dim2, A.data, (int)A.dim2,
      cfac, C.data, (int)C.dim2
    );

#ifdef WB_CLOCK
    if (wbc_zgemm.stop()) {
       wbc_zgeNX.flag+=2;
       wbc_zgeNX.tcpu+= (long unsigned)
       (double(A.dim1*A.dim2*B.dim2)*(4E-9*CLOCKS_PER_SEC)+1.5);
    }
    else {
       wbc_zgeNN.flag+=2;
       wbc_zgeNN.tcpu+= (long unsigned)
       (double(A.dim1*A.dim2*B.dim2)*(4E-9*CLOCKS_PER_SEC)+1.5);
    }
#endif
}


inline void wbMatProd(
    const wbMatrix<double> &A,
    const wbMatrix<wbcomplex> &B, wbMatrix<wbcomplex> &C,
    char aflag='N', char bflag='N',
    const double afac=1., const double cfac=0.,
    const char i00flag=0
){
    wbMatrix<double> R,I, rC, iC;

    B.getReal(R);
    B.getImag(I);

    if (bflag=='C') I.flipSign();

    if (cfac!=0.) { C.getReal(rC); C.getImag(iC); }

    wbMatProd(A, R, rC, aflag, bflag, afac, cfac, i00flag);
    wbMatProd(A, I, iC, aflag, bflag, afac, cfac, i00flag);

    C.set(rC,iC);
}


inline void wbMatProd(
    const wbMatrix<wbcomplex> &A,
    const wbMatrix<double> &B, wbMatrix<wbcomplex> &C,
    char aflag='N', char bflag='N',
    const double afac=1., const double cfac=0.,
    const char i00flag=0
){
    wbMatrix<double> R,I, rC, iC;

    A.getReal(R);
    A.getImag(I);
    
    if (aflag=='C') I.flipSign();

    if (cfac!=0.) { C.getReal(rC); C.getImag(iC); }

    wbMatProd(R, B, rC, aflag, bflag, afac, cfac, i00flag);
    wbMatProd(I, B, iC, aflag, bflag, afac, cfac, i00flag);

    C.set(rC,iC);
}




template<class T>
inline void wbMatProd(
    const wbMatrix<T> &A,
    const wbMatrix<T> &B, wbMatrix<T> &C,
    char aflag, char bflag,
    const T afac, const T cfac,
    const char i00flag
){
    unsigned a1=A.dim1, a2=A.dim2, b1=B.dim1, b2=B.dim2;
    const char flags[]="NTCntc";

#ifdef WB_CLOCK
    wbc_matprod.resume();
#endif

    if (!strchr(flags,aflag) || !strchr(flags,bflag)) wblog(FL,
    "ERR MatProd - invalid flags %c<%d>, %c<%d>", aflag,aflag,bflag,bflag);


    if (aflag!='N') SWAP(a1,a2);
    if (bflag!='N') SWAP(b1,b2);

    if (a2!=b1) {
        wblog(FL, "ERR MatProd - dimension mismatch: (%d,%d) * (%d,%d) ?",
        a1,a2,b1,b2); return;
    }

    if (cfac!=0) {
        if (C.data==NULL && i00flag) {
            if (cfac!=1.) wblog(FL,
            "WRN C = A*B + c*[] with c=%s !?", toStr(cfac).data);

            C.init(a1,b2);
        }
        else
        if (C.dim1!=a1 || C.dim2!=b2) {
            wblog(FL, "ERR MatProd - dimension mismatch: C=(%d,%d) =? (%d,%d).",
            C.dim1, C.dim2, a1, b2);
            return;
        }
    }

    if (a1==0 || b2==0) {
        if (cfac==0) C.init(a1,b2);
        return;
    }

    if (a2==0) {        
        wblog(FL, "ERR MatProd Cannot multiply (%dx%d)*(%dx%d)!",
        a1,a2,b1,b2); return;
    }

    if (&C==&A || &C==&B) { 
        wblog(FL, "ERR I/O spaces must be distinct!\n[%7lX %7lX %7lX]",
        &A, &B, &C); return;
    }

    if (cfac==0)
    C.init (a1,b2);

    if (A.isdiag || B.isdiag)
         MMDIAG(A, B, C,     aflag, bflag, afac, cfac);
    else DZGEMM(A, B, C, a2, aflag, bflag, afac, cfac);

#ifdef WB_CLOCK
    wbc_matprod.stop();
#endif
}


template<class T>
inline void wbMatProd(
    const wbMatrix<T> &A,
    const wbvector<T> &B, wbvector<T> &C,
    char aflag, const T afac, const T cfac,
    const char i00flag=0
){
    unsigned m=A.dim1, n=A.dim2;
    if (aflag!='N') SWAP(m,n);

    if (B.len!=n) wblog(FL,
    "ERR Severe size mismatch Ax=b (%dx%d, %d)",m,n,B.len);

    if (!C.data) {
       if (cfac!=0. && i00flag) wblog(FL,
          "WRN C = A*B + c*[] with c=%s !?",toStr(cfac).data);
       C.init(m);
    }
    else if (C.len!=m) {
       if (cfac==0.) C.init(m); else wblog(FL,
       "ERR Severe size mismatch Ax=b (%d, %d)", m, C.len);
    }

    DZGEMM_VEC(A,B,C,m,n,aflag,afac,cfac);
}


template<>
inline void DZGEMM_VEC(
    const wbMatrix<double> &A,
    const wbvector<double> &B,
    wbvector<double> &C, int m, int n,
    char aflag, const double afac, const double cfac
){
    dgemm (
      'N', aflag, 1, m, n, afac,
       B.data, 1, A.data, A.dim2, cfac,
       C.data, 1
    );
}

template<>
inline void DZGEMM_VEC(
    const wbMatrix<wbcomplex> &A,
    const wbvector<wbcomplex> &B,
    wbvector<wbcomplex> &C, int m, int n,
    char aflag, const wbcomplex afac, const wbcomplex cfac
){
    zgemm (
      'N', aflag, 1, m, n, afac,
       B.data, 1, A.data, A.dim2, cfac,
       C.data, 1
    );
}


inline void DZGEMM_aux(
    const wbMatrix<double> &A,
    const wbMatrix<double> &B, wbMatrix<double> &C,
    const unsigned k, 
    char aflag='N', char bflag='N',
    const double afac=1., const double cfac=0.
){
    if (A.isdiag || B.isdiag)
         MMDIAG(A, B, C,    aflag, bflag, afac, cfac);
    else DZGEMM(A, B, C, k, aflag, bflag, afac, cfac);
}


template<class TA, class TB, class TC>
inline void MMDIAG(
    const wbMatrix<TA> &A,
    const wbMatrix<TB> &B, wbMatrix<TC> &C,
    char aflag, char bflag,
    const TA afac, const TC cfac
){
    unsigned i, s=C.dim1*C.dim2, dmax;

    if (!A.isdiag && !B.isdiag) wblog(FL,
    "ERR Wrong call - elements not diagonal.");
    
    if (A.isdiag && A.dim1>1 && A.dim2>1)  
    if (A(1,0)!=0. || A(0,1)!=0.) wblog(FL,
    "ERR Input A is NOT diag even though isdiag=%d", A.isdiag);

    if (B.isdiag && B.dim1>1 && B.dim2>1)  
    if (B(1,0)!=0. || B(0,1)!=0.) wblog(FL,
    "ERR Input B is NOT diag even though isdiag=%d", B.isdiag);


    if (cfac!=1.) for (i=0; i<s; i++) C.data[i]*=cfac;

    if (A.isdiag && !B.isdiag) {
        MMDIAG(A.getDiag(), B, C, aflag, bflag, afac);
    }
    else
    if ((!A.isdiag) && B.isdiag) {
        MMDIAG(A, B.getDiag(), C, aflag, bflag, afac);
    }
    else { 
        wbvector<TA> adiag;
        wbvector<TB> bdiag;
        A.getDiag(adiag); B.getDiag(bdiag);

        if (aflag=='C')
        for (i=0; i<adiag.len; i++) adiag[i]=conj(adiag[i]);

        if (bflag=='C')
        for (i=0; i<bdiag.len; i++) bdiag[i]=conj(bdiag[i]);

        dmax=MIN( MIN(A.dim1,A.dim2), MIN(B.dim1,B.dim2) );
        for (i=0; i<dmax; i++) C(i,i) += afac * adiag[i] * bdiag[i];
    }

    if (C.isdiag) if (!A.isdiag || !B.isdiag) C.isdiag=0;
}


template<class TA, class TB, class TC>
inline wbMatrix<TC>& MMDIAG(
   const wbvector<TA> &A,
   const wbMatrix<TB> &B, wbMatrix<TC> &C,
   char aflag, char bflag,
   const TA afac, const TC cfac
){
   unsigned i,j, m=B.dim1, n=B.dim2;
   char acc,btr;
   TA dbl;


   if (bflag=='C' && typeid(TB)!=typeid(wbcomplex)) bflag='T';

   if (!strchr("NTC",aflag) || !strchr("NTC",bflag)) wblog(FL,
   "ERR Invalid flag(s) %c<%d>, %c<%d>",aflag,aflag,bflag,bflag);

   acc=(aflag=='C');
   btr=(bflag!='N'); if (btr) SWAP(m,n);

   if (A.len!=m) wblog(FL,
   "ERR Size mismatch (diag(%d) * %dx%d)", A.len,m,n);
   if (&C==&B) wblog(FL,"ERR Output space C equal to input space.");

   if (C.isEmpty() || (cfac==0. && (C.dim1!=m || C.dim2!=n)))
      C.init(m,n);
   else {
      if (C.dim1!=m || C.dim2!=n) wblog(FL,
      "ERR Size mismatch (%dx%d; %dx%d)",C.dim1,C.dim2, m,n);
      if (cfac!=1) {
         if (cfac==0.) C.reset();
         else C*=cfac;
      }
   }

   for (i=0; i<C.dim1; i++) {
      dbl=A[i];   if (acc) dbl=conj(dbl);
      dbl*=afac;  if (dbl==0.) continue;

      switch (bflag) {
        case 'N': for (j=0; j<C.dim2; j++) C(i,j)+=(dbl*B(i,j)); break;
        case 'T': for (j=0; j<C.dim2; j++) C(i,j)+=(dbl*B(j,i)); break;
        case 'C': for (j=0; j<C.dim2; j++) C(i,j)+=(dbl*conj(B(j,i)));
      }
   }

   return C;
}


template<class TA, class TB, class TC>
inline wbMatrix<TC>& MMDIAG(
   const wbMatrix<TA> &A,
   const wbvector<TB> &B, wbMatrix<TC> &C,
   char aflag, char bflag,
   const TA afac, const TC cfac
){
   unsigned i,j, m=A.dim1, n=A.dim2;
   char atr,bcc;
   TB dbl;


   if (aflag=='C' && typeid(TA)!=typeid(wbcomplex)) aflag='T';

   if (!strchr("NTC",aflag) || !strchr("NTC",bflag)) wblog(FL,
   "ERR Invalid flag(s) %c<%d>, %c<%d>",aflag,aflag,bflag,bflag);

   atr=(aflag!='N'); if (atr) SWAP(m,n);
   bcc=(bflag=='C');

   if (B.len!=n) wblog(FL,
      "ERR Size mismatch (%dx%d * diag(%d))", m,n, B.len);
   if (&C==&A) wblog(FL,"ERR Output space C equal to input space.");

   if (C.isEmpty() || (cfac==0. && (C.dim1!=m || C.dim2!=n)))
      C.init(m,n);
   else {
      if (C.dim1!=m || C.dim2!=n) wblog(FL,
      "ERR Size mismatch (%dx%d; %dx%d)",C.dim1,C.dim2, m,n);
      if (cfac!=1) {
         if (cfac==0.) C.reset();
         else C*=cfac;
      }
   }

   for (j=0; j<C.dim2; j++) {
      dbl=B[j];   if (bcc) dbl=conj(dbl);
      dbl*=afac;  if (dbl==0.) continue;

      switch (aflag) {
        case 'N': for (i=0; i<C.dim1; i++) C(i,j)+=(dbl*A(i,j)); break;
        case 'T': for (i=0; i<C.dim1; i++) C(i,j)+=(dbl*A(j,i)); break;
        case 'C': for (i=0; i<C.dim1; i++) C(i,j)+=(dbl*conj(A(j,i)));
      }
   }

   return C;
}


void wbCMatProd(
    const wbCMat &A0,
    const wbCMat &B0,
    wbCMat &C,
    char ia, char ib,
    const wbcomplex &afac, const wbcomplex &cfac,
    const char i00flag
){
    unsigned a1=A0.R.dim1, a2=A0.R.dim2, b1=B0.R.dim1, b2=B0.R.dim2;

    char FLAGS[]="NTCc", aflag, bflag, ta, tb, ca, cb;

    wbCMat XX,CX, *A, *B;

    if (!strchr(FLAGS,ia) || !strchr(FLAGS,ib)) wblog(FL,
    "ERR Invalid flag(s) %c<%d>, %c<%d>", ia,ia,ib,ib);

    ta = (ia=='T' || ia=='C');
    tb = (ib=='T' || ia=='C');

    if (ta) SWAP(a1,a2);
    if (tb) SWAP(b1,b2);

    if (a2!=b1 || a2==0) wblog(FL, 
    "ERR Invalid dimensions / mismatch: (%d,%d) * (%d,%d) ?", a1,a2,b1,b2);

    if (a1==0 || b2==0) {
        if (cfac==wbcomplex(0)) C.init(a1,b2); else
        if (C.R.dim1!=a1 || C.R.dim2!=b2)
            wblog(FL,"ERR Dimension mismatch C=(%d,%d) ? (%d,%d).",
            C.R.dim1, C.R.dim2, a1, b2);
        return;
    }

    if (&A0==&C || &B0==&C) CX=C; 
    
    A = &A0==&C ? &CX : (wbCMat*)&A0; 
    B = &B0==&C ? &CX : (wbCMat*)&B0;

    if (afac!=1.) {
       if ((A->totSize())<(B->totSize())) {
          if (A==&A0 || B==&CX) { XX=A0; A=&XX; }
          if (ia=='C') { A->conj(); ia='T'; } else
          if (ia=='c') { A->conj(); ia='N'; }

          (*A)*=afac;
       }
       else {
          if (B==&B0 || A==&CX) { XX=B0; B=&XX; }
          if (ib=='C') { B->conj(); ib='T'; } else
          if (ib=='c') { B->conj(); ib='N'; }

          (*B)*=afac;
       }
    }

    ca = (ia=='c' || ia=='C');  
    cb = (ib=='c' || ib=='C');  

    if (ia=='c') ia='N';  aflag=ia;
    if (ib=='c') ib='N';  bflag=ib;

    if (cfac!=wbcomplex(0)) {
        if (C.R.data==NULL && i00flag) {
            if (cfac!=0.) wblog(FL,
               "WRN C = A*B + c*[] with c=%g+i%g !?",cfac.r,cfac.i);
            C.init(a1,b2);
        }
        else if (C.R.dim1!=a1 || C.R.dim2!=b2)
            wblog(FL, "ERR Dimension mismatch C=(%d,%d) ? (%d,%d).",
            C.R.dim1, C.R.dim2, a1, b2);

        if (cfac!=1.)
        C*=cfac;
    }
    else
    C.init(a1,b2);

    DZGEMM_aux(A->R, B->R, C.R, a2, aflag, bflag);

    if (A->I.data && B->I.data)
    DZGEMM_aux(A->I, B->I, C.R, a2, aflag, bflag, ca ^ cb ? +1. : -1.);

    if (C.I.data==NULL && (A->I.data || B->I.data))
    C.I.init(C.R.dim1, C.R.dim2);

    if (A->I.data)
    DZGEMM_aux(A->I, B->R, C.I, a2, aflag, bflag, ca ? -1. : +1.);

    if (B->I.data)
    DZGEMM_aux(A->R, B->I, C.I, a2, aflag, bflag, cb ? -1. : +1.);

    return;
}


inline void WbEigenSymmetric (
    const wbMatrix<double> &M, wbMatrix<double> &V, wbvector<double> &E
){

    unsigned n=M.dim1; pINT i;
    wbvector<double> aux;

    if (M.dim1!=M.dim2) wblog(FL,
       "ERR WbEigenSymmetric requires square matrix (%d,%d).",
        M.dim1, M.dim2);

    V=M; E.init(n); if (n==0) return;
    aux.init(3*n-1);

    dsyev (
       'V',
       'L',
        n,
        V.data,
        n,
        E.data,
        aux.data,
        aux.len,
        i
    );


    if (i) wblog(FL,"ERR DSYEV returned %d.", i);
}

inline void WbEigenSymmetric (
    const wbMatrix<wbcomplex> &M, wbMatrix<wbcomplex> &V, wbvector<double> &E
){
    unsigned n=M.dim1; pINT i; double dz=0;
    wbvector<wbcomplex> auxz;
    wbvector<double> auxd;

    if (M.dim1!=M.dim2) wblog(FL,
       "ERR WbEigenSymmetric requires square matrix (%d,%d).",
        M.dim1, M.dim2);
    if (!M.isHConj(1E-12,&dz)) wblog(FL,
       "ERR WbEigenSymmetric got non-hermitian matrix (%.3g)",dz
    );

    V=M; E.init(n); if (n==0) return;
    auxz.init(2*n-1);
    auxd.init(3*n-2);

    zheev (
       'V',
       'L',
        n,
        V.data,
        n,
        E.data,
        auxz.data,
        auxz.len,
        auxd.data,
        i
    );


    if (i) {
       wblog(FL,"ERR ZHEEV returned %d.", i);
    }

    V.Conj();
};


template<class T>
inline void GESVD(
   wbMatrix<T> &A,
   wbMatrix<T> &U,
   wbvector<double> &S,
   wbMatrix<T> &VT
){
   wblog(FL,"ERR GESVD not defined for datatype %s",
   getName(typeid(T)));
};


template<>
inline void GESVD(
   wbMatrix<double> &A,
   wbMatrix<double> &U,
   wbvector<double> &S,
   wbMatrix<double> &VT
){
   unsigned M=A.dim1, N=A.dim2, K=S.len;
   wbvector<double> W;
   double lwrk;
   pINT i;
 

#ifdef WB_CLOCK
    wbc_dgesvd.resume();
#endif

   dgesvd(
     'S','S', N, M, A.data, N, S.data, VT.data, N, U.data, K,
      &lwrk, -1,
      i
   );

   W.init((unsigned)lwrk);

   dgesvd(
     'S','S', N, M, A.data, N, S.data,
      VT.data, N,
      U.data, K,
      W.data, W.len, i
   );

#ifdef WB_CLOCK
    wbc_dgesvd.stop();
#endif

   if (i) wblog(FL,"ERR DGESVD returned %d.", i);
};


template<>
inline void GESVD(
   wbMatrix<wbcomplex> &A,
   wbMatrix<wbcomplex> &U,
   wbvector<double> &S,
   wbMatrix<wbcomplex> &VT
){
   unsigned M=A.dim1, N=A.dim2, K=S.len;
   wbvector<wbcomplex> W, RW(5*MIN(M,N));
   wbcomplex lwrk;
   pINT i;


#ifdef WB_CLOCK
    wbc_zgesvd.resume();
#endif

   zgesvd(
     'S','S', N, M, A.data, N, S.data, VT.data, N, U.data, K,
      &lwrk, -1,
      RW.data, i
   );

   W.init((unsigned)lwrk.r);

   zgesvd(
     'S','S', N, M, A.data, N, S.data,
      VT.data, N, U.data, K,
      W.data, W.len, RW.data, i
   );

#ifdef WB_CLOCK
    wbc_zgesvd.stop();
#endif

   if (i) wblog(FL,"ERR ZGESVD returned %d.", i);
};


template<class T>
inline void wbSVD(
   const wbMatrix<T> &A0,
   wbMatrix<T> &U,
   wbvector<double> &S,
   wbMatrix<T> &VT
){
   wbMatrix<T> A(A0);

   unsigned M=A.dim1, N=A.dim2, K=MIN(M,N);
 
   U.init(M,K); S.init(K); VT.init(K,N);
   if (!M || !N) return;

   GESVD(A,U,S,VT);
};


#endif

